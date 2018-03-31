#include "cv_integrator.h"

#include "bssrdf.h"
#include "camera.h"
#include "scene.h"
#include "interaction.h"
#include "progressreporter.h"
#include "paramset.h"

#include "cv_pixel.h"
#include "cv_film.h"
#include "dualmat.h"

namespace pbrt {

	// Variables for logging.
	STAT_COUNTER("Integrator/Camera rays traced", nCameraRays);
	STAT_PERCENT("Integrator/Zero-radiance paths", zeroRadiancePaths, totalPaths);
	STAT_INT_DISTRIBUTION("Integrator/Path length", pathLength);

	namespace {
	}  // anonymous namespace

	CvPathIntegrator::CvPathIntegrator(int maxDepth,
									   std::shared_ptr<const Camera> camera,
									   std::shared_ptr<Sampler> sampler,
									   const Bounds2i &pixelBounds, Float rrThreshold,
									   const std::string &lightSampleStrategy)
		: PathIntegrator(maxDepth, camera, sampler, pixelBounds,
						 rrThreshold, lightSampleStrategy) {
	}

	void CvPathIntegrator::Render(const Scene &scene) {
		Preprocess(scene, *sampler);
		// Render image tiles in parallel

		// Check Film type and get CvFilm pointer
		//CHECK(typeid(camera->film) == typeid(CvFilm))
		//    << "Film type must be \"CvFilm\"";
		CvFilm *film = reinterpret_cast<CvFilm*>(camera->film);

		// Compute number of tiles, _nTiles_, to use for parallel rendering
		Bounds2i sampleBounds = film->GetSampleBounds();
		Vector2i sampleExtent = sampleBounds.Diagonal();
		const int tileSize = 16;
		Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
					   (sampleExtent.y + tileSize - 1) / tileSize);
		ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
		{
			ParallelFor2D([&](Point2i tile) {
				// Render section of image corresponding to _tile_

				// Allocate _MemoryArena_ for tile
				MemoryArena arena;

				// Get sampler instance for tile
				int seed = tile.y * nTiles.x + tile.x;
				std::unique_ptr<Sampler> tileSampler1 = sampler->Clone(seed);
				std::unique_ptr<Sampler> tileSampler2 = sampler->Clone(seed);

				// Compute sample bounds for tile
				int x0 = sampleBounds.pMin.x + tile.x * tileSize;
				int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
				int y0 = sampleBounds.pMin.y + tile.y * tileSize;
				int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
				Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
				LOG(INFO) << "Starting image tile " << tileBounds;

				// Get _FilmTile_ for tile
				std::unique_ptr<CvFilmTile> filmTile = film->GetCvFilmTile(tileBounds);

				// Loop over pixels in tile to render them
				for (Point2i pixel : tileBounds) {
					{
						ProfilePhase pp(Prof::StartPixel);
						tileSampler1->StartPixel(pixel);
						tileSampler2->StartPixel(pixel);
					}

					// Do this check after the StartPixel() call; this keeps
					// the usage of RNG values from (most) Samplers that use
					// RNGs consistent, which improves reproducability /
					// debugging.
					if (!InsideExclusive(pixel, pixelBounds))
						continue;

					do {
						// Initialize _CameraSample_ for current sample
						CameraSample cameraSample = tileSampler1->GetCameraSample(pixel);

						// Generate camera ray for current sample
						RayDifferential ray;
						Float rayWeight =
							camera->GenerateRayDifferential(cameraSample, &ray);
						ray.ScaleDifferentials(
							1 / std::sqrt((Float)tileSampler1->samplesPerPixel));
						++nCameraRays;

						// Evaluate radiance along camera ray
						CvDualPixel pixel;
						if (rayWeight > 0) {
							pixel = LiControlVariate(ray, scene, *tileSampler1, arena);
						}

						// TODO: Here should be reverted?
						//VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
						//    ray << " -> L = " << L;

						// Add camera ray's contribution to image
						filmTile->AddSample(cameraSample.pFilm, pixel, rayWeight);

						// Free _MemoryArena_ memory from computing image sample
						// value
						arena.Reset();
					} while (tileSampler1->StartNextSample() && tileSampler2->StartNextSample());
				}
				LOG(INFO) << "Finished image tile " << tileBounds;
				// Merge image tile into _Film_
				film->MergeFilmTile(std::move(filmTile));
				reporter.Update();
			}, nTiles);
			reporter.Done();
		}
		LOG(INFO) << "Rendering finished";
		film->WriteImage(1,sampler->samplesPerPixel);    
	}

	CvDualPixel CvPathIntegrator::LiControlVariate(const RayDifferential &r,
												   const Scene &scene, Sampler &sampler,
												   MemoryArena &arena, int depth) const {
		ProfilePhase p(Prof::SamplerIntegratorLi);
		Spectrum L1(0.f), L2(0.f);
		std::vector <Spectrum> betas(2 , 1.f);
		std::vector<Float> reciprocal_pdfs(2 , 1.f);
		RayDifferential ray(r);
		bool specularBounce = false;
		int bounces;
		// Added after book publication: etaScale tracks the accumulated effect
		// of radiance scaling due to rays passing through refractive
		// boundaries (see the derivation on p. 527 of the third edition). We
		// track this value in order to remove it from beta when we apply
		// Russian roulette; this is worthwhile, since it lets us sometimes
		// avoid terminating refracted rays that are about to be refracted back
		// out of a medium and thus have their beta value increased.
		Float etaScale = 1;

		for (bounces = 0;; ++bounces) {
			SurfaceInteraction isect;
			bool foundIntersection = scene.Intersect(ray, &isect); 
      
			if (foundIntersection) {
				Spectrum Le = isect.Le(-ray.d);
				if (!Le.IsBlack()) {
					L1 += betas[0] * Le;
					L2 += betas[1] * Le;
					break;
				}
			} else {
				for (const auto &light : scene.infiniteLights) {
					Spectrum Le = light->Le(ray);
					L1 += betas[0] * Le;
					L2 += betas[1] * Le;
				}
			}

			// Terminate path if ray escaped or _maxDepth_ was reached
			if (!foundIntersection || bounces >= maxDepth) break;

			// Compute scattering functions and skip over medium boundaries
			isect.ComputeScatteringFunctions(ray, arena, true);
			if (!isect.bsdf) {
				VLOG(2) << "Skipping intersection due to null bsdf";
				ray = isect.SpawnRay(ray.d);
				bounces--;
				continue;
			}

			const Distribution1D *distrib = lightDistribution->Lookup(isect.p);

			Vector3f wo = -ray.d;
			Vector3f wi;
			std::vector<Float> pdfs(2 , 1.f);
			BxDFType flag = BxDFType(0);
			std::vector<Spectrum> fs(2);

			//sample a new path direction
			//Compute for two material 0:after(F) 1:befpre(H)
			// Compute scattering functions and skip over medium boundaries
			DualMaterial::SetMaterialId(0);
			isect.ComputeScatteringFunctions(ray, arena, true);

			if (!isect.bsdf) {
				VLOG(2) << "Skipping intersection due to null bsdf";
				ray = isect.SpawnRay(ray.d);
				bounces--;
				continue;
			}

			fs[0] = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdfs[0],
										 BSDF_ALL, &flag);

			DualMaterial::SetMaterialId(1);
			isect.ComputeScatteringFunctions(ray, arena, true);
			fs[1] = isect.bsdf->f_pdf(wo, wi, &pdfs[1], BSDF_ALL);

			DualMaterial::SetMaterialId(0);
			isect.ComputeScatteringFunctions(ray, arena, true);

			VLOG(2) << "Sampled BSDF, f1 = " << fs[0] << ", pdf = " << pdfs[0];
			VLOG(2) << "Sampled BSDF, f2 = " << fs[1] << ", pdf = " << pdfs[1];
			if (fs[0].IsBlack() || fs[1].IsBlack() || pdfs[0] == 0.f || pdfs[1] == 0.f) {
				break;
			}

			Float G = AbsDot(wi, isect.shading.n);
			betas[0] *= fs[0] * G ;
			betas[1] *= fs[1] * G ;
			reciprocal_pdfs[0] /= pdfs[0];
			reciprocal_pdfs[1] /= pdfs[1];

			specularBounce = (flag & BSDF_SPECULAR) != 0;
			if ((flag & BSDF_SPECULAR) && (flag & BSDF_TRANSMISSION)) {
				Float eta = isect.bsdf->eta;
				// Update the term that tracks radiance scaling for refraction
				// depending on whether the ray is entering or leaving the
				// medium.
				etaScale *= (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
			}
		  ray = isect.SpawnRay(wi);

		  if (isect.bssrdf && (flag & BSDF_TRANSMISSION)) {
			  Error("Participating media is not currently supported!!");
		  }

		  // Possibly terminate the path with Russian roulette.
		  // Factor out radiance scaling due to refraction in rrBeta.
		  const Spectrum rrBeta = betas[0] * etaScale * reciprocal_pdfs[0];

		  if (rrBeta.MaxComponentValue() < rrThreshold + 1.f && bounces >  3) {
			  Float q = std::max((Float).05, 1.f - rrBeta.MaxComponentValue());
			  if (sampler.Get1D() < q) {
				  reciprocal_pdfs[0] /=  q;
				  reciprocal_pdfs[1] /=  q;
				  break;
			  }
			  reciprocal_pdfs[0] /= 1 - q;
			  reciprocal_pdfs[1] /= 1 - q;
			  DCHECK(!std::isinf(betas[0].y()));
		  }
		}
		ReportValue(pathLength, bounces);
		return CvDualPixel(L1, L2, reciprocal_pdfs[0]);
	}

	Integrator *CreateCvPathIntegrator(const ParamSet &params,
									   std::shared_ptr<Sampler> sampler,
									   std::shared_ptr<const Camera> camera) {
		int maxDepth = params.FindOneInt("maxdepth", 5);
		int np;
		const int *pb = params.FindInt("pixelbounds", &np);
		Bounds2i pixelBounds = camera->film->GetSampleBounds();
		if (pb) {
			if (np != 4)
				Error("Expected four values for \"pixelbounds\" parameter. Got %d.",np);
			else {
				pixelBounds = Intersect(pixelBounds,
										Bounds2i{{pb[0], pb[2]},
										{pb[1], pb[3]}});
				if (pixelBounds.Area() == 0)
					Error("Degenerate \"pixelbounds\" specified.");
			}
		}
		Float rrThreshold = params.FindOneFloat("rrthreshold", 1.);
		std::string lightStrategy = params.FindOneString("lightsamplestrategy", "spatial");
		return new CvPathIntegrator(maxDepth, camera, sampler, pixelBounds,
									rrThreshold, lightStrategy);
	}
}  // namespace pbrt
