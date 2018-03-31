#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CV_PATH_INTEGRATOR_H
#define PBRT_CV_PATH_INTEGRATOR_H

// integrators/path.h*
#include "integrators/path.h"

#include "pbrt.h"
#include "integrator.h"
#include "lightdistrib.h"

#include "cv_pixel.h"

namespace pbrt {

// PathIntegrator Declarations
class CvPathIntegrator : public PathIntegrator {
public:
    // PathIntegrator Public Methods
    CvPathIntegrator(int maxDepth, 
                     std::shared_ptr<const Camera> camera,
                     std::shared_ptr<Sampler> sampler,
                     const Bounds2i &pixelBounds, Float rrThreshold,
                     const std::string &lightSampleStrategy);

    void Render(const Scene &scene) final override;

    CvDualPixel LiControlVariate(const RayDifferential &ray,
                                 const Scene &scene, Sampler &sampler,
								 MemoryArena &arena,int depth = 0) const;
};

Integrator *CreateCvPathIntegrator(const ParamSet &params,
                                   std::shared_ptr<Sampler> sampler,
                                   std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_CV_PATH_INTEGRATOR_H
