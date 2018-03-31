#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CV_FILM_H
#define PBRT_CV_FILM_H

#include "film.h"
#include "cv_pixel.h"

namespace pbrt {

class CvFilm : public Film {
public:
    CvFilm(const Point2i &resolution, const Bounds2f &cropWindow,
           std::unique_ptr<Filter> filter, Float diagonal,
           const std::string &filename, Float scale,
           Float maxSampleLuminance = Infinity);

    std::unique_ptr<CvFilmTile> GetCvFilmTile(const Bounds2i &sampleBounds);
    void MergeFilmTile(std::unique_ptr<CvFilmTile> tile);
    void WriteImage(Float splatScale = 1, int samplesPerPixel = 0) final override;

    CvDualPixel &GetCvPixel(const Point2i &p);

private:
    std::unique_ptr<CvDualPixel[]> cvPixels;
};

class CvFilmTile {
  public:
    // FilmTile Public Methods
    CvFilmTile(const Bounds2i &pixelBounds, const Vector2f &filterRadius,
             const Float *filterTable, int filterTableSize,
             Float maxSampleLuminance)
        : pixelBounds(pixelBounds),
          filterRadius(filterRadius),
          invFilterRadius(1 / filterRadius.x, 1 / filterRadius.y),
          filterTable(filterTable),
          filterTableSize(filterTableSize),
          maxSampleLuminance(maxSampleLuminance) {
        pixels = std::vector<CvDualPixel>(std::max(0, pixelBounds.Area()));
    }

    void AddSample(const Point2f &pFilm, CvDualPixel splat,
                   Float sampleWeight = 1.) {
        CHECK(sampleWeight == 1.) << "Now the case \"sampleWeight = 1\" is supported!";

        ProfilePhase _(Prof::AddFilmSample);
        //if (L.y() > maxSampleLuminance)
        //    L *= maxSampleLuminance / L.y();
        // Compute sample's raster bounds
        Point2f pFilmDiscrete = pFilm - Vector2f(0.5f, 0.5f);
        Point2i p0 = (Point2i)Ceil(pFilmDiscrete - filterRadius);
        Point2i p1 =
            (Point2i)Floor(pFilmDiscrete + filterRadius) + Point2i(1, 1);
        p0 = Max(p0, pixelBounds.pMin);
        p1 = Min(p1, pixelBounds.pMax);

        // Loop over filter support and add sample to pixel arrays

        // Precompute $x$ and $y$ filter table offsets
        int *ifx = ALLOCA(int, p1.x - p0.x);
        for (int x = p0.x; x < p1.x; ++x) {
            Float fx = std::abs((x - pFilmDiscrete.x) * invFilterRadius.x *
                                filterTableSize);
            ifx[x - p0.x] = std::min((int)std::floor(fx), filterTableSize - 1);
        }
        int *ify = ALLOCA(int, p1.y - p0.y);
        for (int y = p0.y; y < p1.y; ++y) {
            Float fy = std::abs((y - pFilmDiscrete.y) * invFilterRadius.y *
                                filterTableSize);
            ify[y - p0.y] = std::min((int)std::floor(fy), filterTableSize - 1);
        }
        for (int y = p0.y; y < p1.y; ++y) {
            for (int x = p0.x; x < p1.x; ++x) {
                // Evaluate filter value at $(x,y)$ pixel
                int offset = ify[y - p0.y] * filterTableSize + ifx[x - p0.x];
                Float filterWeight = filterTable[offset];

                // Update pixel values with filtered sample contribution
                CvDualPixel &pixel = GetPixel(Point2i(x, y));
                pixel.AddPixel(splat, filterWeight);
           }
        }
    }

    CvDualPixel &GetPixel(const Point2i &p) {
        CHECK(InsideExclusive(p, pixelBounds));
        int width = pixelBounds.pMax.x - pixelBounds.pMin.x;
        int offset =
            (p.x - pixelBounds.pMin.x) + (p.y - pixelBounds.pMin.y) * width;
        return pixels[offset];
    }

    const CvDualPixel &GetPixel(const Point2i &p) const {
        CHECK(InsideExclusive(p, pixelBounds));
        int width = pixelBounds.pMax.x - pixelBounds.pMin.x;
        int offset =
            (p.x - pixelBounds.pMin.x) + (p.y - pixelBounds.pMin.y) * width;
        return pixels[offset];
    }

    Bounds2i GetPixelBounds() const { return pixelBounds; }

  private:
    // FilmTile Private Data
    const Vector2f filterRadius, invFilterRadius;
    const Float *filterTable;
    const int filterTableSize;
    std::vector<CvDualPixel> pixels;
    const Float maxSampleLuminance;
    const Bounds2i pixelBounds;
    friend class CvFilm;
};
Film *CreateCvFilm(const ParamSet &paramSet, std::unique_ptr<Filter> filter);

}  // namespace pbrt

#endif  // PBRT_CV_FILM_H
