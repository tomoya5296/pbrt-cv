#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CV_PIXEL_H
#define PBRT_CV_PIXEL_H

#include "spectrum.h"

namespace pbrt {

class CvDualPixel {
public:
    CvDualPixel();
    CvDualPixel(const Spectrum &L1, const Spectrum &L2, Float reciprocal_pdf);
    CvDualPixel(const CvDualPixel &p);
    /*CvDualPixel operator+(const CvDualPixel &p){
        return CvDualPixel(L1 + p.L1, L2 + p.L2, reciprocal_pdf + p.reciprocal_pdf);
    }
    CvDualPixel operator/(const float p) {
        return CvDualPixel(L1 / p, L2 / p, reciprocal_pdf / p);
    }*/
    CvDualPixel &operator=(const CvDualPixel &p);

    void SetZero();
    void AddPixel(const CvDualPixel &p, Float filterWeight = 1.);

private:
    Spectrum L1, L2, D, L1square, L2square, Dsquare;
    Float reciprocal_pdf;
    Float filterWeightSum;

    friend class CvFilm;
    friend class CvFilmTile;
};

}  // namespace pbrt

#endif  //
