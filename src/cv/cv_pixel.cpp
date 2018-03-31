#include "cv_pixel.h"

namespace pbrt {

CvDualPixel::CvDualPixel() {
    SetZero();
}

CvDualPixel::CvDualPixel(const Spectrum &L1_, const Spectrum &L2_, Float recip_pdf)
    : CvDualPixel() {
    L1 = L1_ * recip_pdf;
    L2 = L2_ * recip_pdf;
    L1square =  L1_ * L1_ * recip_pdf* recip_pdf;
    L2square =  L2_ * L2_ * recip_pdf* recip_pdf;
    reciprocal_pdf = recip_pdf;
    Spectrum d = L1_ - L2_;
    D = d * recip_pdf;
    Dsquare = d * d * recip_pdf * recip_pdf;
    filterWeightSum = 0.f;
}

CvDualPixel::CvDualPixel(const CvDualPixel &p)
    : CvDualPixel() {
    this->operator=(p);
}

CvDualPixel &CvDualPixel::operator=(const CvDualPixel &p) {
    this->L1 = p.L1;
    this->L2 = p.L2;
    this->D = p.D;
    this->L1square = p.L1square;
    this->L2square = p.L2square;
    this->Dsquare = p.Dsquare;
    this->reciprocal_pdf = p.reciprocal_pdf;
    this->filterWeightSum = p.filterWeightSum;
    return *this;
}

void CvDualPixel::SetZero() {
    L1 = Spectrum(0.f);
    L2 = Spectrum(0.f);
    D = Spectrum(0.f);
    L1square = Spectrum(0.f);
    L2square = Spectrum(0.f);
    Dsquare = Spectrum(0.f);
    reciprocal_pdf = Float(0.f);
    filterWeightSum = Float(0.f);
}

void CvDualPixel::AddPixel(const CvDualPixel &p, Float filterWeight) {
    L1 += filterWeight * p.L1;
    L2 += filterWeight * p.L2;
    D += filterWeight * p.D;
    L1square += filterWeight * p.L1square;
    L2square += filterWeight * p.L2square;
    Dsquare += filterWeight * p.Dsquare;
    reciprocal_pdf += filterWeight * p.reciprocal_pdf;
    filterWeightSum += filterWeight;
}

}  // namespace pbrt
