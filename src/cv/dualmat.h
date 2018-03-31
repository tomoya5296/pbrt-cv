#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_DUALMAT_H
#define PBRT_MATERIALS_DUALMAT_H

// materials/mixmat.h*
#include "pbrt.h"
#include "material.h"

namespace pbrt {

// MixMaterial Declarations
class DualMaterial : public Material {
  public:
    // MixMaterial Public Methods
    DualMaterial(const std::shared_ptr<Material> &m1_,
                 const std::shared_ptr<Material> &m2_)
        : m1(m1_), m2(m2_) {}

    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

    static void SetMaterialId(int id) {
        materialId = id;
    }

  private:
    // MixMaterial Private Data
    std::shared_ptr<Material> m1, m2;
    static thread_local int materialId;
};

DualMaterial *CreateDualMaterial(const std::shared_ptr<Material> &m1,
                                 const std::shared_ptr<Material> &m2);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_DUALMAT_H
