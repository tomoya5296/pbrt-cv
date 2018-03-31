// materials/mixmat.cpp*
#include "materials/matte.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"
#include "interaction.h"

#include "dualmat.h"

namespace pbrt {

// MixMaterial Method Definitions
void DualMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
                                              MemoryArena &arena,
                                              TransportMode mode,
                                              bool allowMultipleLobes) const {
    // Compute weights and original _BxDF_s for mix material
    if (materialId == 0) {
        m1->ComputeScatteringFunctions(si, arena, mode, allowMultipleLobes);
    } else {
        m2->ComputeScatteringFunctions(si, arena, mode, allowMultipleLobes);
    }
}

thread_local int DualMaterial::materialId = 0;

DualMaterial *CreateDualMaterial(const std::shared_ptr<Material> &m1,
                                 const std::shared_ptr<Material> &m2) {
    return new DualMaterial(m1, m2);
}

}  // namespace pbrt
