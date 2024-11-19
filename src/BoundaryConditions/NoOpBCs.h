#ifndef DRIFTWAVE_NOOPBCS_H
#define DRIFTWAVE_NOOPBCS_H

#include "CustomBCs.h"

namespace LR = Nektar::LocalRegions;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;
namespace Nektar
{

/**
 * @brief No-op boundary conditions class
 */
class NoOpBCs : public CustomBCs
{
public:
    friend class MemoryManager<NoOpBCs>;

    /// Creates an instance of this class
    static CustomBCsSharedPtr create(
        const LU::SessionReaderSharedPtr &pSession,
        const Array<OneD, MR::ExpListSharedPtr> &pFields,
        const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
        const int field_idx, const int pSpaceDim, const int bcRegion,
        const int cnt, SD::BoundaryConditionShPtr cnd)
    {
        CustomBCsSharedPtr p = MemoryManager<NoOpBCs>::AllocateSharedPtr(
            pSession, pFields, pTraceNormals, field_idx, pSpaceDim, bcRegion,
            cnt, cnd);
        return p;
    }

    /// Name of the class
    static std::string className;

protected:
    virtual void v_Apply(Array<OneD, Array<OneD, NekDouble>> &Fwd,
                         Array<OneD, Array<OneD, NekDouble>> &physarray,
                         const NekDouble &time);

private:
    NoOpBCs(const LU::SessionReaderSharedPtr &pSession,
            const Array<OneD, MR::ExpListSharedPtr> &pFields,
            const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
            const int field_idx, const int pSpaceDim, const int bcRegion,
            const int cnt, SD::BoundaryConditionShPtr cnd);

    virtual ~NoOpBCs(void){};
};

} // namespace Nektar

#endif
