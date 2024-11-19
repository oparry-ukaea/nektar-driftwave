///////////////////////////////////////////////////////////////////////////////
//
// File RogersRicci.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Abstract base class for Rogers & Ricci equation sytems
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_ROGERSRICCI_H
#define NEKTAR_ROGERSRICCI_H

#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

#include "BoundaryConditions/CustomBCs.h"
#include "ImplicitHelper.hpp"

using namespace Nektar::SolverUtils;

namespace Nektar
{

/**
 * @brief An equation system for the drift-wave solver.
 */
class RogersRicci : public AdvectionSystem
{
public:
    // Friend class to allow the memory manager to allocate shared pointers of
    // this class.
    friend class MemoryManager<RogersRicci>;

    /// Default destructor.
    virtual ~RogersRicci() = default;

    virtual void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);
    virtual void GetDriftFluxVector(
        const Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux);

    Array<OneD, NekDouble> &GetDriftNormalVelocity();

protected:
    /// Helper function to define constants.
    NekDouble c(std::string n)
    {
        auto it = m_c.find(n);

        ASSERTL0(it != m_c.end(), "Unknown constant");

        return it->second;
    }

    /// Map of known constants
    std::map<std::string, NekDouble> m_c = {
        {"T_e", 6.0},      {"L_z", 18.0},       {"n_0", 2.0e18},
        {"m_i", 6.67e-27}, {"omega_ci", 9.6e5}, {"lambda", 3.0},
        {"R", 0.5}};

    /// Storage for the (drift) advection velocity. The outer index is
    /// dimension, and inner index the solution nodes (in physical space).
    Array<OneD, Array<OneD, NekDouble>> m_driftVel;
    /// Storage for the dot product of drift velocity with element edge normals,
    /// required for the DG formulation.
    Array<OneD, NekDouble> m_traceVdriftnorm;
    /// A Riemann solver object to solve numerical fluxes arising from DG: in
    /// this case a simple upwind.
    RiemannSolverSharedPtr m_drift_riemannSolver;
    /// Helper object for fully-implicit solve.
    std::shared_ptr<ImplicitHelper> m_implHelper;

    Array<OneD, NekDouble> m_r;

    /// Protected constructor. Since we use a factory pattern, objects should be
    /// constructed via the SolverUtils::EquationSystem factory.
    RogersRicci(const LibUtilities::SessionReaderSharedPtr &session,
                const SpatialDomains::MeshGraphSharedPtr &graph);

    virtual void v_InitObject(bool DeclareField = true) override;

    void InitialiseNonlinSysSolver(void);

    virtual void ExplicitTimeInt(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble time) = 0;

    void set_int_idx(const int &idx, int &int_idx);

    int m_npts;
    int m_ndims;

    // Subclasses set these indices
    int n_idx;
    int Te_idx;
    int w_idx;
    int phi_idx;

    // Indices in integration arrays; also set by subclasses
    int n_int_idx;
    int Te_int_idx;
    int w_int_idx;
    int phi_int_idx;

    // Use this, rather than AdvectionSystem::m_advObject, for drift velocity
    // advection - makes things clearer when multiple advection objects are
    // needed in 3D subclass
    SolverUtils::AdvectionSharedPtr advObj_vdrift;

private:
    /// User defined boundary conditions
    std::vector<CustomBCsSharedPtr> m_custom_BCs;
};

// Helper functions
void check_var_idx(const LibUtilities::SessionReaderSharedPtr session,
                   const int &idx, const std::string var_name);

void check_field_sizes(Array<OneD, MultiRegions::ExpListSharedPtr> fields,
                       const int npts);

void get_flux_vector(const Array<OneD, Array<OneD, NekDouble>> &physfield,
                     Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux,
                     const Array<OneD, Array<OneD, NekDouble>> &adv_vels);

Array<OneD, NekDouble> &get_norm_vels(
    Array<OneD, NekDouble> &trace_norm_vels,
    const Array<OneD, Array<OneD, NekDouble>> &adv_vels,
    const Array<OneD, Array<OneD, NekDouble>> &trace_norms,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields);

void init_advection_obj(LibUtilities::SessionReaderSharedPtr session,
                        RiemannSolverSharedPtr &riemannSolver,
                        SolverUtils::AdvectionSharedPtr &advObj);

} // namespace Nektar

#endif
