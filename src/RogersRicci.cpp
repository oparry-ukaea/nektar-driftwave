///////////////////////////////////////////////////////////////////////////////
//
// File: RogersRicci.cpp
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

#include "RogersRicci.h"

#include <MultiRegions/ContField.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

namespace Nektar
{

class UpwindNeumannSolver : public SolverUtils::RiemannSolver
{

public:
    static RiemannSolverSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession)
    {
        return RiemannSolverSharedPtr(new UpwindNeumannSolver(pSession));
    }

    static std::string solver_name;

    UpwindNeumannSolver(const LibUtilities::SessionReaderSharedPtr &pSession)
        : SolverUtils::RiemannSolver(pSession)
    {
    }

    void SetNeumannIdx(std::set<std::size_t> idx)
    {
        m_neumannIdx = idx;
    }

protected:
    std::set<std::size_t> m_neumannIdx;

    void v_Solve(const int nDim,
                 const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
                 const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
                 Array<OneD, Array<OneD, NekDouble>> &flux) final
    {
        ASSERTL1(CheckScalars("Vn"), "Vn not defined.");
        const Array<OneD, NekDouble> &traceVel = m_scalars["Vn"]();

        for (int j = 0; j < traceVel.size(); ++j)
        {
            const Array<OneD, const Array<OneD, NekDouble>> &tmp =
                traceVel[j] >= 0 ? Fwd : Bwd;
            for (int i = 0; i < Fwd.size(); ++i)
            {
                flux[i][j] = traceVel[j] * tmp[i][j];
            }
        }
    }
};

std::string UpwindNeumannSolver::solver_name =
    SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
        "CustomUpwind", UpwindNeumannSolver::create,
        "Customised upwind solver for 3DRR");

RogersRicci::RogersRicci(const LibUtilities::SessionReaderSharedPtr &session,
                         const SpatialDomains::MeshGraphSharedPtr &graph)
    : AdvectionSystem(session, graph), m_driftVel(3)
{
    // Set up constants
    /*
    m_c["B"] = c("omega_ci") * c("m_i") * c("q_E");
    m_c["c_s0"] = sqrt(c("T_e0") / c("m_i"));
    m_c["rho_s0"] = c("c_s0") / c("omega_ci");
    m_c["S_0n"] = 0.03 * c("n0") * c("c_s0") / c("R");
    m_c["S_0T"] = 0.03 * c("T_e0") * c("c_s0") / c("R");
    m_c["omega"] = 1.5 * c("R") / c("L_z");
    */
}

void RogersRicci::v_InitObject(bool DeclareField)
{
    WARNINGL0(GetNumExpModes() > 1,
              "Something in RogersRicci eqn sys fails with NUMMODES<2...");

    AdvectionSystem::v_InitObject(DeclareField);

    // Store mesh dimension for easy retrieval later.
    m_ndims = m_graph->GetMeshDimension();
    ASSERTL0(m_ndims == 2 || m_ndims == 3,
             "Solver only supports 2D or 3D meshes.");

    // Check fields all have the same number of quad points
    m_npts = m_fields[0]->GetNpoints();
    check_field_sizes(m_fields, m_npts);

    m_fields[phi_idx] =
        MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
            m_session, m_graph, m_session->GetVariable(phi_idx), true, true);

    // Assign storage for drift velocity.
    for (int i = 0; i < m_driftVel.size(); ++i)
    {
        m_driftVel[i] = Array<OneD, NekDouble>(m_npts, 0.0);
    }

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            m_homoInitialFwd = false;

            if (m_fields[0]->GetTrace())
            {
                m_traceVdriftnorm = Array<OneD, NekDouble>(GetTraceNpoints());
            }

            // Advection and Riemann solver for drift velocity
            init_advection_obj(m_session, m_drift_riemannSolver, advObj_vdrift);
            advObj_vdrift->SetFluxVector(&RogersRicci::GetDriftFluxVector,
                                         this);
            m_drift_riemannSolver->SetScalar(
                "Vn", &RogersRicci::GetDriftNormalVelocity, this);
            advObj_vdrift->InitObject(m_session, m_fields);
            break;
        }

        default:
        {
            ASSERTL0(false, "Unsupported projection type: only discontinuous"
                            " projection supported.");
            break;
        }
    }

    m_ode.DefineOdeRhs(&RogersRicci::ExplicitTimeInt, this);
    m_ode.DefineProjection(&RogersRicci::DoOdeProjection, this);

    if (!m_explicitAdvection)
    {
        m_implHelper = std::make_shared<ImplicitHelper>(
            m_session, m_fields, m_ode, m_intVariables.size());
        m_implHelper->InitialiseNonlinSysSolver();
        m_ode.DefineImplicitSolve(&ImplicitHelper::ImplicitTimeInt,
                                  m_implHelper);
    }

    // Store distance of quad points from origin in transverse plane.
    // (used to compute source terms)
    Array<OneD, NekDouble> x = Array<OneD, NekDouble>(m_npts);
    Array<OneD, NekDouble> y = Array<OneD, NekDouble>(m_npts);
    m_r                      = Array<OneD, NekDouble>(m_npts);
    if (m_ndims == 3)
    {
        Array<OneD, NekDouble> z = Array<OneD, NekDouble>(m_npts);
        m_fields[0]->GetCoords(x, y, z);
    }
    else
    {
        m_fields[0]->GetCoords(x, y);
    }
    for (int i = 0; i < m_npts; ++i)
    {
        m_r[i] = sqrt(x[i] * x[i] + y[i] * y[i]);
    }

    // Set up user-defined boundary conditions, if any were specified
    for (int ifld = 0; ifld < m_fields.size(); ifld++)
    {
        int cnt = 0;
        for (int icnd = 0; icnd < m_fields[ifld]->GetBndConditions().size();
             ++icnd)
        {
            SpatialDomains::BoundaryConditionShPtr cnd =
                m_fields[ifld]->GetBndConditions()[icnd];
            if (cnd->GetBoundaryConditionType() != SpatialDomains::ePeriodic)
            {
                std::string type = cnd->GetUserDefined();
                if (!type.empty())
                {
                    CustomBCsSharedPtr BCs_instance =
                        GetCustomBCsFactory().CreateInstance(
                            type, m_session, m_fields, m_traceNormals, ifld,
                            m_spacedim, icnd, cnt, cnd);
                    m_custom_BCs.push_back(BCs_instance);
                }
                cnt +=
                    m_fields[ifld]->GetBndCondExpansions()[icnd]->GetExpSize();
            }
        }
    }
}

/**
 * @brief Perform projection into correct polynomial space.
 *
 * This routine projects the @p inarray input and ensures the @p outarray output
 * lives in the correct space. Since we are hard-coding DG, this corresponds to
 * a simple copy from in to out, since no elemental connectivity is required and
 * the output of the RHS function is polynomial.
 */
void RogersRicci::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvariables = inarray.size(), npoints = GetNpoints();
    SetBoundaryConditions(time);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
    }
}

/**
 * @brief Compute the flux vector for drift velocity advection.
 */
void RogersRicci::GetDriftFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    get_flux_vector(physfield, flux, m_driftVel);
}

/**
 * @brief Compute the normal drift velocities on the
 * trace/skeleton/edges of the mesh.
 */
Array<OneD, NekDouble> &RogersRicci::GetDriftNormalVelocity()
{
    return get_norm_vels(m_traceVdriftnorm, m_driftVel, m_traceNormals,
                         m_fields);
}

void RogersRicci::set_int_idx(const int &idx, int &int_idx)
{
    auto pos = std::find(m_intVariables.begin(), m_intVariables.end(), idx);
    int_idx  = std::distance(m_intVariables.begin(), pos);
}

void check_field_sizes(Array<OneD, MultiRegions::ExpListSharedPtr> fields,
                       const int npts)
{
    for (auto i = 0; i < fields.size(); i++)
    {
        ASSERTL0(fields[i]->GetNpoints() == npts,
                 "Detected fields with different numbers of quadrature points; "
                 "this solver assumes they're all the same");
    }
}

void check_var_idx(const LibUtilities::SessionReaderSharedPtr session,
                   const int &idx, const std::string var_name)
{
    std::stringstream err;
    err << "Expected variable index " << idx << " to correspond to '"
        << var_name << "'. Check your session file.";
    ASSERTL0(session->GetVariable(idx).compare(var_name) == 0, err.str());
}

Array<OneD, NekDouble> &get_norm_vels(
    Array<OneD, NekDouble> &trace_norm_vels,
    const Array<OneD, Array<OneD, NekDouble>> &adv_vels,
    const Array<OneD, Array<OneD, NekDouble>> &trace_norms,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
{
    // Number of trace (interface) points
    int nTracePts = trace_norm_vels.size();

    ASSERTL0(
        adv_vels.size() >= trace_norms.size(),
        "Velocity array dimension is less than trace_norms array dimension");

    ASSERTL0(trace_norms[0].size() == nTracePts,
             "trace_norms array and trace_norm_vels have different number of "
             "points");

    // Auxiliary variable to compute the normal velocity
    Array<OneD, NekDouble> tmp(nTracePts);

    // Reset the normal velocity
    Vmath::Zero(nTracePts, trace_norm_vels, 1);

    // Compute dot product of velocity along trace with trace normals. Store in
    // trace_norm_vels.
    for (int i = 0; i < trace_norms.size(); ++i)
    {
        fields[0]->ExtractTracePhys(adv_vels[i], tmp);

        Vmath::Vvtvp(nTracePts, trace_norms[i], 1, tmp, 1, trace_norm_vels, 1,
                     trace_norm_vels, 1);
    }

    return trace_norm_vels;
}

void get_flux_vector(const Array<OneD, Array<OneD, NekDouble>> &physfield,
                     Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux,
                     const Array<OneD, Array<OneD, NekDouble>> &adv_vels)
{
    ASSERTL0(adv_vels.size() >= flux[0].size(),
             "Velocity array dimension is less than flux array dimension");

    int nq = physfield[0].size();

    for (int i = 0; i < flux.size(); ++i)
    {
        for (int j = 0; j < flux[0].size(); ++j)
        {
            Vmath::Vmul(nq, physfield[i], 1, adv_vels[j], 1, flux[i][j], 1);
        }
    }
}

void init_advection_obj(LibUtilities::SessionReaderSharedPtr session,
                        RiemannSolverSharedPtr &riemannSolver,
                        SolverUtils::AdvectionSharedPtr &advObj)
{
    std::string advName, riemName;
    session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
    session->LoadSolverInfo("UpwindType", riemName, "Upwind");
    advObj =
        SolverUtils::GetAdvectionFactory().CreateInstance(advName, advName);
    riemannSolver = SolverUtils::GetRiemannSolverFactory().CreateInstance(
        riemName, session);
    advObj->SetRiemannSolver(riemannSolver);
}

} // namespace Nektar
