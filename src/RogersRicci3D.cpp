///////////////////////////////////////////////////////////////////////////////
//
// File: RogersRicci3D.cpp
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
// Description: Equation system for the 3D Rogers & Ricci model.
//
///////////////////////////////////////////////////////////////////////////////

#include "RogersRicci3D.h"
#include <MultiRegions/ContField.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

namespace Nektar
{

std::string RogersRicci3D::className =
    GetEquationSystemFactory().RegisterCreatorFunction(
        "RogersRicci3D", RogersRicci3D::create,
        "System for the Rogers-Ricci 3D system of equations.");

RogersRicci3D::RogersRicci3D(
    const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph)
    : RogersRicci(session, graph), UnsteadySystem(session, graph)
{
    n_idx   = 0;
    Te_idx  = 1;
    ue_idx  = 2;
    w_idx   = 3;
    phi_idx = 4;
}

void RogersRicci3D::v_InitObject(bool DeclareField)
{
    // Needs to be set before calling parent member function
    m_intVariables = {n_idx, Te_idx, ue_idx, w_idx};

    RogersRicci::v_InitObject(DeclareField);

    ASSERTL0(m_fields.size() == 5,
             "Incorrect number of variables detected (expected 4): check your "
             "session file.");

    // Check variable order is as expected
    check_var_idx(m_session, n_idx, "n");
    check_var_idx(m_session, Te_idx, "T_e");
    check_var_idx(m_session, w_idx, "w");
    check_var_idx(m_session, ue_idx, "u_e");
    check_var_idx(m_session, phi_idx, "phi");

    // Indices aren't the same in the integration arrays as in m_fields; set the
    // former here
    set_int_idx(n_idx, n_int_idx);
    set_int_idx(Te_idx, Te_int_idx);
    set_int_idx(w_idx, w_int_idx);
    set_int_idx(ue_idx, ue_int_idx);

    m_session->LoadParameter("mu", this->mu);
    m_session->LoadParameter("tau", this->tau);

    // Init some class-level storage
    adv_result = Array<OneD, Array<OneD, NekDouble>>(m_intVariables.size());

    for (auto ii = 0; ii < this->adv_result.size(); ii++)
    {
        this->adv_result[ii] = Array<OneD, NekDouble>(m_npts, 0.0);
    }
    this->ue_vals = Array<OneD, Array<OneD, NekDouble>>(this->m_ndims);
    for (auto ii = 0; ii < this->ue_vals.size(); ii++)
    {
        this->ue_vals[ii] = Array<OneD, NekDouble>(m_npts, 0.0);
    }

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {

            if (m_fields[0]->GetTrace())
            {
                m_traceUenorm = Array<OneD, NekDouble>(GetTraceNpoints());
            }

            // Advection object and Riemann solver for electron parallel
            // velocity
            init_advection_obj(m_session, ue_riemannSolver, advObj_ue);
            advObj_ue->SetFluxVector(&RogersRicci3D::GetUeFluxVector, this);
            ue_riemannSolver->SetScalar(
                "Vn", &RogersRicci3D::GetUeNormalVelocity, this);
            advObj_ue->InitObject(m_session, m_fields);
            break;
        }

        default:
            break;
            // Valid projection checked elsewhere; no action here
    }
}

/**
 * @brief Evaluate the right-hand side of the ODE system used to integrate in
 * time.
 *
 * @param inarray    Array containing each field's current state.
 * @param outarray   The result of the right-hand side operator for each field
 *                   being time integrated.
 * @param time       Current value of time.
 */
void RogersRicci3D::ExplicitTimeInt(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    Array<OneD, NekDouble> n       = inarray[n_int_idx];
    Array<OneD, NekDouble> T_e     = inarray[Te_int_idx];
    Array<OneD, NekDouble> w       = inarray[w_int_idx];
    Array<OneD, NekDouble> ue      = inarray[ue_int_idx];
    Array<OneD, NekDouble> n_out   = outarray[n_int_idx];
    Array<OneD, NekDouble> T_e_out = outarray[Te_int_idx];
    Array<OneD, NekDouble> w_out   = outarray[w_int_idx];
    Array<OneD, NekDouble> ue_out  = outarray[ue_int_idx];
    Array<OneD, NekDouble> phi     = m_fields[phi_idx]->UpdatePhys();

    // Set up factors for electrostatic potential solve. We support a generic
    // Helmholtz solve of the form (\nabla^2 - \lambda) u = f, so this sets
    // \lambda to zero.
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda]   = 0.0;
    factors[StdRegions::eFactorCoeffD22] = 0.0;

    Vmath::Zero(m_fields[phi_idx]->GetNcoeffs(),
                m_fields[phi_idx]->UpdateCoeffs(), 1);

    // Solve for phi. Output of this routine is in coefficient (spectral) space,
    // so backwards transform to physical space since we'll need that for the
    // advection step & computing drift velocity.
    m_fields[phi_idx]->HelmSolve(w, m_fields[phi_idx]->UpdateCoeffs(), factors);
    m_fields[phi_idx]->BwdTrans(m_fields[phi_idx]->GetCoeffs(), phi);

    // Calculate drift velocity v_E: PhysDeriv takes input and computes spatial
    // derivatives.
    Array<OneD, NekDouble> dummy = Array<OneD, NekDouble>(m_npts);
    m_fields[phi_idx]->PhysDeriv(m_fields[phi_idx]->GetPhys(), m_driftVel[1],
                                 m_driftVel[0], dummy);

    // We frequently use vector math (Vmath) routines for one-line operations
    // like negating entries in a vector.
    Vmath::Neg(m_npts, m_driftVel[1], 1);

    // // Add parallel velocity to advection
    // Vmath::Vcopy(m_npts, ue, 1, m_driftVel[2], 1);

    // First compute drift advection terms
    advObj_vdrift->Advect(inarray.size(), m_fields, m_driftVel, inarray,
                          outarray, time);
    Vmath::Smul(m_npts, -40.0, n_out, 1, n_out, 1);
    Vmath::Smul(m_npts, -40.0, T_e_out, 1, T_e_out, 1);
    Vmath::Smul(m_npts, -40.0 / 2.0, ue_out, 1, ue_out, 1);
    Vmath::Smul(m_npts, -40.0, w_out, 1, w_out, 1);

    // Copy ue phys vals into class-level storage to make callback funcs work
    Vmath::Vcopy(m_npts, ue, 1, this->ue_vals[2], 1);
    // Compute ue advection terms for n, T, ue itself; add to outarray
    advObj_ue->Advect(3, m_fields, this->ue_vals, inarray, adv_result, time);
    Vmath::Vadd(m_npts, n_out, 1, adv_result[0], 1, n_out, 1);
    Vmath::Vadd(m_npts, T_e_out, 1, adv_result[1], 1, T_e_out, 1);
    Vmath::Vadd(m_npts, ue_out, 1, adv_result[2], 1, ue_out, 1);

    // Put advection term on the right hand side.
    const NekDouble rho_s0 = 1.2e-2;
    const NekDouble r_s = 20 * rho_s0; // can't find this in list of constants,
                                       // stolen from rr.py... fingers crossed

    // stolen from Ed/Owen's code, rr.py
    const NekDouble Ls_boost = 2.0;
    const NekDouble L_s      = 0.5 * rho_s0 * Ls_boost; // maybe wrong

    for (auto i = 0; i < m_npts; ++i)
    {
        NekDouble et = exp(3 - phi[i] / sqrt(T_e[i] * T_e[i] + 1e-4));
        NekDouble st = 0.03 * (1.0 - tanh((rho_s0 * m_r[i] - r_s) / L_s));
        n_out[i]     = n_out[i] - 1.0 / 24.0 * et * n[i] + st;
        T_e_out[i] = T_e_out[i] - 1.0 / 36.0 * (1.71 * et - 0.71) * T_e[i] + st;
        w_out[i]   = w_out[i] + 1.0 / 24.0 * (1 - et);
    }
}

/**
 * @brief Compute the flux vector for parallel electron velocity advection.
 */
void RogersRicci3D::GetUeFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    get_flux_vector(physfield, flux, this->ue_vals);
}

/**
 * @brief Compute the normal parallel electron velocities on the
 * trace/skeleton/edges of the mesh.
 */
Array<OneD, NekDouble> &RogersRicci3D::GetUeNormalVelocity()
{
    return get_norm_vels(m_traceUenorm, this->ue_vals, m_traceNormals,
                         m_fields);
}

} // namespace Nektar
