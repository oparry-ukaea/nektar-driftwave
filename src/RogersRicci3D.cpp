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
    w_idx   = 2;
    ue_idx  = 3;
    phi_idx = 4;
}

void RogersRicci3D::v_InitObject(bool DeclareField)
{
    // Needs to be set before calling parent member function
    m_intVariables = {n_idx, Te_idx, w_idx, ue_idx};

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
}

/**
 * @brief Evaluate the right-hand side of the ODE system used to integrate in
 * time.
 *
 * This routine performs the bulk of the work in this class, and essentially
 * computes the right hand side term of the generalised ODE system
 *
 * \f\[ \frac{\partial \mathbf{u}}{\partial t} = \mathbf{R}(\mathbf{u}) \f\]
 *
 * The order of operations is as follows:
 *
 * - First, compute the electrostatic potential \f$ \phi \f$, given the
 * - Using this, compute the drift velocity \f$ (\partial_y\phi,
 *   -\partial_x\phi).
 * - Then evaluate the \f$ \nabla\cdot\mathbf{F} \f$ operator using the
 *   advection object #m_advObject.
 * - Finally put this on the right hand side and evaluate the source terms for
 *   each field.
 *
 * The assumption here is that fields are ordered inside `m_fields` so that
 * field 0 is vorticity \f$ \zeta \f$, field 1 is number density \f$ n \f$, and
 * field 2 is electrostatic potential. Only \f$ \zeta \f$ and \f$ n \f$ are time
 * integrated.
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
    m_fields[phi_idx]->PhysDeriv(m_fields[phi_idx]->GetPhys(), m_advVel[1],
                                 m_advVel[0], dummy);

    // We frequently use vector math (Vmath) routines for one-line operations
    // like negating entries in a vector.
    Vmath::Neg(m_npts, m_advVel[1], 1);

    // Do advection for zeta, n. The hard-coded '3' here indicates that we
    // should only advect the first two components of inarray.
    m_advObject->Advect(inarray.size(), m_fields, m_advVel, inarray, outarray,
                        time);

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
        n_out[i]     = -40 * n_out[i] - 1.0 / 24.0 * et * n[i] + st;
        T_e_out[i] =
            -40 * T_e_out[i] - 1.0 / 36.0 * (1.71 * et - 0.71) * T_e[i] + st;
        w_out[i] = -40 * w_out[i] + 1.0 / 24.0 * (1 - et);
    }
}
} // namespace Nektar
