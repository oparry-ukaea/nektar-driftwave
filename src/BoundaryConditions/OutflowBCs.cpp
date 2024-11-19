///////////////////////////////////////////////////////////////////////////////
//
// File: OutflowBCs.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Sheath boundary conditions class
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include "OutflowBCs.h"

namespace Nektar
{

std::string OutflowBCs::className =
    GetCustomBCsFactory().RegisterCreatorFunction("Outflow", OutflowBCs::create,
                                                  "OutflowBCs.");

OutflowBCs::OutflowBCs(const LU::SessionReaderSharedPtr &pSession,
                       const Array<OneD, MR::ExpListSharedPtr> &pFields,
                       const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
                       const int field_idx, const int pSpaceDim,
                       const int bcRegion, const int cnt,
                       SD::BoundaryConditionShPtr cnd)
    : CustomBCs(pSession, pFields, pTraceNormals, field_idx, pSpaceDim,
                bcRegion, cnt, cnd),
      m_mom_idx(field_idx)
{

    // m_dens_idx = m_field_indices["ne"];
    // ASSERTL0(m_dens_idx >= 0, "Expected a density field called ne - not
    // found")
}

void OutflowBCs::v_Apply(Array<OneD, Array<OneD, NekDouble>> &Fwd,
                         Array<OneD, Array<OneD, NekDouble>> &physarray,
                         const NekDouble &time)
{
    // boost::ignore_unused(time);

    // const Array<OneD, const int> &traceBndMap =
    //     m_fields[m_dens_idx]->GetTraceBndMap();

    // // Loop over all explists in this region
    // auto explists = m_fields[m_dens_idx]->GetBndCondExpansions()[m_bcRegion];
    // for (int e = 0; e < explists->GetExpSize(); ++e)
    // {
    //     // Current explist
    //     LR::ExpansionSharedPtr explist = explists->GetExp(e);
    //     // Offset in the field arrays for this explist
    //     int explist_offset = explists->GetPhys_Offset(e);
    //     // Offset in the trace map for this explist
    //     int trace_offset = m_fields[m_dens_idx]->GetTrace()->GetPhys_Offset(
    //         traceBndMap[m_offset + e]);

    //     Array<OneD, NekDouble> mass_times_vel(explist->GetTotPoints());
    //     this->evaluate_expression(explist, mass_times_vel);

    //     // Loop over points in this explist
    //     for (int idx_in_explist = 0; idx_in_explist <
    //     explist->GetTotPoints();
    //          idx_in_explist++)
    //     {
    //         int pt_idx = trace_offset + idx_in_explist;

    //         // Set momentum = density * const_mass * const_velocity
    //         // Assumes density, momentum explists have the same sizes
    //         (m_fields[m_mom_idx]
    //              ->GetBndCondExpansions()[m_bcRegion]
    //              ->UpdatePhys())[explist_offset + idx_in_explist] =
    //             m_fields[m_dens_idx]->GetPhys()[pt_idx] *
    //             mass_times_vel[idx_in_explist];
    //     }
    // }
}
} // namespace Nektar
