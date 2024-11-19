#include "CustomBCs.h"

namespace Nektar
{
// Static factory singleton
CustomBCsFactory &GetCustomBCsFactory()
{
    static CustomBCsFactory instance;
    return instance;
}

/**
 * Abstract base class for custom boundary conditions.
 */
CustomBCs::CustomBCs(const LU::SessionReaderSharedPtr &session,
                     const Array<OneD, MR::ExpListSharedPtr> &fields,
                     const Array<OneD, Array<OneD, NekDouble>> &trace_normals,
                     const int field_idx, const int space_dim,
                     const int bc_region, const int cnt,
                     SD::BoundaryConditionShPtr cnd)
    : m_session(session), m_fields(fields), m_traceNormals(trace_normals),
      m_spacedim(space_dim), m_bcRegion(bc_region), m_offset(cnt), m_cnd(cnd),
      m_field_indices()
{
    auto vars = session->GetVariables();
    for (auto ii = 0; ii < vars.size(); ii++)
    {
        m_field_indices[vars[ii]] = ii;
    }

    ASSERTL0(space_dim == 3, "CustomBCs only set up for space_dim=3");
    SD::BoundaryConditionType BCtype = m_cnd->GetBoundaryConditionType();
    ASSERTL0(BCtype == SD::eDirichlet || BCtype == SD::eNeumann,
             "CustomBCs must be either Dirichlet or Neumann type");
}

void CustomBCs::evaluate_expression(LR::ExpansionSharedPtr explist,
                                    Array<OneD, NekDouble> result)
{

    int npts = explist->GetTotPoints();
    // Coords of quad points in this explist
    Array<OneD, NekDouble> tmp_x(npts), tmp_y(npts), tmp_z(npts);
    explist->GetCoords(tmp_x, tmp_y, tmp_y);

    switch (m_cnd->GetBoundaryConditionType())
    {
        case SD::eDirichlet:
        {
            auto dcond =
                std::dynamic_pointer_cast<SD::DirichletBoundaryCondition>(
                    m_cnd);
            dcond->m_dirichletCondition.Evaluate(tmp_x, tmp_y, tmp_y, result);
        }
        break;
        case SD::eNeumann:
        {
            auto ncond =
                std::dynamic_pointer_cast<SD::NeumannBoundaryCondition>(m_cnd);
            ncond->m_neumannCondition.Evaluate(tmp_x, tmp_y, tmp_y, result);
        }
        break;
        default:
        {
        }
        break;
    }
}

/**
 * @param   bcRegion      id of the boundary region
 * @param   cnt
 * @param   Fwd
 * @param   physarray
 * @param   time
 */
void CustomBCs::Apply(Array<OneD, Array<OneD, NekDouble>> &Fwd,
                      Array<OneD, Array<OneD, NekDouble>> &physarray,
                      const NekDouble &time)
{
    v_Apply(Fwd, physarray, time);
}

/**
 * @ brief Newly added bc should specify this virtual function
 * if the Bwd/value in m_bndCondExpansions is the target value, like Dirichlet,
 * bc weight should be 1.0.
 * if some average Fwd and Bwd/value in m_bndCondExpansions
 * is the target value like WallViscousBC weight should be 0.5.
 */
void CustomBCs::v_ApplyBwdWeight()
{
    for (int i = 0; i < m_fields.size(); ++i)
    {
        m_fields[i]->SetBndCondBwdWeight(m_bcRegion, m_weight);
    }
}

} // namespace Nektar
