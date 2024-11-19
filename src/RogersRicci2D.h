///////////////////////////////////////////////////////////////////////////////
//
// File RogersRicci2D.h
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
// Description: Equation system for the 2D Rogers & Ricci model.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_ROGERSRICCI2D_H
#define NEKTAR_ROGERSRICCI2D_H

#include "RogersRicci.h"

namespace Nektar
{

/**
 * @brief An equation system for the drift-wave solver.
 */
class RogersRicci2D : public RogersRicci
{
public:
    // Friend class to allow the memory manager to allocate shared pointers of
    // this class.
    friend class MemoryManager<RogersRicci2D>;

    /// Creates an instance of this class. This static method is registered with
    /// a factory.
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &session,
        const SpatialDomains::MeshGraphSharedPtr &graph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<RogersRicci2D>::AllocateSharedPtr(session, graph);
        p->InitObject();
        return p;
    }

    /// Name of class, used to statically initialise a function pointer for the
    /// create method above.
    static std::string className;

protected:
    RogersRicci2D(const LibUtilities::SessionReaderSharedPtr &session,
                  const SpatialDomains::MeshGraphSharedPtr &graph);

    virtual void ExplicitTimeInt(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble time) override;

    void v_InitObject(bool DeclareField) override;
};

} // namespace Nektar

#endif
