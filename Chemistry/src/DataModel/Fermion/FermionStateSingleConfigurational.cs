// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
    // This class is for Fermion terms that are not grouped into Hermitian bunches.
    // Maybe need a stype for quantum state?

    /// <summary>
    /// Class representing a sequence of fermionic raising and lowering operators, subject to the additional constraints: 
    /// 1) Normal-ordered, where all raising operators are to the left of all lowering operators.
    /// 2) Index-ordered, where are raising(lowering) operators are in ascending(descending) order.
    /// 3) Contains only creation operators.
    /// </summary>
    public class FermionStateSingleConfigurational : FermionTermSingle
    {

        public SingleFermionTerm term;
        

    }


    // An indexing convention is important here.
    public class StateMultiConfigurational
    {


        public Dictionary<SingleFermionTerm, double> terms;
        

    }
}



