// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
    // This class is for Fermion terms that are not grouped into Hermitian bunches.
    // Maybe need a stype for quantum state?

    // An indexing convention is important here.
    public class SingleConfigurational
    {

        public SingleFermionTerm term;
        

    }


    // An indexing convention is important here.
    public class StateMultiConfigurational
    {


        public Dictionary<SingleFermionTerm, double> terms;
        

    }
}



