// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
    // This class is for Fermion terms that are not grouped into Hermitian bunches.
    // Maybe need a stype for quantum state?

    public class FermionStateMultiConfigurational : FermionHamiltonian
    {
        public FermionStateMultiConfigurational() : base() { }
        public object referenceState;
    }

}



