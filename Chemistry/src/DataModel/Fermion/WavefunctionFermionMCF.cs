// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.LadderOperators;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    // This class is for Fermion terms that are not grouped into Hermitian bunches.
    // Maybe need a stype for quantum state?

    /*public class WavefunctionFermionMCF
    {
        public List<((double, double), WavefunctionFermionSCF)> Superposition = new List<((double, double), WavefunctionFermionSCF)>();
        public WavefunctionFermionSCF reference = new WavefunctionFermionSCF();
    }*/

    // To do: work in progress
    public class WavefunctionFermionMCFandUCCPlaceholder
    {
        public List<((double, double), IndexOrderedLadderSequence)> Superposition = new List<((double, double), IndexOrderedLadderSequence)>();
        public WavefunctionFermionSCF reference = new WavefunctionFermionSCF();
    }

}



