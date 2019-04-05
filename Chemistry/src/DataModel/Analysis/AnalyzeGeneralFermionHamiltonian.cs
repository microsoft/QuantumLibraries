// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using Microsoft.Extensions.Logging;
using static System.Math;

namespace Microsoft.Quantum.Chemistry
{

    
    /// <summary>
    /// Representation of a general Fermion Hamiltonian. This Hamiltonian is
    /// assumed to be a linear combination of sequences of creation
    /// and annihilation operators. 
    /// </summary>
    public partial class FermionHamiltonian
    {
        public Dictionary<FermionTermType, Double> ComputeOneNorms()
        {
            return FermionTerms.ToDictionary(
                (termTypePair) => termTypePair
                    .Key,
                (termTypePair) => termTypePair
                    .Value
                    .AsParallel()
                    .Select(o => Abs(o.coeff))
                    .Sum()
            );
        }
       
    }
    
    
}