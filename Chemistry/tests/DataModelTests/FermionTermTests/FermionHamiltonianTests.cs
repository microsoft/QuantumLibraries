// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Xunit;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry.Tests
{
    using FT = FermionTermHermitian;

    public class FermionHamiltonianTests
    {
        public FermionHamiltonian GenerateTestHamiltonian()
        {
            var hamiltonian = new FermionHamiltonian();

            List<(FT, double)> fermionTerms = new List<(int[], double)>()
            {
                (new int[] {}, 10.0),
                (new[] {0,0}, 1.0),
                (new[] {1,1}, 1.0),
                (new[] {2,2}, 1.0),
                (new[] {0,2}, 1.0),
                (new[] {1,3}, 1.0),
                (new[] {2,6}, 1.0),
                (new[] {0,2,2,0}, 1.0),
                (new[] {1,3,3,1}, 1.0),
                (new[] {2,6,6,2}, 1.0),
                (new[] {0,2,2,1}, 1.0),
                (new[] {1,3,3,2}, 1.0),
                (new[] {2,6,6,5}, 1.0),
                (new[] {0,2,4,3}, 1.0),
                (new[] {1,4,3,2}, 1.0),
                (new[] {2,4,5,3}, 1.0)
            }.Select(o => (new FT(o.Item1.ToLadderSequence()), o.Item2)).ToList();

            hamiltonian.AddTerms(fermionTerms);
            return hamiltonian;
        }


        [Fact]
        void CountTerms()
        {
            var hamiltonian = GenerateTestHamiltonian();
            var nTerms = 16;
            Assert.Equal(nTerms, hamiltonian.CountTerms());

            hamiltonian.AddHamiltonian(GenerateTestHamiltonian());
            Assert.Equal(nTerms, hamiltonian.CountTerms());

            hamiltonian.AddHamiltonian(GenerateTestHamiltonian());
            Assert.Equal(nTerms, hamiltonian.CountTerms());
        }

        [Fact]
        void NormTerms()
        {
            var hamiltonian = GenerateTestHamiltonian();
            var oneNorm = hamiltonian.Norm();

            hamiltonian.AddHamiltonian(GenerateTestHamiltonian());
            Assert.Equal(oneNorm * 2.0, hamiltonian.Norm(), 5);

            hamiltonian.AddHamiltonian(GenerateTestHamiltonian());
            Assert.Equal(oneNorm * 3.0, hamiltonian.Norm(), 5);

        }

    }
}