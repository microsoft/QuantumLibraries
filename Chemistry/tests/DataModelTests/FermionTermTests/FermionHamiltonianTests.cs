// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Xunit;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;


namespace Microsoft.Quantum.Chemistry.Tests
{
    using FT = FermionTermHermitian;

    public class FermionHamiltonianTests
    {
        public FermionHamiltonian GenerateTestHamiltonian()
        {
            var hamiltonian = new FermionHamiltonian();

            List<(FT, DoubleCoeff)> fermionTerms = new List<(int[], double)>()
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
            }.Select(o => (new FT(o.Item1.ToLadderSequence()), o.Item2.ToDouble())).ToList();

            hamiltonian.AddTerms(fermionTerms);
            return hamiltonian;
        }

        [Fact]
        void CheckKeyPresence()
        {
            var hamiltonian = GenerateTestHamiltonian();
            var check = new FermionTermHermitian(new[] { 0, 2, 2, 0 }.ToLadderSequence());

            hamiltonian.AddTerm(check, 100.0.ToDouble());

            var sourceDict = hamiltonian.terms[TermType.Fermion.PQQP];

            var check2 = new FermionTermHermitian(new[] { 0, 2, 2, 0 }.ToLadderSequence());
            var coeff = sourceDict[check2];

            Assert.Equal(101.0, coeff.Value);
            Assert.Equal(101.0, hamiltonian.terms[TermType.Fermion.PQQP][check].Value);

            // Check using GetTerm method.
            Assert.Equal(101.0, hamiltonian.GetTerm(new FermionTermHermitian(check)).Value);
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

        [Theory]
        [MemberData(nameof(OrbitalsData))]
        public void BuildFromOrbitalIntegrals(
            OrbitalIntegral orbitalIntegral, 
            TermType.Fermion termType, 
            (int, int[], double)[] termsRaw)
        {
            var terms = termsRaw.Select(o => (new FermionTermHermitian(o.Item2.ToLadderSequence()),o.Item3));
            var orbitalhamiltonian = new OrbitalIntegralHamiltonian();
            orbitalhamiltonian.AddTerm(orbitalIntegral);
            var hamiltonian = orbitalhamiltonian.ToFermionHamiltonian(SpinOrbital.IndexConvention.HalfUp);

            Assert.True(hamiltonian.terms.ContainsKey(termType));
            // Check that expected terms are found
            foreach (var term in terms)
            {
                Assert.Contains(term.Item1, hamiltonian.terms[termType].Keys);
                Assert.Equal(term.Item2, hamiltonian.terms[termType][term.Item1].Value);
            }
            // Check only expected terms are found
            foreach (var term in hamiltonian.terms[termType])
            {
                Assert.Contains(term.Key, terms.Select(o => o.Item1));
            }
        }

        public static IEnumerable<object[]> OrbitalsData =>
            new List<object[]>
            {
                new object[] { new OrbitalIntegral(new[] {0,0 },1.0, OrbitalIntegral.Convention.Dirac), TermType.Fermion.PP,
                    new (int, int[], double)[] {
                        (1, new[] { 0, 0 }, 1.0 ),
                        (1, new[] { 1, 1 }, 1.0 )}},
                new object[] { new OrbitalIntegral(new[] {0,1 },1.0, OrbitalIntegral.Convention.Dirac),  TermType.Fermion.PQ,
                    new (int, int[], double)[] {
                        (2, new[] { 0, 1 }, 2.0 ),
                        (2, new[] { 2, 3 }, 2.0 )}},
                new object[] { new OrbitalIntegral(new[] {0,1,1,0 },1.0, OrbitalIntegral.Convention.Dirac),  TermType.Fermion.PQQP,
                    new (int, int[], double)[] {
                        (2, new[] { 0, 1, 1, 0 }, 1.0 ),
                        (2, new[] { 2, 3, 3, 2 }, 1.0 ),
                        (2, new[] { 0, 3, 3, 0 }, 1.0 ),
                        (2, new[] { 1, 2, 2, 1 }, 1.0 )}},
                new object[] { new OrbitalIntegral(new[] {0,1,0,1 },1.0, OrbitalIntegral.Convention.Dirac),  TermType.Fermion.PQQP,
                    new (int, int[], double)[] {
                        (2, new[] { 0, 1, 1, 0 }, -1.0 ),
                        (2, new[] { 2, 3, 3, 2 }, -1.0 ),
                    } },
                new object[] { new OrbitalIntegral(new[] {0,1,0,1 },1.0, OrbitalIntegral.Convention.Dirac),  TermType.Fermion.PQRS,
                    new (int, int[], double)[] {
                        (2, new[] { 0, 3, 2, 1 }, 2.0 ),
                        (2,  new[] { 0, 2, 3, 1 }, 2.0 )
                    } },
                new object[] { new OrbitalIntegral(new[] {0,1,0,0 },1.0, OrbitalIntegral.Convention.Dirac),  TermType.Fermion.PQQR,
                    new (int, int[], double)[] {
                        (2, new[] { 0, 2, 2, 1 }, 2.0 ),
                        (2, new[] { 0, 2, 3, 0 }, 2.0 ),
                    } },
                new object[] {new OrbitalIntegral(new[] {0,0,1,2 },1.0, OrbitalIntegral.Convention.Dirac),  TermType.Fermion.PQQR,
                    new (int, int[], double)[] {
                        (3, new[] { 0, 1, 2, 0 }, -2.0 ),
                        (3, new[] { 3, 4, 5, 3 }, -2.0 ),
                    } },
                new object[] { new OrbitalIntegral(new[] {0,0,1,2 },1.0, OrbitalIntegral.Convention.Dirac),  TermType.Fermion.PQRS,
                    new (int, int[], double)[] {
                        (3, new[] { 0, 3, 4, 2 }, 2.0 ),
                        (3, new[] { 0, 3, 5, 1 }, 2.0 ),
                        (3, new[] { 0, 4, 3, 2 }, 2.0 ),
                        (3, new[] { 0, 5, 3, 1 }, 2.0 ),
                    } },
                new object[] { new OrbitalIntegral(new[] {0,1,2,0 },1.0, OrbitalIntegral.Convention.Dirac),  TermType.Fermion.PQQR,
                    new (int, int[], double)[] {
                        (3, new[] { 0, 1, 2, 0 }, 2.0 ),
                        (3, new[] { 1, 3, 3, 2 }, 2.0 ),
                        (3, new[] { 0, 4, 5, 0 }, 2.0 ),
                        (3, new[] { 3, 4, 5, 3 }, 2.0 ),
                    } },
                new object[] { new OrbitalIntegral(new[] {0,1,2,3 },1.0, OrbitalIntegral.Convention.Dirac),  TermType.Fermion.PQRS,
                    new (int, int[], double)[] {
                        (4, new[] { 0, 5, 6, 3 }, 2.0 ),
                        (4, new[] { 0, 6, 5, 3 }, 2.0 ),
                        (4, new[] { 4, 5, 7, 6 }, -2.0 ),
                        (4, new[] { 0, 1, 3, 2 }, -2.0 ),
                        (4, new[] { 1, 4, 7, 2 }, 2.0 ),
                        (4, new[] { 4, 6, 7, 5 }, -2.0 ),
                        (4, new[] { 0, 2, 3, 1 }, -2.0 ),
                        (4, new[] { 1, 7, 4, 2 }, 2.0 ),
                    } },
                new object[] { new OrbitalIntegral(new[] {3,1,0,2 },1.0, OrbitalIntegral.Convention.Dirac),  TermType.Fermion.PQRS,
                    new (int, int[], double)[] {
                        (4, new[] { 2, 4, 5, 3 }, 2.0 ),
                        (4, new[] { 0, 6, 7, 1 }, 2.0 ),
                        (4, new[] { 4, 6, 7, 5 }, 2.0 ),
                        (4, new[] { 0, 2, 3, 1 }, 2.0 ),
                        (4, new[] { 2, 5, 4, 3 }, 2.0 ),
                        (4, new[] { 0, 7, 6, 1 }, 2.0 ),
                        (4, new[] { 4, 7, 6, 5 }, 2.0 ),
                        (4, new[] { 0, 3, 2, 1 }, 2.0 ),
                    } },
            };


    }
}