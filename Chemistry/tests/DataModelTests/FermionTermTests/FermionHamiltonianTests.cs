// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;
using System.IO;
using System.Linq;

using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

using Newtonsoft.Json;

using Xunit;

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

            var sourceDict = hamiltonian.Terms[TermType.Fermion.PQQP];

            var check2 = new FermionTermHermitian(new[] { 0, 2, 2, 0 }.ToLadderSequence());
            var coeff = sourceDict[check2];

            Assert.Equal(101.0, coeff.Value);
            Assert.Equal(101.0, hamiltonian.Terms[TermType.Fermion.PQQP][check].Value);

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


        [Fact]
        public void JsonEncoding()
        {
            var filename = "Broombridge/broombridge_v0.2.yaml";
            CurrentVersion.Data broombridge = Deserializers.DeserializeBroombridge(filename);
            CurrentVersion.ProblemDescription problemData = broombridge.ProblemDescriptions.First();

            OrbitalIntegralHamiltonian orbitalIntegralHamiltonian = problemData.CreateOrbitalIntegralHamiltonian();
            FermionHamiltonian original = orbitalIntegralHamiltonian.ToFermionHamiltonian(SpinOrbital.IndexConvention.HalfUp);

            var json = JsonConvert.SerializeObject(original);
            File.WriteAllText("fermion.original.json", json);

            var serialized = JsonConvert.DeserializeObject<FermionHamiltonian>(json);
            File.WriteAllText("fermion.serialized.json", JsonConvert.SerializeObject(serialized));

            Assert.Equal(original.SystemIndices.Count, serialized.SystemIndices.Count);
            Assert.Equal(original.Terms.Count, serialized.Terms.Count);
            Assert.Equal(original.Norm(), serialized.Norm());
            Assert.Equal(original.ToString(), serialized.ToString());
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

            Assert.True(hamiltonian.Terms.ContainsKey(termType));
            // Check that expected terms are found
            foreach (var term in terms)
            {
                Assert.Contains(term.Item1, hamiltonian.Terms[termType].Keys);
                Assert.Equal(term.Item2, hamiltonian.Terms[termType][term.Item1].Value);
            }
            // Check only expected terms are found
            foreach (var term in hamiltonian.Terms[termType])
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