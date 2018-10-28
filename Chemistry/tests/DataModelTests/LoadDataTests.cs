// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.


using Xunit;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;
using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;
using System.IO;

using YamlDotNet.Core;
using YamlDotNet.Serialization;

namespace Microsoft.Quantum.Chemistry.Tests
{
    using static FermionTermType.Common;
    using FermionTerm = FermionTerm;
    using FermionTermType = FermionTermType;
    using SpinOrbital = SpinOrbital;

    public class LoadFromYAMLTests
    {
        /*
        [Fact]
        public void LoadFromYAMLTest()
        {
            var filename = "resources/yamltest.yaml";

            using (var reader = File.OpenText(filename))
            {
                var deserializer = new DeserializerBuilder().Build();
                var yamlData = deserializer.Deserialize<IntegralDataSchema>(reader);

                //yamlData.IntegralSets
            }

            

            //var hamiltonian = FermionHamiltonian.LoadFromYAML(filename);
        }*/
    }

    public class LoadFromLiquidTests
    {
        [Theory]
        [MemberData(nameof(LiquidOrbitalsData))]
        public void LoadFromLiquidTest(string line, FermionTermType termType, FermionTerm[] terms)
        {
            //var test = terms.Item1;
            //string[] lines, FermionTermType termType, FermionTerm[] terms
            var hamiltonian = LoadData.LoadFromLiquid(line);
            Assert.True(hamiltonian.FermionTerms.ContainsKey(termType));
            // Check that expected terms are found
            foreach(var term in terms)
            {
                Assert.Contains(term, hamiltonian.FermionTerms[termType], new Comparers.FermionTermComparer());
            }
            // Check only expected terms are found
            foreach(var term in hamiltonian.FermionTerms[termType])
            {
                Assert.Contains(term, terms, new Comparers.FermionTermComparer());
            }
            // Verify that each term is in the canonical order.
            Assert.True(hamiltonian.VerifyFermionTerms());
        }

        public static IEnumerable<object[]> LiquidOrbitalsData =>
            new List<object[]>
            {
                new object[] { "0,0=1.0", PPTermType,
                    new FermionTerm[] {
                        new FermionTerm(1, new Int64[] {1,0}, new Int64[] { 0, 0 }, 1.0 ),
                        new FermionTerm(1, new Int64[] {1,0}, new Int64[] { 1, 1 }, 1.0 )}},
                new object[] { "0,1=1.0", PQTermType,
                    new FermionTerm[] {
                        new FermionTerm(2, new Int64[] {1,0}, new Int64[] { 0, 1 }, 2.0 ),
                        new FermionTerm(2, new Int64[] {1,0}, new Int64[] { 2, 3 }, 2.0 )}},
                new object[] { "0,1,1,0=1.0", PQQPTermType,
                    new FermionTerm[] {
                        new FermionTerm(2, new Int64[] {1,1,0,0}, new Int64[] { 0, 1, 1, 0 }, 1.0 ),
                        new FermionTerm(2, new Int64[] {1,1,0,0}, new Int64[] { 2, 3, 3, 2 }, 1.0 ),
                        new FermionTerm(2, new Int64[] {1,1,0,0}, new Int64[] { 0, 3, 3, 0 }, 1.0 ),
                        new FermionTerm(2, new Int64[] {1,1,0,0}, new Int64[] { 1, 2, 2, 1 }, 1.0 )}},
                new object[] { "0,1,0,1=1.0", PQQPTermType,
                    new FermionTerm[] {
                        new FermionTerm(2, new Int64[] {1,1,0,0}, new Int64[] { 0, 1, 1, 0 }, -1.0 ),
                        new FermionTerm(2, new Int64[] {1,1,0,0}, new Int64[] { 2, 3, 3, 2 }, -1.0 ),
                    } },
                new object[] { "0,1,0,1=1.0", PQRSTermType,
                    new FermionTerm[] {
                        new FermionTerm(2, new Int64[] {1,1,0,0}, new Int64[] { 0, 3, 2, 1 }, 2.0 ),
                        new FermionTerm(2, new Int64[] {1,1,0,0}, new Int64[] { 0, 2, 3, 1 }, 2.0 )
                    } },
                new object[] { "0,1,0,0=1.0", PQQRTermType,
                    new FermionTerm[] {
                        new FermionTerm(2, new Int64[] {1,1,0,0}, new Int64[] { 0, 2, 2, 1 }, 2.0 ),
                        new FermionTerm(2, new Int64[] {1,1,0,0}, new Int64[] { 0, 2, 3, 0 }, 2.0 ),
                    } },
                new object[] { "0,0,1,2=1.0", PQQRTermType,
                    new FermionTerm[] {
                        new FermionTerm(3, new Int64[] {1,1,0,0}, new Int64[] { 0, 1, 2, 0 }, -2.0 ),
                        new FermionTerm(3, new Int64[] {1,1,0,0}, new Int64[] { 3, 4, 5, 3 }, -2.0 ),
                    } },
                new object[] { "0,0,1,2=1.0", PQRSTermType,
                    new FermionTerm[] {
                        new FermionTerm(3, new Int64[] {1,1,0,0}, new Int64[] { 0, 3, 4, 2 }, 2.0 ),
                        new FermionTerm(3, new Int64[] {1,1,0,0}, new Int64[] { 0, 3, 5, 1 }, 2.0 ),
                        new FermionTerm(3, new Int64[] {1,1,0,0}, new Int64[] { 0, 4, 3, 2 }, 2.0 ),
                        new FermionTerm(3, new Int64[] {1,1,0,0}, new Int64[] { 0, 5, 3, 1 }, 2.0 ),
                    } },
                new object[] { "0,1,2,0=1.0", PQQRTermType,
                    new FermionTerm[] {
                        new FermionTerm(3, new Int64[] {1,1,0,0}, new Int64[] { 0, 1, 2, 0 }, 2.0 ),
                        new FermionTerm(3, new Int64[] {1,1,0,0}, new Int64[] { 1, 3, 3, 2 }, 2.0 ),
                        new FermionTerm(3, new Int64[] {1,1,0,0}, new Int64[] { 0, 4, 5, 0 }, 2.0 ),
                        new FermionTerm(3, new Int64[] {1,1,0,0}, new Int64[] { 3, 4, 5, 3 }, 2.0 ),
                    } },
                new object[] { "0,1,2,3=1.0", PQRSTermType,
                    new FermionTerm[] {
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 0, 5, 6, 3 }, 2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 0, 6, 5, 3 }, 2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 4, 5, 7, 6 }, -2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 0, 1, 3, 2 }, -2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 1, 4, 7, 2 }, 2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 4, 6, 7, 5 }, -2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 0, 2, 3, 1 }, -2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 1, 7, 4, 2 }, 2.0 ),
                    } },
                new object[] { "3,1,0,2=1.0", PQRSTermType,
                    new FermionTerm[] {
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 2, 4, 5, 3 }, 2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 0, 6, 7, 1 }, 2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 4, 6, 7, 5 }, 2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 0, 2, 3, 1 }, 2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 2, 5, 4, 3 }, 2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 0, 7, 6, 1 }, 2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 4, 7, 6, 5 }, 2.0 ),
                        new FermionTerm(4, new Int64[] {1,1,0,0}, new Int64[] { 0, 3, 2, 1 }, 2.0 ),
                    } },
            };
        

    }



}