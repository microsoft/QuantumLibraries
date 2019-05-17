// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;
using System.IO;
using System.Linq;
using Microsoft.Jupyter.Core;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Magic;
using Newtonsoft.Json;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Xunit;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.JordanWigner;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Tests
{
    public class FermionHamiltonianLoadMagicTests
    {
        public (FermionHamiltonianLoadMagic, MockChannel) Init() =>
            (new FermionHamiltonianLoadMagic(), new MockChannel());

        [Fact]
        public void LoadNoInput()
        {
            var (magic, channel) = Init();

            Assert.Equal("%chemistry.fh.load", magic.Name);
            var result = magic.Run("", channel);

            Assert.Equal(ExecuteStatus.Error, result.Status);
        }

        [Fact]
        public void LoadInvalidFile()
        {
            var (magic, channel) = Init();
            var args = JsonConvert.SerializeObject(new FermionHamiltonianLoadMagic.Arguments
            {
                FileName = "foo_bar.yaml"
            });

            Assert.Throws<FileNotFoundException>(() => magic.Run(args, channel));
        }

        [Fact]
        public void LoadFromBroombridgeFile()
        {
            var (magic, channel) = Init();
            var args = JsonConvert.SerializeObject(new FermionHamiltonianLoadMagic.Arguments
            {
                FileName = "broombridge_v0.2.yaml"
            });

            var result = magic.Run(args, channel);
            var hamiltonian = result.Output as FermionHamiltonian;
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            Assert.Equal(12, hamiltonian.SystemIndices.Count);
            Assert.Equal(6, hamiltonian.Terms.Count);
            Assert.Equal(64.5730917943, hamiltonian.Norm());
        }

        /*
        [Fact]
        public void LoadFromProblemDescription()
        {
            var (magic, channel) = Init();
            var broombridgeMagic = new BroombridgeMagic();
            var broombridge = (V0_2.Data)broombridgeMagic.Run("broombridge_v0.2.yaml", channel).Output;

            var args = JsonConvert.SerializeObject(new FermionHamiltonianLoadMagic.Arguments
            {
                ProblemDescription = broombridge.ProblemDescriptions.First()
            });

            var result = magic.Run(args, channel);
            var hamiltonian = result.Output as FermionHamiltonian;
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            Assert.Equal(12, hamiltonian.SystemIndices.Count);
            Assert.Equal(6, hamiltonian.Terms.Count);
            Assert.Equal(64.5730917943, hamiltonian.Norm());
        }
        */
    }

    public class FermionHamiltonianAddTermsMagicTests
    {
        public (FermionHamiltonianAddTermsMagic, MockChannel) Init() =>
            (new FermionHamiltonianAddTermsMagic(), new MockChannel());

        [Fact]
        public void AddHamiltonianTerms()
        {
            var (magic, channel) = Init();
            Assert.Equal("%chemistry.fh.add_terms", magic.Name);

            var original = new FermionHamiltonian();
            var fermionTerms = new List<(int[], double)>
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
                };

            var args = JsonConvert.SerializeObject(new FermionHamiltonianAddTermsMagic.Arguments
            {
                Hamiltonian = original,
                FermionTerms = fermionTerms
            });

            var result = magic.Run(args, channel);
            var hamiltonian = result.Output as FermionHamiltonian;
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            Assert.Equal(16, hamiltonian.CountTerms());
        }

        [Fact]
        public void AddTermsToExistingHamiltonian()
        {
            var (magic, channel) = Init();
            Assert.Equal("%chemistry.fh.add_terms", magic.Name);

            var loadMagic = new FermionHamiltonianLoadMagic();
            var args = JsonConvert.SerializeObject(new FermionHamiltonianLoadMagic.Arguments { FileName = "broombridge_v0.2.yaml" });
            var original = loadMagic.Run(args, channel).Output as FermionHamiltonian;
            var fermionTerms = new List<(int[], double)>
            {
                (new int[] {}, 10.0),   
                (new[] {0,6}, 1.0),
                (new[] {0,2,2,0}, 1.0)
            };

            args = JsonConvert.SerializeObject(new FermionHamiltonianAddTermsMagic.Arguments
            {
                Hamiltonian = original,
                FermionTerms = fermionTerms
            });
            File.WriteAllText("add_terms.args.json", args);

            var result = magic.Run(args, channel);
            var hamiltonian = result.Output as FermionHamiltonian;
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            Assert.True(hamiltonian.CountTerms() > original.CountTerms());
        }
    }
}