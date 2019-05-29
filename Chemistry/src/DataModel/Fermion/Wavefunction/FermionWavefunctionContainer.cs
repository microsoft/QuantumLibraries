// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Numerics;
using System.Collections.Generic;

using System.Runtime.Serialization.Formatters.Binary;
using System.IO;
using System.IO.Compression;
using YamlDotNet.Serialization;
using Microsoft.Extensions.Logging;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    /// <summary>
    /// Class for storing any fermion wavefunction type.
    /// </summary>
    /// <typeparam name="TIndex">Index type used for all fermion operators.</typeparam>
    public class FermionWavefunction<TIndex>
    where TIndex : IEquatable<TIndex>, IComparable<TIndex>
    {
        /// <summary>
        /// Type of wavefunction.
        /// </summary>
        public StateType Method { get; set; } = StateType.NotRecognized;
        /// <summary>
        /// Energy of wavefunction relative to a Hamiltonian containing this object.
        /// </summary>
        public double Energy { get; set; } = 0;
        /// <summary>
        /// Single configurational wavefunction data.
        /// </summary>
        public SingleCFWavefunction<TIndex> SCFData { get; set; } = new SingleCFWavefunction<TIndex>();
        /// <summary>
        /// Sparse multi configurational wavefunction data.
        /// </summary>
        public SparseMultiCFWavefunction<TIndex> MCFData { get; set; } = new SparseMultiCFWavefunction<TIndex>();
        /// <summary>
        /// Unitary coupled-cluster wavefunction data.
        /// </summary>
        public UnitaryCCWavefunction<TIndex> UCCData { get; set; } = new UnitaryCCWavefunction<TIndex>();

        /// <summary>
        /// Constructor for empty fermion wavefunction object.
        /// </summary>
        internal FermionWavefunction(){}


        /// <summary>
        /// Construct a single-reference wave function
        /// </summary>
        /// <param name="term">
        /// Sequence of indices of creation operators acting
        /// on the vacuum state.
        /// </param>
        public FermionWavefunction(IEnumerable<TIndex> term)
        {
            // This is deliberately set to SparseMultiConfigurational
            // instead of SingleConfigurational.
            Method = StateType.SparseMultiConfigurational;

            var singleReferenceWavefunction = new SingleCFWavefunction<TIndex>(term);

            // We use Single() to throw an exception if multiple terms
            // are specified for this type of state.
            MCFData.Set(singleReferenceWavefunction, new Complex(1.0, 0.0));
        }

        /// <summary>
        /// Construct a sparse multi-reference wave function
        /// </summary>
        /// <param name="terms">
        /// List of tuples specifying an unnormalized superposition of
        /// basis states. The first item of each tuple is a list of 
        /// indices to the creation operator sequence acting on the vacuum state.
        /// The second item of each tuple is the unnormalized amplitude of the
        /// specified basis state.
        /// </param>
        public FermionWavefunction(IEnumerable<(TIndex[], double)> terms) {
            Method = StateType.SparseMultiConfigurational;
            foreach (var term in terms) {
                MCFData.Set(new SingleCFWavefunction<TIndex>(term.Item1), new Complex(term.Item2, 0.0));    
            }
        }

        /// <summary>
        /// Construct a unitary coupled-cluster wave function represented
        /// by a unitary coupled-cluster operator acting on a single-reference
        /// state.
        /// </summary>
        /// <param name="reference"> Sequence of indices of creation operators acting
        /// on the vacuum state.</param>
        /// <param name="excitations"></param>
        public FermionWavefunction(
            IEnumerable<TIndex> reference,
            IEnumerable<(TIndex[], double)> excitations
            )
        {
            Method = StateType.UnitaryCoupledCluster;

            UCCData.Reference = new SingleCFWavefunction<TIndex>(reference);

            foreach (var term in excitations)
            {
                UCCData.Set(new IndexOrderedSequence<TIndex>(
                    new LadderSequence<TIndex>(term.Item1)), 
                    new Complex(term.Item2, 0.0));
            }
        }

    }

    public static partial class Extensions
    {
        /// <summary>
        /// Convert spin-orbital indices to integer indices
        /// </summary>
        /// <param name="wavefunction">Change indexing scheme of this wavefunction.</param>
        /// <param name="indexConvention">Convention for mapping spin-orbitals to indices.</param>
        /// <returns>A fermion wavefunction where spin-orbitals are indexed by integers
        /// according to the chosen indexing scheme.
        /// </returns>
        public static FermionWavefunction<int> ToIndexing(this FermionWavefunction<SpinOrbital> wavefunction, IndexConvention indexConvention)
         => new FermionWavefunction<int>()
         {
             Method = wavefunction.Method,
             Energy = wavefunction.Energy,
             SCFData = wavefunction.SCFData.SelectIndex<int>((x) => x.ToInt(indexConvention)),
             MCFData = wavefunction.MCFData.SelectIndex<int>((x) => x.ToInt(indexConvention)),
             UCCData = wavefunction.UCCData.SelectIndex<int>((x) => x.ToInt(indexConvention))
         };

    }



}