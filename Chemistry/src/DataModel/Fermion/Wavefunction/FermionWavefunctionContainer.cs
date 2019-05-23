// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;

using System.Runtime.Serialization.Formatters.Binary;
using System.IO;
using System.IO.Compression;
using YamlDotNet.Serialization;
using Microsoft.Extensions.Logging;
using Microsoft.Quantum.Chemistry;
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
    }

    public static partial class Extensions
    {

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