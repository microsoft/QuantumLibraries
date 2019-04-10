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
    /// Container containing references to all fermion wavefunction types.
    /// </summary>
    /// <typeparam name="TIndex">Index type used for all fermion operators.</typeparam>
    public class FermionWavefunction<TIndex>
    where TIndex : IEquatable<TIndex>, IComparable<TIndex>
    {
        public StateType Method { get; set; }
        public double Energy { get; set; }
        public SingleCFWavefunction<TIndex> SCFData { get; set; } = new SingleCFWavefunction<TIndex>();
        public SparseMultiCFWavefunction<TIndex> MCFData { get; set; } = new SparseMultiCFWavefunction<TIndex>();
        public UnitaryCCWavefunction<TIndex> UCCData { get; set; } = new UnitaryCCWavefunction<TIndex>();
    }

    public static partial class Extensions
    {

        public static FermionWavefunction<int> ToIndexing(this FermionWavefunction<SpinOrbital> wavefunction, IndexConvention indexConvention)
         => new FermionWavefunction<int>()
         {
             Method = wavefunction.Method,
             Energy = wavefunction.Energy,
             SCFData = wavefunction.SCFData.ToNewIndex<int>((x) => x.ToInt(indexConvention)),
             MCFData = wavefunction.MCFData.ToNewIndex<int>((x) => x.ToInt(indexConvention))
         };
    }
}