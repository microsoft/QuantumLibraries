// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using static System.Math;


namespace Microsoft.Quantum.Chemistry
{
    using HTermArray = QArray<HTerm>;
    using HTermList = List<HTerm>;
    using Microsoft.Quantum.Chemistry.JordanWigner;
    using Microsoft.Extensions.Logging;

    /// <summary>
    /// <para>
    /// Jordan-Wigner representation of a general Fermion Hamiltonian <see cref="FermionHamiltonian"/>.
    /// This representation may only be created from instances of <see cref="FermionHamiltonian"/>,
    /// and stores term data in a format suitable for consumption by Q#,
    /// and optimized for a product formula Hamiltonian simulation algorithm.
    /// </para>
    /// <para>
    /// This supports the following <see cref="FermionTermType"/>: 
    /// <see cref="IdentityTermType"/>,
    /// <see cref="PPTermType"/>,
    /// <see cref="PQTermType"/>,
    /// <see cref="PQQPTermType"/>,
    /// <see cref="PQQRTermType"/>,
    /// <see cref="PQRSTermType"/>.
    /// </para>
    /// <para>
    /// Some optimizations are performed: 
    /// <list type="bullet">
    /// <item>PQQP and PP terms are merged where needed.</item>
    /// <item>PQQR and PR terms are merged where possible.</item>
    /// <item>
    /// All PQRS terms with the same set of spin-orbital indices are performed simultaneously,
    /// and only the XXXX, XXYY, XYXY, YXXY, YYYY, YYXX, YXYX, XYYX terms are only performed as
    /// needed.
    /// </item>
    /// <item>Terms in each group of term types are applied in lexicographic ordering.</item>
    /// </list>
    /// </para>
    /// </summary>
    public partial class JordanWignerEncoding
    {
        // This has contribution from PP and PQQP terms.
        public Double energyOffset = 0.0;
        // This has contribution from PP and PQQP terms.
        private HTermList hZTerms = new HTermList();
        // This has contribution from PQQP terms.
        private HTermList hZZTerms = new HTermList();
        // This has contributions from PQ and PQQR terms.
        private HTermList hPQandPQQRTerms = new HTermList();
        // This has contributions from PQRS terms.
        private HTermList h0123Terms = new HTermList();
        // Number of orbitals
        /// <summary>
        /// <see cref="NOrbitals"/> is the number of distinct orbitals.
        /// <see cref="NSpinOrbitals"/> is the number of distinct spin orbitals.
        /// </summary>
        public Int64 NOrbitals, NSpinOrbitals;
        /// <summary>
        /// This indexes the spin-orbitals to be occupied by electrons.
        /// </summary>
        public QArray<Int64> statePrep = new QArray<Int64>();

        // Supported Fermion term types
        /// <summary>
        /// Container all supported <see cref="FermionTermType"/>s that
        /// may be processed from a <see cref="FermionHamiltonian"/>.
        /// </summary>
        private static HashSet<FermionTermType> SupportedFermionTermTypes = new HashSet<FermionTermType>(new Comparers.FermionTermTypeComparer());

        /// <summary>
        /// Returns term data for consumption by Q#.
        /// </summary>
        public JWOptimizedHTerms Terms
        {
            get => new JWOptimizedHTerms((new HTermArray(hZTerms), new HTermArray(hZZTerms), new HTermArray(hPQandPQQRTerms), new HTermArray(h0123Terms)));
        }

        /// <summary>
        /// Returns data for consumption by Q#.
        /// </summary>
        public JordanWignerEncodingData QSharpData(string selectInputState = "Greedy")
        {
            var inputState = new QArray<JordanWignerInputState>();
            if(selectInputState == "Greedy")
            {
                inputState = new QArray<JordanWignerInputState>(new[] { InputStateFromGreedyAlgorithm });
            }
            else
            {
                inputState = InputStateFromFile[selectInputState];
            }
            return new JordanWignerEncodingData(
                (NSpinOrbitals
                , Terms
                , inputState
                , energyOffset));
        }

        /// <summary>
        /// Constructor for <see cref="JordanWignerEncoding"/>.
        /// </summary>
        /// <param name="nOrbitals">Number of distinct orbitals.</param>
        /// <param name="nSpinOrbitals">Number of distinct spin-orbitals. This usually 2x <paramref name="nOrbitals"/>.</param>
        /// <param name="statePrep">Indices of occupied spin-orbitals.</param>
        /// <param name="globalPhase">Residual global phase after collecting terms.</param>
        /// <param name="termData">Term data for consumption by Q#.</param>
        public JordanWignerEncoding(
            Int64 nOrbitals, Int64 nSpinOrbitals, QArray<Int64> statePrep, Double globalPhase,
            (HTermList, HTermList, HTermList, HTermList) termData
        )
        {
            this.NOrbitals = nOrbitals;
            this.NSpinOrbitals = nSpinOrbitals;
            this.statePrep = statePrep;
            this.energyOffset = globalPhase;
            this.hZTerms = termData.Item1;
            this.hZZTerms = termData.Item2;
            this.hPQandPQQRTerms = termData.Item3;
            this.h0123Terms = termData.Item4;
            
        }

        
        /// <summary>
        /// Counts the total number of distinct spin-orbitals.
        /// </summary>
        /// <param name="currentCount">Current count for the number of distinct spin-orbitals.</param>
        /// <param name="data">Term data proceesed to update count.</param>
        /// <returns>Updated number of distint spin-orbitals.</returns>
        public static Int64 CountNSpinOrbitals(Int64 currentCount, HTermList data)
        {
            if (data.Any()) { 
                var dataMax = data.Select(dataidx => dataidx.Item1.Max() + 1).Max();
                return Max(currentCount, dataMax);
            }
            else
            {
                return currentCount;
            }

        }

        /// <summary>
        /// Constructor of <see cref="JordanWignerEncoding"/> by processing
        /// a Fermion Hamiltonian represnted by <see cref="FermionHamiltonian"/>.
        /// </summary>
        /// <param name="hamiltonian">Fermion Hamiltonian to be processed.</param>
        /// <returns>Representation of <paramref name="hamiltonian"/> suitable for consumption
        /// by Q#, and optimized for a product formula Hamiltonian simulation algorithm.</returns>
        public static JordanWignerEncoding Create(FermionHamiltonian hamiltonian)
        {
            var logger = Logging.LoggerFactory.CreateLogger<JordanWignerEncoding>();

            // These are currently supported Fermion term types.
            HashSet<FermionTermType> SupportedFermionTermTypes = new HashSet<FermionTermType>(new Comparers.FermionTermTypeComparer());
            SupportedFermionTermTypes.Add(FermionTermType.Common.PPTermType);
            SupportedFermionTermTypes.Add(FermionTermType.Common.PQTermType);
            SupportedFermionTermTypes.Add(FermionTermType.Common.PQQPTermType);
            SupportedFermionTermTypes.Add(FermionTermType.Common.PQQRTermType);
            SupportedFermionTermTypes.Add(FermionTermType.Common.PQRSTermType);

            // data from file general Hamiltonian format
            var fermionTermTypes = new HashSet<FermionTermType>(hamiltonian.FermionTerms.Keys);
            var generalFermionTerms = hamiltonian.FermionTerms;
            Int64 NOrbitals = hamiltonian.NOrbitals;
            Int64 NSpinOrbitals = 2 * NOrbitals;
            Int64 nElectrons = hamiltonian.NElectrons;

            // Clear data
            var globalPhase = hamiltonian.EnergyOffset;
            var hZTerms = new HTermList();
            var hZZTerms = new HTermList();
            var hPQandPQQRTerms = new HTermList();
            var h0123Terms = new HTermList();
            var statePrep = new QArray<Int64>();

            logger.LogInformation(
                "Supported term types:\n" +
                String.Join(", ",
                    SupportedFermionTermTypes
                )
            );
            logger.LogInformation(
                "Input term types:\n" +
                String.Join(", ",
                    fermionTermTypes
                )
            );
            // Map Hamiltonian to Jordan-Wigner encoding

            if (fermionTermTypes.IsSubsetOf(SupportedFermionTermTypes))
            {
                if (generalFermionTerms.ContainsKey(FermionTermType.Common.PPTermType))
                {
                    foreach (var generalPPTerm in generalFermionTerms[FermionTermType.Common.PPTermType])
                    {
                        globalPhase += 0.5 * generalPPTerm.coeff;
                        var term = generalPPTerm.SpinOrbitalIndices.Take(1).ToInts(nOrbitals: NOrbitals);
                        hZTerms.Add(new HTerm((new QArray<Int64>(term), new QArray<Double>(new[] { -0.5 * generalPPTerm.coeff }))));
                    }
                }
                
                if (generalFermionTerms.ContainsKey(FermionTermType.Common.PQQPTermType))
                {
                    foreach (var generalPQQPTerm in generalFermionTerms[FermionTermType.Common.PQQPTermType])
                    {
                        var ZZterm = generalPQQPTerm.SpinOrbitalIndices.Take(2).ToInts(nOrbitals: NOrbitals);
                        var Zterm0 = generalPQQPTerm.SpinOrbitalIndices.Take(1).ToInts(nOrbitals: NOrbitals);
                        var Zterm1 = generalPQQPTerm.SpinOrbitalIndices.Skip(1).Take(1).ToInts(nOrbitals: NOrbitals);
                        globalPhase += 0.25 * generalPQQPTerm.coeff;
                        hZTerms.Add(new HTerm((new QArray<Int64>(Zterm0), new QArray<Double> (new []{ -0.25 * generalPQQPTerm.coeff }))));
                        hZTerms.Add(new HTerm((new QArray<Int64>(Zterm1), new QArray<Double> (new []{ -0.25 * generalPQQPTerm.coeff }))));
                        hZZTerms.Add(new HTerm((new QArray<Int64>(ZZterm), new QArray<Double> (new []{ 0.25 * generalPQQPTerm.coeff }))));
                    }
                }
                
                if (generalFermionTerms.ContainsKey(FermionTermType.Common.PQTermType))
                {
                    foreach (var generalPQTerm in generalFermionTerms[FermionTermType.Common.PQTermType])
                    {
                        var pqterm = generalPQTerm.SpinOrbitalIndices.ToInts(nOrbitals: NOrbitals);
                        hPQandPQQRTerms.Add(new HTerm((new QArray<Int64>(pqterm), new QArray<Double> (new[] { 0.25 * generalPQTerm.coeff }))));
                    }
                }
                
                if (generalFermionTerms.ContainsKey(FermionTermType.Common.PQQRTermType))
                {
                    foreach (var generalPQQRTerm in generalFermionTerms[FermionTermType.Common.PQQRTermType])
                    {
                        var pqqrterm = new QArray<Int64>(generalPQQRTerm.SpinOrbitalIndices.ToInts(nOrbitals: NOrbitals));
                        var multiplier = 1.0;
                        
                        if (pqqrterm.First() == pqqrterm.Last())
                        {
                            // This means terms are ordered like QPRQ. So we reorder to PQQR
                            pqqrterm[0] = pqqrterm[1];
                            pqqrterm[1] = pqqrterm[3];
                            pqqrterm[3] = pqqrterm[2];
                            pqqrterm[2] = pqqrterm[1];
                        }
                        else if(pqqrterm[1] == pqqrterm[3])
                        {
                            // This means terms are ordered like PQRQ. So we reorder to PQQR
                            pqqrterm[3] = pqqrterm[2];
                            pqqrterm[2] = pqqrterm[1];
                            multiplier = -1.0;
                        }
                        hPQandPQQRTerms.Add(new HTerm((pqqrterm, new QArray<Double> (new []{ -0.125 * multiplier* generalPQQRTerm.coeff }))));
                        // PQ term
                        var pqterm = new List<Int64>(pqqrterm);
                        pqterm.RemoveRange(1, 2);
                        hPQandPQQRTerms.Add(new HTerm((new QArray<Int64>(pqterm), new QArray<Double>( new[] { 0.125 * multiplier * generalPQQRTerm.coeff }))));
                    }
                }
                
                if (generalFermionTerms.ContainsKey(FermionTermType.Common.PQRSTermType))
                {
                    foreach (var generalPQRSTerm in generalFermionTerms[FermionTermType.Common.PQRSTermType])
                    {
                        var so = generalPQRSTerm.SpinOrbitalIndices.ToInts(nOrbitals: NOrbitals);
                        var pqrsSorted = new List<Int64> { so[0], so[1], so[2], so[3] };
                        pqrsSorted.Sort(new Extensions.IntIComparer());

                        h0123Terms.Add(new HTerm(IdentifyHpqrsPermutation((pqrsSorted, so, generalPQRSTerm.coeff))));
                    }
                }
                

                using (logger.BeginScope("Sorting terms..."))
                {
                    hZTerms.Sort(new Extensions.HTermIndexIComparer());
                    hZTerms = hZTerms.AccumulateTermArray();

                    hZZTerms.Sort(new Extensions.HTermIndexIComparer());
                    hZZTerms = hZZTerms.AccumulateTermArray();

                    hPQandPQQRTerms.Sort(new HPQandPQQRIComparer());
                }
                logger.LogInformation("Sorted terms");


                hPQandPQQRTerms = hPQandPQQRTerms.AccumulateTermArray();

                h0123Terms.Sort(new Extensions.HTermIndexIComparer());
                h0123Terms = h0123Terms.AccumulateTermArray();
                
                NSpinOrbitals = CountNSpinOrbitals(NSpinOrbitals, hZTerms);
                NSpinOrbitals = CountNSpinOrbitals(NSpinOrbitals, hZZTerms);
                NSpinOrbitals = CountNSpinOrbitals(NSpinOrbitals, hPQandPQQRTerms);
                NSpinOrbitals = CountNSpinOrbitals(NSpinOrbitals, h0123Terms);
            }
            else
            {
                logger.LogError("FermionHamiltonian `hamiltonian` contains unsupported terms");
            }


            // Map Hamiltonian input trial states to Jordan-Wigner encoding.
            var jordanWignerEncoding = new JordanWignerEncoding(NOrbitals, NSpinOrbitals, statePrep, globalPhase, (hZTerms, hZZTerms, hPQandPQQRTerms, h0123Terms));
            var greedyState = hamiltonian.GreedyStatePreparation().Superposition.First();
            jordanWignerEncoding.InputStateFromGreedyAlgorithm = jordanWignerEncoding.InitialStatePrep(greedyState.complexCoeff, greedyState.term);

            var inputStateFromFile = new Dictionary<string, QArray<JordanWignerInputState>>();
            foreach(var superposition in hamiltonian.InputStates)
            {
                inputStateFromFile.Add(superposition.Label, jordanWignerEncoding.InitialStatePrep(superposition.Superposition));
            }
            jordanWignerEncoding.InputStateFromFile = inputStateFromFile;

            return jordanWignerEncoding;
        }
        
        /// <summary>
        /// Function for classifying PQRS terms with the same set of
        /// spin-orbital indices.
        /// </summary>
        static (QArray<Int64>, QArray<Double>) IdentifyHpqrsPermutation((List<Int64>, Int64[], Double) term)
        {
            //We only consider permutations pqrs || psqr || prsq || qprs || spqr || prqs
            var (pqrsSorted, pqrsPermuted, coeff) = term;
            coeff = coeff * 0.5 * 0.125;
            var h123 = new QArray<Double> (new[] { .0, .0, .0 });
            var v0123 = new QArray<Double> (new[] { .0, .0, .0, .0 });

            //Console.WriteLine($"{pqrsSorted}, {pqrsPermuted}");

            var prsq = new QArray<Int64> (new[] { pqrsSorted[0], pqrsSorted[2], pqrsSorted[3], pqrsSorted[1] });
            var pqsr = new QArray<Int64> (new[] { pqrsSorted[0], pqrsSorted[1], pqrsSorted[3], pqrsSorted[2] });
            var psrq = new QArray<Int64> (new[] { pqrsSorted[0], pqrsSorted[3], pqrsSorted[2], pqrsSorted[1] });

            //Console.WriteLine($"{pqrs}, {psqr}, {prsq}, {qprs}, {spqr}, {prqs}");

            if (Enumerable.SequenceEqual(pqrsPermuted, prsq))
            {
                h123 = new QArray<Double> (new[] { 0.0, 0.0, coeff });
            }
            else if (Enumerable.SequenceEqual(pqrsPermuted, pqsr))
            {
                h123 = new QArray<Double> (new[] { -coeff, 0.0, 0.0 });
            }
            else if (Enumerable.SequenceEqual(pqrsPermuted, psrq))
            {
                h123 = new QArray<Double> (new[] { 0.0, -coeff, 0.0 });
            }
            else
            {
                h123 = new QArray<Double> (new[] { 0.0, 0.0, 0.0 });
            }

            v0123 = new QArray<Double> (new[] { -h123[0] - h123[1] + h123[2],
                                        h123[0] - h123[1] + h123[2],
                                        -h123[0] - h123[1] - h123[2],
                                        -h123[0] + h123[1] + h123[2] });

            // DEBUG for output h123
            //v0123 = h123;
            //v0123.Add(0.0);

            return (new QArray<Int64>(pqrsSorted), v0123);
        }

        #region Comparer
        /// <summary>
        /// IComparer for sorting PQ and PQQR terms.
        /// </summary>
        public class HPQandPQQRIComparer : IComparer<HTerm>
        {
            public int Compare(HTerm x, HTerm y)
            {

                var xArr = x.Item1.Length == 2 ? new Int64[] { x.Item1[0], x.Item1[1], -1, -1 } : new Int64[] { x.Item1[0], x.Item1[3], x.Item1[2], x.Item1[1] };
                var yArr = y.Item1.Length == 2 ? new Int64[] { y.Item1[0], y.Item1[1], -1, -1 } : new Int64[] { y.Item1[0], y.Item1[3], y.Item1[2], y.Item1[1] };
                return Extensions.CompareIntArray(xArr, yArr);
            }
        }
        #endregion

        /// <summary>
        /// Prints all processed spin-orbitals in Jordan-Wigner
        /// representation.
        /// </summary>
        public void LogSpinOrbitals(LogLevel logLevel = LogLevel.Debug)
        {
            var logger = Logging.LoggerFactory.CreateLogger<JordanWignerEncoding>();

            logger.Log(logLevel, $@"Terms processed. {NOrbitals} Orbitals, {NSpinOrbitals} SpinOrbitals, {energyOffset} phase, {hZTerms.Count()} Z, {hZZTerms.Count()} ZZ, {hPQandPQQRTerms.Count()} PQandPQQR, {h0123Terms.Count()} 0123.");

            logger.Log(logLevel, statePrep.ToString());
           

            foreach (var item in hZTerms)
            {
                logger.Log(logLevel, $"Z: {item.Data}");
            }
            foreach (var item in hZZTerms)
            {
                logger.Log(logLevel, $"ZZ: {item.Data}");
            }
            foreach (var item in hPQandPQQRTerms)
            {
                logger.Log(logLevel, $"PQ and PQQR: {item.Data}");
            }
            foreach (var item in h0123Terms)
            {
                logger.Log(logLevel, $"PQRS: {item.Data}");
            }
        }

    }
    
}