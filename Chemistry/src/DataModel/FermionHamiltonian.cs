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

namespace Microsoft.Quantum.Chemistry
{
    /// <summary>
    /// Representation of a general Fermion Hamiltonian. This Hamiltonian is
    /// assumed to be a linear combination of sequences of creation
    /// and annihilation operators. 
    /// </summary>
    public partial class FermionHamiltonian
    {
        public readonly Configuration Configuration = new Configuration();

        // This will sort terms and accumulate terms in a canonical format.
        // This format is sorted by:
        // 1) List {1,1,...,0,...} of ones for each `a^\dag` followed by zeros for each `a` 
        // 2) Number of unique spin-orbital indices
        // 3) spin-orbital index for a^\dag terms in ascending order, then spin-orbital index for a term in descending order.

        // A Fermionic term has orbital and spin indices.
        // A Fermionic term e.g. a^{\dag}_{4,0} a^{\dag}_{7,1} a_{2,0} a^{\dag}_{9,1} is
        // represented by {{1,1,0,1},{{4,0},{7,1},{2,0},{9,1}}}
        // This is indexes by QArray<(QArray<Int64>)>

        // Every Fermionic term is classified into a term type based on:
        // 1) The list of ones and zeroes representing creation and annihilation operators.
        // 2) The number of distinct spin-orbital indices.
        
        /// <summary>
        ///     A label for this particular Hamiltonian.
        ///     Can be used to identify the Hamiltonian out of set
        ///     loaded from the same file.
        /// </summary>
        public string Name { get; set; } = "<unknown>";

        /// <summary>
        /// Container for all terms in a Fermion Hamiltonian.
        /// </summary>
        public Dictionary<FermionTermType, List<FermionTerm>> FermionTerms = new Dictionary<FermionTermType, List<FermionTerm>>(new Comparers.FermionTermTypeComparer());
        public Int64 NOrbitals = 0;
        public Double EnergyOffset = 0.0;
        public Int64 NElectrons = 0;
        public string MiscellaneousInformation;
        
        /// <summary>
        /// Empty constructor for <see cref="FermionHamiltonian"/>.
        /// </summary>
        public FermionHamiltonian(Int64 nOrbitals = 0, Int64 nElectrons = 0)
        {
            NOrbitals = nOrbitals;
            NElectrons = nElectrons;
        }

        /// <summary>
        /// Constructor for <see cref="FermionHamiltonian"/>.
        /// </summary>
        /// <param name="fermionTerms">Dictionary of all fermion terms.</param>
        /// <param name="nOrbitals">Total number of distinct orbital indices.</param>
        public FermionHamiltonian(Dictionary<FermionTermType, List<FermionTerm>> fermionTerms, Int64 nOrbitals, Int64 nElectrons = 0, Double energyOffset = 0.0)
        {
            FermionTerms = fermionTerms;
            NOrbitals = nOrbitals;
            NElectrons = nElectrons;
            EnergyOffset = energyOffset;
            SortAndAccumulate();
        }
        

        /// <summary>
        /// Method for adding a <see cref="FermionTerm"/> of type <see cref="FermionTermType"/>
        /// to a <see cref="FermionHamiltonian"/>.
        /// </summary>
        /// <param name="termType">Type of fermion term to be added.</param>
        /// <param name="term">Fermion term to be added.</param>
        public void AddFermionTerm(FermionTermType termType, FermionTerm term)
        {
            if (!FermionTerms.ContainsKey(termType))
            {
                FermionTerms.Add(termType, new List<FermionTerm>());
            }
            FermionTerms[termType].Add(term);
        }

        /// <summary>
        /// Method for adding a <see cref="FermionTerm"/> of type <see cref="FermionTermType"/>
        /// to a <see cref="FermionHamiltonian"/>.
        /// </summary>
        /// <param name="termType">Type of fermion term to be added.</param>
        /// <param name="fermionIdxArray">indices of Fermion term to be added.</param>
        public void AddFermionTerm(FermionTermType termType, Int64[] fermionIdxArray, Double coefficient)
        {
            if (!FermionTerms.ContainsKey(termType))
            {
                FermionTerms.Add(termType, new List<FermionTerm>());
            }
            var conjugateArray = termType.GetConjugateSequence();
            var fermionTerm = new FermionTerm(
                nOrbitals: NOrbitals,
                caArray: conjugateArray,
                fermionIdxArray: fermionIdxArray,
                coeffIn: coefficient
                );
            AddFermionTerm(termType, fermionTerm);
        }

        /// <summary>
        /// Method for adding a <see cref="FermionTerm"/>
        /// to a <see cref="FermionHamiltonian"/>.
        /// </summary>
        /// <param name="term">Fermion term to be added.</param>
        public void AddFermionTerm(FermionTerm term)
        {
            AddFermionTerm(term.GetFermionTermType(), term);
        }

        /// <summary>
        /// Method for adding a <see cref="FermionTerm"/>
        /// to a <see cref="FermionHamiltonian"/>.
        /// </summary>
        /// <param name="term">Fermion term to be added.</param>
        public void AddFermionTerm(IEnumerable<SpinOrbital> spinOrbital, Double coefficient)
        {
            AddFermionTerm(new FermionTerm(spinOrbital, coefficient));
        }

        /// <summary>
        /// Method for adding to a <see cref="FermionHamiltonian"/>. 
        /// All <see cref="FermionTerm"/> generated by all symmetries of
        /// <see cref="OrbitalIntegral.EnumerateOrbitalSymmetries"/> and
        /// <see cref="OrbitalIntegral.EnumerateSpinOrbitals"/>.
        /// </summary>
        /// <param name="termType"><see cref="OrbitalIntegral"/> representing terms to be added.</param>
        public void AddFermionTerm(OrbitalIntegral orbitalIntgral)
        {
            if(orbitalIntgral.Length() == 2)
            {
                CreateTwoBodySpinOrbitalTerms(orbitalIntgral);
            }
            else if(orbitalIntgral.Length() == 4)
            {
                CreateFourBodySpinOrbitalTerms(orbitalIntgral);
            }
            else
            {
                throw new System.NotImplementedException();
            }
        }

        /// <summary>
        /// Method for adding to a <see cref="FermionHamiltonian"/>. 
        /// All <see cref="FermionTerm"/> generated by all symmetries of
        /// <see cref="OrbitalIntegral.EnumerateOrbitalSymmetries"/> and
        /// <see cref="OrbitalIntegral.EnumerateSpinOrbitals"/>.
        /// </summary>
        /// <param name="termType"><see cref="OrbitalIntegral"/> representing terms to be added.</param>
        public void AddFermionTerm(IEnumerable<OrbitalIntegral> orbitalIntgrals)
        {
            foreach(var orbitalIntegral in orbitalIntgrals)
            {
                AddFermionTerm(orbitalIntegral);
            }
        }

        #region Create canonical fermion terms from orbitals
        /// <summary>
        /// Updates an instance of <see cref="FermionHamiltonian"/>
        /// with all spin-orbitals from described by a sequence of two-body orbital integrals.
        /// </summary>
        /// <param name="nOrbitals">Total number of distinct orbitals.</param>
        /// <param name="hpqTerms">Sequence of two-body orbital integrals.</param>
        /// <param name="hamiltonian">Fermion Hamiltonian to be updated.</param>
        public void CreateTwoBodySpinOrbitalTerms(OrbitalIntegral orbitalIntegral)
        {
            // One-electron orbital integral symmetries
            // ij = ji
            var pqSpinOrbitals = orbitalIntegral.EnumerateOrbitalSymmetries().EnumerateSpinOrbitals();

            var coefficient = orbitalIntegral.Coefficient;
            
            foreach (var pq in pqSpinOrbitals)
            {
                var p = pq[0];
                var q = pq[1];

                var pInt = p.ToInt();
                var qInt = q.ToInt();
                if (pInt == qInt)
                {
                    AddFermionTerm(FermionTermType.Common.PPTermType, new FermionTerm { CreationAnnihilationIndices = new Int64[] { 1, 0 }, SpinOrbitalIndices = pq, coeff = coefficient });
                }
                else if(pInt < qInt)
                {
                    AddFermionTerm(FermionTermType.Common.PQTermType, new FermionTerm { CreationAnnihilationIndices = new Int64[] { 1, 0 }, SpinOrbitalIndices = pq, coeff = 2.0 * coefficient });
                }
            }
        }

        /// <summary>
        /// Updates an instance of <see cref="FermionHamiltonian"/>
        /// with all spin-orbitals from described by a sequence of four-body orbital integrals.
        /// </summary>
        /// <param name="nOrbitals">Total number of distinct orbitals.</param>
        /// <param name="rawPQRSTerms">Sequence of four-body orbital integrals.</param>
        /// <param name="hamiltonian">Fermion Hamiltonian to be updated.</param>
        public void CreateFourBodySpinOrbitalTerms(OrbitalIntegral orbitalIntegral)
        {
            // Two-electron orbital integral symmetries
            // ijkl = lkji = jilk = klij = ikjl = ljki = kilj = jlik.
            var pqrsSpinOrbitals = orbitalIntegral.EnumerateOrbitalSymmetries().EnumerateSpinOrbitals();
            var coefficient = orbitalIntegral.Coefficient;


            // We only need to see one of these.
            // Now iterate over pqrsArray
            foreach (var pqrs in pqrsSpinOrbitals)
            {
                var p = pqrs[0];
                var q = pqrs[1];
                var r = pqrs[2];
                var s = pqrs[3];

                var pInt = p.ToInt();
                var qInt = q.ToInt();
                var rInt = r.ToInt();
                var sInt = s.ToInt();

                // Only consider terms on the lower diagonal due to Hermitian symmetry.

                // For terms with two different orbital indices, possibilities are
                // PPQQ (QQ = 0), PQPQ, QPPQ (p<q), PQQP, QPQP (p<q), QQPP (PP=0)
                // Hence, if we only count PQQP, and PQPQ, we need to double the coefficient.
                // iU jU jU iU | iU jD jD iD | iD jU jU iD | iD jD jD iD
                if (pInt == sInt && qInt == rInt && pInt < qInt)
                {   // PQQP
                    AddFermionTerm(FermionTermType.Common.PQQPTermType, new FermionTerm { CreationAnnihilationIndices = new Int64[] { 1, 1, 0, 0 }, SpinOrbitalIndices = new SpinOrbital[] { p, q, r, s }, coeff = 1.0 * coefficient });
                }
                else if (pInt == rInt && qInt == sInt && pInt < qInt)
                {
                    // iU jU iU jU | iD jD iD jD
                    // PQPQ
                    AddFermionTerm(FermionTermType.Common.PQQPTermType, new FermionTerm { CreationAnnihilationIndices = new Int64[] { 1, 1, 0, 0 }, SpinOrbitalIndices = new SpinOrbital[] { p, q, s, r }, coeff = -1.0 * coefficient });
                }
                else if (qInt == rInt && pInt < sInt && rInt != sInt && pInt != qInt)
                {
                    // PQQR
                    // For any distinct pqr, [i;j;j;k] generates PQQR ~ RQQP ~ QPRQ ~ QRPQ. We only need to record one.
                    if (rInt < sInt)
                    {
                        if (pInt < qInt)
                        {
                            AddFermionTerm(FermionTermType.Common.PQQRTermType, new FermionTerm { CreationAnnihilationIndices = new Int64[] { 1, 1, 0, 0 }, SpinOrbitalIndices = new SpinOrbital[] { p, q, s, r }, coeff = -2.0 * coefficient });
                        }
                        else
                        {
                            AddFermionTerm(FermionTermType.Common.PQQRTermType, new FermionTerm { CreationAnnihilationIndices = new Int64[] { 1, 1, 0, 0 }, SpinOrbitalIndices = new SpinOrbital[] { q, p, s, r }, coeff = 2.0 * coefficient });
                        }

                    }
                    else
                    {
                        if (pInt < qInt)
                        {
                            AddFermionTerm(FermionTermType.Common.PQQRTermType, new FermionTerm { CreationAnnihilationIndices = new Int64[] { 1, 1, 0, 0 }, SpinOrbitalIndices = new SpinOrbital[] { p, q, r, s }, coeff = 2.0 * coefficient });
                        }
                        else
                        {
                            AddFermionTerm(FermionTermType.Common.PQQRTermType, new FermionTerm { CreationAnnihilationIndices = new Int64[] { 1, 1, 0, 0 }, SpinOrbitalIndices = new SpinOrbital[] { q, p, r, s }, coeff = -2.0 * coefficient });
                        }
                    }
                }
                else if (qInt == sInt && pInt < rInt && rInt != sInt && pInt != sInt)
                {
                    // PQRQ
                    // For any distinct pqr, [i;j;k;j] generates {p, q, r, q}, {q, r, q, p}, {q, p, q, r}, {r, q, p, q}. We only need to record one.
                    if (pInt < qInt)
                    {
                        if (rInt > qInt)
                        {
                            AddFermionTerm(FermionTermType.Common.PQQRTermType, new FermionTerm { CreationAnnihilationIndices = new Int64[] { 1, 1, 0, 0 }, SpinOrbitalIndices = new SpinOrbital[] { p, q, r, s }, coeff = 2.0 * coefficient });
                        }
                        else
                        {
                            AddFermionTerm(FermionTermType.Common.PQQRTermType, new FermionTerm { CreationAnnihilationIndices = new Int64[] { 1, 1, 0, 0 }, SpinOrbitalIndices = new SpinOrbital[] { p, q, s, r }, coeff = -2.0 * coefficient });
                        }
                    }
                    else
                    {
                        AddFermionTerm(FermionTermType.Common.PQQRTermType, new FermionTerm { CreationAnnihilationIndices = new Int64[] { 1, 1, 0, 0 }, SpinOrbitalIndices = new SpinOrbital[] { q, p, r, s }, coeff = -2.0 * coefficient });
                    }
                }
                else if (pInt < qInt && pInt < rInt && pInt < sInt && qInt != rInt && qInt != sInt && rInt != sInt)
                {
                    // PQRS
                    // For any distinct pqrs, [i;j;k;l] generates 
                    // {{p, q, r, s}<->{s, r, q, p}<->{q, p, s, r}<->{r, s, p, q}, 
                    // {1,2,3,4}<->{4,3,2,1}<->{2,1,4,3}<->{3,4,1,2}
                    // {p, r, q, s}<->{s, q, r, p}<->{r, p, s, q}<->{q, s, p, r}}
                    // 1324, 4231, 3142, 2413
                    if (rInt < sInt)
                    {
                        AddFermionTerm(FermionTermType.Common.PQRSTermType, new FermionTerm { CreationAnnihilationIndices = new Int64[] { 1, 1, 0, 0 }, SpinOrbitalIndices = new SpinOrbital[] { p, q, s, r }, coeff = -2.0 * coefficient });
                    }
                    else
                    {
                        AddFermionTerm(FermionTermType.Common.PQRSTermType, new FermionTerm { CreationAnnihilationIndices = new Int64[] { 1, 1, 0, 0 }, SpinOrbitalIndices = new SpinOrbital[] { p, q, r, s }, coeff = 2.0 * coefficient });
                    }
                }
            }
        }

        #endregion

        /// <summary>
        /// Checks that all terms in a <see cref="FermionHamiltonian"/>
        /// are in canonical order. 
        /// </summary>
        /// <returns>
        /// Returns <c>true</c> if all terms are in canonical order, and <c>false</c><
        /// otherwise.
        /// /returns>
        public bool VerifyFermionTerms()
        {
            // Check that only valid term types are present.
            foreach (var termType in FermionTerms.Keys)
            {
                if (!termType.IsInCanonicalOrder())
                {
                    return false;
                }
            }
            // Then check every element.
            foreach (var termTypeKV in FermionTerms)
            {
                foreach(var term in termTypeKV.Value)
                {
                    if (!term.CreationAnnihilationIndices.SequenceEqual(termTypeKV.Key.type.Item2))
                    {
                        return false;
                    }
                    
                    if (!term.IsInCanonicalOrder())
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        /// <summary>
        ///  Converts a <see cref="FermionHamiltonian"/> to canonical order. This generates
        ///  new terms and modifies the coefficient as needed.
        /// </summary>
        public void ToCanonicalOrder()
        {
            // For now, this assumes that the creation and annihilation operators are in canonical order.
            // Thus anticommutation here only introduces a minus sign.
            // This will loop over FermionTerm.ToCanonicalOrder();
            // This will then call SortAndAccumulate();
            throw new System.NotImplementedException();
        }

        /// <summary>
        /// Assuming that all <see cref="FermionTerm"/>s of a <see cref="FermionHamiltonian"/> are
        /// in canonical order, this sorts all terms according to <see cref="Comparers.FermionTermIComparer"/>
        /// and combines duplicates term types.
        /// </summary>
        public void SortAndAccumulate()
        {
            foreach (var termType in FermionTerms)
            {
                // This assumes that Fermion terms are in canonical order.
                var NOnes = termType.Key.type.Item2.Where(o => o == 1).ToArray().Length;
                termType.Value.Sort(new Comparers.FermionTermIComparer(nOnesIn: NOnes));
            }
            foreach (var termType in new HashSet<FermionTermType>(FermionTerms.Keys))
            {
                FermionTerms[termType] = FermionTerms[termType].AccumulateFermionTerm();
            }
        }

        /// <summary>
        /// String representation of Fermion Hamiltonian.
        /// </summary>
        /// <returns>String representation of Fermion Hamiltonian.</returns>
        public override string ToString()
        {
            var output = "";
            foreach (var termType in FermionTerms)
            {
                output += $"(FermionTermType, Number of entries): ({termType.Key}, {termType.Value.Count()}).\n";
                foreach(var term in termType.Value)
                {
                    output += $"{term}\n";
                }
            }
            return output;
        }


        /// <summary>
        ///      Prints all Hamiltonian terms.
        /// </summary>
        public void LogSpinOrbitals(LogLevel logLevel = LogLevel.Debug)
        {
            var logger = Logging.LoggerFactory.CreateLogger<FermionHamiltonian>();
            using (logger.BeginScope("Logging spin orbitals"))
            {
                foreach (var termType in FermionTerms)
                {
                    logger.Log(logLevel, $"FermionTermType: {termType.Key} has {termType.Value.Count()} entries");
                    foreach (var term in termType.Value)
                    {
                        logger.Log(logLevel, $"{term}");
                    }
                }
            }
        }
        
        public void SaveToBinary(string filename)
        {
            using (FileStream write = new FileStream(filename, FileMode.Create, FileAccess.Write))
            {
                using (GZipStream zip = new GZipStream(write, CompressionMode.Compress))
                {
                    BinaryFormatter binary = new BinaryFormatter();
                    binary.Serialize(zip, (FermionTerms, NOrbitals));
                }
            }
        }

        public static FermionHamiltonian LoadFromBinary(string filename)
        {
            using (FileStream write = new FileStream(filename, FileMode.Open, FileAccess.Read))
            {
                using (GZipStream zip = new GZipStream(write, CompressionMode.Decompress))
                {
                    BinaryFormatter binary = new BinaryFormatter();
                    var (FermionTerms, NOrbitals) = ((Dictionary<FermionTermType, List<FermionTerm>>, Int64))(binary.Deserialize(zip));
                   return new FermionHamiltonian(FermionTerms, NOrbitals);
                }
            }
        }


    }

}