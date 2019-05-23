﻿// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;
using System.Numerics;

using Microsoft.Quantum.Chemistry.LadderOperators;

using Microsoft.Quantum.Chemistry.Fermion;

namespace Microsoft.Quantum.Chemistry.OrbitalIntegrals
{
    /// <summary>
    /// Extensions for converting orbital integrals to fermion terms.
    /// </summary>
    public static partial class Extensions
    {

        /// <summary>
        /// Method for constructing a fermion Hamiltonian from an orbital integral Hamiltonina.
        /// </summary>
        /// <param name="sourceHamiltonian">Input orbital integral Hamiltonian.</param>
        /// <param name="indexConvention">Indexing scheme from spin-orbitals to integers.</param>
        /// <returns>Fermion Hamiltonian constructed from orbital integrals.</returns>
        public static FermionHamiltonian ToFermionHamiltonian(
            this OrbitalIntegralHamiltonian sourceHamiltonian, 
            IndexConvention indexConvention)
        {
            var nOrbitals = sourceHamiltonian.SystemIndices.Max() + 1;
            var hamiltonian = new FermionHamiltonian();
            Func<OrbitalIntegral, double, IEnumerable<(HermitianFermionTerm, DoubleCoeff)>> conversion = 
                (orb, coeff) => new OrbitalIntegral(orb.OrbitalIndices, coeff).ToHermitianFermionTerms(nOrbitals, indexConvention)
                .Select(o => (o.Item1, o.Item2.ToDoubleCoeff()));

            foreach (var termType in sourceHamiltonian.Terms)
            {
                foreach(var term in termType.Value)
                {
                    hamiltonian.AddRange(conversion(term.Key, term.Value.Value));
                }
            }
            // Number of fermions is twice the number of orbitals.
            hamiltonian.SystemIndices = new HashSet<int>(Enumerable.Range(0, 2 * nOrbitals));
            return hamiltonian;
        }


        /// <summary>
        /// Creates all fermion terms generated by all symmetries of an orbital integral.
        /// </summary>
        /// <param name="nOrbitals">Total number of distinct orbitals.</param>
        /// <param name="orbitalIntegral">Input orbital integral.</param>
        /// <param name="indexConvention">Indexing scheme from spin-orbitals to integers.</param>
        /// <returns>List of fermion terms generated by all symmetries of an orbital integral.</returns>
        public static IEnumerable<(HermitianFermionTerm, double)> ToHermitianFermionTerms(
            this OrbitalIntegral orbitalIntegral,
            int nOrbitals,
            IndexConvention indexConvention)
        {
            var termType = orbitalIntegral.TermType;
            if (termType == TermType.OrbitalIntegral.OneBody)
            {
                return orbitalIntegral.ToOneBodySpinOrbitalTerms(nOrbitals, indexConvention);
            }
            else if (termType == TermType.OrbitalIntegral.TwoBody)
            {
                return orbitalIntegral.ToTwoBodySpinOrbitalTerms(nOrbitals, indexConvention);
            }
            else if(termType == TermType.OrbitalIntegral.Identity)
            {
                return new List<(HermitianFermionTerm, double)>() { (new HermitianFermionTerm() , orbitalIntegral.Coefficient) }; 
            }
            else
            {
                throw new System.NotImplementedException();
            }
        }

        #region Create canonical fermion terms from orbitals
        /// <summary>
        /// Creates all fermion terms generated by all symmetries of a one-body orbital integral.
        /// </summary>
        /// <param name="nOrbitals">Total number of distinct orbitals.</param>
        /// <param name="orbitalIntegral">Input orbital integral.</param>
        /// <param name="indexConvention">Indexing scheme from spin-orbitals to integers.</param>
        private static IEnumerable<(HermitianFermionTerm, double)> ToOneBodySpinOrbitalTerms(
            this OrbitalIntegral orbitalIntegral,
            int nOrbitals,
            IndexConvention indexConvention)
        {
            // One-electron orbital integral symmetries
            // ij = ji
            var pqSpinOrbitals = orbitalIntegral.EnumerateOrbitalSymmetries().EnumerateSpinOrbitals();

            var coefficient = orbitalIntegral.Coefficient;

            foreach (var pq in pqSpinOrbitals)
            {
                var pInt = pq[0].ToInt(indexConvention, nOrbitals);
                var qInt = pq[1].ToInt(indexConvention, nOrbitals);
                if (pInt == qInt)
                {
                    yield return (new HermitianFermionTerm(new[] { pInt, qInt }.ToLadderSequence()), orbitalIntegral.Coefficient);
                }
                else if (pInt < qInt)
                {
                    yield return (new HermitianFermionTerm(new[] { pInt, qInt }.ToLadderSequence()), 2.0 * orbitalIntegral.Coefficient);
                }
            }
        }

        /// <summary>
        /// Updates an instance of <see cref="FermionHamiltonian"/>
        /// with all spin-orbitals from described by a sequence of four-body orbital integrals.
        /// </summary>
        /// <param name="nOrbitals">Total number of distinct orbitals.</param>
        /// <param name="orbitalIntegral">Sequence of four-body orbital integrals.</param>
        /// <param name="indexConvention">Indexing scheme from spin-orbitals to integers.</param>
        private static IEnumerable<(HermitianFermionTerm, double)> ToTwoBodySpinOrbitalTerms(
            this OrbitalIntegral orbitalIntegral,
            int nOrbitals,
            IndexConvention indexConvention)
        {
            List<(HermitianFermionTerm, double)> fermionTerms = new List<(HermitianFermionTerm, double)>();
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

                var pInt = p.ToInt(indexConvention, nOrbitals);
                var qInt = q.ToInt(indexConvention, nOrbitals);
                var rInt = r.ToInt(indexConvention, nOrbitals);
                var sInt = s.ToInt(indexConvention, nOrbitals);

                // Only consider terms on the lower diagonal due to Hermitian symmetry.

                // For terms with two different orbital indices, possibilities are
                // PPQQ (QQ = 0), PQPQ, QPPQ (p<q), PQQP, QPQP (p<q), QQPP (PP=0)
                // Hence, if we only count PQQP, and PQPQ, we need to double the coefficient.
                // iU jU jU iU | iU jD jD iD | iD jU jU iD | iD jD jD iD
                if (pInt == sInt && qInt == rInt && pInt < qInt)
                {   // PQQP
                    yield return (new HermitianFermionTerm(new[] { pInt, qInt, rInt, sInt }.ToLadderSequence()), 1.0 * coefficient );
                }
                else if (pInt == rInt && qInt == sInt && pInt < qInt)
                {
                    // iU jU iU jU | iD jD iD jD
                    // PQPQ
                    yield return (new HermitianFermionTerm(new[] { pInt, qInt, sInt, rInt }.ToLadderSequence()), -1.0 * coefficient );
                }
                else if (qInt == rInt && pInt < sInt && rInt != sInt && pInt != qInt)
                {
                    // PQQR
                    // For any distinct pqr, [i;j;j;k] generates PQQR ~ RQQP ~ QPRQ ~ QRPQ. We only need to record one.
                    if (rInt < sInt)
                    {
                        if (pInt < qInt)
                        {
                            yield return (new HermitianFermionTerm(new[]  { pInt, qInt, sInt, rInt }.ToLadderSequence()), -2.0 * coefficient );
                        }
                        else
                        {
                            yield return (new HermitianFermionTerm(new[]  { qInt, pInt, sInt, rInt }.ToLadderSequence()), 2.0 * coefficient );
                        }

                    }
                    else
                    {
                        if (pInt < qInt)
                        {
                            yield return (new HermitianFermionTerm(new[]  { pInt, qInt, rInt, sInt }.ToLadderSequence()), 2.0 * coefficient );
                        }
                        else
                        {
                            yield return (new HermitianFermionTerm(new[]  { qInt, pInt, rInt, sInt }.ToLadderSequence()), -2.0 * coefficient );
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
                            yield return (new HermitianFermionTerm(new[]  { pInt, qInt, rInt, sInt }.ToLadderSequence()), 2.0 * coefficient );
                        }
                        else
                        {
                            yield return (new HermitianFermionTerm(new[]  { pInt, qInt, sInt, rInt }.ToLadderSequence()), -2.0 * coefficient );
                        }
                    }
                    else
                    {
                        yield return (new HermitianFermionTerm(new[]  { qInt, pInt, rInt, sInt }.ToLadderSequence()), -2.0 * coefficient );
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
                        yield return (new HermitianFermionTerm(new[]  { pInt, qInt, sInt, rInt }.ToLadderSequence()), -2.0 * coefficient );
                    }
                    else
                    {
                        yield return (new HermitianFermionTerm(new[]  { pInt, qInt, rInt, sInt }.ToLadderSequence()), 2.0 * coefficient );
                    }
                }
            }
        }

        #endregion

    }
    

    /// <summary>
    /// Extensions for converting spin orbitals to integers
    /// </summary>
    public static partial class Extensions
    {
        /// <summary>
        /// Converts an array of spin-orbitals into an array of integers representing each spin orbital.
        /// </summary>
        public static IEnumerable<int> ToInts(this IEnumerable<SpinOrbital> spinOrbitals, IndexConvention indexConvention, int nOrbitals) =>
            spinOrbitals.Select(x => x.ToInt(indexConvention, nOrbitals));

        /// <summary>
        /// Converts an array of spin-orbitals into an array of integers representing each spin orbital.
        /// </summary>
        public static IEnumerable<int> ToInts(this IEnumerable<SpinOrbital> spinOrbitals) =>
            spinOrbitals.Select(x => x.ToInt());

        /// <summary>
        /// Converts an array of orbital index and spin index tuples into an array of spin-orbitals.
        /// </summary>
        public static SpinOrbital[] ToSpinOrbitals(this IEnumerable<(int, Spin)> spinOrbitalIndices) => spinOrbitalIndices.Select(o => new SpinOrbital(o)).ToArray();

        /// <summary>
        /// Converts an array of orbital index and spin index tuples into an array of spin-orbitals.
        /// </summary>
        public static SpinOrbital[] ToSpinOrbitals(this IEnumerable<(int, int)> spinOrbitalIndices) =>
            spinOrbitalIndices.Select(o => new SpinOrbital(o)).ToArray();

        /// <summary>
        /// Enumerates all spin-orbitals described by an array of orbital integrals by
        /// applying <see cref="OrbitalIntegral.EnumerateSpinOrbitals"/> to each.
        /// </summary>
        /// <param name="orbitalIntegrals">Array of orbital integrals.</param>
        /// <returns>Array of Array of spin-orbitals.</returns>
        public static SpinOrbital[][] EnumerateSpinOrbitals(this IEnumerable<OrbitalIntegral> orbitalIntegrals) => orbitalIntegrals.SelectMany(o => o.EnumerateSpinOrbitals()).ToArray();
    }
}



