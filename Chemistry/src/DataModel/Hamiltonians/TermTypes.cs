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

    public class Hamiltonian<TermClassification, TermIndexing>
        //where TermClassification: IEquatable<TermClassification>
        where TermIndexing: IEquatable<TermIndexing>, HamiltonianTerm<TermClassification>
    {
        /// <summary>
        /// Container for all terms in a Hamiltonian.
        /// </summary>
        Dictionary<TermClassification, Dictionary<TermIndexing, double>> terms;

        public Hamiltonian()
        {
            terms = new Dictionary<TermClassification, Dictionary<TermIndexing, double>>();
        }

        public Hamiltonian(Hamiltonian<TermClassification, TermIndexing> hamiltonian)
        {
            terms = hamiltonian.terms;
        }

        public void AddTerm(TermClassification type, TermIndexing index, double coefficient)
        {
            if (!terms.ContainsKey(type))
            {
                terms.Add(type, new Dictionary<TermIndexing, double>());
            }
            if (terms[type].ContainsKey(index))
            {
                terms[type][index] = terms[type][index] + coefficient;
            }
            else
            {
                terms[type].Add(index, coefficient);
            }
        }

        public void AddTerm(TermIndexing index, double coefficient)
        {
            AddTerm(index.GetTermType(), index, coefficient);
        }

        public void AddTerms(IEnumerable<(TermIndexing, double)> terms)
        {
            foreach (var term in terms)
            {
                AddTerm(term.Item1, term.Item2);
            }
        }

        public int CountTerms()
        {
            return terms.Select(o => o.Value.Count()).Sum();
        }

        public double Norm(double power = 1.0)
        {
            return Norm(terms.Keys, power);
        }
        public double Norm(IEnumerable<TermClassification> termTypes, double power = 1.0)
        {
            return Math.Pow(terms.Where(o => termTypes.Contains(o.Key)).Select(termType => termType.Value.Select(termValue => Math.Pow(Math.Abs(termValue.Value), power)).Sum()).Sum(),1.0/power);
        }
    }


    public interface HamiltonianTerm <Classification>
    {
        Classification GetTermType();
    }

    public static class TermType
    {
        public enum OrbitalIntegral
        {
            Identity, OneBody, TwoBody
        }

        public enum Fermion
        {
            Identity, PP, PQ, PQQP, PQQR, PQRS
        }


        /// <summary>
        /// Represents the type of a term in a Fermion Hamiltonian.
        /// </summary>
        [Serializable]
        public struct FermionTermType
        {

            /// <summary>
            /// First parameter <c>Int64</c> is the number of different spin-orbitals.
            /// Second parameter <c>QArray<Int64></c> represents a sequence of creation and annihilation operators.
            /// Third parameter <c>string</c> annotates the term type.
            /// </summary>
            public readonly (Int64, Int64[]) type;

            /// <summary>
            /// Constructor of <c>FermionTermTyp</c>
            /// </summary>
            /// <param name="uniqueTerms">Number of distint indices.</param>
            /// <param name="creationAnnihilation">Sequence of creation and annihilation operators.</param>
            /// <param name="annotation">Text annotation of term type.</param>
            public FermionTermType(Int64 uniqueTerms, Int64[] creationAnnihilation)
            {
                type = (uniqueTerms, creationAnnihilation);
            }

            /// <summary>
            /// Returns the sequence of creation and annihilation operators.
            /// </summary>
            /// <returns>Sequence of creation and annihilation operators of <see cref="FermionTermType></returns>
            public Int64[] GetConjugateSequence()
            {
                return type.Item2;
            }

            /// <summary>
            /// Commonly used term types in a Fermion Hamiltonian.
            /// </summary>
            public sealed class Common
            {
                #region Common term types
                public readonly static FermionTermType IdentityTermType = new FermionTermType(0, new Int64[] { });
                public readonly static FermionTermType PPTermType = new FermionTermType(1, new Int64[] { 1, 0 });
                public readonly static FermionTermType PQTermType = new FermionTermType(2, new Int64[] { 1, 0 });
                public readonly static FermionTermType PQQPTermType = new FermionTermType(2, new Int64[] { 1, 1, 0, 0 });
                public readonly static FermionTermType PQQRTermType = new FermionTermType(3, new Int64[] { 1, 1, 0, 0 });
                public readonly static FermionTermType PQRSTermType = new FermionTermType(4, new Int64[] { 1, 1, 0, 0 });
                #endregion

            }

            /// <summary>
            /// Checks whether a <c>FermionTermType</c> is in canonical order. This means
            /// that <c>CreationAnnihilation</c> is a sequence of all ones followed by a sequence of all zeroes.
            /// </summary>
            /// <returns>Returns <c>true</c> if in canonical oder, and <c>false</c> otherwise.</returns>
            public bool IsInCanonicalOrder()
            {
                if (type.Item1 > type.Item2.Length)
                {
                    return false;
                }
                for (int i = 0; i < type.Item2.Length - 1; i++)
                {
                    if (type.Item2[i] < type.Item2[i + 1])
                    {
                        return false;
                    }
                }
                return true;
            }


            /// <summary>
            /// String representation of <see cref="FermionTermType"/>.
            /// </summary>
            /// <returns>String representation of <see cref="FermionTermType"/>.</returns>
            public override string ToString()
            {
                return $"({type.Item1}, [{String.Join(",", type.Item2)}])";
            }

            /// <summary>
            /// Boolean tests.
            /// </summary>
            public static bool operator ==(FermionTermType x, FermionTermType y)
            {
                return x.type.Item1 == y.type.Item1 && x.type.Item2.SequenceEqual(y.type.Item2);
            }
            public static bool operator !=(FermionTermType x, FermionTermType y)
            {
                return !(x == y);
            }

            public override bool Equals(object x)
            {
                return (x is FermionTermType ft) ? Equals(ft) : false;
            }

            public bool Equals(FermionTermType x)
            {
                if (ReferenceEquals(null, x))
                {
                    return false;
                }
                else if (ReferenceEquals(this, x))
                {
                    return true;
                }
                else if (GetType() != x.GetType())
                {
                    return false;
                }
                else
                    return this == x;
            }

            public override int GetHashCode()
            {
                int h = 19;
                foreach (var i in type.Item2)
                {
                    h = h * 53 + (i * type.Item1).GetHashCode();
                }
                return h;
            }
        }
    }
}