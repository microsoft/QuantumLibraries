// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;
using Newtonsoft.Json;

namespace Microsoft.Quantum.Chemistry.Generic
{
    /// <summary>
    /// Generic Hamiltonian class. This is the base class for any Hamiltonians,
    /// which are collections of categorized terms.
    /// </summary>
    /// <typeparam name="TermClassification">Index to categories of terms.</typeparam>
    /// <typeparam name="TermIndexing">Index to individual terms.</typeparam>
    /// <typeparam name="TermValue">The type of the value of each Term</typeparam>
    public partial class Hamiltonian<TermClassification, TermIndexing, TermValue>
        //TODO: Restore `where TermClassification: IEquatable<TermClassification>`
        // in the future if we want more complicated term classifications.
        where TermIndexing : ITermIndex<TermClassification>
        where TermValue : ITermValue<TermValue>
        // TODO: Restore `IEquatable<TermIndexing>` in the future if we expand to more types of terms.
    {
        /// <summary>
        /// Represents a single Terms in the Hamiltonian.
        /// </summary>
        [JsonConverter(typeof(HamiltonianTermsJsonConverter))]
        public class HamiltonianTerm : Dictionary<TermIndexing, TermValue> { }

        /// <summary>
        /// Represents the collection of all Terms in the Hamiltonian.
        /// </summary>
        [JsonConverter(typeof(HamiltonianTermsJsonConverter))]
        public class HamiltonianTerms : Dictionary<TermClassification, HamiltonianTerm> { }

        /// <summary>
        /// Container for all terms in a Hamiltonian.
        /// </summary>
        public HamiltonianTerms Terms = new HamiltonianTerms();

        /// <summary>
        /// Indices to systems (e.g. fermions, qubits, or orbitals) the Hamiltonian acts on.
        /// </summary>
        public HashSet<int> SystemIndices = new HashSet<int>();

        /// <summary>
        /// Constructor for empty Hamiltonian.
        /// </summary>
        public Hamiltonian()
        {
        }

        /// <summary>
        /// Constructor for copying a Hamiltonian.
        /// </summary>
        public Hamiltonian(Hamiltonian<TermClassification, TermIndexing, TermValue> hamiltonian)
        {
            Terms = hamiltonian.Terms;
        }

        /// <summary>
        /// Method for adding a term to a Hamiltonian.
        /// </summary>
        /// <param name="type">Category of term.</param>
        /// <param name="index">Index to term.</param>
        /// <param name="coefficient">Coefficient of term.</param>
        public void AddTerm(TermClassification type, TermIndexing index, TermValue coefficient)
        {
            if (!Terms.ContainsKey(type))
            {
                Terms.Add(type, new HamiltonianTerm());
            }
            if (Terms[type].ContainsKey(index))
            {
                Terms[type][index] = Terms[type][index].AddValue(coefficient);
            }
            else
            {
                Terms[type].Add(index, coefficient);
                AddToSystemIndices(index);
            }
        }
        

        /// <summary>
        /// Add multiple terms to a Hamiltonian.
        /// </summary>
        /// <param name="type">Category of terms.</param>
        /// <param name="terms">Enumerable sequence of terms and coefficients.</param>
        public void AddTerms(TermClassification type, IEnumerable<(TermIndexing, TermValue)> terms)
        {
            foreach (var term in terms)
            {
                AddTerm(term.Item1, term.Item2);
            }
        }

        /// <summary>
        /// Method for adding a term to a Hamiltonian. This method 
        /// infers the term category from the term index if possible.
        /// </summary>
        /// <param name="index">Index to term.</param>
        /// <param name="coefficient">Coefficient of term.</param>
        public void AddTerm(TermIndexing index, TermValue coefficient)
        {
            AddTerm(index.GetTermType(), index, coefficient);
        }

        /// <summary>
        /// Add multiple terms to a Hamiltonian. This method 
        /// infers the term category from the term index if possible.
        /// </summary>
        /// <param name="terms">
        /// Enumerable sequence of terms and coefficients.
        /// </param>
        public void AddTerms(IEnumerable<(TermIndexing, TermValue)> terms)
        {
            foreach (var term in terms)
            {
                AddTerm(term.Item1, term.Item2);
            }
        }

        /// <summary>
        /// Method for retrieving a term to a Hamiltonian. This method 
        /// infers the term category from the term index if possible.
        /// </summary>
        /// <param name="index">Index to term.</param>
        public TermValue GetTerm(TermIndexing index)
        {
            var type = index.GetTermType();
            if (!Terms.ContainsKey(type))
            {
                return default;
            }
            if (!Terms[type].ContainsKey(index))
            {
                return default;
            }
            return Terms[type][index];
        }

        /// <summary>
        /// Method for add all terms from a source Hamiltonian into this Hamiltonian.
        /// </summary>
        /// <param name="sourceHamiltonian">Source Hamiltonian.</param>
        public void AddHamiltonian(Hamiltonian<TermClassification, TermIndexing, TermValue> sourceHamiltonian)
        {
            foreach(var termType in sourceHamiltonian.Terms)
            {
                AddTerms(termType.Key, termType.Value.Select(o => (o.Key, o.Value)));
            }
        }

        /// <summary>
        /// Counts the number of terms in a Hamiltonian.
        /// </summary>
        /// <returns>Number of terms in a Hamiltonian.</returns>
        public int CountTerms()
        {
            return Terms.Select(o => o.Value.Count()).Sum();
        }

        /// <summary>
        /// Counts the number of systems (fermions) in a Hamiltonian.
        /// </summary>
        /// <returns>Number of systems in a Hamiltonian.</returns>
        public int CountUniqueSystemIndices()
        {
            return SystemIndices.Count();
        }

        /// <summary>
        /// Computes the L_p norm of coefficicients of all terms in a Hamiltonian.
        /// </summary>
        /// <param name="power">Selects type of norm.</param>
        /// <returns>L_p norm of Hamiltonian coefficients.</returns>
        public double Norm(double power = 1.0)
        {
            return Norm(Terms.Keys, power);
        }

        /// <summary>
        /// Method that add system indices to the systemIndices hashset.
        /// </summary>
        /// <param name="index"></param>
        public virtual void AddToSystemIndices(TermIndexing index) {
        }

        /// <summary>
        /// Computes the L_p norm of coefficicients of categories of terms in a Hamiltonian.
        /// </summary>
        /// <param name="termTypes">Selects the categories of Hamiltonian terms.</param>
        /// <param name="power">Selects type of norm.</param>
        /// <returns>L_p norm of Hamiltonian coefficients.</returns>
        public double Norm(IEnumerable<TermClassification> termTypes, double power = 1.0)
        {
            var typesEnum = Terms.Where(o => termTypes.Contains(o.Key));
            return typesEnum
                .Select(termType => termType.Value
                .Select(termIndex => termIndex.Value.Norm(power)).Norm(power)).Norm(power);
        }
        
        /// <summary>
        /// String representation of Hamiltonian.
        /// </summary>
        /// <returns>String representation of Hamiltonian.</returns>
        public override string ToString()
        {
            var output = "";
            foreach (var termType in Terms)
            {
                output += $"{termType.Key} has {termType.Value.Count()} entries).\n";
                foreach (var term in termType.Value)
                {
                    output += $"    {term.ToString()}\n";
                }
            }
            return output;
        }
    }

}