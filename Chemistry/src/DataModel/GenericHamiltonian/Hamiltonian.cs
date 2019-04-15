// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;
using Newtonsoft.Json;
using Microsoft.Quantum.Chemistry;

namespace Microsoft.Quantum.Chemistry.Generic
{
    /// <summary>
    /// Generic Hamiltonian class. This is the base class for any Hamiltonians,
    /// which are collections of categorized terms.
    /// </summary>
    /// <typeparam name="TTermClassification">Index to categories of terms.</typeparam>
    /// <typeparam name="TTermIndexing">Index to individual terms.</typeparam>
    /// <typeparam name="TTermValue">The type of the value of each Term</typeparam>
    public class Hamiltonian<TTermClassification, TTermIndexing, TTermValue>
        //TODO: Restore `where TTermClassification: IEquatable<TTermClassification>`
        // in the future if we want more complicated term classifications.
        where TTermIndexing : ITermIndex<TTermClassification, TTermIndexing>
        where TTermValue: ITermValue<TTermValue>
        // TODO: Restore `IEquatable<TTermIndexing>` in the future if we expand to more types of terms.
    {
        /// <summary>
        /// Represents a single Terms in the Hamiltonian.
        /// </summary>
        [JsonConverter(typeof(Json.HamiltonianTermsJsonConverter))]
        public class HamiltonianTerm : Dictionary<TTermIndexing, TTermValue> { }

        /// <summary>
        /// Represents the collection of all Terms in the Hamiltonian.
        /// </summary>
        [JsonConverter(typeof(Json.HamiltonianTermsJsonConverter))]
        public class HamiltonianTerms : Dictionary<TTermClassification, HamiltonianTerm> { }

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
        public Hamiltonian(Hamiltonian<TTermClassification, TTermIndexing, TTermValue> hamiltonian)
        {
            Terms = hamiltonian.Terms;
        }

        /// <summary>
        /// Adds a term to a Hamiltonian. 
        /// </summary>
        /// <param name="type">Category of term.</param>
        /// <param name="index">Index to term.</param>
        /// <param name="coefficient">Coefficient of term.</param>
        // Warning be careful of case where coefficient is pass by reference.
        public void Add(TTermClassification type, TTermIndexing index, TTermValue coefficient)
        {
            
            var newIndex = index.Clone();
            var newCoeff = coefficient.Clone();
            var sign = newIndex.GetSign();
            newIndex.ResetSign();
            // Some terms have an internal +- sign that multiples the coefficient.
            // We reset this sign once it has been accounted for.
            
            if (!Terms.ContainsKey(type))
            {
                Terms.Add(type, new HamiltonianTerm());
            }
            if (Terms[type].ContainsKey(index))
            {
                // The index can contain a sign coefficient e.g. +-1. 
                Terms[type][newIndex] = Terms[type][newIndex].AddValue(newCoeff, sign);
                //Terms[type][index] = Terms[type][index].AddValue(coefficient, 1);
            }
            else
            {

                Terms[type][newIndex] = newCoeff.SetValue(newCoeff, sign);
                //Terms[type][index] = coefficient;
                //index.ResetSign();
                AddToSystemIndices(index);
            }
        }
        
        /// <summary>
        /// Adds multiple term to a Hamiltonian. 
        /// </summary>
        /// <param name="type">Category of terms.</param>
        /// <param name="terms">Enumerable sequence of terms and coefficients.</param>
        public void AddRange(TTermClassification type, IEnumerable<(TTermIndexing, TTermValue)> terms)
        {
            foreach (var term in terms)
            {
                Add(term.Item1, term.Item2);
            }
        }

        /// <summary>
        /// Adds a term to a Hamiltonian. This method 
        /// infers the term category from the term index if possible.
        /// </summary>
        /// <param name="index">Index to term.</param>
        /// <param name="coefficient">Coefficient of term.</param>
        public void Add(TTermIndexing index, TTermValue coefficient)
        {
            Add(index.GetTermType(), index, coefficient);
        }

        /// <summary>
        /// Add multiple terms to a Hamiltonian. This method 
        /// infers the term category from the term index if possible.
        /// </summary>
        /// <param name="terms">
        /// Enumerable sequence of terms and coefficients.
        /// </param>
        public void AddRange(IEnumerable<(TTermIndexing, TTermValue)> terms)
        {
            foreach (var term in terms)
            {
                Add(term.Item1, term.Item2);
            }
        }

        /// <summary>
        /// Method for retrieving a term to a Hamiltonian. This method 
        /// infers the term category from the term index if possible.
        /// </summary>
        /// <param name="index">Index to term.</param>
        public TTermValue GetTerm(TTermIndexing index)
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
        public void AddHamiltonian(Hamiltonian<TTermClassification, TTermIndexing, TTermValue> sourceHamiltonian)
        {
            foreach(var termType in sourceHamiltonian.Terms)
            {
                AddRange(termType.Key, termType.Value.Select(o => (o.Key, o.Value)));
            }
        }

        /// <summary>
        /// Counts the number of terms in a Hamiltonian.
        /// </summary>
        /// <returns>Number of terms in a Hamiltonian.</returns>
        public int CountTerms() => Terms.Select(o => o.Value.Count()).Sum();

        /// <summary>
        /// Counts the number of systems (fermions) in a Hamiltonian.
        /// </summary>
        /// <returns>Number of systems in a Hamiltonian.</returns>
        public int CountUniqueSystemIndices() =>SystemIndices.Count();

        /// <summary>
        /// Computes the L_p norm of coefficicients of all terms in a Hamiltonian.
        /// </summary>
        /// <param name="power">Selects type of norm.</param>
        /// <returns>L_p norm of Hamiltonian coefficients.</returns>
        public double Norm(double power = 1.0) => Norm(Terms.Keys, power);

        /// <summary>
        /// Method that add system indices to the systemIndices hashset.
        /// </summary>
        /// <param name="index"></param>
        public virtual void AddToSystemIndices(TTermIndexing index) {
        }

        /// <summary>
        /// Computes the L_p norm of coefficicients of categories of terms in a Hamiltonian.
        /// </summary>
        /// <param name="termTypes">Selects the categories of Hamiltonian terms.</param>
        /// <param name="power">Selects type of norm.</param>
        /// <returns>L_p norm of Hamiltonian coefficients.</returns>
        public double Norm(IEnumerable<TTermClassification> termTypes, double power = 1.0)
        {
            var typesEnum = Terms.Where(o => termTypes.Contains(o.Key));
            return typesEnum
                .Select(termType => termType.Value
                    .Select(termIndex => termIndex.Value
                        .Norm(power))
                    .Norm(power))
                .Norm(power);
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