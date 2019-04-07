// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry.Hamiltonian
{
    /// <summary>
    /// Generic Hamiltonian class. This is the base class for any Hamiltonians,
    /// which are collections of categorized terms.
    /// </summary>
    /// <typeparam name="TermClassification">Index to categories of terms.</typeparam>
    /// <typeparam name="TermIndexing">Index to individual terms.</typeparam>
    public class GenericHamiltonian<TermClassification, TermIndexing, TermValue>
        //TODO: Restore `where TermClassification: IEquatable<TermClassification>`
        // in the future if we want more complicated term classifications.
        where TermIndexing : ITermIndex<TermClassification>
        where TermValue: ITermValue<TermValue>
        // TODO: Restore `IEquatable<TermIndexing>` in the future if we expand to more types of terms.
    {
        /// <summary>
        /// Container for all terms in a Hamiltonian.
        /// </summary>
        public Dictionary<TermClassification, Dictionary<TermIndexing, TermValue>> terms = new Dictionary<TermClassification, Dictionary<TermIndexing, TermValue>>();

        /// <summary>
        /// Indices to systems (e.g. fermions, qubits, or orbitals) the Hamiltonian acts on.
        /// </summary>
        public HashSet<int> systemIndices = new HashSet<int>();

        /// <summary>
        /// Constructor for empty Hamiltonian.
        /// </summary>
        public GenericHamiltonian()
        {
        }

        /// <summary>
        /// Constructor for copying a Hamiltonian.
        /// </summary>
        public GenericHamiltonian(GenericHamiltonian<TermClassification, TermIndexing, TermValue> hamiltonian)
        {
            terms = hamiltonian.terms;
        }

        /// <summary>
        /// Method for adding a term to a Hamiltonian.
        /// </summary>
        /// <param name="type">Category of term.</param>
        /// <param name="index">Index to term.</param>
        /// <param name="coefficient">Coefficient of term.</param>
        public void AddTerm(TermClassification type, TermIndexing index, TermValue coefficient)
        {
            if (!terms.ContainsKey(type))
            {
                terms.Add(type, new Dictionary<TermIndexing, TermValue>());
            }
            if (terms[type].ContainsKey(index))
            {
                terms[type][index] = terms[type][index].AddValue(coefficient);
            }
            else
            {
                terms[type].Add(index, coefficient);
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
            if (!terms.ContainsKey(type))
            {
                return default;
            }
            if (!terms[type].ContainsKey(index))
            {
                return default;
            }
            return terms[type][index];
        }

        /// <summary>
        /// Method for add all terms from a source Hamiltonian into this Hamiltonian.
        /// </summary>
        /// <param name="sourceHamiltonian">Source Hamiltonian.</param>
        public void AddHamiltonian(GenericHamiltonian<TermClassification, TermIndexing, TermValue> sourceHamiltonian)
        {
            foreach(var termType in sourceHamiltonian.terms)
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
            return terms.Select(o => o.Value.Count()).Sum();
        }

        /// <summary>
        /// Computes the L_p norm of coefficicients of all terms in a Hamiltonian.
        /// </summary>
        /// <param name="power">Selects type of norm.</param>
        /// <returns>L_p norm of Hamiltonian coefficients.</returns>
        public double Norm(double power = 1.0)
        {
            return Norm(terms.Keys, power);
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
            var typesEnum = terms.Where(o => termTypes.Contains(o.Key));
            return typesEnum
                .Select(termType => termType.Value
                .Select(termIndex => termIndex.Value.Norm(power)).Norm(power)).Norm(power);
        }
        
    }
    
}