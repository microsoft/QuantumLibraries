// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;

namespace Microsoft.Quantum.Chemistry
{
    /// <summary>
    /// Generic Hamiltonian class. This is the base class for any Hamiltonians,
    /// which are collections of categorized terms.
    /// </summary>
    /// <typeparam name="TermClassification">Index to categories of terms.</typeparam>
    /// <typeparam name="TermIndexing">Index to individual terms.</typeparam>
    public class Hamiltonian<TermClassification, TermIndexing> 
        //where TermClassification: IEquatable<TermClassification>
        where TermIndexing: HamiltonianTerm<TermClassification>//, IEquatable<TermIndexing>
    {
        /// <summary>
        /// Container for all terms in a Hamiltonian.
        /// </summary>
        public Dictionary<TermClassification, Dictionary<TermIndexing, double>> terms;

        /// <summary>
        /// Indices to systems (e.g. fermions, qubits, or orbitals) the Hamiltonian acts on.
        /// </summary>
        public HashSet<int> systemIndices = new HashSet<int>();

        /// <summary>
        /// Constructor for empty Hamiltonian.
        /// </summary>
        public Hamiltonian()
        {
            terms = new Dictionary<TermClassification, Dictionary<TermIndexing, double>>();
        }

        /// <summary>
        /// Constructor for copying a Hamiltonian.
        /// </summary>
        public Hamiltonian(Hamiltonian<TermClassification, TermIndexing> hamiltonian)
        {
            terms = hamiltonian.terms;
        }

        /// <summary>
        /// Method for adding a term to a Hamiltonian.
        /// </summary>
        /// <param name="type">Category of term.</param>
        /// <param name="index">Index to term.</param>
        /// <param name="coefficient">Coefficient of term.</param>
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
                AddToSystemIndices(index);
            }
        }

        /// <summary>
        /// Add multiple terms to a Hamiltonian.
        /// </summary>
        /// <param name="type">Category of terms.</param>
        /// <param name="terms">Enumerable sequence of terms and coefficients.</param>
        public void AddTerms(TermClassification type, IEnumerable<(TermIndexing, double)> terms)
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
        public void AddTerm(TermIndexing index, double coefficient)
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
        public void AddTerms(IEnumerable<(TermIndexing, double)> terms)
        {
            foreach (var term in terms)
            {
                AddTerm(term.Item1, term.Item2);
            }
        }

        /// <summary>
        /// Method for add all terms from a source Hamiltonian into this Hamiltonian.
        /// </summary>
        /// <param name="sourceHamiltonian">Source Hamiltonian.</param>
        public void AddHamiltonian(Hamiltonian<TermClassification, TermIndexing> sourceHamiltonian)
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
            return terms.Select(o => o.Value.Count()).AsParallel().Sum();
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
            return Math.Pow(terms.Where(o => termTypes.Contains(o.Key)).Select(termType => termType.Value.AsParallel().Select(termValue => Math.Pow(Math.Abs(termValue.Value), power)).Sum()).Sum(),1.0/power);
        }
        
    }
    
}