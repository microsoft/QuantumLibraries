// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

namespace Microsoft.Quantum.Qaoa.QaoaHybrid

{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Accord.Math.Optimization;
    using Microsoft.Quantum.Qaoa;
    using Microsoft.Quantum.Simulation.Core;
    using Microsoft.Quantum.Simulation.Simulators;
    using Microsoft.Quantum.Qaoa.QaoaHybrid.Helpers;

    /// <summary>
    /// This class runs a hybrid (quantum-classical) QAOA given an instance of a combinatorial optimization problem encoded into a 2-local Hamiltonian. The classical part is used for optimizing QAOA input parameters and is implemented using the Cobyla optimizer. QAOA input parameters can be optionally specified by a user and they will be treated as a starting point for optimization. Otherwise, input parameters are initialized randomly (possibly many times, as specified by the numberOfRandomStartingPoints variable).
    /// </summary>
    public class HybridQaoa
    {
        private readonly int numberOfIterations;
        private readonly int p;
        private readonly ProblemInstance problemInstance;
        private readonly bool shouldLog;
        private QaoaSolution solution;
        private QaoaLogger? logger;

        /// <summary>
        /// Initializes a new instance of the <see cref="HybridQaoa"/> class.
        /// </summary>
        /// <param name="numberOfIterations">
        /// The number of iterations for sampling the average value of the objective function Hamiltonian.
        /// </param>
        /// <param name="p">
        /// A parameter related to the depth of a QAOA circuit. It corresponds to the number of times evolution operators are applied.
        /// </param>
        /// <param name="problemInstance">
        /// An object that describes a combinatorial optimization problem to be solved by the algorithm.
        /// </param>
        /// <param name="shouldLog">
        /// A flag that specifies whether a log from the hybrid QAOA should be saved in a text file.
        /// </param>
        public HybridQaoa(int numberOfIterations, int p, ProblemInstance problemInstance, bool shouldLog = false)
        {
            this.numberOfIterations = numberOfIterations;
            this.p = p;
            this.problemInstance = problemInstance;
            this.shouldLog = shouldLog;
        }

        /// <summary>
        /// Calculates the value of the cost function based on costs provided.
        /// </summary>
        /// <param name="result">
        /// A binary string. In this context it is a result that we get after measuring the QAOA state.
        /// </param>
        /// <param name="costs">
        /// A list of costs for the cost function.
        /// </param>
        /// <returns>
        /// The value of the costfunction.
        /// </returns>
        public double EvaluateCostFunction(bool[] result, double[] costs)
        {
            double costFunctionValue = 0;
            for (var i = 0; i < this.problemInstance.ProblemSizeInBits; i++)
            {
                costFunctionValue += costs[i] * Convert.ToInt32(result[i]);
            }

            return costFunctionValue;
        }

        /// <summary>
        /// Uses a quantum function to get a solution string from the QAOA that relies on the current values of beta and gamma vectors.
        /// To get a reasonable estimate for the expectation value of a Hamiltonian that encodes the problem, we run the QAOA many times and calculate the expectation based on solutions obtained.
        /// If the expectation of the Hamiltonian is smaller than our current best, we update our best solution to the current solution. The solution vector for the current best solution is the mode of boolean strings that we obtained from the QAOA.
        /// </summary>
        /// <param name="concatenatedQaoaParameters">
        /// Betas and gamma vectors concatenated.
        /// </param>
        /// <returns>
        /// The expected value of a Hamiltonian that we calculated in this run.
        /// </returns>
        private double CalculateObjectiveFunction(double[] concatenatedQaoaParameters)
        {
            var qaoaParameters = new QaoaParameters(concatenatedQaoaParameters);
            var hamiltonianExpectationValue = 0.0;
            var allSolutionVectors = new List<bool[]>();

            var betas = new QArray<double>(qaoaParameters.Betas);
            var gammas = new QArray<double>(qaoaParameters.Gammas);

            var oneLocalHamiltonianCoefficients = new QArray<double>(this.problemInstance.OneLocalHamiltonianCoefficients);
            var twoLocalHamiltonianCoefficients = new QArray<double>(this.problemInstance.TwoLocalHamiltonianCoefficients);

            using (var qsim = new QuantumSimulator())
            {

                for (var i = 0; i < this.numberOfIterations; i++)
                {
                    var result = RunQaoa.Run(qsim, this.problemInstance.ProblemSizeInBits, betas, gammas, oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients, this.p).Result;
                    allSolutionVectors.Add(result.ToArray());
                    var hamiltonianValue = this.problemInstance.EvaluateHamiltonian(result.ToArray());
                    hamiltonianExpectationValue += hamiltonianValue / this.numberOfIterations;

                }
            }

            this.UpdateBestSolution(hamiltonianExpectationValue, allSolutionVectors, qaoaParameters);

            if (this.shouldLog)
            {
                this.logger.LogCurrentBestSolution(betas, gammas, this.solution.SolutionHamiltonianValue, this.solution.SolutionVector);
            }

            return hamiltonianExpectationValue;
        }

        /// <summary>
        /// Updates the currently best solution if a new solution is better.
        /// </summary>
        /// <param name="hamiltonianExpectationValue">
        /// Expectation value of a Hamiltonian.
        /// </param>
        /// <param name="allSolutionVectors">
        /// A vector of all binary solutions that were found by a QAOA.
        /// </param>
        /// <param name="qaoaParameters">
        /// Betas and gamma coefficients.
        /// </param>
        private void UpdateBestSolution(double hamiltonianExpectationValue, List<bool[]> allSolutionVectors, QaoaParameters qaoaParameters)
        {
            if (hamiltonianExpectationValue < this.solution.SolutionHamiltonianValue)
            {
                var mostProbableSolutionVectorTemp = ModeFinder.FindModeInBoolList(allSolutionVectors);
                this.solution.SolutionHamiltonianValue = hamiltonianExpectationValue;
                this.solution.SolutionVector = mostProbableSolutionVectorTemp;
                this.solution.SolutionQaoaParameters = qaoaParameters;
            }
        }

        /// <summary>
        /// Generates constraints for elements in beta and gamma vectors.
        /// </summary>
        /// <returns>
        /// Generated constraints.
        /// </returns>
        /// <remarks>
        /// For the canonical choice of the mixing Hamiltonian (i.e. the sum of X operators acting on single qubits), the range of values in the beta vector is 0 <= beta_i <= PI.
        /// For the objective function Hamiltonian based on Z operators, the range of values in the gamma vector is 0 <= beta_i <= 2PI.
        /// </remarks>
        private NonlinearConstraint[] GenerateConstraints()
        {
            var constraints = new NonlinearConstraint[4*this.p];
            foreach (var i in Enumerable.Range(0, this.p).Select(x => x * 2))
            {
                int gammaIndex = (2 * this.p) + i;
                constraints[i] = new NonlinearConstraint(2 * this.p, x => x[i / 2] >= 0);
                constraints[i + 1] = new NonlinearConstraint(2 * this.p, x => x[i / 2] <= Math.PI);
                constraints[gammaIndex] = new NonlinearConstraint(2 * this.p, x => x[gammaIndex / 2] >= 0);
                constraints[gammaIndex + 1] = new NonlinearConstraint(2 * this.p, x => x[gammaIndex / 2] <= 2 * Math.PI);
            }

            return constraints;
        }

        /// <summary>
        /// Uses a classical optimizer to change beta and gamma parameters so that the objective function is minimized. The optimization is performed some number of times to decrease the chance of getting stuck in a local minimum.
        /// </summary>
        /// <param name="numberOfRandomStartingPoints">
        /// A number of times the hybrid QAOA will be ran with randomly initiated values of beta and gamma. The bigger the number, the lower the chance of getting a solution that is a local minimum.
        /// </param>
        /// <returns>
        /// Optimal solution to the optimization problem input by the user.
        /// </returns>
        /// <remarks>
        /// Currently used optimizer is Cobyla which is a gradient-free optimization technique.
        /// The objective function Hamiltonian is based on Z operators, the range of values in the beta vector is 0 <= beta_i <= PI, the range of values in the gamma vector is 0 <= beta_i <= 2PI.
        /// </remarks>
        public QaoaSolution RunOptimization(int numberOfRandomStartingPoints)
        {
            if (this.shouldLog)
            {
                this.logger = new QaoaLogger();
            }

            this.solution = new QaoaSolution(null, double.MaxValue, null);

            Func<double[], double> objectiveFunction = this.CalculateObjectiveFunction;

            var optimizerObjectiveFunction = new NonlinearObjectiveFunction(2 * this.p, objectiveFunction);
            var constraints = this.GenerateConstraints();
            var cobyla = new Cobyla(optimizerObjectiveFunction, constraints);

            for (int i = 0; i < numberOfRandomStartingPoints; i++)
            {
                var concatenatedQaoaParameters = new QaoaParameters(p).getConcatenatedQaoaParameters();
                var success = cobyla.Minimize(concatenatedQaoaParameters);

                this.logger?.LogSuccess(success);
            }

            if (this.shouldLog)
            {
                this.logger.Dispose();
            }

            return this.solution;
        }

        /// <summary>
        /// Uses a classical optimizer to change beta and gamma parameters so that the objective function is minimized. Initial beta and gamma parameters are provided by a user.
        /// </summary>
        /// <param name="qaoaParameters">
        /// User-defined initial values of beta and gamma parameters.
        /// <returns>
        /// Optimal solution to the optimization problem input by the user.
        /// </returns>
        /// <remarks>
        /// Currently used optimizer is Cobyla which is a gradient-free optimization technique.
        /// The objective function Hamiltonian is based on Z operators, the range of values in the beta vector is 0 <= beta_i <= PI, the range of values in the gamma vector is 0 <= beta_i <= 2PI.
        /// </remarks>
        public QaoaSolution RunOptimization(QaoaParameters qaoaParameters)
        {
            if (this.shouldLog)
            {
                this.logger = new QaoaLogger();
            }

            this.solution = new QaoaSolution(null, double.MaxValue, null);

            Func<double[], double> objectiveFunction = this.CalculateObjectiveFunction;
            var optimizerObjectiveFunction = new NonlinearObjectiveFunction(2 * this.p, objectiveFunction);
            var constraints = this.GenerateConstraints();

            var cobyla = new Cobyla(optimizerObjectiveFunction, constraints);
            var concatenatedQaoaParameters = qaoaParameters.getConcatenatedQaoaParameters();
            var success = cobyla.Minimize(concatenatedQaoaParameters);

            if (this.shouldLog)
            {
                this.logger.LogSuccess(success);
                this.logger.Dispose();
            }

            return this.solution;
        }
    }
}
