// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using System;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;



/// <summary>
/// This code space provides suggested interoperability classes for running the
/// Q# Qccc (quantum circuit-centric classifier) code on Microsoft quantum simulator
/// </summary>
namespace Microsoft.Quantum.MachineLearning.Interop
{
    using Microsoft.Quantum.Logical;
    using Microsoft.Quantum.MachineLearning;
    using System.Runtime.CompilerServices;
    using System.Runtime.ExceptionServices;
    using System.Runtime.InteropServices.ComTypes;
    using System.Xml;


    /// <summary>
    /// Quick conversions into IQArray format
    /// </summary>
    public class Qonvert
    {

        public static long ToC(char pauli)
        {
            if (pauli.Equals('I'))
            {
                return 0L;
            }
            if (pauli.Equals('X'))
            {
                return 1L;
            }
            if (pauli.Equals('Y'))
            {
                return 2L;
            }
            if (pauli.Equals('Z'))
            {
                return 3L;
            }
            return -1L;
        }

        public static IQArray<IQArray<long>> ToQ(List<long[]> src)
        {
            List<IQArray<long>> tmp = new List<IQArray<long>>(src.Count);
            for (int ix = 0; ix < src.Count; ix++)
            {
                tmp.Add(new QArray<long>(src[ix]));
            }
            return new QArray<IQArray<long>>(tmp.ToArray());
        }

        public static IQArray<long> ToQ(List<long> src)
        {
            return new QArray<long>(src.ToArray());
        }

        public static IQArray<IQArray<double>> ToQ(List<double[]> src)
        {
            List<IQArray<double>> tmp = new List<IQArray<double>>(src.Count);
            for (int ix = 0; ix < src.Count; ix++)
            {
                tmp.Add(new QArray<double>(src[ix]));
            }
            return new QArray<IQArray<double>>(tmp.ToArray());
        }
    } //Qonvert

    public class ClassificationModel
    {
        long _nQubits;
        IQArray<IQArray<long>> _structure;
        IQArray<double> _cachedParameters;
        double _bias;

        public ClassificationModel(long nQubits)
        {
            _nQubits = nQubits;
            _structure = null;
            _cachedParameters = null;
            _bias = -2.0;
        }

        public ClassificationModel(long nQubits, List<long[]> structure)
        {
            _nQubits = nQubits;
            _structure = Qonvert.ToQ(structure);
            _cachedParameters = null;
            _bias = -2.0;
        }

        public ClassificationModel(long nQubits, List<long[]> structure,double[] parameters)
        {
            _nQubits = nQubits;
            _structure = Qonvert.ToQ(structure);
            _cachedParameters = new QArray<double>(parameters);
            _bias = -2.0;
        }

        public ClassificationModel(long nQubits, List<long[]> structure, double[] parameters, double bias)
        {
            _nQubits = nQubits;
            _structure = Qonvert.ToQ(structure);
            _cachedParameters = new QArray<double>(parameters);
            _bias = bias;
        }

        public bool isTrained
        {
            get { return (_bias > -1.5) && (_structure != null) && (_cachedParameters != null); }
        }

        public IQArray<IQArray<long>> CircuitStructure
        {
            get { return _structure; }
            set { _structure = value;  }
        }

        public IQArray<double> CachedParameters
        {
            get { return _cachedParameters; }
        }

        public double Bias
        {
            get { return _bias; }
        }

        /// <summary>
        /// Creates a layer of nQubits Pauli rotations
        /// </summary>
        /// <param name="nQubits">Number of qubits to rotate</param>
        /// <param name="pauli">Type of Pauli gate</param>
        /// <returns>Sequence of nQubits rotation templates</returns>
        public static List<long[]> LocalRotationsLayer(long nQubits, char pauli)
        {
            List<long[]> ret = new List<long[]>((int)nQubits);
            for (long iq = 0; iq < nQubits; iq++)
            {
                long[] localRp = new long[] { -1, Qonvert.ToC(pauli), iq };
                ret.Add(localRp);
            }
            return ret;
        }

        /// <summary>
        /// Creates a layer of nQubits Pauli rotations
        /// </summary>
        /// <param name="nQubits">Number of qubits to rotate</param>
        /// <param name="pauli">Type of Pauli gate</param>
        /// <returns>Sequence of nQubits rotation templates</returns>
        public static List<long[]> PartialLocalLayer(long[] indices, char pauli)
        {
            List<long[]> ret = new List<long[]>(indices.Length);
            foreach (long iq in indices)
            {
                long[] localRp = new long[] { -1, Qonvert.ToC(pauli), iq };
                ret.Add(localRp);
            }
            return ret;
        }

        /// <summary>
        /// Creates a cyclic block of nQubits controlled rotations that starts
        /// with contol qubit (nQubits-1), target qubit (cspan-1) % n , followed by the
        /// ladder of entanglers with control qubit iq and target qubit (iq+cspan) % n
        /// </summary>
        /// <param name="nQubits">Number of qubits to entangle</param>
        /// <param name="pauli"></param>
        /// <param name="cspan">index offset between control and target qubits</param>
        /// <returns></returns>
        public static List<long[]> CyclicEntanglerLayer(long nQubits, char pauli, long cspan)
        {
            List<long[]> ret = new List<long[]>((int)nQubits);
            ret.Add(new long[] { -1, Qonvert.ToC(pauli), (long) ((cspan-1) % nQubits), nQubits - 1 });
            for (long iq = 1; iq < nQubits; iq++)
            {
                long[] entRp = new long[] { -1, Qonvert.ToC(pauli), (long)(((iq + 1) * cspan - 1) % nQubits), (long)((iq*cspan - 1) % nQubits) };
                ret.Add(entRp);
            }
            return ret;
        }

        public static List<long[]> CyclicEntanglerLayer(long nQubits, char pauli)
        {
            return CyclicEntanglerLayer(nQubits, pauli, 1L);
        }

        public static List<long[]> JoinLayers(List<List<long[]>> layers)
        {
            List<long[]> structure = new List<long[]>(layers.Count * layers[0].Count);
            for (int ila = 0; ila < layers.Count; ila++)
            {
                structure.AddRange(layers[ila]);
            }
            return structure;
        }

        public static void reindex(List<long[]> struc)
        {
            for (int ix=0; ix < struc.Count; ix++)
            {
                long[] gt = struc[ix];
                gt[0] = ix;
                struc[ix] = gt;
            }
        }

        public void QcccTrainSequential(IQArray<IQArray<double>> parameterSource, IQArray<IQArray<double>> trainingSet, IQArray<long> trainingLabels, IQArray<IQArray<long>> trainingSchedule,
            IQArray<IQArray<long>> validationSchedule, double learningRate, double tolerance, long miniBatchSize, long maxEpochs, long nMeasurements, uint randomizationSeed)
        {
            var sim = new QuantumSimulator(false, randomizationSeed);
            (this._cachedParameters, this._bias) =
            TrainQcccSequential.Run(sim, this._nQubits, this._structure, parameterSource, trainingSet, trainingLabels, trainingSchedule, validationSchedule, learningRate, tolerance, miniBatchSize, maxEpochs, nMeasurements).Result;
        }
        public void QcccTrainSequential(List<double[]> parameterSource, List<double[]> trainingSet, List<long> trainingLabels, List<long[]> trainingSchedule,
             List<long[]> validationSchedule, double learningRate, double tolerance, long miniBatchSize, long maxEpochs, long nMeasurements, uint randomizationSeed)
        {
            QcccTrainSequential(Qonvert.ToQ(parameterSource), Qonvert.ToQ(trainingSet), Qonvert.ToQ(trainingLabels), Qonvert.ToQ(trainingSchedule),
               Qonvert.ToQ(validationSchedule), learningRate, tolerance, miniBatchSize, maxEpochs, nMeasurements, randomizationSeed);
        }

        public void QcccTrainParallel(IQArray<IQArray<double>> parameterSource, IQArray<IQArray<double>> trainingSet, IQArray<long> trainingLabels, IQArray<IQArray<long>> trainingSchedule,
            IQArray<IQArray<long>> validationSchedule, double learningRate, double tolerance, long miniBatchSize, long maxEpochs, long nMeasurements, uint randomizationSeed)
        {
            var simAll = new List<QuantumSimulator>(parameterSource.Count);
            var resultsAll = new List<(IQArray<double>, double)>(parameterSource.Count);
            var parameterComb = new List<IQArray<IQArray<double>>>(parameterSource.Count);

            var indices = new int[parameterSource.Count];
            for (int j = 0; j < parameterSource.Count; j++)
            {
                indices[j] = j;
                simAll.Add(new QuantumSimulator(false, randomizationSeed));
                resultsAll.Add((new QArray<double>(),0.0));
                parameterComb.Add(new QArray<IQArray<double>>(new IQArray<double>[] { parameterSource[j] })); //Isolating parameter starts - one per thread
            }
            Parallel.ForEach(indices,
                (j) =>
                {

                    var rslt =
                    TrainQcccSequential.Run(simAll[j], this._nQubits, this._structure, parameterComb[j], trainingSet, trainingLabels, trainingSchedule, validationSchedule, learningRate, tolerance, miniBatchSize, maxEpochs, nMeasurements).Result;
                    resultsAll[j] = rslt;
                }
                );
            //Estimated parameters and biases for each proposed parameter start. Now postprocess
            long bestValidation = long.MaxValue;
            int bestJ = -1;
            var sim = new QuantumSimulator(false, randomizationSeed);
            for (int j = 0; j < parameterSource.Count; j++)
            {
                var (pars, bias) = resultsAll[j];
                long misses = CountValidationMisses.Run(sim, tolerance, this._nQubits, trainingSet, trainingLabels, validationSchedule, this._structure, pars, bias, nMeasurements).Result;
                if (bestValidation > misses)
                {
                    bestValidation = misses;
                    bestJ = j;
                }
            }
            (this._cachedParameters, this._bias) = resultsAll[bestJ];
        } //QcccTrainParallel

        public void QcccTrainParallel(List<double[]> parameterSource, List<double[]> trainingSet, List<long> trainingLabels, List<long[]> trainingSchedule,
             List<long[]> validationSchedule, double learningRate, double tolerance, long miniBatchSize, long maxEpochs, long nMeasurements, uint randomizationSeed)
        {
            QcccTrainParallel(Qonvert.ToQ(parameterSource), Qonvert.ToQ(trainingSet), Qonvert.ToQ(trainingLabels), Qonvert.ToQ(trainingSchedule),
               Qonvert.ToQ(validationSchedule), learningRate, tolerance, miniBatchSize, maxEpochs, nMeasurements, randomizationSeed);
        }

        public long CountMisclassifications(double tolerance, IQArray<IQArray<double>> samples, IQArray<long> knownLabels, IQArray<IQArray<long>> validationSchedule, long nMeasurements, uint randomizationSeed)
        {
            if (this.isTrained)
            {
                var sim = new QuantumSimulator(false, randomizationSeed);
                return CountValidationMisses.Run(sim, tolerance, this._nQubits, samples, knownLabels, validationSchedule, this._structure, this.CachedParameters, this.Bias, nMeasurements).Result;
            }
            return long.MaxValue;
        }

        public long CountMisclassifications(double tolerance, List<double[]> samples, List<long> knownLabels, List<long[]> validationSchedule,  long nMeasurements, uint randomizationSeed)
        {
            return CountMisclassifications(tolerance, Qonvert.ToQ(samples), Qonvert.ToQ(knownLabels), Qonvert.ToQ(validationSchedule), nMeasurements, randomizationSeed);
        }

        public long CountMisclassifications(double tolerance, List<double[]> samples, List<long> knownLabels, long nMeasurements, uint randomizationSeed)
        {
            var validationSchedule = new List<long[]>(1);
            validationSchedule.Add(new long[] { 0L, 1L, ((long)(samples.Count - 1)) });
            return CountMisclassifications(tolerance, Qonvert.ToQ(samples), Qonvert.ToQ(knownLabels), Qonvert.ToQ(validationSchedule), nMeasurements, randomizationSeed);
        }

    } //class ClassificationModel

}
