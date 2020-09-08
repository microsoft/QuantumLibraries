// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Runtime.InteropServices;
using System.Runtime.ExceptionServices;

using Microsoft.Quantum.Intrinsic;
using Microsoft.Quantum.Simulation;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using Microsoft.Quantum.Simulation.Simulators.Exceptions;
using Microsoft.Quantum.Standard.Emulation;

namespace Microsoft.Quantum.Characterization
{
    public partial class EstimateFrequencyA
    {
        /// <summary>
        ///  Provides a native emulation of the EstimateFrequency operation for adjointable operations when
        ///  the operation is executed using the full-state QuantumSimulator and the given
        ///  state preparation function does not contain any captured qubits via partial application.
        ///  
        /// The way the emulation works is to invoke the state-preparation only once, and then look 
        /// into the resulting QuantumSimulator's state to get the JointProbability and then
        /// use a classical binomial sampling to get a sample for the resulting probability.
        /// This is typically faster compared to run the state-preparation operation n-times and
        /// calculate the binomial estimation from it.
        /// </summary>
        public class Native : EstimateFrequencyA
        {
            [DllImport(QuantumSimulator.QSIM_DLL_NAME, ExactSpelling = true, CallingConvention = CallingConvention.Cdecl, EntryPoint = "JointEnsembleProbability")]
            private static extern double JointEnsembleProbability(uint id, uint n, Pauli[] b, uint[] q);

            private System.Random _random = new System.Random();

            private QuantumSimulator Simulator { get; }

            protected Allocate Allocate { get; set; }
            protected Release Release { get; set; }
            protected ResetAll ResetAll { get; set; }

            public Native(IOperationFactory m) : base(m)
            {
                this.Simulator = m as QuantumSimulator;
            }

            public override void __Init__()
            {
                base.__Init__();

                this.Allocate = this.__Factory__.Get<Allocate>(typeof(Microsoft.Quantum.Intrinsic.Allocate));
                this.Release = this.__Factory__.Get<Release>(typeof(Microsoft.Quantum.Intrinsic.Release));
                this.ResetAll = this.__Factory__.Get<ResetAll>(typeof(Microsoft.Quantum.Intrinsic.ResetAll));
            }

            /// <summary>
            /// Overrides the body to do the emulation when possible. If emulation is not possible, then
            /// it just invokes the default Q# implementation.
            /// </summary>
            public override Func<(IAdjointable, ICallable, long, long), double> __Body__ => (_args) =>
            {
                var (preparation, measure, count, samples) = _args;

                if (!CanEmulate(preparation, measure))
                {
                    return base.__Body__(_args);
                }

                // Find the basis used for measurement from the captured Paulis:
                var paulis = FindPaulis(measure, count);
                if (paulis.Length != count) throw new InvalidOperationException("The number of paulis must match the number of qubits.");

                var qubits = this.Allocate.Apply(count);
                Exception innerException = null;
                double result = 0.0;
                try
                {
                    preparation.Apply(qubits);
                    var p = 1.0 - JointEnsembleProbability(Simulator.Id, (uint)count, paulis, qubits.GetIds());

                    var random = this.Simulator.Seed == 0 ? new System.Random() : new System.Random((int)this.Simulator.Seed);
                    var dist = new BinomialDistribution(samples, p, random);
                    result = (double)dist.NextSample() / (double)samples;
                    return result;
                }
                // If releasing fails due to not being in the |0⟩ state
                // (as commonly happens when an ExecutionFailException
                // is caught above), we'll get an exception in the finally block
                // below.
                //
                // To prevent that, we need to first catch
                // the ExecutionFailException, since finally blocks are only
                // guaranteed to work for handled exceptions. Next, we'll need
                // to check for an exception caused by releasing qubits that
                // weren't in the |0⟩ state and discard it.
                catch (ExecutionFailException ex)
                {
                    innerException = ex;
                    return result;
                }
                finally
                {
                    try
                    {
                        ResetAll.Apply(qubits);
                        Release.Apply(qubits);
                        // If we got to this finally block by handling an
                        // exception, and didn't hit a second exception
                        // while resetting and releasing, then we can
                        // go on and throw our original exception, being
                        // careful to preserve its stack trace.
                        if (innerException != null)
                        {
                            ExceptionDispatchInfo.Capture(innerException).Throw();
                        }
                    }
                    catch (Exception ex)
                    {
                        // If we were already handling an exception, make a
                        // new aggregate exception and throw it.
                        if (innerException != null)
                        {
                            throw new AggregateException(ex, innerException);
                        }
                        // Otherwise, rethrow the exception that happened
                        // during resetting and releasing, being
                        // careful to preserve its stack trace.
                        else
                        {
                            ExceptionDispatchInfo.Capture(ex).Throw();
                        }
                    }
                }
            };

            /// <summary>
            ///  Helper method to extract the array of Paulis. It requires the measurement operation
            ///  is a Partial application of Primitive.Measure
            /// </summary>
            private static Pauli[] FindPaulis(ICallable measure, long count)
            {
                if (measure.FullName == typeof(MeasureAllZ).FullName)
                {
                    return Enumerable.Repeat<Pauli>(Pauli.PauliZ, (int)count).ToArray();
                }
                else
                {
                    var p = measure as OperationPartial<IQArray<Qubit>, (IQArray<Pauli>, IQArray<Qubit>), Result>;
                    return p.Mapper(null).Item1.ToArray();
                }
            }

            /// <summary>
            /// Determines whether we can do classical emulation for the given preparation and measure operations.
            /// Emulation is only possible if:
            /// 1. If we're running this operation on the full state QuantumSimulator.
            /// 2. The preparation operation has no captured qubits
            /// 3. We're using the Primitive.Measure operation for measurement.
            /// 
            /// If all conditions are met, this method returns true.
            /// </summary>
            public virtual bool CanEmulate(IAdjointable preparation, ICallable measure) =>
                    this.Simulator != null &&
                    (preparation.Qubits == null || !preparation.Qubits.Where(q => q != null).Any()) &&
                    // NB: the "Primitive" deprecation stub has been removed, so we check against the exact
                    //     full name instead. That check should be removed in the future.
                    (measure.FullName == "Microsoft.Quantum.Primitive.Measure" || measure.FullName == typeof(Intrinsic.Measure).FullName || measure.FullName == typeof(MeasureAllZ).FullName);
        }
    }
}
