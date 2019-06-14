using System;
using System.Linq;
using System.Runtime.InteropServices;

using Microsoft.Quantum.Intrinsic;
using Microsoft.Quantum.Simulation;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
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

            public Native(IOperationFactory m) : base(m)
            {
                this.Simulator = m as QuantumSimulator;
            }

            public override void Init()
            {
                base.Init();

                this.Allocate = this.Factory.Get<Allocate>(typeof(Microsoft.Quantum.Intrinsic.Allocate));
                this.Release = this.Factory.Get<Release>(typeof(Microsoft.Quantum.Intrinsic.Release));
            }

            /// <summary>
            /// Overrides the body to do the emulation when possible. If emulation is not possible, then
            /// it just invokes the default Q# implementation.
            /// </summary>
            public override Func<(IAdjointable, ICallable, long, long), double> Body => (_args) =>
            {
                var (preparation, measure, count, samples) = _args;

                if (!CanEmulate(preparation, measure))
                {
                    return base.Body(_args);
                }

                // Find the basis used for measurement from the captured Paulis:
                var paulis = FindPaulis(measure, count);
                if (paulis.Length != count) throw new InvalidOperationException("The number of paulis must match the number of qubits.");

                var qubits = this.Allocate.Apply(count);
                try
                {
                    preparation.Apply(qubits);
                    var p = 1.0 - JointEnsembleProbability(Simulator.Id, (uint)count, paulis, qubits.GetIds());
                    preparation.Adjoint.Apply(qubits);

                    var dist = new BinomialDistribution(samples, p);
                    return (double)dist.NextSample() / (double)samples;
                }
                finally
                {
                    Release.Apply(qubits);
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
                    (measure.FullName == typeof(Primitive.Measure).FullName || measure.FullName == typeof(Intrinsic.Measure).FullName || measure.FullName == typeof(MeasureAllZ).FullName);
        }
    }
}
