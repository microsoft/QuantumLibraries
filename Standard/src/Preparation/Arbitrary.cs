// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using System.Runtime.InteropServices;
using Microsoft.Quantum.Simulation;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;


namespace Microsoft.Quantum.Preparation
{
    public partial class ApproximatelyPrepareArbitraryState
    {
        /// <summary>
        ///  Provides a native emulation of the ApproximatelyPrepareArbitraryState operation when
        ///  the operation is executed using the full-state QuantumSimulator.
        /// </summary>
        public class Native : ApproximatelyPrepareArbitraryState
        {
            [DllImport(QuantumSimulator.QSIM_DLL_NAME, ExactSpelling = true,
                CallingConvention = CallingConvention.Cdecl, EntryPoint = "InjectState")]
            private static extern double InjectState(uint sid, uint n, uint[] q, double[] re, double[] im);

            private QuantumSimulator? Simulator { get; }

            public Native(IOperationFactory m) : base(m)
            {
                this.Simulator = m as QuantumSimulator;
            }

            /// <summary>
            /// Overrides the body to do the emulation when possible. If emulation is not possible, then
            /// it just invokes the default Q# implementation.
            /// </summary>
            public override Func<(double, IQArray<Math.ComplexPolar>, Arithmetic.LittleEndian), QVoid>__Body__ => (_args) =>
            {
                var (tolerance, polarAmplitudes, qubits) = _args;

                // TODO: benchmark for small `qubits` arrays to find out in which cases emulation is actually
                // benefitial.
                if (this.Simulator == null)
                {
                    return base.__Body__(_args);
                }

                // calculate the norm as we might need to normalize the state
                var norm = 0.0;
                foreach (var pa in polarAmplitudes) { norm += pa.Magnitude * pa.Magnitude; }
                norm = System.Math.Sqrt(norm);

                var state_size = polarAmplitudes.Length;
                var re = new double[state_size];
                var im = new double[state_size];
                for (int i = 0; i < state_size; i++)
                {
                    var pa = polarAmplitudes[i];
                    re[i] = (System.Math.Abs(pa.Magnitude) * System.Math.Cos(pa.Argument))/norm;
                    im[i] = (System.Math.Abs(pa.Magnitude) * System.Math.Sin(pa.Argument))/norm;
                }

                // Replace implementation in the library to stop throwing and return success/failure result instead.
                try
                {
                    // State injection is precise so tolerance isn't used.
                    InjectState(Simulator.Id, (uint)qubits.Data.Length, qubits.Data.GetIds(), re, im);
                }
                catch
                {
                    return base.__Body__(_args);
                }
                return QVoid.Instance;
            };
        }
    }
}
