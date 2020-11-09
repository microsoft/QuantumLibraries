// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using System;
using System.Collections;
using System.Runtime.InteropServices;
using Microsoft.Quantum.Simulation;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;


namespace Microsoft.Quantum.Preparation
{
    public partial class PrepareArbitraryState
    {
        /// <summary>
        ///  Provides a native emulation of the PrepareArbitraryState operation when
        ///  the operation is executed using the full-state QuantumSimulator.
        /// </summary>
        public class Native : PrepareArbitraryState
        {
            [DllImport(QuantumSimulator.QSIM_DLL_NAME, ExactSpelling = true,
                CallingConvention = CallingConvention.Cdecl, EntryPoint = "InjectState")]
            private static extern double InjectState(uint sid, uint n, uint[] q, double[] re, double[] im);

            private QuantumSimulator Simulator { get; }

            public Native(IOperationFactory m) : base(m)
            {
                this.Simulator = m as QuantumSimulator;
            }

            /// <summary>
            /// Overrides the body to do the emulation when possible. If emulation is not possible, then
            /// it just invokes the default Q# implementation.
            /// </summary>
            public override Func<(IQArray<Math.ComplexPolar>, Arithmetic.LittleEndian), QVoid>__Body__ => (_args) =>
            {
                var (polar_amplitudes, qubits) = _args;

                if (this.Simulator == null)
                {
                    return base.__Body__(_args);
                }

                // calculate the norm as we might need to normalize the state
                var norm = 0.0;
                foreach (var pa in polar_amplitudes) { norm += pa.Magnitude * pa.Magnitude; }
                norm = System.Math.Sqrt(norm);

                var state_size = polar_amplitudes.Length;
                var re = new double[state_size];
                var im = new double[state_size];
                for (int i = 0; i < state_size; i++)
                {
                    var pa = polar_amplitudes[i];
                    re[i] = (System.Math.Abs(pa.Magnitude) * System.Math.Cos(pa.Argument))/norm;
                    im[i] = (System.Math.Abs(pa.Magnitude) * System.Math.Sin(pa.Argument))/norm;
                }

                InjectState(Simulator.Id, (uint)qubits.Data.Length, qubits.Data.GetIds(), re, im);
                return QVoid.Instance;
            };
        }
    }
}
