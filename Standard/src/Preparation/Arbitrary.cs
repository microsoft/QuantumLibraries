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
    public partial class _PrepareAmplitudesFromZeroState
    {
        /// <summary>
        ///  Provides a native emulation of the ApproximatelyPrepareArbitraryState operation when
        ///  the operation is executed using the full-state QuantumSimulator.
        /// </summary>
        public class Native : _PrepareAmplitudesFromZeroState
        {
            [DllImport(QuantumSimulator.QSIM_DLL_NAME, ExactSpelling = true,
                CallingConvention = CallingConvention.Cdecl, EntryPoint = "InjectState")]
            private static extern bool InjectState(uint sid, uint n, uint[] q, double[] re, double[] im);

            private QuantumSimulator? Simulator { get; }

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
                var (polarAmplitudes, qubits) = _args;

                // TODO: benchmark for small `qubits` arrays to find out in which cases emulation is actually
                // benefitial.
                if (this.Simulator == null)
                {
                    return base.__Body__(_args);
                }

                // Calculate the norm as we might need to normalize the requsted state.
                var norm = 0.0;
                foreach (var pa in polarAmplitudes) { norm += pa.Magnitude * pa.Magnitude; }
                norm = System.Math.Sqrt(norm);

                // Setup the amplitudes arrays for the call to native (it needs to translate from polar to cartesian and
                // might need to pad the tail of an incomplete amplitudes' array with zeros).
                var stateSize = (long)1 << (int)qubits.Data.Length;
                var re = new double[stateSize];
                var im = new double[stateSize];
                for (long i = 0; i < polarAmplitudes.Length; i++)
                {
                    var pa = polarAmplitudes[i];
                    re[i] = (System.Math.Abs(pa.Magnitude) * System.Math.Cos(pa.Argument))/norm;
                    im[i] = (System.Math.Abs(pa.Magnitude) * System.Math.Sin(pa.Argument))/norm;
                }
                for (long i = polarAmplitudes.Length; i < stateSize; i++)
                {
                    re[i] = 0.0;
                    im[i] = 0.0;
                }

                // Emulation might fail if the target qubits are entangled or not all in state |0>. In this case
                // we should fallback to the quantum state preparation as it guarantees the operation to be a proper
                // unitary no matter the state of the qubits.
                if (!InjectState(Simulator.Id, (uint)qubits.Data.Length, qubits.Data.GetIds(), re, im))
                {
                    return base.__Body__(_args);
                }
                return QVoid.Instance;
            };
        }
    }
}
