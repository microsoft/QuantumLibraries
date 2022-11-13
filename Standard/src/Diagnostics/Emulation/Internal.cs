// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using System.Linq;
using System.Runtime.InteropServices;
using Microsoft.Quantum.Diagnostics.Emulation;
using Microsoft.Quantum.Simulation.Common;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;

namespace Microsoft.Quantum.Diagnostics
{
    internal class ArrayDumper : QuantumSimulator.StateDumper
    {
        internal ComplexSquareMatrix? Data = null;
        private double scaleFactor = 1.0;
        public ArrayDumper(QuantumSimulator sim) : base(sim)
        {
        }

        public override bool Callback([MarshalAs(UnmanagedType.LPStr)] string idx, double real, double img)
        {
            if (Data == null) throw new Exception("Expected data buffer to be initialized before callback, but it was null.");
            var index = (int)(CommonNativeSimulator.DisplayableState.BasisStateLabelToBigInt(idx));
            var row = index / Data.Dimension;
            var col = index % Data.Dimension;

            Data[row, col, 0] = scaleFactor * real;
            Data[row, col, 1] = scaleFactor * img;

            return true;
        }

        public override bool Dump(IQArray<Qubit>? qubits = null)
        {
            var count = qubits?.Length ?? Simulator.QubitManager!.AllocatedQubitsCount;
            var nQubitsPerRegister = ((int)count / 2);
            Data = new ComplexSquareMatrix(1 << nQubitsPerRegister);
            scaleFactor = System.Math.Sqrt(1 << nQubitsPerRegister);
            return base.Dump(qubits);
        }
    }

    internal partial class DumpReferenceAndTarget
    {
        public class Native : DumpReferenceAndTarget
        {
            private SimulatorBase? Simulator;

            public Native(IOperationFactory m) : base(m)
            {
                Simulator = m as SimulatorBase;
            }

            private QVoid DumpUnitaryFromChoiState(QuantumSimulator simulator, IQArray<Qubit> reference, IQArray<Qubit> target)
            {
                var arrayDumper = new ArrayDumper(simulator);
                arrayDumper.Dump(new QArray<Qubit>(reference.Concat(target)));
                Simulator?.MaybeDisplayDiagnostic(
                    new DisplayableUnitaryOperator
                    {
                        Data = arrayDumper.Data,
                        Qubits = target.ToList()
                    }
                );
                return QVoid.Instance;
            }

            public override Func<(IQArray<Qubit>, IQArray<Qubit>), QVoid> __Body__ => __in__ =>
            {
                var (reference, target) = __in__;
                return Simulator switch
                {
                    QuantumSimulator sim => this.DumpUnitaryFromChoiState(sim, reference, target),
                    // TODO: Add other simulators here as appropriate.
                    _ => base.__Body__(__in__)
                };
            };

        }
    }

}
