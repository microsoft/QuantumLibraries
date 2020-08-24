// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using System.Linq;
using Microsoft.Quantum.Diagnostics.Emulation;
using Microsoft.Quantum.Simulation.Common;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using NumSharp;
using static NumSharp.Slice;

namespace Microsoft.Quantum.Diagnostics
{
    internal class ArrayDumper : QuantumSimulator.StateDumper
    {
        // NB: NumSharp does not yet support complex numbers, so we store data
        //     as an array with a trailing index of length 2.
        internal NumSharp.NDArray? Data = null;
        public ArrayDumper(QuantumSimulator sim) : base(sim)
        {
        }

        public override bool Callback(uint idx, double real, double img)
        {
            if (Data as object == null) throw new Exception("Expected data buffer to be initialized before callback, but it was null.");
            Data[(int)idx, 0] = real;
            Data[(int)idx, 1] = img;
            return true;
        }

        public override bool Dump(IQArray<Qubit>? qubits = null)
        {
            var count = qubits?.Length ?? Simulator.QubitManager!.GetAllocatedQubitsCount();
            var nQubitsPerRegister = ((int)count / 2);
            Data = np.empty(new Shape(1 << ((int)count), 2));
            var result = base.Dump(qubits);

            // At this point, _data should be filled with the full state
            // vector, so let's display it, counting on the right display
            // encoder to be there to pack it into a table.
            var scaleFactor = System.Math.Sqrt(1 << nQubitsPerRegister);
            Data = scaleFactor * Data.reshape(1 << nQubitsPerRegister, 1 << nQubitsPerRegister, 2);

            return result;
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
