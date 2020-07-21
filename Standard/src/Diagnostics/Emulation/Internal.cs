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

    internal partial class FormattedFailure<__T__>
    {

        public class Native : FormattedFailure<__T__>
        {
            private SimulatorBase? Simulator;

            public Native(IOperationFactory m) : base(m)
            {
                Simulator = m as SimulatorBase;
            }

            public override Func<(__T__,__T__,String), QVoid> Body => (__in__) =>
            {
                if (Simulator != null)
                {
                    var (actual, expected, message) = __in__;
                    Simulator.MaybeDisplayDiagnostic(
                        new FailureRecord<__T__>
                        {
                            Message = message,
                            Actual = actual,
                            Expected = expected
                        }
                    );
                }
                return base.Body(__in__);
            };
        }

    }

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
            var count = qubits == null
                        ? this.Simulator.QubitManager!.GetAllocatedQubitsCount()
                        : qubits.Length;
            var nQubitsPerRegister = ((int) count / 2);
            Data = np.empty(new Shape(1 << ((int)count), 2));
            var result = base.Dump(qubits);

            // At this point, _data should be filled with the full state
            // vector, so let's display it, counting on the right display
            // encoder to be there to pack it into a table.
            var scaleFactor = System.Math.Sqrt(1 << nQubitsPerRegister);
            Data = scaleFactor * Data.reshape(1 << nQubitsPerRegister, 1 << nQubitsPerRegister, 2);

            // Clean up the state vector buffer.
            Data = null;

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

            private QVoid DumpUnitary(QuantumSimulator simulator, IQArray<Qubit> reference, IQArray<Qubit> target)
            {
                var arrayDumper = new ArrayDumper(simulator);
                arrayDumper.Dump(new QArray<Qubit>(reference.Concat(target)));
                // TODO: do something to get state of the two registers out from simulator.
                // TODO: write operation out as unitary.
                return QVoid.Instance;
            }

            public override Func<(IQArray<Qubit>, IQArray<Qubit>), QVoid> Body => (__in__) =>
            {
                var (reference, target) = __in__;
                return Simulator switch
                {
                    QuantumSimulator sim => DumpUnitary(sim, reference, target),
                    // TODO: Add Toffoli simulator here.
                    _ => base.Body(__in__)
                };
            };

        }
    }

}