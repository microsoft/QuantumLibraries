// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Diagnostics.CodeAnalysis;
using Microsoft.Quantum.Diagnostics.Emulation;
using Microsoft.Quantum.Simulation.Common;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using Newtonsoft.Json;
using Newtonsoft.Json.Converters;
using Newtonsoft.Json.Linq;

namespace Microsoft.Quantum.Diagnostics
{

    public partial class AllowAtMostNQubits
    {
        public class Native : AllowAtMostNQubits
        {
            private SimulatorBase? Simulator;
            
            private static Stack<(IDisposable, IDisposable)> Handlers =
                new Stack<(IDisposable, IDisposable)>();

            public Native(IOperationFactory m) : base(m)
            {
                Simulator = m as SimulatorBase;
            }

            public override Func<(long, string), QVoid> Body => _args =>
            {
                if (Simulator == null) return QVoid.Instance;

                var (nQubitsAllowed, message) = _args;
                var nQubitsAllocated = 0L;

                Handlers.Push((
                    new ActionDisposer<Action<long>>(
                        nQubits =>
                        {
                            nQubitsAllocated += nQubits;
                            if (nQubitsAllocated > nQubitsAllowed)
                            {
                                throw new ExecutionFailException(
                                    $"{nQubitsAllocated} were allocated, but at most {nQubitsAllowed} are allowed:\n{message}."
                                );
                            }
                        },
                        setup: handler => Simulator.OnAllocateQubits += handler,
                        cleanup: handler => Simulator.OnAllocateQubits -= handler
                    ),

                    new ActionDisposer<Action<IQArray<Qubit>>>(
                        register =>
                        {
                            nQubitsAllocated -= register.Length;
                        },
                        setup: handler => Simulator.OnReleaseQubits += handler,
                        cleanup: handler => Simulator.OnReleaseQubits -= handler
                    )
                ));

                return QVoid.Instance;
            };

            public override Func<(long, string), QVoid> AdjointBody => _args =>
            {
                if (Simulator == null) return QVoid.Instance;

                var (start, end) = Handlers.Pop();
                start.Dispose();
                end.Dispose();

                return QVoid.Instance;
            };
        }

    }

}
