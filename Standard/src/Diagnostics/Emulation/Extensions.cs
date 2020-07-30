// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using Microsoft.Quantum.Simulation.Common;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;

namespace Microsoft.Quantum.Diagnostics.Emulation
{

    internal class SimulatorEventDisposer : IDisposable
    {
        private SimulatorBase Simulator;
        private Action<ICallable, IApplyData> StartOperation, EndOperation;

        public SimulatorEventDisposer(
                SimulatorBase simulator,
                Action<ICallable, IApplyData> startOperation,
                Action<ICallable, IApplyData> endOperation
        )
        {
            Simulator = simulator;
            StartOperation = startOperation;
            EndOperation = endOperation;

            Simulator.OnOperationStart += startOperation;
            Simulator.OnOperationEnd += endOperation;
        }

        public void Dispose()
        {
            System.Console.WriteLine("Removing events.");
            Simulator.OnOperationStart -= StartOperation;
            Simulator.OnOperationEnd -= EndOperation;
        }
    }

    internal class ActionDisposer<TData> : IDisposable
    {
        private Action<TData> cleanup;
        private TData data;
        public ActionDisposer(TData data, Action<TData> setup, Action<TData> cleanup)
        {
            this.data = data;
            this.cleanup = cleanup;
            setup(data);
        }

        public void Dispose() =>
            cleanup(data);
    }

    internal static class Extensions
    {
        internal static SimulatorEventDisposer RegisterOperationHandlers(
                this SimulatorBase simulator,
                Action<ICallable, IApplyData> startOperation,
                Action<ICallable, IApplyData>? endOperation = null
        ) => new SimulatorEventDisposer(
            simulator,
            startOperation,
            endOperation == null
                ? (callable, data) => {}
                : endOperation
        );
    }
}
