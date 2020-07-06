// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using Microsoft.Quantum.Intrinsic;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;

namespace Microsoft.Quantum.Canon
{
    /// <summary>
    /// Uses CCNOT implementation when executed using ToffoliSimulator.
    /// </summary>
    public partial class ApplyAnd
    {
        public class Native : ApplyAnd
        {
            private ToffoliSimulator? simulator = null;
            protected CCNOT? CCNOT { get; set; } = null;

            public Native(IOperationFactory m) : base(m)
            {
                simulator = m as ToffoliSimulator;
            }

            public override void Init()
            {
                base.Init();

                if (simulator != null)
                {
                    CCNOT = Factory.Get<CCNOT>(typeof(CCNOT));
                }
            }

            public override Func<(Qubit, Qubit, Qubit), QVoid> Body => CCNOT?.Body ?? base.Body;

            public override Func<(Qubit, Qubit, Qubit), QVoid> AdjointBody => CCNOT?.AdjointBody ?? base.AdjointBody;

            public override Func<(IQArray<Qubit>, (Qubit, Qubit, Qubit)), QVoid> ControlledBody => CCNOT?.ControlledBody ?? base.ControlledBody;

            public override Func<(IQArray<Qubit>, (Qubit, Qubit, Qubit)), QVoid> ControlledAdjointBody => CCNOT?.ControlledAdjointBody ?? base.ControlledAdjointBody;
        }
    }

    /// <summary>
    /// Uses CCNOT implementation when executed using ToffoliSimulator.
    /// </summary>
    public partial class ApplyLowDepthAnd
    {
        public class Native : ApplyLowDepthAnd
        {
            private ToffoliSimulator simulator = null;
            protected CCNOT CCNOT { get; set; } = null;

            public Native(IOperationFactory m) : base(m)
            {
                simulator = m as ToffoliSimulator;
            }

            public override void Init()
            {
                base.Init();

                if (simulator != null)
                {
                    CCNOT = Factory.Get<CCNOT>(typeof(CCNOT));
                }
            }

            public override Func<(Qubit, Qubit, Qubit), QVoid> Body => CCNOT?.Body ?? base.Body;

            public override Func<(Qubit, Qubit, Qubit), QVoid> AdjointBody => CCNOT?.AdjointBody ?? base.AdjointBody;

            public override Func<(IQArray<Qubit>, (Qubit, Qubit, Qubit)), QVoid> ControlledBody => CCNOT?.ControlledBody ?? base.ControlledBody;

            public override Func<(IQArray<Qubit>, (Qubit, Qubit, Qubit)), QVoid> ControlledAdjointBody => CCNOT?.ControlledAdjointBody ?? base.ControlledAdjointBody;
        }
    }
}
