// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using Microsoft.Quantum.Intrinsic;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;

namespace Microsoft.Quantum.Canon
{
    public partial class ApplyAnd
    {
        /// <summary>
        /// Uses CCNOT implementation when executed using ToffoliSimulator.
        /// </summary>
        public class Native : ApplyAnd
        {
            private ToffoliSimulator? simulator = null;
            protected CCNOT? CCNOT { get; set; } = null;

            public Native(IOperationFactory m) : base(m)
            {
                simulator = m as ToffoliSimulator;
            }

            public override void __Init__()
            {
                base.__Init__();

                if (simulator != null)
                {
                    CCNOT = __Factory__.Get<CCNOT>(typeof(CCNOT));
                }
            }

            public override Func<(Qubit, Qubit, Qubit), QVoid> __Body__ => CCNOT?.__Body__ ?? base.__Body__;

            public override Func<(Qubit, Qubit, Qubit), QVoid> __AdjointBody__ => CCNOT?.__AdjointBody__ ?? base.__AdjointBody__;

            public override Func<(IQArray<Qubit>, (Qubit, Qubit, Qubit)), QVoid> __ControlledBody__ => CCNOT?.__ControlledBody__ ?? base.__ControlledBody__;

            public override Func<(IQArray<Qubit>, (Qubit, Qubit, Qubit)), QVoid> __ControlledAdjointBody__ => CCNOT?.__ControlledAdjointBody__ ?? base.__ControlledAdjointBody__;
        }
    }

    public partial class ApplyLowDepthAnd
    {
        /// <summary>
        /// Uses CCNOT implementation when executed using ToffoliSimulator.
        /// </summary>
        public class Native : ApplyLowDepthAnd
        {
            private ToffoliSimulator? simulator = null;
            protected CCNOT? CCNOT { get; set; } = null;

            public Native(IOperationFactory m) : base(m)
            {
                simulator = m as ToffoliSimulator;
            }

            public override void __Init__()
            {
                base.__Init__();

                if (simulator != null)
                {
                    CCNOT = __Factory__.Get<CCNOT>(typeof(CCNOT));
                }
            }

            public override Func<(Qubit, Qubit, Qubit), QVoid> __Body__ => CCNOT?.__Body__ ?? base.__Body__;

            public override Func<(Qubit, Qubit, Qubit), QVoid> __AdjointBody__ => CCNOT?.__AdjointBody__ ?? base.__AdjointBody__;

            public override Func<(IQArray<Qubit>, (Qubit, Qubit, Qubit)), QVoid> __ControlledBody__ => CCNOT?.__ControlledBody__ ?? base.__ControlledBody__;

            public override Func<(IQArray<Qubit>, (Qubit, Qubit, Qubit)), QVoid> __ControlledAdjointBody__ => CCNOT?.__ControlledAdjointBody__ ?? base.__ControlledAdjointBody__;
        }
    }
}
