// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using Microsoft.Quantum.Simulation.Common;
using Microsoft.Quantum.Simulation.Core;

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

}