// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.QAOA.JupyterTests
{
    using System;
    using System.Collections.Generic;
    using Microsoft.Jupyter.Core;

    public class MockChannel : IChannel
    {
        public List<string> errors = new List<string>();
        public List<string> msgs = new List<string>();

        public void Display(object displayable)
        {
            throw new NotImplementedException();
        }

        public void Stderr(string message) => errors.Add(message);

        public void Stdout(string message) => msgs.Add(message);
    }
}
