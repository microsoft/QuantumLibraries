// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.



using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;

namespace Microsoft.Quantum.Chemistry
{
    public static class Settings
    {
        public enum IndexConventionType
        {
            UpDown, HalfUp
        }
        public static IndexConventionType SetIndexConvention = IndexConventionType.HalfUp;
        
    }
}

