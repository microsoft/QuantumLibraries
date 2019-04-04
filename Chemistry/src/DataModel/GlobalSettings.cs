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
        public static class IndexConvention {
            public enum Type
            {
                UpDown, HalfUp
            }
            public static Type Default()
            {
                return Type.HalfUp;
            }
            public static Type Current = Default();
                
        }        
    }
}

