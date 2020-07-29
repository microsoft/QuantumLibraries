// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System.Text;
using Microsoft.Jupyter.Core;
using Microsoft.Extensions.Logging;
using Newtonsoft.Json;
using NumSharp;
using static NumSharp.Slice;
using Microsoft.Quantum.IQSharp.Jupyter;
using System.Linq;
using System;
using System.Collections.Generic;

namespace Microsoft.Quantum.Diagnostics.Emulation
{
    
    internal static class Extensions
    {
        internal static IEnumerable<object> EnumerateOverAxis(this NumSharp.NDArray array, int axis = 0)
        {
            var prefix = Enumerable.Repeat(All, axis).ToArray();
            foreach (var idx in Enumerable.Range(0, array.Shape[axis]))
            {
                Slice[] slice = prefix.Concat(new Slice[] { idx, Ellipsis }).ToArray();
                yield return array[slice];
            }
        }


    }

}
