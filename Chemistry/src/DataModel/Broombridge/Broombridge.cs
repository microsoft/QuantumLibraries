// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
using System.Linq;
using System.Text.RegularExpressions;
using YamlDotNet.Core;
using YamlDotNet.Serialization;

namespace Microsoft.Quantum.Chemistry
{
    /// <summary>
    /// Functionality to validate and parse Broombridge
    /// </summary>
    public static partial class Broombridge
    {
        /// <summary>
        /// Alias for Broombridge latest supported format.
        /// </summary>
        public class Current : V0_2
        {

        }

        /// <summary>
        /// Class containing Broombridge version number management.
        /// </summary>
        public static class Version
        {
            /// <summary>
            /// Enumerable item for Broombridge version numbers.
            /// </summary>
            public enum Number
            {
                v0_1 = 0, v0_2 = 1
            }

            /// <summary>
            /// Dictionary from version number strings to version number types.
            /// </summary>
            public static Dictionary<string, Number> VersionNumberDict = new Dictionary<string, Number>()
            {
                {"0.1", Number.v0_1 },
                {"0.2", Number.v0_2 }
            };

            /// <summary>
            /// Parse version number string
            /// </summary>
            /// <param name="versionNumber">Version number string</param>
            /// <returns>Version number in enum `Number` format.</returns>
            public static Number ParseVersionNumber(string versionNumber)
            {
                return VersionNumberDict[versionNumber];
            }
        }

    }
    // Parts of this might be merged intro Broombridge parsing due to overlap.
    
}