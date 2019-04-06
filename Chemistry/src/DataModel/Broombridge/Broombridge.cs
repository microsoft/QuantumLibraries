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
        /// Enumerable item for Broombridge version numbers.
        /// </summary>
        public enum VersionNumber
        {
            NotRecognized, v0_1 = 0, v0_2 = 1
        }
        
        /// <summary>
        /// Dictionary from version number strings to version number types.
        /// </summary>
        internal static Dictionary<string, VersionNumber> VersionNumberDict = new Dictionary<string, VersionNumber>()
        {
            {"0.1", VersionNumber.v0_1 },
            {"0.2", VersionNumber.v0_2 }
        };
        
        /// <summary>
        /// Parse version number string.
        /// </summary>
        /// <param name="versionNumber">Version number string.</param>
        /// <param name="parsedVersionNumber">Parsed version number in enum `Number` format.</param>
        /// <returns>Returns <c>true</c> if version number string parsed successfully. Returns <c>false</c> otherwise.</returns>
        public static bool TryParseVersionNumber(string versionNumber, out VersionNumber parsedVersionNumber)
        {
            if (VersionNumberDict.ContainsKey(versionNumber))
            {
                parsedVersionNumber = VersionNumberDict[versionNumber];
                return true;
            }
            else
            {
                parsedVersionNumber = VersionNumber.NotRecognized;
                return false;
            }
            return false;
        }

    }
    // Parts of this might be merged intro Broombridge parsing due to overlap.
    
}