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
        /// Enumerable item for Broombridge version numbers.
        /// </summary>
        public enum VersionNumber
        {
            NotRecognized = -1, v0_1 = 0, v0_2 = 1
        }

        /// <summary>
        /// Alias for Broombridge latest supported format.
        /// </summary>
        public class Current : V0_2 { }

        /// <summary>
        /// Returns Broombridge deserialized into the current version data structure.
        /// Data structure is automatically updated to the current Broombridge version.
        /// </summary>
        /// <param name="filename">Broombridge file address.</param>
        /// <returns>Deserializer Broombridge data strauture.</returns>
        public static Current.Data DeserializeBroombridge(string filename)
        {
            VersionNumber versionNumber = Deserializers.GetVersionNumber(filename);

            if (versionNumber == VersionNumber.v0_1)
            {
                return Updater.Data(Deserializers.DeserializeBroombridgev0_1(filename));
            }
            else if (versionNumber == VersionNumber.v0_2)
            {
                return Deserializers.DeserializeBroombridgev0_2(filename);
            }
            else
            {
                throw new System.InvalidOperationException("Unrecognized Broombridge version number.");
            }
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
        }

    }
    
}