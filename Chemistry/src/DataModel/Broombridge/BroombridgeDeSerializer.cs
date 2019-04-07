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

namespace Microsoft.Quantum.Chemistry.Broombridge
{
    /// <summary>
    /// Broombridge deserializers
    /// </summary>
    public static class Deserializers
    {
        /// <summary>
        /// Returns version number of Broombridge file.
        /// </summary>
        /// <param name="filename">Broombridge file address.</param>
        /// <returns>Version number of Broombridge file</returns>
        public static VersionNumber GetVersionNumber(string filename)
        {
            using (var reader = File.OpenText(filename))
            {
                var deserializer = new DeserializerBuilder().Build();
                var data = deserializer.Deserialize<Dictionary<string, object>>(reader);
                var format = (Dictionary<object, object>)data["format"];
                var version = format["version"] as string ?? "";
                VersionNumber versionNumberFound = VersionNumber.NotRecognized;
                TryParseVersionNumber(version, out versionNumberFound);
                return versionNumberFound;
            }
        }

        /// <summary>
        /// Enumerable item for Broombridge version numbers.
        /// </summary>
        public enum VersionNumber
        {
            NotRecognized = -1, v0_1 = 0, v0_2 = 1
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
        /// Returns Broombridge deserialized into the current version data structure.
        /// Data structure is automatically updated to the current Broombridge version.
        /// </summary>
        /// <param name="filename">Broombridge file address.</param>
        /// <returns>Deserializer Broombridge data strauture.</returns>
        public static CurrentVersion.Data DeserializeBroombridge(string filename)
        {
            VersionNumber versionNumber = Deserializers.GetVersionNumber(filename);

            if (versionNumber == VersionNumber.v0_1)
            {
                return DataStructures.Update(Deserializers.DeserializeBroombridgev0_1(filename));
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

        /// <summary>
        /// Deserialize Broombridge v0.1 from a file into the Broombridge v0.1 data structure.
        /// </summary>
        /// <param name="filename">Broombridge filename to deserialize</param>
        /// <returns>Deserialized Broombridge v0.1 data.</returns>
        public static DataStructures.V0_1.Data DeserializeBroombridgev0_1(string filename)
        {
            using (var reader = File.OpenText(filename))
            {
                var deserializer = new DeserializerBuilder().Build();
                return deserializer.Deserialize<DataStructures.V0_1.Data>(reader);
            }
        }

        /// <summary>
        /// Deserialize Broombridge v0.2 from a file into the Broombridge v0.2 data structure.
        /// </summary>
        /// <param name="filename">Broombridge filename to deserialize</param>
        /// <returns>Deserialized Broombridge v0.2 data.</returns>
        public static DataStructures.V0_2.Data DeserializeBroombridgev0_2(string filename)
        {
            using (var reader = File.OpenText(filename))
            {
                var deserializer = new DeserializerBuilder().Build();
                return deserializer.Deserialize<DataStructures.V0_2.Data>(reader);
            }
        }
    }

    /// <summary>
    /// Broombridge serializers
    /// </summary>
    public static class Serializers
    {
        /// <summary>
        /// Broombridge serializer
        /// </summary>
        /// <param name="filename">Broombridge filename to serialize</param>
        /// <returns>Serialized Broombridge</returns>
        public static void SerializeBroombridgev0_2(DataStructures.V0_2.Data data, string filename)
        {
            var stringBuilder = new StringBuilder();
            var serializer = new Serializer();
            stringBuilder.AppendLine(serializer.Serialize(data));
            Console.WriteLine(stringBuilder);
            Console.WriteLine("");
        }
    } 
}