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
    /// Enumerable item for Broombridge version numbers.
    /// </summary>
    public enum VersionNumber
    {
        NotRecognized = -1, v0_1 = 0, v0_2 = 1
    }

    /// <summary>
    /// Broombridge deserializers
    /// </summary>
    public static class Deserializers
    {
        /// <summary>
        /// Dictionary from version number strings to version number types.
        /// </summary>
        internal static Dictionary<string, VersionNumber> VersionNumberDict = new Dictionary<string, VersionNumber>()
        {
            // https://github.com/Microsoft/Quantum/blob/master/Chemistry/Schema/broombridge-0.1.schema.json
            {"0.1", VersionNumber.v0_1 },
            {"broombridge-0.1.schema", VersionNumber.v0_1 },
            // TODO: URL of 0.2 schema.
            {"0.2", VersionNumber.v0_2 },
            {"broombridge-0.2.schema", VersionNumber.v0_2 }
        };

        
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
                var schema = data["$schema"] as string;
                VersionNumber versionNumber = VersionNumber.NotRecognized;
                if(schema != null)
                {
                    foreach (var kv in VersionNumberDict)
                    {
                        if (schema.Contains(kv.Key))
                        {
                            versionNumber = kv.Value;
                            break;
                        }
                    }
                }
                return versionNumber;
            }
        }


        /// <summary>
        /// Returns Broombridge deserialized into the current version data structure.
        /// Data structure is automatically updated to the current Broombridge version.
        /// </summary>
        /// <param name="filename">Broombridge file address.</param>
        /// <returns>Deserializer Broombridge data strauture.</returns>
        public static Data DeserializeBroombridge(string filename)
        {
            VersionNumber versionNumber = GetVersionNumber(filename);
            var output = new V0_2.Data();
            if (versionNumber == VersionNumber.v0_1)
            {
                output = DataStructures.Update(Deserialize<V0_1.Data>(filename));
            }
            else if (versionNumber == VersionNumber.v0_2)
            {
                output = Deserialize<V0_2.Data>(filename);
            }
            else
            {
                throw new System.InvalidOperationException("Unrecognized Broombridge version number.");
            }
            return new Data(output);
        }

        /// <summary>
        /// Generic deserializer from a file into a data structure of type `TData`.
        /// </summary>
        /// <typeparam name="TData">Type of data to be deserialized.</typeparam>
        /// <param name="filename">Path to data to be deserialized.</param>
        /// <returns></returns>
        public static TData Deserialize<TData>(string filename)
        {
            using (var reader = File.OpenText(filename))
            {
                return new DeserializerBuilder()
                    .Build()
                    .Deserialize<TData>(reader);
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
        public static void SerializeBroombridgev0_2(V0_2.Data data, string filename)
        {
            var stringBuilder = new StringBuilder();
            var serializer = new Serializer();
            stringBuilder.AppendLine(serializer.Serialize(data));
            Console.WriteLine(stringBuilder);
            Console.WriteLine("");
        }
    } 
}