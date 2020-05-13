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
        ///     Returns version number of a Broombridge file.
        /// </summary>
        /// <param name="filename">Path to a Broombridge file.</param>
        /// <returns>Version number of Broombridge file</returns>
        public static VersionNumber GetVersionNumber(string filename)
        {
            using var reader = File.OpenText(filename);
            return GetVersionNumber(reader);
        }

        /// <summary>
        ///     Returns version number of a Broombridge file.
        /// </summary>
        /// <param name="reader">Stream for reading Broombridge data.</param>
        /// <returns>Version number of Broombridge file</returns>
        public static VersionNumber GetVersionNumber(TextReader reader)
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


        /// <summary>
        /// Returns Broombridge deserialized into the current version data structure.
        /// Data structure is automatically updated to the current Broombridge version.
        /// </summary>
        /// <param name="reader">Stream for reading Broombridge data.</param>
        /// <returns>Deserializer Broombridge data structure.</returns>
        public static Data DeserializeBroombridge(TextReader reader)
        {
            // We'll need the stream twice: once to get the version of
            // Broombridge used to serialize the data,
            // and again to actually deserialize once we know the right
            // version. Since we can't actually duplicate the stream,
            // we'll read into memory first, then use a string reader
            // in both cases.
            var rawData = reader.ReadToEnd();
            var versionNumber = GetVersionNumber(new StringReader(rawData));
            var stringReader = new StringReader(rawData);
            return new Data(
                versionNumber switch
                {
                    VersionNumber.v0_1 => DataStructures.Update(
                        Deserialize<V0_1.Data>(stringReader)
                    ),
                    VersionNumber.v0_2 => Deserialize<V0_2.Data>(stringReader),
                    _ => throw new System.InvalidOperationException(
                        "Unrecognized Broombridge version number."
                    )
                }
            );
        }

        /// <summary>
        /// Returns Broombridge deserialized into the current version data structure.
        /// Data structure is automatically updated to the current Broombridge version.
        /// </summary>
        /// <param name="filename">Path to a Broombridge file.</param>
        /// <returns>Deserializer Broombridge data structure.</returns>
        public static Data DeserializeBroombridge(string filename)
        {
            using var reader = File.OpenText(filename);
            return DeserializeBroombridge(filename);
        }

        /// <summary>
        /// Generic deserializer from a file into a data structure of type `TData`.
        /// </summary>
        /// <typeparam name="TData">Type of data to be deserialized.</typeparam>
        /// <param name="reader">Stream for reading Broombridge data.</param>
        /// <returns></returns>
        public static TData Deserialize<TData>(TextReader reader) =>
            new DeserializerBuilder()
                .Build()
                .Deserialize<TData>(reader);

        /// <summary>
        /// Generic deserializer from a file into a data structure of type `TData`.
        /// </summary>
        /// <typeparam name="TData">Type of data to be deserialized.</typeparam>
        /// <param name="filename">Path to data to be deserialized.</param>
        /// <returns></returns>
        public static TData Deserialize<TData>(string filename)
        {
            using var reader = File.OpenText(filename);
            return Deserialize<TData>(reader);
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