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
        /// Broombridge deserializers
        /// </summary>
        public static class Deserializers
        {
            /// <summary>
            /// Returns version number of Broombridge file.
            /// </summary>
            /// <param name="filename">Broombridge file address.</param>
            /// <returns>Version number of Broombridge file</returns>
            public static Version.Number GetVersionNumber(string filename)
            {
                using (var reader = File.OpenText(filename))
                {
                    var deserializer = new DeserializerBuilder().Build();
                    var data = deserializer.Deserialize<Dictionary<string, object>>(reader);
                    var format = (Dictionary<object, object>)data["format"];
                    var version = (string)format["version"];
                    return Version.ParseVersionNumber(version);
                }
            }

            /// <summary>
            /// Returns deserializer Broombridge data strauture.
            /// Data structure is automatically updated to the current Broombridge version.
            /// </summary>
            /// <param name="filename">Broombridge file address.</param>
            /// <returns>Deserializer Broombridge data strauture.</returns>
            public static Current.Data Source(string filename)
            {
                Version.Number versionNumber = GetVersionNumber(filename);

                if (versionNumber == Version.Number.v0_1)
                {
                    return Updater.Data(v0_1(filename));
                }
                else if (versionNumber == Version.Number.v0_2)
                {
                    return v0_2(filename);
                }
                else
                {
                    throw new System.InvalidOperationException("Unrecognized Broombridge version number.");
                }
            }

            /// <summary>
            /// Broombridge v0.1 deserializer
            /// </summary>
            /// <param name="filename">Broombridge filename to deserialize</param>
            /// <returns>Deserialized Broombridge v0.1 data.</returns>
            public static V0_1.Data v0_1(string filename)
            {
                using (var reader = File.OpenText(filename))
                {
                    var deserializer = new DeserializerBuilder().Build();
                    return deserializer.Deserialize<V0_1.Data>(reader);
                }
            }

            /// <summary>
            /// Broombridge deserializer v0.2
            /// </summary>
            /// <param name="filename">Broombridge filename to deserialize</param>
            /// <returns>Deserialized Broombridge v0.2 data.</returns>
            public static V0_2.Data v0_2(string filename)
            {
                using (var reader = File.OpenText(filename))
                {
                    var deserializer = new DeserializerBuilder().Build();
                    return deserializer.Deserialize<V0_2.Data>(reader);
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
            public static void v0_2(V0_2.Data data, string filename)
            {
                var stringBuilder = new StringBuilder();
                var serializer = new Serializer();
                stringBuilder.AppendLine(serializer.Serialize(data));
                Console.WriteLine(stringBuilder);
                Console.WriteLine("");
            }
        }
    }       
}