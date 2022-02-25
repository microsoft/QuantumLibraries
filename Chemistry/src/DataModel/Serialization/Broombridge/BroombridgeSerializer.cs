// Copyright (c) Microsoft Corporation.
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
using YamlDotNet.Core.Events;
using System.Reflection;
using System.Runtime.Serialization;

namespace Microsoft.Quantum.Chemistry.Broombridge
{
    /// <summary>
    /// Enumerable item for Broombridge version numbers.
    /// </summary>
    public enum VersionNumber
    {
        NotRecognized = -1, v0_1 = 0, v0_2 = 1, v0_3
    }

    internal class YamlStringEnumConverter : IYamlTypeConverter
    {
        public bool Accepts(Type type) => type.IsEnum;

        public object ReadYaml(IParser parser, Type type)
        {
            var parsedEnum = parser.Consume<Scalar>();
            var serializableValues = type.GetMembers()
                .Select(m => new KeyValuePair<string, MemberInfo>(m.GetCustomAttributes<EnumMemberAttribute>(true).Select(ema => ema.Value).FirstOrDefault(), m))
                .Where(pa => !string.IsNullOrEmpty(pa.Key))
                .ToDictionary(pa => pa.Key, pa => pa.Value);

            if (!serializableValues.ContainsKey(parsedEnum.Value))
            {
                throw new YamlException(parsedEnum.Start, parsedEnum.End, $"Value '{parsedEnum.Value}' not found in enum '{type.Name}'");
            }

            return Enum.Parse(type, serializableValues[parsedEnum.Value].Name);
        }

        public void WriteYaml(IEmitter emitter, object value, Type type)
        {
            var enumMember = type.GetMember(value.ToString()).FirstOrDefault();
            var yamlValue = enumMember?.GetCustomAttributes<EnumMemberAttribute>(true).Select(ema => ema.Value).FirstOrDefault() ?? value.ToString();
            emitter.Emit(new Scalar(yamlValue));
        }
    }

    public static class BroombridgeSerializer
    {
        public static IEnumerable<ElectronicStructureProblem> Deserialize(TextReader reader)
        {
            var data = Deserializers.DeserializeBroombridge(reader);
            return data
                .Raw
                .ProblemDescriptions
                .Select(
                    problem => new ElectronicStructureProblem
                    {
                        BasisSet = problem.BasisSet?.FromBroombridgeV0_1(),
                        CoulombRepulsion = problem.CoulombRepulsion.FromBroombridgeV0_1(),
                        Geometry = problem.Geometry?.FromBroombridgeV0_1(),
                        EnergyOffset = problem.EnergyOffset.FromBroombridgeV0_1(),
                        FciEnergy = problem.FciEnergy?.FromBroombridgeV0_1(),
                        InitialStates = problem.InitialStates?.FromBroombridgeV0_2(),
                        Metadata = problem.Metadata,
                        NElectrons = problem.NElectrons,
                        NOrbitals = problem.NOrbitals,
                        OrbitalIntegralHamiltonian = V0_3.ToOrbitalIntegralHamiltonian(problem),
                        ScfEnergy = problem.ScfEnergy?.FromBroombridgeV0_1(),
                        ScfEnergyOffset = problem.ScfEnergyOffset?.FromBroombridgeV0_1()
                    }
                );
        }

        public static void Serialize(TextWriter writer, IEnumerable<ElectronicStructureProblem> problems)
        {
            Serializers.SerializeData(
                new Broombridge.V0_3.Data
                {
                    // TODO: fix additional properties by converting IEnumerable<ESP> to
                    //       new problem collection class. See https://github.com/microsoft/QuantumLibraries/issues/287.
                    Bibliography = null,
                    Format = new V0_1.Format
                    {
                        Version = "0.2"
                    },
                    Generator = new V0_1.Generator
                    {
                        Source = "qdk-chem",
                        Version = typeof(BroombridgeSerializer).Assembly.GetName().Version.ToString()
                    },
                    Schema = V0_2.SchemaUrl,
                    ProblemDescriptions = problems
                        .Select(
                            problem => problem.ToBroombridgeV0_3()
                        )
                        .ToList()
                },
                writer
            );
        }
    }

    /// <summary>
    /// Broombridge deserializers
    /// </summary>
    [Obsolete(
        "Please use BroombridgeSerializer instead.",
        false
    )]
    public static class Deserializers
    {
        /// <summary>
        /// Dictionary from version number strings to version number types.
        /// </summary>
        internal static Dictionary<string, VersionNumber> VersionNumberDict = new Dictionary<string, VersionNumber>()
        {
            // https://github.com/microsoft/Quantum/blob/main/Chemistry/Schema/broombridge-0.1.schema.json
            ["0.1"] = VersionNumber.v0_1,
            ["broombridge-0.1.schema"] = VersionNumber.v0_1,
            // TODO: URL of 0.2 schema.
            ["0.2"] = VersionNumber.v0_2,
            ["broombridge-0.2.schema"] = VersionNumber.v0_2,
            // TODO: Update actual schema for 0.3.
            ["0.3"] = VersionNumber.v0_3,
            ["broombridge-0.3.schema"] = VersionNumber.v0_3
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
            var versionNumber = VersionNumber.NotRecognized;
            if (schema != null)
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
                    VersionNumber.v0_2 => DataStructures.Update(
                        Deserialize<V0_2.Data>(stringReader)
                    ),
                    VersionNumber.v0_3 => Deserialize<V0_3.Data>(stringReader),
                    VersionNumber.NotRecognized => throw new System.InvalidOperationException(
                        "Unrecognized Broombridge version number."
                    ),
                    _ => throw new System.Exception($"Internal error occurred; version {versionNumber} is valid but was not deserialized.")
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
            return DeserializeBroombridge(reader);
        }

        /// <summary>
        /// Generic deserializer from a file into a data structure of type `TData`.
        /// </summary>
        /// <typeparam name="TData">Type of data to be deserialized.</typeparam>
        /// <param name="reader">Stream for reading Broombridge data.</param>
        /// <returns></returns>
        public static TData Deserialize<TData>(TextReader reader) =>
            new DeserializerBuilder()
                .WithTypeConverter(new YamlStringEnumConverter())
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
        /// <param name="data">Broombridge v0.2 data to be serialized.</param>
        /// <param name="filename">Name of the file to write serialized data to.</param>
        internal static void SerializeData<TData>(TData data, string filename)
        {
            using var writer = new StreamWriter(File.OpenWrite(filename));
            var stringBuilder = new StringBuilder();
            var serializer = new Serializer();
            stringBuilder.AppendLine(serializer.Serialize(data));
        }

        /// <summary>
        /// Broombridge serializer
        /// </summary>
        /// <param name="data">Broombridge v0.2 data to be serialized.</param>
        /// <param name="writer">Text writer to write serialized Broombridge data to.</param>
        internal static void SerializeData<TData>(TData data, TextWriter writer)
        {
            var stringBuilder = new StringBuilder();
            var serializer =
                new SerializerBuilder()
                .ConfigureDefaultValuesHandling(DefaultValuesHandling.OmitNull)
                .WithTypeConverter(new DoubleYamlConverter())
                .Build();
            writer.WriteLine(serializer.Serialize(data));
            Console.WriteLine("");
        }
    }

    /// <summary>
    ///     YAML type converter used to ensure that the double value <c>0</c>
    ///     in <c>"G17"</c> converts to <c>0.0</c> to help assist compatibility
    ///     with dynamically typed languages.
    /// </summary>
    internal class DoubleYamlConverter : IYamlTypeConverter
    {
        /// <inheritdoc />
        public bool Accepts(Type type) => type == typeof(double);

        /// <inheritdoc />
        public object ReadYaml(IParser parser, Type type) =>
            Double.TryParse(parser.Consume<Scalar>()?.Value, out var value)
            ? value as object
            : null;

        /// <inheritdoc />
        public void WriteYaml(IEmitter emitter, object value, Type type) =>
            emitter.Emit(new Scalar(
                null, null,
                ((double)value) == 0
                ? "0.0"
                : ((double)value).ToString("G17"),
                ScalarStyle.Any, true, false
            ));
    }
}
