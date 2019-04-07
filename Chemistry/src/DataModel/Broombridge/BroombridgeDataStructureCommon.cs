// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

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
    using System;

    /// <summary>
    /// Alias for Broombridge latest supported format.
    /// </summary>
    public class CurrentVersion : DataStructures.V0_2 { }
    /// <summary>
    /// Functionality to validate and parse Broombridge
    /// </summary>
    public static partial class DataStructures
    {

        public static class DataStructure
        {
            public struct Format
            {
                [YamlMember(Alias = "version", ApplyNamingConventions = false)]
                public string Version { get; set; }
            }

            public struct Generator
            {
                [YamlMember(Alias = "source", ApplyNamingConventions = false)]
                public string Source { get; set; }

                [YamlMember(Alias = "version", ApplyNamingConventions = false)]
                public string Version { get; set; }
            }


            public enum BibliographyKind
            {
                arXiv, DOI, URL
            }
            public struct BibliographyItem : IYamlConvertible
            {
                public string Value { get; set; }
                public BibliographyKind Kind { get; set; }

                public void Read(IParser parser, Type expectedType, ObjectDeserializer nestedObjectDeserializer)
                {
                    var dict = (Dictionary<string, string>)nestedObjectDeserializer(typeof(Dictionary<string, string>));
                    // There should only be one item.
                    Debug.Assert(dict.Count == 1);
                    var kind = dict.Keys.ElementAt(0);

                    Value = dict[kind];
                    if (kind.ToLowerInvariant() == "arxiv")
                    {
                        Kind = BibliographyKind.arXiv;
                    }
                    else if (kind.ToLowerInvariant() == "doi")
                    {
                        Kind = BibliographyKind.DOI;
                    }
                    else
                    {
                        Kind = BibliographyKind.URL;
                    }
                }

                public void Write(IEmitter emitter, ObjectSerializer nestedObjectSerializer)
                {
                    nestedObjectSerializer(new Dictionary<string, string>
                    {
                        [Kind.Map(
                            (BibliographyKind.arXiv, () => "arXiv"),
                            (BibliographyKind.DOI, () => "doi"),
                            (BibliographyKind.URL, () => "url")
                        )] = Value
                    });
                }
            }

            public struct BasisSet
            {
                [YamlMember(Alias = "type", ApplyNamingConventions = false)]
                public string Type { get; set; }

                [YamlMember(Alias = "name", ApplyNamingConventions = false)]
                public string Name { get; set; }
            }

            public class HasUnits
            {
                [YamlMember(Alias = "units", ApplyNamingConventions = false)]
                // FIXME: make this an enum of allowed units.
                public string Units { get; set; }
            }

            public class Geometry : HasUnits
            {
                [YamlMember(Alias = "coordinate_system", ApplyNamingConventions = false)]
                public string CoordinateSystem { get; set; }

                [YamlMember(Alias = "symmetry", ApplyNamingConventions = false)]
                public string Symmetry { get; set; }

                [YamlMember(Alias = "atoms", ApplyNamingConventions = false)]
                // FIXME: new struct for atoms
                public List<Dictionary<string, object>> Atoms { get; set; }
            }

            public struct HamiltonianData
            {
                [YamlMember(Alias = "particle_hole_representation", ApplyNamingConventions = false)]
                // TODO: make this not object
                // FIXME: currently strips off the last element as the "value", but the
                //        present schema requires us to pull off the last two as a (double, string).
                public DataStructure.ArrayQuantity<object, object> ParticleHoleRepresentation { get; set; }

                [YamlMember(Alias = "one_electron_integrals", ApplyNamingConventions = false)]
                public DataStructure.ArrayQuantity<long, double> OneElectronIntegrals { get; set; }

                [YamlMember(Alias = "two_electron_integrals", ApplyNamingConventions = false)]
                public DataStructure.ArrayQuantity<long, double> TwoElectronIntegrals { get; set; }

            }

            public class SimpleQuantity : HasUnits
            {
                [YamlMember(Alias = "value", ApplyNamingConventions = false)]
                public double Value { get; set; }
            }

            public class BoundedQuantity : HasUnits
            {
                [YamlMember(Alias = "value", ApplyNamingConventions = false)]
                public double? Value { get; set; }

                [YamlMember(Alias = "upper", ApplyNamingConventions = false)]
                public double Upper { get; set; }

                [YamlMember(Alias = "lower", ApplyNamingConventions = false)]
                public double Lower { get; set; }
            }

            public class ArrayQuantity<TIndex, TValue> : HasUnits, IYamlConvertible
            {
                // TODO: make this an enum.
                public string Format { get; set; }
                public List<(TIndex[], TValue)> Values { get; set; }

                public void Read(IParser parser, Type expectedType, ObjectDeserializer nestedObjectDeserializer)
                {
                    // We need to read units ourselves since we've overriden IYamlConvertible interface.
                    var data = (Dictionary<string, object>)nestedObjectDeserializer(typeof(Dictionary<string, object>));
                    Units = (string)data["units"];
                    Format = (string)data["format"];
                    Values = new List<(TIndex[], TValue)>();
                    foreach (var value in ((IEnumerable<object>)data["values"]))
                    {
                        var entries = (IEnumerable<object>)value;
                        var a = entries.Take(entries.Count() - 1).Select(e => (TIndex)Convert.ChangeType(e, typeof(TIndex))).ToArray();
                        var q = (TValue)Convert.ChangeType(entries.Last(), typeof(TValue));
                        Values.Add((a, q));
                    }
                }

                public void Write(IEmitter emitter, ObjectSerializer nestedObjectSerializer)
                {
                    nestedObjectSerializer(new Dictionary<string, object>
                    {
                        ["units"] = Units,
                        ["format"] = Format,
                        ["values"] = Values
                            .Select(entry =>
                                entry.Item1.Select((idx) => (object)idx).Concat(new object[] { entry.Item2 })
                            )
                            .ToList()
                    });
                }
            }

            /// <summary>
            /// Coupled-cluster operator
            /// </summary>
            public struct ClusterOperator
            {
                [YamlMember(Alias = "reference_state", ApplyNamingConventions = false)]
                public List<string> Reference { get; set; }

                [YamlMember(Alias = "one_body_amplitudes", ApplyNamingConventions = false)]
                public List<List<string>> OneBodyAmplitudes { get; set; }

                [YamlMember(Alias = "two_body_amplitudes", ApplyNamingConventions = false)]
                public List<List<string>> TwoBodyAmplitudes { get; set; }
            }

        }



    }
}