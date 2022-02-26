// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
using System.Linq;
using System.Text.RegularExpressions;
using YamlDotNet.Core;
using YamlDotNet.Serialization;
using System.Numerics;

using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Newtonsoft.Json;
using System.Runtime.Serialization;

using static Microsoft.Quantum.Chemistry.OrbitalIntegrals.IndexConventionConversions;
using YamlDotNet.Core.Events;
using System.Reflection;
using System.Collections.Immutable;

namespace Microsoft.Quantum.Chemistry.Broombridge
{
    // What data structures are unmodified from previous versions?
    using Format = V0_1.Format;
    using Generator = V0_1.Generator;
    using BibliographyItem = V0_1.BibliographyItem;
    using Geometry = V0_1.Geometry;
    using HasUnits = V0_1.HasUnits;
    using BasisSet = V0_1.BasisSet;
    using SimpleQuantity = V0_1.SimpleQuantity;
    using BoundedQuantity = V0_1.BoundedQuantity;
    using State = V0_2.State;
    using ClusterOperator = V0_2.ClusterOperator;

    internal static class BroombridgeExtensionsV0_3
    {
        internal static V0_3.ProblemDescription ToBroombridgeV0_3(
                this ElectronicStructureProblem problem
        ) => new V0_3.ProblemDescription
            {
                BasisSet = problem.BasisSet != null
                            ? new V0_1.BasisSet
                            {
                                Name = problem.BasisSet?.Name,
                                Type = problem.BasisSet?.Type
                            }
                            : null,
                CoulombRepulsion = problem.CoulombRepulsion.ToBroombridgeV0_2(),
                EnergyOffset = problem.EnergyOffset.ToBroombridgeV0_2(),
                FciEnergy = problem.FciEnergy?.ToBroombridgeV0_2(),
                Geometry = problem.Geometry?.ToBroombridgeV0_2(),
                Hamiltonian = problem.OrbitalIntegralHamiltonian.ToBroombridgeV0_3(),
                InitialStates = problem.InitialStates?.ToBroombridgeV0_2(),
                Metadata = problem.Metadata,
                NElectrons = problem.NElectrons,
                NOrbitals = problem.NOrbitals,
                ScfEnergy = problem.ScfEnergy?.ToBroombridgeV0_2(),
                ScfEnergyOffset = problem.ScfEnergyOffset?.ToBroombridgeV0_2()
            };

        internal static V0_3.ArrayQuantityWithSymmetry<TOut, TValue> TransformKeys<TIn, TOut, TValue>(
            this V0_3.ArrayQuantityWithSymmetry<TIn, TValue> arrayQuantity,
            Func<TIn, TOut> transform
        ) =>
            new V0_3.ArrayQuantityWithSymmetry<TOut, TValue>
            {
                Format = arrayQuantity.Format,
                IndexConvention = arrayQuantity.IndexConvention,
                Symmetry = arrayQuantity.Symmetry,
                Units = arrayQuantity.Units,
                Values = arrayQuantity
                    .Values
                    .Select(item => new V0_3.ArrayQuantityWithSymmetry<TOut, TValue>.Item
                    {
                        Key = transform(item.Key),
                        Value = item.Value
                    })
                    .ToList()
            };

        internal static V0_3.ArrayQuantityWithSymmetry<TKey, TValue> WithSymmetry<TKey, TValue>(
            this V0_3.ArrayQuantity<TKey, TValue> arrayQuantity,
            V0_3.Symmetry symmetry
        ) =>
            new V0_3.ArrayQuantityWithSymmetry<TKey, TValue>
            {
                Format = arrayQuantity.Format,
                IndexConvention = arrayQuantity.IndexConvention,
                Symmetry = symmetry,
                Units = arrayQuantity.Units,
                Values = arrayQuantity.Values
            };

        internal static V0_3.HamiltonianData ToBroombridgeV0_3(this OrbitalIntegralHamiltonian hamiltonian)
        {
            var twoElectronIntegrals = hamiltonian
                    .Terms[TermType.OrbitalIntegral.TwoBody]
                .ToBroombridgeV0_3(new V0_3.Symmetry
                {
                    // List terms explicitly with no compression.
                    Permutation = V0_3.PermutationSymmetry.Trivial
                });
            twoElectronIntegrals.IndexConvention = OrbitalIntegral.Convention.Mulliken;
            return new V0_3.HamiltonianData
            {
                OneElectronIntegrals = hamiltonian
                    .Terms[TermType.OrbitalIntegral.OneBody]
                    .ToBroombridgeV0_3(new V0_3.Symmetry
                    {
                        // List terms explicitly with no compression.
                        Permutation = V0_3.PermutationSymmetry.Trivial
                    })
                    .TransformKeys(idxs => (idxs[0], idxs[1])),
                TwoElectronIntegrals = twoElectronIntegrals
                    .TransformKeys(idxs => (idxs[0], idxs[1], idxs[2], idxs[3]))
            };
        }

        internal static V0_3.ArrayQuantityWithSymmetry<long[], double> ToBroombridgeV0_3(
            this Dictionary<OrbitalIntegrals.OrbitalIntegral, DoubleCoeff> terms,
            V0_3.Symmetry symmetry
        ) =>
            new V0_3.ArrayQuantityWithSymmetry<long[], double>()
            {
                Format = V0_3.ArrayFormat.Sparse,
                Units = "hartree",
                Values = terms.Select(term =>
                {
                    var idxs = ConvertIndices(
                                    term
                                    .Key
                                    .ToCanonicalForm(symmetry.Permutation.Value.FromBroombridgeV0_3())
                                    .OrbitalIndices,
                                    OrbitalIntegral.Convention.Dirac,
                                    OrbitalIntegral.Convention.Mulliken
                                )
                                .ToOneBasedIndices()
                                .Select(idx => (long)idx)
                                .ToArray();
                    return new V0_3.ArrayQuantity<long[], double>.Item
                    {
                        Key = idxs,
                        Value = term.Value.Value
                    };
                }).ToList(),
                Symmetry = symmetry
            };

        internal static V0_3.HamiltonianData ToBroombridgeV0_3(this V0_1.HamiltonianData hamiltonianData) =>
            new V0_3.HamiltonianData
            {
                ParticleHoleRepresentation = hamiltonianData.ParticleHoleRepresentation?.ToBroombridgeV0_3(),
                OneElectronIntegrals = hamiltonianData
                    .OneElectronIntegrals
                    .ToBroombridgeV0_3()
                    .WithSymmetry(new V0_3.Symmetry
                    {
                        Permutation = V0_3.PermutationSymmetry.Eightfold
                    })
                    .TransformKeys(key => (key[0], key[1])),
                TwoElectronIntegrals = hamiltonianData
                    .TwoElectronIntegrals
                    .ToBroombridgeV0_3()
                    .WithSymmetry(new V0_3.Symmetry
                    {
                        Permutation = V0_3.PermutationSymmetry.Eightfold
                    })
                    .TransformKeys(key => (key[0], key[1], key[2], key[3])),
            };

        internal static V0_3.ArrayQuantity<TKey[], TValue> ToBroombridgeV0_3<TKey, TValue>(this V0_1.ArrayQuantity<TKey, TValue> array) =>
            new V0_3.ArrayQuantity<TKey[], TValue>
            {
                Format = Enum.TryParse<V0_3.ArrayFormat>(array.Format, true, out var format)
                         ? format
                         : throw new Exception($"Invalid array format {array.Format} when converting 0.1 array quantity to 0.3 array quantity."),
                IndexConvention = array.IndexConvention,
                Units = array.Units,
                Values = array
                    .Values
                    .Select(item => new V0_3.ArrayQuantity<TKey[], TValue>.Item
                    {
                        Key = item.Item1,
                        Value = item.Item2
                    })
                    .ToList()
            };

        internal static (O, O) Select<I, O>(this (I, I) value, Func<I, O> func) =>
            (func(value.Item1), func(value.Item2));
        internal static (O, O, O) Select<I, O>(this (I, I, I) value, Func<I, O> func) =>
            (func(value.Item1), func(value.Item2), func(value.Item3));
        internal static (O, O, O, O) Select<I, O>(this (I, I, I, I) value, Func<I, O> func) =>
            (func(value.Item1), func(value.Item2), func(value.Item3), func(value.Item4));

        internal static T[] ToArray<T>(this (T, T) value) =>
            new T[] { value.Item1, value.Item2 };
        internal static T[] ToArray<T>(this (T, T, T) value) =>
            new T[] { value.Item1, value.Item2, value.Item3 };
        internal static T[] ToArray<T>(this (T, T, T, T) value) =>
            new T[] { value.Item1, value.Item2, value.Item3, value.Item4 };
    }

    /// <summary>
    ///     Broombridge v0.3 format.
    ///
    ///     Changes from v0.2:
    ///     <list type="bullet">
    ///         <item>Addition of new symmetry key in HamiltonianData.</item>
    ///         <item>Sparse-format arrays are now represented with separate keys and values to simplify parsing logic.</item>
    ///     </list>
    /// </summary>
    #region Broombridge v0.3 format
    public static class V0_3
    {
        // TODO: This URL is not yet valid!
        public static string SchemaUrl = "https://raw.githubusercontent.com/microsoft/Quantum/main/Chemistry/Schema/broombridge-0.3.schema.json";

        // Root of Broombridge data structure
        public struct Data
        {
            public static readonly Format DefaultFormat = new Broombridge.V0_1.Format
            {
                Version = "0.3"
            };

            [YamlMember(Alias = "$schema", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "$schema")]
            public string Schema { get; set; }

            [YamlMember(Alias = "format", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "format")]
            public Format Format { get; set; }

            [YamlMember(Alias = "generator", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "generator")]
            public Generator Generator { get; set; }

            [YamlMember(Alias = "bibliography", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "bibliography")]
            public List<BibliographyItem> Bibliography { get; set; }

            [YamlMember(Alias = "problem_description", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "problem_description")]
            public List<ProblemDescription> ProblemDescriptions { get; set; }

        }

        public struct ProblemDescription
        {
            [YamlMember(Alias = "metadata", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "metadata")]
            public Dictionary<string, object> Metadata { get; set; }

            [YamlMember(Alias = "basis_set", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "basis_set")]
            public BasisSet? BasisSet { get; set; }

            [YamlMember(Alias = "geometry", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "geometry")]
            public Geometry? Geometry { get; set; }

            [YamlMember(Alias = "coulomb_repulsion", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "coulomb_repulsion")]
            public SimpleQuantity CoulombRepulsion { get; set; }

            [YamlMember(Alias = "scf_energy", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "scf_energy")]
            public SimpleQuantity? ScfEnergy { get; set; }

            [YamlMember(Alias = "scf_energy_offset", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "scf_energy_offset")]
            public SimpleQuantity? ScfEnergyOffset { get; set; }

            [YamlMember(Alias = "fci_energy", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "fci_energy")]
            public BoundedQuantity? FciEnergy { get; set; }

            [YamlMember(Alias = "n_orbitals", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "n_orbitals")]
            public int NOrbitals { get; set; }

            [YamlMember(Alias = "n_electrons", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "n_electrons")]
            public int NElectrons { get; set; }

            [YamlMember(Alias = "energy_offset", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "energy_offset")]
            public SimpleQuantity EnergyOffset { get; set; }

            [YamlMember(Alias = "hamiltonian", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "hamiltonian")]
            public HamiltonianData Hamiltonian { get; set; }

            [YamlMember(Alias = "initial_state_suggestions", ApplyNamingConventions = false)]
            [JsonProperty(PropertyName = "initial_state_suggestions")]
            public List<State>? InitialStates { get; set; }
        }

        public struct HamiltonianData
        {
            [YamlMember(Alias = "particle_hole_representation", ApplyNamingConventions = false)]
            // TODO: Placeholder object for ParticleHoleRepresentation, which we do not
            // yet support.
            // NB: Array quantity no longer implicitly adds [], so the type declaration is different
            //     from in 0.1.
            public ArrayQuantity<object[], object>? ParticleHoleRepresentation { get; set; }

            [YamlMember(
                Alias = "one_electron_integrals",
                ApplyNamingConventions = false,
                // Don't allow subclasses here.
                SerializeAs = typeof(ArrayQuantity<(long, long), double>)
            )]
            public ArrayQuantity<(long, long), double> OneElectronIntegrals { get; set; }

            [YamlMember(Alias = "two_electron_integrals", ApplyNamingConventions = false)]
            public ArrayQuantityWithSymmetry<(long, long, long, long), double> TwoElectronIntegrals { get; set; }

        }

        // TODO: move out of v0.3.
        // TODO: finish.
        public struct StringEnum<T> : IYamlConvertible
        where T: struct, Enum
        {
            public T Value;
            public StringEnum(T value)
            {
                Value = value;
            }
            private static ImmutableDictionary<string, T> EnumValues;

            static StringEnum()
            {
                EnumValues = typeof(T).GetMembers()
                    .Select(m => new KeyValuePair<string, MemberInfo>(
                        m
                            .GetCustomAttributes<EnumMemberAttribute>(true)
                            .Select(ema => ema.Value)
                            .FirstOrDefault(),
                        m
                    ))
                    .Where(pa => !string.IsNullOrEmpty(pa.Key))
                    .ToImmutableDictionary(pa => pa.Key, pa => Enum.Parse<T>(pa.Value.Name));
            }

            public void Read(IParser parser, Type expectedType, ObjectDeserializer nestedObjectDeserializer)
            {
                var field = (string)nestedObjectDeserializer(typeof(string));
                Value = EnumValues[field];
            }

            public void Write(IEmitter emitter, ObjectSerializer nestedObjectSerializer)
            {
                var value = Value;
                var name = EnumValues.Where(pair => pair.Value.Equals(value)).Single().Key;
                nestedObjectSerializer(name);
            }

            public static implicit operator T(StringEnum<T> e) => e.Value;
            public static implicit operator StringEnum<T>(T e) => new StringEnum<T>(e);
        }

        public enum ArrayFormat
        {
            [EnumMember(Value = "sparse")]
            Sparse
        }

        public class ArrayQuantity<TIndex, TValue> : HasUnits
        {
            public struct Item : IYamlConvertible
            {
                [YamlMember(Alias = "key", ApplyNamingConventions = false)]
                public TIndex Key { get; set; }

                [YamlMember(Alias = "value", ApplyNamingConventions = false)]
                public TValue Value { get; set; }

                public void Read(IParser parser, Type expectedType, ObjectDeserializer nestedObjectDeserializer)
                {
                    parser.Consume<MappingStart>();
                    var readKey = false;
                    var readValue = false;
                    while (!readKey || !readValue)
                    {
                        var name = nestedObjectDeserializer(typeof(string));
                        switch (name)
                        {
                            case "key":
                                this.Key = ReadKey(parser, nestedObjectDeserializer);
                                readKey = true;
                                break;

                            case "value":
                                this.Value = (TValue)nestedObjectDeserializer(typeof(TValue));
                                readValue = true;
                                break;
                        };
                    }
                    parser.Consume<MappingEnd>();
                }

                private TIndex ReadKey(IParser parser, ObjectDeserializer nestedObjectDeserializer)
                {
                    if (typeof(TIndex).FullName.StartsWith("System.ValueTuple`"))
                    {
                        // TODO [perf]: Cache the create method.
                        var createMethod = typeof(ValueTuple).GetMethods(
                                BindingFlags.Static | BindingFlags.Public
                            )
                            .Where(meth => meth.Name == "Create")
                            .Where(meth => meth.GetParameters().Length == typeof(TIndex).GenericTypeArguments.Length)
                            .Single();
                        var args = new List<object>();
                        parser.Consume<SequenceStart>();
                        foreach (var elementType in typeof(TIndex).GenericTypeArguments)
                        {
                            var element = nestedObjectDeserializer(elementType);
                            args.Add(element);
                        }
                        parser.Consume<SequenceEnd>();
                        return (TIndex)createMethod.MakeGenericMethod(typeof(TIndex).GenericTypeArguments).Invoke(null, args.ToArray());
                    }
                    else return (TIndex)nestedObjectDeserializer(typeof(TIndex));
                }

                public void Write(IEmitter emitter, ObjectSerializer nestedObjectSerializer)
                {
                    emitter.Emit(new MappingStart(null, null, true, MappingStyle.Flow));
                    nestedObjectSerializer("key");
                    if (typeof(TIndex).FullName.StartsWith("System.ValueTuple`"))
                    {
                        emitter.Emit(new SequenceStart(null, null, true, SequenceStyle.Flow));
                        foreach (var (_, idx) in typeof(TIndex).GenericTypeArguments.Select((_, idx) => (_, idx)))
                        {
                            var field = typeof(TIndex).GetField($"Item{idx + 1}");
                            nestedObjectSerializer(field.GetValue(Key));
                        }
                        emitter.Emit(new SequenceEnd());
                    }
                    else
                    {
                        nestedObjectSerializer(Key);
                    }
                    nestedObjectSerializer("value");
                    nestedObjectSerializer(Value);
                    emitter.Emit(new MappingEnd());
                }
            }

            [YamlMember(Alias = "format")]
            public StringEnum<ArrayFormat> Format { get; set; }
            
            [YamlMember(Alias = "values")]
            public List<Item> Values { get; set; }

            [YamlMember(Alias = "index_convention")]
            public OrbitalIntegral.Convention? IndexConvention { get; set; } = null;
        }

        public enum PermutationSymmetry
        {
            [EnumMember(Value = "eightfold")]
            Eightfold,
            [EnumMember(Value = "fourfold")]
            Fourfold,
            [EnumMember(Value = "trivial")]
            Trivial
        }

        public struct Symmetry
        {
            [YamlMember(Alias = "permutation")]
            public StringEnum<PermutationSymmetry> Permutation { get; set; }
        }

        public class ArrayQuantityWithSymmetry<TIndex, TValue> : ArrayQuantity<TIndex, TValue>
        {
            [YamlMember(Alias = "symmetry")]
            public Symmetry Symmetry { get; set; }
        }

        /// <summary>
        /// Builds Hamiltonian from Broombridge if data is available.
        /// </summary>
        internal static OrbitalIntegralHamiltonian ToOrbitalIntegralHamiltonian(ProblemDescription broombridge)
        {
            // Add the identity terms
            var identityterm = broombridge.CoulombRepulsion.Value + broombridge.EnergyOffset.Value;
            var hamiltonian =  new OrbitalIntegralHamiltonian();
            var hamiltonianData = broombridge.Hamiltonian;

            // This will convert from Broombridge 1-indexing to 0-indexing.
            hamiltonian.Add
                (hamiltonianData.OneElectronIntegrals.Values
                .Select(o =>
                    new OrbitalIntegral(
                        o.Key.Select(k => (int)(k - 1)).ToArray(),
                        o.Value,
                        OrbitalIntegral.Convention.Mulliken
                    )
                .ToCanonicalForm(PermutationSymmetry.Eightfold.FromBroombridgeV0_3()))
                .Distinct());

            // This will convert from Broombridge 1-indexing to 0-indexing.
            // This will convert to Dirac-indexing.
            hamiltonian.Add
                (hamiltonianData.TwoElectronIntegrals.Values
                .Select(o =>
                    new OrbitalIntegral(
                        o.Key.Select(k => (int)(k - 1)).ToArray(),
                        o.Value,
                        OrbitalIntegral.Convention.Mulliken
                    )
                .ToCanonicalForm(hamiltonianData.TwoElectronIntegrals.Symmetry.Permutation.Value.FromBroombridgeV0_3()))
                .Distinct());

            hamiltonian.Add(new OrbitalIntegral(), identityterm);
            return hamiltonian;
        }

    }
    #endregion
    
}
