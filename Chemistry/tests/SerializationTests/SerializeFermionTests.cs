using System;
using Xunit;
using Newtonsoft.Json;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;

using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Pauli;
using Microsoft.Quantum.Chemistry.QSharpFormat;
using Microsoft.Quantum.Chemistry.Generic;

using Microsoft.Quantum.Chemistry;

using Newtonsoft.Json.Linq;

namespace SerializationTests
{
    public class SerializeFermionTests
    {
        [Fact]

        public void SerializeHermitianFermionTerm()
        {
            HermitianFermionTerm term0 = new HermitianFermionTerm(new[] { 1, 2, 3, 4 });

            string json = JsonConvert.SerializeObject(term0, Formatting.None);

            HermitianFermionTerm term1 = new HermitianFermionTerm( JsonConvert.DeserializeObject<LadderSequence<int>>(json));

            Debug.WriteLine(@json);
        }

        // Lithium Hydride filename.
        static string filename = "LiH_0.1.yaml";

        [Fact]
        public void SerializeEmptyHamiltonian()
        {
            FermionHamiltonian hamiltonian = new FermionHamiltonian();

            string json = JsonConvert.SerializeObject(hamiltonian, Formatting.Indented);

            Console.WriteLine(json);

            FermionHamiltonian deserializedHamiltonian = JsonConvert.DeserializeObject<FermionHamiltonian>(json);
        }

        [Fact]
        public void SerializeLithiumHydrideHamiltonian()
        {
            var broombridge = Deserializers.DeserializeBroombridge(filename).ProblemDescriptions.First();

            var orbHam = broombridge.OrbitalIntegralHamiltonian;

            var ferHam = orbHam.ToFermionHamiltonian(IndexConvention.UpDown);

            string json = JsonConvert.SerializeObject(ferHam);

            Debug.WriteLine(@json);

            var deserializedHamiltonian = JsonConvert.DeserializeObject<FermionHamiltonian>(json);
        }
    }

    

    public class ExampleSerializationTest
    {
        [Fact]
        public void Test1()
        {
            User user = new User
            {
                UserName = @"domain\username"
            };

            string json = JsonConvert.SerializeObject(user, Formatting.Indented);

            Console.WriteLine(json);
        }
    }

    public class UserConverter : JsonConverter
    {
        public override void WriteJson(JsonWriter writer, object value, JsonSerializer serializer)
        {
            User user = (User)value;

            writer.WriteValue(user.UserName);
        }

        public override object ReadJson(JsonReader reader, Type objectType, object existingValue, JsonSerializer serializer)
        {
            User user = new User();
            user.UserName = (string)reader.Value;

            return user;
        }

        public override bool CanConvert(Type objectType)
        {
            return objectType == typeof(User);
        }
    }

    [JsonConverter(typeof(UserConverter))]
    public class User
    {
        public string UserName { get; set; }
    }
}
