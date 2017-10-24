using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace Microsoft.Quantum.Canon.Tests
{
    public class TestBase
    {

        public static IOperationFactory GetSimulator()
        {
            var sim = new QuantumSimulator();
            RegisterNativeCanonOperations(sim);
            return sim;
        }

        private static Dictionary<String, TypeInfo> GetOperationTypesByNamespace(String namespaceName) {
            var matchingTypes = (
                from definedType in
                    Assembly
                    .GetExecutingAssembly()
                    .DefinedTypes
                where definedType.Namespace == namespaceName
                where typeof(Operation).IsAssignableFrom(definedType)
                select definedType
            );
            return matchingTypes.ToDictionary(
                definedType => definedType.Name,
                definedType => definedType
            );
        }

        private static Dictionary<TKey, (TValueLeft, TValueRight)>
            DictionaryIntersection<TKey, TValueLeft, TValueRight>(
                Dictionary<TKey, TValueLeft> left,
                Dictionary<TKey, TValueRight> right
            ) => (
                    from key in left.Keys
                    where right.Keys.Contains(key)
                    select key
                ).ToDictionary(
                    key => key,
                    key => (left[key], right[key])
                );

        private static void RegisterNativeCanonOperations(AbstractFactory<Operation> simulator)
        {
            // Find anything in the Native namespace.
            var nativeOperationTypes = GetOperationTypesByNamespace("Microsoft.Quantum.Canon.Native");
            var possibleStubOperationTypes = GetOperationTypesByNamespace("Microsoft.Quantum.Canon");

            // If a name is in both, then it's a stub.
            var stubOperationTypes = DictionaryIntersection(
                possibleStubOperationTypes,
                nativeOperationTypes
            );

            foreach (var replacement in stubOperationTypes.Values)
            {
                simulator.Register(replacement.Item1, replacement.Item2);
            }
        }
    }
}
