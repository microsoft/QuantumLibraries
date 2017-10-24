using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using Xunit;
using Xunit.Abstractions;
using Xunit.Sdk;
using Microsoft.Quantum.Simulation.Core;

namespace Microsoft.Quantum.Testing
{
    public class QuantumTestCaseException : Exception
    {
        public QuantumTestCaseException(string message) : base(message)
        {
        }
    }

    public class QuantumTestCaseOperationData : IXunitSerializable
    {
        public string assemblyName;
        public string className;
        public string fullClassName;
        public string skipReason;

        /// <summary>
        /// Synchronous operation runner for the simulator. This way we can actually catch expections without 
        /// them being wrapped into AggregateException. 
        /// </summary>
        public Action<IOperationFactory> OperationRunner
        {
            get
            {
                return (IOperationFactory sim) =>
                {
                    if (skipReason == null)
                    {
                        System.Reflection.Assembly assembly = System.Reflection.Assembly.Load(assemblyName);
                        Type operationType = assembly.GetType(fullClassName);
                        ICallable<QVoid, QVoid> op = (typeof(IOperationFactory)
                               .GetMethod("Get")
                               .MakeGenericMethod(typeof(ICallable<QVoid, QVoid>), operationType)
                               .Invoke(sim, null)
                           as ICallable<QVoid, QVoid>);
                        if (op == null)
                        {
                            throw new QuantumTestCaseException($"Operation with name: {fullClassName} was not found in assembly {assemblyName}");
                        }
                        op.Apply(QVoid.Instance);
                    }
                };
            }
        }

        public void Deserialize(IXunitSerializationInfo info)
        {
            assemblyName = info.GetValue<string>("assemblyName");
            className = info.GetValue<string>("className");
            fullClassName = info.GetValue<string>("fullClassName");
            skipReason = info.GetValue<string>("skipReason");
        }

        public void Serialize(IXunitSerializationInfo info)
        {
            info.AddValue("assemblyName", assemblyName);
            info.AddValue("className", className);
            info.AddValue("fullClassName", fullClassName);
            info.AddValue("skipReason", skipReason);
        }
    }

    [XunitTestCaseDiscoverer("Microsoft.Quantum.Testing.QuantumTestTargetDiscoverer", "Microsoft.Quantum.Testing")]
    public class QuantumTestTargetAttribute : FactAttribute
    {
        /// <summary>
        /// Suffix of the operation name that signifies that given operation is a test. Defaults to Test.
        /// </summary>
        public string Suffix { get; set; } = "Test";

        /// <summary>
        /// Assembly where look for the test cases. Defaults to the assembly in which test method is defined.
        /// </summary>
        public string AssemblyName { get; set; }

        /// <summary>
        /// Namespace where we look for the test cases. Defaults to the namespace in which test method is defined.
        /// </summary>
        public string TestNamespace { get; set; }
    }

    public class QuantumTestTargetDiscoverer : IXunitTestCaseDiscoverer
    {
        private QuantumTestCase SkipErrorTestCase(ITestFrameworkDiscoveryOptions discoveryOptions, ITestMethod testMethod, string message )
        {
            return new QuantumTestCase(
                DiagnosticMessageSink,
                discoveryOptions.MethodDisplay() ?? TestMethodDisplay.Method,
                testMethod,
                new QuantumTestCaseOperationData()
                {
                    assemblyName = "",
                    className = "",
                    fullClassName = "",
                    skipReason = message
                });
        }

        public IEnumerable<IXunitTestCase> Discover(ITestFrameworkDiscoveryOptions discoveryOptions, ITestMethod testMethod, IAttributeInfo factAttribute)
        {
            List<IXunitTestCase> result = new List<IXunitTestCase>();

            string assemblyName = factAttribute.GetNamedArgument<string>("AssemblyName");
            System.Reflection.Assembly testAssembly = null;
            if (assemblyName == null)
            {
                testAssembly = testMethod.TestClass.Class.ToRuntimeType().Assembly;
            }
            else
            {
                try
                {
                    testAssembly = System.Reflection.Assembly.Load(assemblyName);
                }
                catch( Exception e )
                {
                    result.Add(SkipErrorTestCase(discoveryOptions, testMethod, "Assembly issue: " + e.Message));
                    return result;
                }
            }

            string testNamespace = factAttribute.GetNamedArgument<string>("TestNamespace") ?? testMethod.TestClass.Class.ToRuntimeType().Namespace;
            if (testNamespace == null)
            {
                result.Add(SkipErrorTestCase(discoveryOptions, testMethod, "Could not find namespase with tests"));
                return result;
            }

            string Suffix = factAttribute.GetNamedArgument<string>("Suffix");
            string skipReason = factAttribute.GetNamedArgument<string>("Skip");

            IEnumerable<Type> ourTypes =
                from definedType in testAssembly.DefinedTypes
                where definedType.Name.EndsWith(Suffix)
                where definedType.Namespace == testNamespace
                where typeof(Operation).IsAssignableFrom(definedType)
                where typeof(ICallable<QVoid, QVoid>).IsAssignableFrom(definedType)
                select definedType;

            foreach (Type operationType in ourTypes)
            {
                result.Add(
                    new QuantumTestCase(DiagnosticMessageSink, discoveryOptions.MethodDisplay() ?? TestMethodDisplay.Method, testMethod,
                        new QuantumTestCaseOperationData()
                        {
                            assemblyName = testAssembly.FullName,
                            className = operationType.Name,
                            fullClassName = operationType.FullName,
                            skipReason = skipReason
                        })
                    );
            }

            if (result.Count == 0)
            {
                result.Add(SkipErrorTestCase(discoveryOptions, testMethod, "No tests found."));
            }

            return result;
        }

        /// <summary>
        /// Gets the message sink to be used to send diagnostic messages.
        /// </summary>
        protected IMessageSink DiagnosticMessageSink { get; }

        /// <summary>
        /// Initializes a new instance of the <see cref="QuantumFactDiscoverer"/> class.
        /// </summary>
        /// <param name="diagnosticMessageSink">The message sink used to send diagnostic messages</param>
        public QuantumTestTargetDiscoverer(IMessageSink diagnosticMessageSink)
        {
            DiagnosticMessageSink = diagnosticMessageSink;
        }
    }

    [Serializable]
    public class QuantumTestCase : XunitTestCase
    {
        QuantumTestCaseOperationData operationData;
        /// <summary/>
        [EditorBrowsable(EditorBrowsableState.Never)]
        [Obsolete("Called by the de-serializer; should only be called by deriving classes for de-serialization purposes")]
        public QuantumTestCase() : base()
        {
        }

        public QuantumTestCase(IMessageSink messageSink, TestMethodDisplay defaultTestMethodDisplay, ITestMethod testMethod, QuantumTestCaseOperationData operationData)
        : base(messageSink, defaultTestMethodDisplay, testMethod, new[] { operationData })
        {
            this.operationData = operationData;
        }

        protected override string GetDisplayName(IAttributeInfo factAttribute, string displayName)
        {
            if (operationData.className == "")
            {
                return displayName;
            }
            else
            {
                return operationData.className;
            }
        }

        protected override string GetSkipReason(IAttributeInfo factAttribute)
        {
            return operationData.skipReason;
        }

        public override void Serialize(IXunitSerializationInfo data)
        {
            data.AddValue("operationData", operationData);
            base.Serialize(data);
        }

        public override void Deserialize(IXunitSerializationInfo data)
        {
            base.Deserialize(data);
            operationData = data.GetValue<QuantumTestCaseOperationData>("operationData");
            this.TestMethodArguments = new[] { operationData };
        }
    }
}
