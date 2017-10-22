using Microsoft.Quantum.Simulation.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Xunit;

namespace Microsoft.Quantum.Canon.Tests
{
    public class NativeOperationTests : TestBase
    {
        [Theory(DisplayName ="Native trig")]
        [RandomDoubleData(nSamples: 5, min: -1, max: 1)]
        public void TestNativeTrig(double theta)
        {
            var sim = GetSimulator();

            var sin = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.Sin>();
            Assert.Equal(sin.Apply(theta), Math.Sin(theta));

            var cos = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.Cos>();
            Assert.Equal(cos.Apply(theta), Math.Cos(theta));

            var asin = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.ArcSin>();
            Assert.Equal(asin.Apply(theta), Math.Asin(theta));

            var acos = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.ArcCos>();
            Assert.Equal(acos.Apply(theta), Math.Acos(theta));
        }
    }
}
