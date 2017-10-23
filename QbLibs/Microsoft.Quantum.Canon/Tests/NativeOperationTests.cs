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
        [Theory(DisplayName ="Native math")]
        [RandomDoubleData(nSamples: 5, min: -1, max: 1)]
        public void TestNativeMath(double theta)
        {
            var sim = GetSimulator();

            var sin = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.Sin>();
            Assert.Equal(sin.Apply(theta), Math.Sin(theta));

            var cos = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.Cos>();
            Assert.Equal(cos.Apply(theta), Math.Cos(theta));

            var tan = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.Tan>();
            Assert.Equal(tan.Apply(theta), Math.Tan(theta));

            var asin = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.ArcSin>();
            Assert.Equal(asin.Apply(theta), Math.Asin(theta));

            var acos = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.ArcCos>();
            Assert.Equal(acos.Apply(theta), Math.Acos(theta));
            
            var atan = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.ArcTan>();
            Assert.Equal(atan.Apply(theta), Math.Atan(theta));

            var sqrt = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.Sqrt>();
            Assert.Equal(sqrt.Apply(theta), Math.Sqrt(theta));

            var log = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.Log>();
            Assert.Equal(log.Apply(theta), Math.Log(theta));

            var sinh = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.Sinh>();
            Assert.Equal(sinh.Apply(theta), Math.Sinh(theta));

            var cosh = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.Cosh>();
            Assert.Equal(cosh.Apply(theta), Math.Cosh(theta));

            var tanh = sim.Get<ICallable<double, double>, Microsoft.Quantum.Canon.Tanh>();
            Assert.Equal(tanh.Apply(theta), Math.Tanh(theta));

            var ceil = sim.Get<ICallable<double, long>, Microsoft.Quantum.Canon.Ceiling>();
            Assert.Equal(ceil.Apply(theta), Math.Ceiling(theta));

            var floor = sim.Get<ICallable<double, long>, Microsoft.Quantum.Canon.Floor>();
            Assert.Equal(floor.Apply(theta), Math.Floor(theta));
        }
    }
}
