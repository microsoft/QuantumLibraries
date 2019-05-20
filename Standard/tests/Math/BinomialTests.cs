using Xunit;

namespace Microsoft.Quantum.Standard.Emulation
{
    public class BinomialTests
    {
        [Fact]
        public void TestExtrema()
        {
            var dist = new BinomialDistribution(100, 0.0);
            Assert.Equal(0, dist.NextSample());
            dist = new BinomialDistribution(100, 1.0);
            Assert.Equal(100, dist.NextSample());
        }
    }
}
