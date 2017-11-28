namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Samples.BitFlipCode;

    operation BitFlipSampleParityTest() : ()  {
        body {
            CheckBitFlipCodeStateParity();
        }
    }

    operation BitFlipSampleWt1CorrectionTest() : ()  {
        body {
            CheckBitFlipCodeCorrectsBitFlipErrors();
        }
    }

    operation BitFlipSampleWCanonTest() : () {
        body {
            CheckCanonBitFlipCodeCorrectsBitFlipErrors();
        }
    }
}
