namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;

    function FeatureRegisterSize(sample : Double[]) : Int {
        return Ceiling(Lg(IntAsDouble(Length(sample))));
    }

}
