namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;

    /// Sample container access method
    @Deprecated("")
    function getSample(samples: LabeledSampleContainer, ix: Int): LabeledSample {
        return (samples!)[ix];
    }

    /// Access the raw data in a labeled sample
    @Deprecated("")
    function getData(samp: LabeledSample): Double[] {
        return Fst(samp!);
    }

    /// Access the label in a labeled sample
    @Deprecated("")
    function getLabel(samp:LabeledSample) : Int
    {
        return Snd(samp!);
    }


    /// Abstraction for a container of labeled samples
    @Deprecated("")
    newtype LabeledSampleContainer = LabeledSample[];

    @Deprecated("Microsoft.Quantum.Diagnostics.DumpRegister")
    function dumpRegisterToConsole ( qs: Qubit[]) : Unit
    {}
    //{DumpRegister((),qs);} //Swap for empty body when some dumping of registers is needed

    @Deprecated("Microsoft.Quantum.MachineLearning.NQubitsRequired")
    function qubitSpan(seq : GateSequence) : Int {
        return NQubitsRequired(seq);
    }

    /// Set force a qubit into a desired basis state
    @Deprecated("Microsoft.Quantum.Measurement.SetToBasisState")
    operation Set (desired: Result, q1: Qubit) : Unit
    {
        //body
        //{
            let current = M(q1);
            if (desired != current)
            {
                X(q1);
            }
        //}
    }

    @Deprecated("Microsoft.Quantum.Math.SquaredNorm")
    function squareNorm(v:Double[]):Double
    {
        mutable ret = 0.0;
        for (u in v)
        {
            set ret = ret + u*u;
        }
        return ret;
    }

    @Deprecated("") // replace with ForEach.
    operation randomizeArray(src:Double[], relativeFuzz: Double) : Double[]
    {
        mutable ret = new Double[Length(src)];
        for (ix in 0..(Length(src)-1))
        {
            set ret w/=ix <- _RandomlyRescale(src[ix], relativeFuzz);
        }
        return ret;
    }

    @Deprecated("Microsoft.Quantum.Math.NearlyEqualD")
    function nearIdenticalDoubles(x:Double,y:Double):Bool {
        return NearlyEqualD(x, y); //Note key tolerance constant here
    }


}
