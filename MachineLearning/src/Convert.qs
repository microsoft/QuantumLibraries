// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;

    function unFlattenSchedule(sc : Int[][]) : SamplingSchedule
    {
        mutable ret = new Range[0];
        for (flattenedRange in sc) {
            set ret += [flattenedRange[0]..flattenedRange[1]..flattenedRange[2]];
        }
        return SamplingSchedule(ret);
    }

    function unFlattenLabeledSamples(dat:Double[][], labs:Int[]) : LabeledSample[] {
        mutable cnt = MinI(Length(dat), Length(labs));
        mutable ret = new LabeledSample[cnt];
        for (j in 0..(cnt - 1)) {
            set ret w/= j <- LabeledSample(dat[j], labs[j]);
        }
        return ret;
    }

    /// Debugging prop
    operation unFlattenPauli(p:Int): Pauli
    {
        if (p==1)
        {
            return PauliX;
        }
        if (p==2)
        {
            return PauliY;
        }
        if (p==3)
        {
            return PauliZ;
        }
        return PauliI;
    }

    /// Debugging prop
    /// upcasting controlled rotation in flat representation (paramIx,pauliIx,gateSpan)
    operation unFlattenControlledRotation(cod:Int[]): ControlledRotation {
        return ControlledRotation(
            GateSpan(
                cod[2], cod[3...]
            ),
            unFlattenPauli(cod[1]),
            cod[0]
        );
    }

    /// Debugging prop
    operation unFlattenGateSequence(seq: Int[][]) : GateSequence {
        mutable tmp = new ControlledRotation[Length(seq)];
        for (icr in 0..(Length(seq) - 1)) {
            set tmp w/= icr <- unFlattenControlledRotation(seq[icr]);
        }
        return GateSequence(tmp);
    }

}