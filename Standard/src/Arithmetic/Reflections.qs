namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

	///Flip the sign of just one amplitude
	operation ReflectAboutInteger(index : Int, reg : LittleEndian): Unit is Adj + Ctl {
        let nQubits = Length(reg!);
        let bitstring = IntAsBoolArray(index, nQubits);
        if (nQubits < 2) {
            within {
                if (not bitstring[0]) {
                    X(reg![0]);
                }
            } apply {
                Z(reg![0]);
            }
        } else {
            within {
                ApplyToEachCA(CControlledCA(X), Zip(bitstring, reg!));
            } apply {
                (Controlled Z)(Most(reg!), Tail(reg!)); //The true complexity of this operation is in O(nQubits)
            }
        }
	} //_amplitudeSignFlip


}
