// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    // FIXME: remove once generics code gen.
    // NB: copied from AmplitudeAmplification.qb as a temporary measure.
    function QubitSlice(indices : Int[], qubits : Qubit[]) : Qubit[] 
    {
        let nSliced = Length(indices);
        mutable sliced = new Qubit[nSliced];
        for( idx in 0..nSliced - 1 ) {
            set sliced[idx] = qubits[indices[idx]];
        }

        return sliced;
    }

	//QubitExclude([idxQubitFlag], qubits)
	//TODO replace with Exclude in generics
	function QubitExclude(remove : Int[], qubits : Qubit[]) : Qubit[] 
	{

		let nSliced = Length(remove);
		let nQubits = Length(qubits);
		//Would be better with sort function
		//Or way to add elements to array

		mutable arrayKeep = new Int[nQubits];
		mutable sliced = new Qubit[nQubits - nSliced];
		mutable counter = 0;

		for ( idx in 0..nQubits - 1) {
			set arrayKeep[idx] = idx;
		}
		for ( idx in 0..nSliced - 1 ) {
			set arrayKeep[remove[idx]] = -1;
		}
		for ( idx in 0..nQubits - 1 ) {
			if(arrayKeep[idx] >= 0){
				set sliced[counter] = qubits[arrayKeep[idx]];
				set counter = counter + 1;
			}
		}

		return sliced;
	}

    operation ApplyToSubregister(op : (Qubit[] => ()), idxs : Int[], target : Qubit[]) : () {
        body {
            // FIXME: change to Slice<Qubit[]> once generics code gen.
            let subregister = QubitSlice(idxs, target);
            op(subregister);
        }
    }

    operation ApplyToSubregisterA(op : (Qubit[] => () : Adjoint), idxs : Int[], target : Qubit[]) : () {
        body {
            ApplyToSubregister(op, idxs, target);
        }
        adjoint {
            ApplyToSubregister(Adjoint op, idxs, target);
        }
    }

    operation ApplyToSubregisterC(op : (Qubit[] => () : Controlled), idxs : Int[], target : Qubit[]) : () {
        body {
            ApplyToSubregister(op, idxs, target);
        }
        controlled (controls) {
            let cop = (Controlled op);
            ApplyToSubregister(cop(controls, _), idxs, target);
        }
    }

    operation ApplyToSubregisterCA(op : (Qubit[] => () : Controlled, Adjoint), idxs : Int[], target : Qubit[]) : () {
        body {
            ApplyToSubregister(op, idxs, target);
        }
        adjoint {
            ApplyToSubregister(Adjoint op, idxs, target);
        }
        controlled (controls) {
            let cop = (Controlled op);
            ApplyToSubregister(cop(controls, _), idxs, target);
        }
        controlled adjoint (controls) {
            let cop = (Controlled Adjoint op);
            ApplyToSubregister(cop(controls, _), idxs, target);
        }
    }

    function RestrictToSubregister(op : (Qubit[] => ()), idxs : Int[]) : (Qubit[] => ()) {
        return ApplyToSubregister(op, idxs, _);
    }
    function RestrictToSubregisterA(op : (Qubit[] => () : Adjoint), idxs : Int[]) : (Qubit[] => () : Adjoint) {
        return ApplyToSubregisterA(op, idxs, _);
    }
    function RestrictToSubregisterC(op : (Qubit[] => () : Controlled), idxs : Int[]) : (Qubit[] => () : Controlled) {
        return ApplyToSubregisterC(op, idxs, _);
    }
    function RestrictToSubregisterCA(op : (Qubit[] => () : Adjoint, Controlled), idxs : Int[]) : (Qubit[] => () : Adjoint, Controlled) {
        return ApplyToSubregisterCA(op, idxs, _);
    }

}
