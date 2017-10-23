namespace Microsoft.Quantum.Canon {

    // FIXME: the names of these functions are not compliant with the style guide.
    // FIXME: define recursively instead of special casing to first, second, and fourth.
    //TODO revert to generics

    operation Trotter1ImplCA<'T>(evolutionGenerator : (Int, ((Int, Double, Qubit[]) => () : Adjoint, Controlled)), stepSize : Double, target : Qubit[]) : () {
        Body {
            let (nSteps, op) = evolutionGenerator
            for(idx in 0..nSteps-1){
                op(idx, stepSize, target)
            }
        }
        Adjoint auto
        Controlled auto
        Controlled Adjoint auto
    }

    operation Trotter2ImplCA<'T>(evolutionGenerator : (Int, ((Int, Double, Qubit[]) => () : Adjoint, Controlled)), stepSize : Double, target : Qubit[]) : () {
        Body {
            let (nSteps, op) = evolutionGenerator
            for(idx in 0..nSteps-1){
                op(idx, stepSize * 0.5, target)
            }
            for(idx in (nSteps-1)..-1..0){
                op(idx, stepSize * 0.5, target)
            }
        }
        Adjoint auto
        Controlled auto
        Controlled Adjoint auto
    }

    function DecomposeIntoTimeStepsCA<'T>(evolutionGenerator : (Int, ((Int, Double, Qubit[]) => () : Adjoint, Controlled)), order : Int) : ((Double, Qubit[]) => () : Adjoint, Controlled) {
        if (order == 1) {
            return Trotter1ImplCA(evolutionGenerator, _, _)
        } elif (order == 2) {
            return Trotter2ImplCA(evolutionGenerator, _, _)
        } else {
            fail "Order $order not yet supported."
        }

        // FIXME: needed so we have a return value of the right type in all cases, but
        //        this line is unreachable and should be removed.
        return Trotter1ImplCA(evolutionGenerator, _, _)
    }

    // TODO: write variants of the above that do not assume adjoint and controlled.

}