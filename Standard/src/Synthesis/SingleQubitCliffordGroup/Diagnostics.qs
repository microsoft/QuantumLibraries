namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    // TODO
    function IdentityFact1C(op : SingleQubitClifford, message : String) : Unit {
        let c = CanonicalForm1C(op);
        if c::Omega != 0 or c::E != 0 or c::X != 0 or c::S != 0 {
            FormattedFailure(op, Identity1C(), "Single-qubit Clifford operator was not identity.");
        }
    }

}
