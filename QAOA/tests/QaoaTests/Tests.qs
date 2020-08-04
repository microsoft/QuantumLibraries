namespace Quantum.QaoaTests {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.QAOA;

    
    @Test("QuantumSimulator")
    operation AllocateQubit () : Unit {

        using (q = Qubit()) {
            Assert([PauliZ], [q], Zero, "Newly allocated qubit must be in |0> state.");
        }
        
        Message("Test passed.");
    }
}
