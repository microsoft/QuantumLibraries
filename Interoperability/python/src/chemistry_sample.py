import logging
logging.basicConfig(level=logging.INFO)

import qsharp.chemistry
from qsharp.chemistry import load_broombridge, load_fermion_hamiltonian, load_input_state, encode

# Load a fermion hamiltonian:
fh1 = load_fermion_hamiltonian("broombridge.yaml")
logging.info("fh1 ready.")
logging.debug(fh1)

# optionally, load first the broombridge data, and from there
# load the fermion hamiltonian from a problem description:
broombridge = load_broombridge("broombridge.yaml")
fh2 = broombridge.problem_description[0].load_fermion_hamiltonian()
logging.info("fh2 ready.")
logging.debug(fh2.terms)

# Here I show how to add a couple of completely made-up terms to the fh2 Hamiltonian:
terms = [ ([], 10.0), ([0,6], 1.0), ([0,2,2,0], 1.0)]
fh2.add_terms(terms)
logging.info("terms added successfully")
logging.debug(fh2)



# Simularly, you can load an input state either directly from a broombridge.yaml file,
is1 = load_input_state("broombridge.yaml", "|E1>")
logging.info("is1 ready.")
logging.debug(is1)

# or from the problem description:
is2 = broombridge.problem_description[0].load_input_state("|E1>")
logging.info("is2 ready.")
logging.debug(is2)

# Once we have the hamiltonian and input state, we can call 'encode' to generate a Jordan-Wigner
# representation, suitable for quantum simulation:
qsharp_encoding = encode(fh1, is1)

# Compile a Q# function on the fly as an entry point for our simulation.
# For now, it is better if you compile the entry-point operation from Python (instead of using a .qs file):
TrotterEstimateEnergy = qsharp.compile("""
    open Microsoft.Quantum.Chemistry.JordanWigner;
    
    operation TrotterEstimateEnergy (qSharpData: JordanWignerEncodingData, nBitsPrecision : Int, trotterStepSize : Double) : (Double, Double) {
        
        let (nSpinOrbitals, data, statePrepData, energyShift) = qSharpData!;
        
        // Order of integrator
        let trotterOrder = 1;
        let (nQubits, (rescaleFactor, oracle)) = TrotterStepOracle(qSharpData, trotterStepSize, trotterOrder);
        
        // Prepare ProductState
        let statePrep =  PrepareTrialState(statePrepData, _);
        let phaseEstAlgorithm = RobustPhaseEstimation(nBitsPrecision, _, _);
        let estPhase = EstimateEnergy(nQubits, statePrep, oracle, phaseEstAlgorithm);
        let estEnergy = estPhase * rescaleFactor + energyShift;
        return (estPhase, estEnergy);
    }
""")


logging.info('Starting simulation.')
result = TrotterEstimateEnergy.simulate(qSharpData=qsharp_encoding, nBitsPrecision=10, trotterStepSize=.4)
print(f'Trotter simulation complete. (phase, energy): {result}')
