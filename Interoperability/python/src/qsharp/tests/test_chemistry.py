import logging
logging.basicConfig(level=logging.INFO)

import qsharp.chemistry
from qsharp.chemistry import load_broombridge, load_fermion_hamiltonian, load_input_state, encode

print ( qsharp.component_versions() )

def test_load_broombridge():
    """
    Checks that we can load a broombridget schema file.
    """
    broombridge = load_broombridge("broombridge.yaml")
    assert(len(broombridge.problem_description) == 1)
    assert(len(broombridge.bibliography) == 3)


def test_load_fermion_hamiltonian():
    """
    Checks that loading a fermion hamiltonian from file or from broombridget gives the same result.
    """
    broombridge = load_broombridge("broombridge.yaml")

    fh1 = broombridge.problem_description[0].load_fermion_hamiltonian()
    fh2 = load_fermion_hamiltonian("broombridge.yaml")

    assert(len(fh1.terms) == 6)
    assert(fh1 == fh2)
    
    fh3 = broombridge.problem_description[0].load_fermion_hamiltonian("HalfUp")
    fh4 = load_fermion_hamiltonian("broombridge.yaml", "HalfUp")

    assert(len(fh3.terms) == 6)
    assert(fh3 == fh4)
    
    # should this be true? assert(fh3 != fh1)


def test_load_input_state():
    """
    Checks that loading an input state from file or from broombridget gives the same result.
    """
    broombridge = qsharp.chemistry.load_broombridge("broombridge.yaml")

    is1 = broombridge.problem_description[0].load_input_state("UCCSD |G>")
    is2 = load_input_state("broombridge.yaml", "UCCSD |G>")

    assert(is1.label == "UCCSD |G>")
    assert(is1 == is2)
    
    is3 = broombridge.problem_description[0].load_input_state("UCCSD |G>", "HalfUp")
    is4 = load_input_state("broombridge.yaml", "UCCSD |G>", "HalfUp")

    assert(is3.label == "UCCSD |G>")
    assert(is3 == is4)
    assert(is1 != is4)


def test_load_greedy_state():
    """
    Checks that not passing a wavefunction label generates the greedy input state.
    """
    broombridge = qsharp.chemistry.load_broombridge("broombridge.yaml")

    is1 = broombridge.problem_description[0].load_input_state("", "HalfUp")
    is2 = load_input_state("broombridge.yaml", "", "HalfUp")

    assert(is1.label == "Greedy")
    assert(is1 == is2)

    
def test_jw_encode():
    """
    Checks that we can encode a hamiltonian + input state
    """
    broombridge = qsharp.chemistry.load_broombridge("broombridge.yaml")

    fh1 = broombridge.problem_description[0].load_fermion_hamiltonian()
    is1 = broombridge.problem_description[0].load_input_state()

    jw = encode(fh1, is1)
    assert(len(jw.data) == 4)