import logging
logging.basicConfig(level=logging.INFO)

from qsharp_chemistry import load_broombridge, load_fermion_hamiltonian, load_input_state, encode
from qsharp_chemistry.problem_description import IndexConvention
from qsharp_chemistry.fermion_hamiltonian import InputState

import qsharp

print ( qsharp.component_versions() )

def test_load_broombridge():
    """
    Checks that we can load a broombridge schema file.
    """
    broombridge = load_broombridge("broombridge.yaml")
    assert(len(broombridge.problem_description) == 1)
    assert(len(broombridge.bibliography) == 3)


def test_load_fermion_hamiltonian():
    """
    Checks that loading a fermion hamiltonian from file or from broombridge gives the same result.
    """
    broombridge = load_broombridge("broombridge.yaml")

    fh1 = broombridge.problem_description[0].load_fermion_hamiltonian()
    fh2 = load_fermion_hamiltonian("broombridge.yaml")

    assert(len(fh1.terms) == 6)
    assert(fh1 == fh2)
    
    fh3 = broombridge.problem_description[0].load_fermion_hamiltonian(IndexConvention.HalfUp)
    fh4 = load_fermion_hamiltonian("broombridge.yaml", IndexConvention.HalfUp)

    assert(len(fh3.terms) == 6)
    assert(fh3 == fh4)
    
    # should this be true? assert(fh3 != fh1)


def test_load_input_state():
    """
    Checks that loading an input state from file or from broombridge gives the same result.
    """
    broombridge = load_broombridge("broombridge.yaml")

    is1 = broombridge.problem_description[0].load_input_state("UCCSD |G>")
    is2 = load_input_state("broombridge.yaml", "UCCSD |G>")

    assert(is1.Method == "UnitaryCoupledCluster")
    assert(is1 == is2)
    
    is3 = broombridge.problem_description[0].load_input_state("UCCSD |G>", IndexConvention.HalfUp)
    is4 = load_input_state("broombridge.yaml", "UCCSD |G>", IndexConvention.HalfUp)

    assert(is3.Method == "UnitaryCoupledCluster")
    assert(is3 == is4)
    assert(is1 != is4)


def test_load_greedy_state():
    """
    Checks that not passing a wavefunction label generates the greedy input state.
    """
    broombridge = load_broombridge("broombridge.yaml")

    is1 = broombridge.problem_description[0].load_input_state("", IndexConvention.HalfUp)
    is2 = load_input_state("broombridge.yaml", "", IndexConvention.HalfUp)

    assert(is1.Method == "SparseMultiConfigurational")
    assert(is1 == is2)

    
def test_jw_encode():
    """
    Checks that we can encode a hamiltonian + input state
    """
    broombridge = load_broombridge("broombridge.yaml")

    fh1 = broombridge.problem_description[0].load_fermion_hamiltonian()
    is1 = broombridge.problem_description[0].load_input_state()

    jw = encode(fh1, is1)
    assert(len(jw) == 4)
