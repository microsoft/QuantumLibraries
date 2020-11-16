import pytest
import os
import logging

import numpy as np

from pathlib import Path

from qsharp_chemistry import (
    load_broombridge, 
    load_fermion_hamiltonian, 
    load_input_state, 
    encode, 
    load_and_encode
)
from qsharp_chemistry.problem_description import IndexConvention, InputState

import qsharp

logging.basicConfig(level=logging.INFO)
print ( qsharp.component_versions() )

@pytest.fixture()
def base_path():
    return Path(__file__).parent.absolute()


def test_load_broombridge(base_path):
    """
    Checks that we can load a broombridge schema file.
    """
    broombridge = load_broombridge(os.path.join(base_path, "broombridge.yaml"))
    assert(len(broombridge.problem_description) == 1)
    assert(len(broombridge.bibliography) == 3)


def test_load_fermion_hamiltonian(base_path):
    """
    Checks that loading a fermion hamiltonian from file or from broombridge gives the same result.
    """
    broombridge = load_broombridge(os.path.join(base_path, "broombridge.yaml"))

    fh1 = broombridge.problem_description[0].load_fermion_hamiltonian()
    fh2 = load_fermion_hamiltonian(os.path.join(base_path, "broombridge.yaml"))

    assert(len(fh1.terms) == 6)
    assert(fh1 == fh2)
    
    fh3 = broombridge.problem_description[0].load_fermion_hamiltonian(IndexConvention.HalfUp)
    fh4 = load_fermion_hamiltonian(os.path.join(base_path, "broombridge.yaml"), IndexConvention.HalfUp)

    assert(len(fh3.terms) == 6)
    assert(fh3 == fh4)
    
    # should this be true? assert(fh3 != fh1)


def test_load_input_state(base_path):
    """
    Checks that loading an input state from file or from broombridge gives the same result.
    """
    broombridge = load_broombridge(os.path.join(base_path, "broombridge.yaml"))

    is1 = broombridge.problem_description[0].load_input_state("UCCSD |G>")
    is2 = load_input_state(os.path.join(base_path, "broombridge.yaml"), "UCCSD |G>")

    assert(is1.Method == "UnitaryCoupledCluster")
    assert(is1 == is2)
    
    is3 = broombridge.problem_description[0].load_input_state("UCCSD |G>", IndexConvention.HalfUp)
    is4 = load_input_state(os.path.join(base_path, "broombridge.yaml"), "UCCSD |G>", IndexConvention.HalfUp)

    assert(is3.Method == "UnitaryCoupledCluster")
    assert(is3 == is4)
    assert(is1 != is4)


def test_load_greedy_state(base_path):
    """
    Checks that not passing a wavefunction label generates the greedy input state.
    """
    broombridge = load_broombridge(os.path.join(base_path, "broombridge.yaml"))

    is1 = broombridge.problem_description[0].load_input_state("", IndexConvention.HalfUp)
    is2 = load_input_state(os.path.join(base_path, "broombridge.yaml"), "", IndexConvention.HalfUp)

    assert(is1.Method == "SparseMultiConfigurational")
    assert(is1 == is2)

    
def test_jw_encode(base_path):
    """
    Checks that we can encode a hamiltonian + input state
    """
    broombridge = load_broombridge(os.path.join(base_path, "broombridge.yaml"))

    fh1 = broombridge.problem_description[0].load_fermion_hamiltonian()
    is1 = broombridge.problem_description[0].load_input_state()

    jw = encode(fh1, is1)
    assert(len(jw) == 4)


def test_load_and_encode(base_path):
    jw = load_and_encode(
        file_name=os.path.join(base_path, "broombridge.yaml"),
        problem_description_index=0,
        initial_state_label="UCCSD |G>"
    )
    (num_qubits, hamiltonian_term_list, input_state_terms, energy_offset) = jw
    assert num_qubits == 12
    assert np.isclose(energy_offset, -3.7893, 3)
    assert len(hamiltonian_term_list) == 4
    assert len(input_state_terms[1]) == 5
