#!/bin/env python
# -*- coding: utf-8 -*-
##
# tomography.py: Single qubit tomography of Q# operations.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##

## IMPORTS ##

import qsharp
import numpy as np

try:
    import qinfer as qi
except ImportError:
    qi = None

try:
    import qutip as qt
except ImportError:
    qt = None

def projector(P):
    if qt is None:
        raise ImportError("This function requires QuTiP.")

    return (qt.qeye(2) + P) / 2.0

def single_qubit_process_tomography(
    operation,
    n_measurements=2000,
    n_particles=4000
):
    """
    :param operation: A Q# operation of type ((Pauli, Pauli) => Result) whose
        inputs are named `prep` and `meas` (respectively), representing
        a state preparation, evolution, and measurement.
    """
    if qt is None:
        raise ImportError("This function requires QuTiP.")
    if qi is None:
        raise ImportError("This function requires QInfer.")

    print("Preparing tomography model...")
    state_basis = qi.tomography.pauli_basis(1)
    prior = qi.tomography.BCSZChoiDistribution(state_basis)
    model = qi.tomography.TomographyModel(prior.basis)

    updater = qi.SMCUpdater(model, n_particles, prior)

    print("Performing tomography...")
    for idx_experiment in range(n_measurements):
        prep = qsharp.Pauli.sample()
        meas = qsharp.Pauli.sample()

        # Convert into a QuTiP object by using the standard transformation
        # between state and process tomography.
        qobj = 2.0 * qt.tensor(
            projector(prep.as_qobj()).trans(), projector(meas.as_qobj())
        )
        expparams = np.array(
            [(model.basis.state_to_modelparams(qobj),)],
            dtype=model.expparams_dtype
        )

        datum = 1 - operation.simulate(prep=prep, meas=meas)

        updater.update(datum, expparams)

    return {
        # We multiply by 2 to turn into a Choi–Jamiłkowski operator instead
        # of a Choi–Jamiłkowski state.
        'est_channel': 2.0 * model.basis.modelparams_to_state(updater.est_mean()),
        # Returning the updater allows for exploring properties not extracted
        # elsewhere.
        'posterior': updater
    }
