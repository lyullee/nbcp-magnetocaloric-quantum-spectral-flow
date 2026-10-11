"""Frozen M0--M3 bond models and finite-cluster operators."""
from ewald_cluster import EwaldSettings, m3_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator, weak_model_bonds
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ


MODEL_NAMES={
    'M0':'2D nearest-neighbour XXZ',
    'M1':'M0 plus additional Jpp',
    'M2':'M1 plus nearest-neighbour dipole tensor',
    'M3':'M1 plus full 3D-periodic Ewald dipoles',
}


def model_bonds(graph,model,boundary='tin_foil'):
    if model=='M0':
        return weak_model_bonds(graph,'XXZ_2026_parameters')
    if model=='M1':
        return weak_model_bonds(graph,'XXZ_plus_additional_Jpp')
    if model=='M2':
        return weak_model_bonds(graph,'XXZ_plus_nn_dipoles_plus_additional_Jpp')
    if model=='M3':
        return m3_bonds(graph,EwaldSettings(boundary=boundary))
    raise ValueError(f'unknown model {model}')


def model_operators(graph,model,field_T,boundary='tin_foil'):
    bonds=model_bonds(graph,model,boundary)
    zeeman=GZ*MU_B_OVER_KB*field_T
    return [SpinHamiltonianOperator(graph['N'],bonds,zeeman,parity=p) for p in (0,1)]
