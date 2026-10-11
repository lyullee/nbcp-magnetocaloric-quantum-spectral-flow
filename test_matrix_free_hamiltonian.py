import numpy as np

from matrix_free_hamiltonian import SpinHamiltonianOperator, weak_model_bonds
from screen_model import lattice
from screen_weak_terms import matrices


def run():
    rng=np.random.default_rng(1741)
    for graph in [lattice((3,0),(0,3)),lattice((2,2),(-2,4))]:
        for model in ['XXZ_2026_parameters','XXZ_plus_additional_Jpp',
                      'XXZ_plus_nn_dipoles_plus_additional_Jpp']:
            reference=matrices(graph,model)
            bonds=weak_model_bonds(graph,model)
            for parity,(dense,mz) in enumerate(reference):
                op=SpinHamiltonianOperator(graph['N'],bonds,parity=parity)
                assert op.shape==dense.shape
                # N=12 dense multiplication is expensive enough that one
                # independent vector per block is sufficient here.
                repeats=3 if graph['N']==9 else 1
                for _ in range(repeats):
                    x=rng.normal(size=op.shape[0])+1j*rng.normal(size=op.shape[0])
                    assert np.max(np.abs((op@x)-dense@x))<3e-12
                    y=rng.normal(size=op.shape[0])+1j*rng.normal(size=op.shape[0])
                    assert abs(np.vdot(x,op@y)-np.vdot(op@x,y))<5e-9
    print('matrix-free Hamiltonian validation passed')


if __name__=='__main__':run()
