"""Map the 3D Ewald dipolar sum onto finite triangular quantum clusters."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from dipole_ewald import (ewald_tensor, nbcpo_primitive_cell,
                          spin_coupling_microeV)
from matrix_free_hamiltonian import SpinBond, weak_model_bonds
from screen_weak_terms import MEV_TO_K


MICROEV_TO_K=MEV_TO_K/1000.


@dataclass(frozen=True)
class EwaldSettings:
    alpha_per_angstrom: float=.22
    real_cutoff_angstrom: float=42.
    reciprocal_cutoff_per_angstrom: float=8.
    boundary: str='tin_foil'


def cluster_cell(graph,a_angstrom=5.312,c_angstrom=6.983):
    primitive=nbcpo_primitive_cell(a_angstrom,c_angstrom)
    a,b=graph['translations']
    return np.column_stack((primitive[:,:2]@np.asarray(a,dtype=float),
                            primitive[:,:2]@np.asarray(b,dtype=float),
                            primitive[:,2]))


def _quotient_key(graph,x,y):
    a,b=graph['translations']
    n=graph['N']
    return ((b[1]*x-b[0]*y)%n,(-a[1]*x+a[0]*y)%n)


def ewald_dipole_bonds(graph,settings=EwaldSettings(),g=(4.2,4.2,4.716)):
    """Return one Ewald-summed dipolar tensor for every unordered site pair.

    Interactions with periodic images of the same finite-cluster site reduce
    to a spin-1/2 constant for the real symmetric self tensor and are omitted.
    The constant does not affect entropy, heat capacity, or field optima.
    """
    primitive=nbcpo_primitive_cell()
    cell=cluster_cell(graph)
    coordinates=graph['coordinates']
    cache={}
    bonds=[]
    for i in range(graph['N']):
        xi,yi=coordinates[i]
        for j in range(i+1,graph['N']):
            xj,yj=coordinates[j]
            dx,dy=xj-xi,yj-yi
            forward=_quotient_key(graph,dx,dy)
            reverse=_quotient_key(graph,-dx,-dy)
            if reverse<forward:
                dx,dy=-dx,-dy
                key=reverse
            else:
                key=forward
            if key not in cache:
                displacement=primitive[:,:2]@np.asarray([dx,dy],dtype=float)
                geometric=ewald_tensor(
                    displacement,cell,settings.alpha_per_angstrom,
                    settings.real_cutoff_angstrom,
                    settings.reciprocal_cutoff_per_angstrom,
                    settings.boundary)
                cache[key]=spin_coupling_microeV(geometric,g)*MICROEV_TO_K
            bonds.append(SpinBond(i,j,cache[key]))
    return bonds


def m3_bonds(graph,settings=EwaldSettings()):
    """M1 exchange/Jpp plus the complete periodic long-range dipolar sum."""
    return (weak_model_bonds(graph,'XXZ_plus_additional_Jpp')+
            ewald_dipole_bonds(graph,settings))
