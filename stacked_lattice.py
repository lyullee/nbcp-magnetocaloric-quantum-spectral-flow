"""Periodic AA-stacked triangular multigraph for the NBCP M4 model.

The interlayer convention follows the inset of Fig. 9(b) in arXiv:2607.23173:
one vertical J2 bond and two symmetry-inequivalent sets of three diagonal
bonds (Ja and Jb) connect a site to the next layer.  Choosing which oriented
set is called Ja is a coordinate convention; reversing c swaps Ja and Jb.

Bond multiplicities are retained.  This is essential for small periodic
cells, where different bonds in the infinite lattice can map to the same pair
of finite-cell basis states.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Tuple


@dataclass(frozen=True)
class Bond:
    i: int
    j: int
    kind: str
    displacement_primitive: Tuple[int, int, int]


def stacked_triangular_multigraph(lx: int, ly: int, lz: int):
    if min(lx,ly,lz)<2:
        raise ValueError('all periodic lengths must be at least two')
    sites=[(x,y,z) for z in range(lz) for y in range(ly) for x in range(lx)]
    index={s:i for i,s in enumerate(sites)}
    def at(x,y,z):return index[(x%lx,y%ly,z%lz)]
    bonds=[]
    inplane=[(1,0,0),(0,1,0),(1,-1,0)]
    plus=[(1,0,1),(0,1,1),(1,-1,1)]
    minus=[(-1,0,1),(0,-1,1),(-1,1,1)]
    for i,(x,y,z) in enumerate(sites):
        for d in inplane:
            bonds.append(Bond(i,at(x+d[0],y+d[1],z),'J1',d))
        bonds.append(Bond(i,at(x,y,z+1),'J2',(0,0,1)))
        for d in plus:
            bonds.append(Bond(i,at(x+d[0],y+d[1],z+1),'Ja',d))
        for d in minus:
            bonds.append(Bond(i,at(x+d[0],y+d[1],z+1),'Jb',d))
    return {'shape':(lx,ly,lz),'N':len(sites),'sites':sites,'bonds':bonds,
            'convention':'Ja uses +a1,+a2,+a3 into layer z+1; Jb uses the opposite three. Reversing c swaps labels.'}


def degree_by_kind(graph):
    degree={kind:[0]*graph['N'] for kind in ['J1','J2','Ja','Jb']}
    for b in graph['bonds']:
        degree[b.kind][b.i]+=1
        degree[b.kind][b.j]+=1
    return degree


def multiplicities(graph):
    return Counter((min(b.i,b.j),max(b.i,b.j),b.kind) for b in graph['bonds'])

