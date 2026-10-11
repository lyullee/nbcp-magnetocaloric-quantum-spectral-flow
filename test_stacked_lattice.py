from stacked_lattice import degree_by_kind, multiplicities, stacked_triangular_multigraph


def run():
    g=stacked_triangular_multigraph(3,3,3)
    assert g['N']==27
    counts={k:sum(b.kind==k for b in g['bonds']) for k in ['J1','J2','Ja','Jb']}
    assert counts=={'J1':81,'J2':27,'Ja':81,'Jb':81}
    degree=degree_by_kind(g)
    assert set(degree['J1'])=={6}
    assert set(degree['J2'])=={2}
    assert set(degree['Ja'])=={6}
    assert set(degree['Jb'])=={6}
    assert max(multiplicities(g).values())==1
    # N=18 is the largest planned first 3D quantum cluster.  It cannot
    # represent both 3-sublattice in-plane order and qz=1/3 simultaneously;
    # retaining multiplicities still gives a defined thermodynamic sensitivity
    # cluster rather than a phase-structure claim.
    g18=stacked_triangular_multigraph(3,3,2)
    assert g18['N']==18
    assert max(multiplicities(g18).values())>=2
    assert set(degree_by_kind(g18)['J2'])=={2}
    print('stacked lattice validation passed')


if __name__=='__main__':run()
