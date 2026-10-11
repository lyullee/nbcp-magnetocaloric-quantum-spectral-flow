from finite_cluster_stability_alignment import cluster_critical_mode
from screen_model import lattice


def run():
    # Every three-colour-compatible M0 cluster contains K exactly, so its
    # transverse instability equals the analytic bulk M0 value.
    for graph in (lattice((3,0),(0,3)),lattice((4,1),(1,4)),
                  lattice((4,1),(-2,4))):
        r=cluster_critical_mode(graph,'M0')
        assert abs(r['critical_field_T']-1.7743602002000802)<1e-10
    print('finite-cluster stability validation passed')


if __name__=='__main__':run()
