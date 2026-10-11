"""Add the two half-grid N=18 fields required for a 0.05 T Gate-1 test."""
from n18_gate1_pilot import CACHE, calculate
from screen_model import lattice


def main():
    CACHE.mkdir(parents=True,exist_ok=True)
    graph=lattice((3,-3),(3,3))
    for model in ('M0','M3'):
        for field in (1.15,1.25):
            calculate(graph,model,field)


if __name__=='__main__':
    main()
