import numpy as np

from bulk_to_engineering_surrogate import inverse_field, transform_field


def run():
    x=np.array([0.,1.5,2.])
    for method in ('reduced_field','critical_offset'):
        y=transform_field(x,1.75,1.87,method)
        assert np.max(np.abs(inverse_field(y,1.75,1.87,method)-x))<1e-14
    print('bulk-to-engineering bridge validation passed')


if __name__=='__main__':run()
