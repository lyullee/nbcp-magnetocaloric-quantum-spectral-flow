"""Low-energy-corrected finite-temperature Lanczos thermodynamics.

The lowest eigenpairs are treated exactly.  Random-vector Lanczos quadrature
then estimates only the orthogonal complement.  This removes the dominant
low-temperature variance of plain FTLM while retaining matrix-free scaling.
Energies and temperatures are in kelvin (k_B = 1).
"""
from __future__ import annotations

from dataclasses import dataclass
import time

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.sparse.linalg import eigsh


@dataclass
class LanczosRule:
    energies: np.ndarray
    weights: np.ndarray
    residual: float


@dataclass
class SectorRules:
    dimension: int
    low_energies: np.ndarray
    random_rules: list[LanczosRule]


def lanczos_rule(operator, start, krylov_dim=60, tolerance=1e-12,
                 full_reorthogonalization=True):
    """Return Gauss quadrature nodes and weights for a normalized start."""
    q=np.asarray(start,dtype=np.complex128).copy()
    norm=np.linalg.norm(q)
    if norm==0:
        raise ValueError('Lanczos start vector is zero')
    q/=norm
    dimension=len(q)
    limit=min(int(krylov_dim),dimension)
    basis=np.empty((limit,dimension),dtype=np.complex128)
    alpha=[]
    beta=[]
    previous=np.zeros_like(q)
    beta_previous=0.
    residual=0.
    for step in range(limit):
        basis[step]=q
        z=operator@q-beta_previous*previous
        a=float(np.vdot(q,z).real)
        z-=a*q
        if full_reorthogonalization:
            # Two passes are deliberate: the low-T trace is sensitive to
            # duplicated Ritz values caused by loss of orthogonality.
            for _ in range(2):
                coefficients=basis[:step+1].conj()@z
                z-=coefficients@basis[:step+1]
        b=float(np.linalg.norm(z))
        alpha.append(a)
        residual=b
        if b<tolerance or step+1==limit:
            break
        beta.append(b)
        previous,q=q,z/b
        beta_previous=b
    eigenvalues,eigenvectors=eigh_tridiagonal(np.asarray(alpha),np.asarray(beta))
    weights=np.abs(eigenvectors[0])**2
    return LanczosRule(eigenvalues,weights,float(residual))


def lowest_eigenpairs(operator,count,tolerance=1e-10,seed=0):
    """Compute sorted low eigenpairs; dense fallback is for validation only."""
    dimension=operator.shape[0]
    count=min(int(count),max(0,dimension-1))
    if count==0:
        return np.empty(0),np.empty((dimension,0),dtype=np.complex128)
    if dimension<=64 or count>=dimension-1:
        eye=np.eye(dimension,dtype=np.complex128)
        dense=np.column_stack([operator@eye[:,i] for i in range(dimension)])
        values,vectors=np.linalg.eigh(dense)
        return values[:count],vectors[:,:count]
    rng=np.random.default_rng(seed)
    v0=rng.normal(size=dimension)+1j*rng.normal(size=dimension)
    values,vectors=eigsh(operator,k=count,which='SA',tol=tolerance,v0=v0)
    order=np.argsort(values)
    return np.real(values[order]),vectors[:,order]


def build_sector_rules(operator,random_vectors=8,krylov_dim=60,
                       exact_low_states=12,seed=0):
    """Build reusable quadrature rules for one symmetry sector."""
    dimension=operator.shape[0]
    low_values,low_vectors=lowest_eigenpairs(operator,exact_low_states,seed=seed)
    complement_dimension=dimension-len(low_values)
    rng=np.random.default_rng(seed)
    rules=[]
    if complement_dimension:
        for _ in range(int(random_vectors)):
            start=rng.normal(size=dimension)+1j*rng.normal(size=dimension)
            if low_vectors.size:
                start-=low_vectors@(low_vectors.conj().T@start)
            start_norm=np.linalg.norm(start)
            if start_norm<1e-12:
                raise RuntimeError('random start vanished after low-state projection')
            rules.append(lanczos_rule(operator,start/start_norm,krylov_dim))
    return SectorRules(dimension,np.asarray(low_values),rules)


def thermodynamics(sectors,temperatures,spins):
    """Combine symmetry-sector rules into Z, E, C and S per spin.

    The trace of the orthogonal complement is estimated by multiplying each
    normalized random-vector quadrature by (D-k)/R.  Exact low states carry
    unit trace weight.
    """
    temperatures=np.asarray(temperatures,dtype=float)
    if np.any(temperatures<=0):
        raise ValueError('temperatures must be positive')
    all_min=[]
    for sector in sectors:
        if len(sector.low_energies):
            all_min.append(float(np.min(sector.low_energies)))
        all_min.extend(float(np.min(rule.energies)) for rule in sector.random_rules)
    if not all_min:
        raise ValueError('no trace contributions')
    reference=min(all_min)
    rows=[]
    for temperature in temperatures:
        z=0.
        e1=0.
        e2=0.
        for sector in sectors:
            if len(sector.low_energies):
                weights=np.exp(-(sector.low_energies-reference)/temperature)
                z+=float(np.sum(weights))
                e1+=float(np.dot(weights,sector.low_energies))
                e2+=float(np.dot(weights,sector.low_energies**2))
            complement=sector.dimension-len(sector.low_energies)
            if sector.random_rules:
                scale=complement/len(sector.random_rules)
                for rule in sector.random_rules:
                    weights=rule.weights*np.exp(-(rule.energies-reference)/temperature)
                    z+=scale*float(np.sum(weights))
                    e1+=scale*float(np.dot(weights,rule.energies))
                    e2+=scale*float(np.dot(weights,rule.energies**2))
        mean=e1/z
        variance=max(0.,e2/z-mean**2)
        entropy_total=np.log(z)+(mean-reference)/temperature
        rows.append({
            'temperature_K':float(temperature),
            'log_partition':float(np.log(z)-reference/temperature),
            'energy_K_total':float(mean),
            'heat_capacity_kB_total':float(variance/temperature**2),
            'entropy_kB_total':float(entropy_total),
            'entropy_kB_per_spin':float(entropy_total/spins),
        })
    return rows


def run_ftlm(operators,temperatures,spins,random_vectors=8,krylov_dim=60,
             exact_low_states=12,seed=0):
    start=time.perf_counter()
    sectors=[]
    for index,operator in enumerate(operators):
        sectors.append(build_sector_rules(
            operator,random_vectors=random_vectors,krylov_dim=krylov_dim,
            exact_low_states=exact_low_states,seed=seed+104729*index))
    rows=thermodynamics(sectors,temperatures,spins)
    return {
        'method':'low-energy-corrected FTLM',
        'spins':int(spins),
        'sector_dimensions':[sector.dimension for sector in sectors],
        'random_vectors_per_sector':int(random_vectors),
        'krylov_dimension':int(krylov_dim),
        'exact_low_states_per_sector':int(exact_low_states),
        'seed':int(seed),
        'max_last_lanczos_beta_K':float(max(
            (rule.residual for sector in sectors for rule in sector.random_rules),
            default=0.)),
        'low_energies_K_by_sector':[sector.low_energies.tolist() for sector in sectors],
        'elapsed_seconds':float(time.perf_counter()-start),
        'rows':rows,
    }
