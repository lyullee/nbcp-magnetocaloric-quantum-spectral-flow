"""Create the compact figure used in the first feasibility report."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'
FIG=ROOT/'figures'
FIG.mkdir(exist_ok=True)

fine=pd.read_csv(OUT/'n12_fine_field_profiles.csv')
decision=pd.read_csv(OUT/'n12_fine_field_decisions.csv')
fixed=pd.read_csv(OUT/'finite_size_stability_check.csv')

plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,3,figsize=(10.8,3.35),constrained_layout=True)

for model,label,color in [('base_XXZ','Base XXZ','#3478bf'),
                          ('representative_Jpp_plus_NN_dipole','Realistic weak terms','#d1495b')]:
    x=fine[(fine.model==model)&(fine.T_cold_K==.1)].sort_values('B_end_T')
    ax[0].plot(x.B_end_T,x.Q_model_mJ_per_cm3,'o-',ms=2.8,lw=1.4,label=label,color=color)
ax[0].set(xlabel='Endpoint field (T)',ylabel=r'Conditional $Q$ (mJ cm$^{-3}$)',title='(a) 100 mK field scan, N=12')
ax[0].legend(frameon=False,fontsize=8)

for field,color in [(1.3,'#2a9d8f'),(1.5,'#f4a261'),(1.66,'#7b2cbf')]:
    x=fixed[fixed.B_end_T==field].sort_values('T_cold_K')
    # The stored percentage is absolute, so restore the correction sign.
    signed=np.sign(x.delta_Q_N12_mJ_per_cm3)*x.N12_abs_percent_of_public_best_Q
    ax[1].plot(x.T_cold_K,signed,'o-',ms=3,lw=1.3,label=f'{field:g} T',color=color)
ax[1].axhline(10,color='0.6',ls='--',lw=.8);ax[1].axhline(-10,color='0.6',ls='--',lw=.8)
ax[1].set(xlabel='Cold temperature (K)',ylabel=r'$\Delta Q$ / public best $Q$ (%)',title='(b) Fixed-field correction, N=12')
ax[1].legend(frameon=False,fontsize=8)

ax[2].plot(decision.T_cold_K,decision.base_optimal_B_on_fine_grid_T,'o-',label='Base XXZ',color='#3478bf')
ax[2].plot(decision.T_cold_K,decision.perturbed_optimal_B_on_fine_grid_T,'o-',label='Realistic weak terms',color='#d1495b')
ax[2].set(xlabel='Cold temperature (K)',ylabel='Best endpoint field (T)',title='(c) Re-optimized field, N=12')
ax[2].legend(frameon=False,fontsize=8)

for a in ax:a.grid(alpha=.2)
fig.suptitle('First finite-cluster feasibility screen (classical exact diagonalization)',fontsize=11)
fig.savefig(FIG/'first_feasibility_screen.png',dpi=300)
fig.savefig(FIG/'first_feasibility_screen.pdf')
print(FIG/'first_feasibility_screen.png')
