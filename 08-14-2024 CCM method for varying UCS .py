# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:08:51 2024

@author: vikra
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%%
#input properties
sp_wt = 2.7e-2 #MN/m3
H_over = 150 #m overbuden mass
po = sp_wt*H_over #Mpa vertical stress
R = 2.5 #m raduis of the tunnel
L = 1.5 #m Face effect

#Rock mass properties
Sig_ci = np.linspace(10,25,8) #Mpa
GSI = 15
mi = 7
nu = 0.3
fric = math.radians(25)

#%%
#ground reaction curve
# Initialize the figure
f, axs = plt.subplots(1, 2, figsize=(12, 6))

sns.set_palette("colorblind")
# Loop through each GSI value and plot
for si_ci in Sig_ci:
    mb = mi * math.exp((GSI - 100) / 28)
    if GSI >= 25:
        s = math.exp((GSI - 100) / 9)
        a = 0.5
    else:
        s = 0
        a = 0.65 - GSI / 200

    Erm = math.sqrt(si_ci / 100) * 10**((GSI - 10) / 40)   # GPa
    Grm = Erm / (2 * (1 + nu))  # GPa
    K = (1 + math.sin(fric)) / (1 - math.sin(fric))
    So = (po / (mb * si_ci)) + s / mb**2  # MPa far-field stress
    Pi_cr = 1/16 * (1 - math.sqrt(1 + 16 * So))**2  # Mpa scaled internal critical pressure
    pi_cr = (Pi_cr - s / mb**2) * mb * si_ci  # MPa actual internal critical pressure

    # Dataframe creation for GRC and LDP
    GRC = pd.DataFrame({"Point": np.linspace(1, 20, num=20)})
    GRC['pi, Mpa'] = pi_cr * (20 - GRC["Point"]) / 19  # pip_grc
    GRC["Pi, Mpa"] = (GRC['pi, Mpa'] / (mb * si_ci)) + s / mb**2  # pips_grc
    GRC["Rpl/R"] = np.exp(2 * Pi_cr**0.5 - 2 * GRC["Pi, Mpa"]**0.5)  # xi_grc
    GRC["pi_r"] = GRC['pi, Mpa'] + sp_wt * (GRC["Rpl/R"] * R - R)
    GRC["pi_f"] = GRC['pi, Mpa'] - sp_wt * (GRC["Rpl/R"] * R - R)
    a = ((K - 1) / (K + 1)) + (2 / (K + 1) * GRC["Rpl/R"]**(K + 1))
    b = ((1 - 2 * nu) / (4 * (So - Pi_cr))) * np.log(GRC["Rpl/R"])**2
    c = ((1 - 2 * nu) / (1 + K)) * (math.sqrt(Pi_cr)) / (So - Pi_cr) + \
        ((1 - nu) / 2) * ((K - 1) / (K + 1)**2) * (1 / (So - Pi_cr))
    d = (K + 1) * np.log(GRC["Rpl/R"]) - (GRC["Rpl/R"])**(K + 1) + 1
    e = (R * (po - pi_cr) / (2 * Grm * 1000)) * 1000
    GRC["Urm, mm"] = (a + b - c * d) * e

    # For LDP
    u_max = GRC['Urm, mm'].max()
    GRC['Lr/R'] = -4 + (GRC['Point'] - 1) * 12 / 11
    GRC['Lf, m'] = GRC['Lr/R'] * R
    GRC['U, mm'] = u_max * (1 + np.exp(- GRC['Lf, m'] / 1.1 / R))**-1.7

    # Plotting on the same figure
    sns.lineplot(data=GRC, x='Urm, mm', y="pi, Mpa", ax=axs[0], markers=True, label=f'UCS={si_ci}')
    sns.lineplot(data=GRC, x='Lf, m', y='U, mm', ax=axs[1], markers=True, label=f'UCS={si_ci}')

# Customizing plots
axs[0].grid(True)
axs[0].set_xlabel("Radial displacement, mm")
axs[0].set_ylabel("Insitu stress, Mpa")
axs[0].legend(title="UCS, MPa")
axs[0].set_title("GRC for different UCS values")
# axs[0].set_xticks(np.arange(0, axs[0].get_xlim()[1] + 2, 2))

axs[1].set_xlabel("Distance from Face, m")
axs[1].set_ylabel("Displacement, mm")
axs[1].grid(True)
axs[1].legend(title="UCS, MPa")
axs[1].set_title("LDP for different UCS values")
axs[1].axvline(x=0, color='black', linestyle='--')
axs[1].axvline(x=1.5, color='red', linestyle='--') 
# axs[1].set_yticks(np.arange(0, axs[1].get_ylim()[1] + 5, 5))
# axs[1].set_xticks(np.arange(-10, axs[1].get_xlim()[1] + 5, 5))  # Set y-axis grid interval to 5 for second plot
  # Set y-axis grid interval to 5 for second plot

# Find the intersection point at x=1.5

  # Add vertical line at x=0
# axs[1].set_xlim(left=-axs[1].get_xlim()[1]) 
# plt.suptitle("GRC and LDP for various GSI value")
sns.set_style("dark")
f.tight_layout()
plt.savefig("GRC and LDP for varying UCS value.png", dpi=300)
plt.show()