#!/usr/bin/env python
# coding: utf-8

# In[41]:


import warnings
warnings.filterwarnings('ignore')
import flowermd
import hoomd
import gsd
import matplotlib.pyplot as plt
import numpy as np
import gsd.hoomd
from flowermd.base import Pack, Simulation
from flowermd.library import KremerGrestBeadSpring, LJChain
from flowermd.utils import get_target_box_number_density
from flowermd.base import Molecule
from flowermd.library import Graphene
from flowermd.library.forcefields import BeadSpring
from flowermd.base.forcefield import BaseHOOMDForcefield
import unyt as u
import hoomd
import mbuild as mb
import numpy as np
from cmeutils.visualize import FresnelGSD
cpu = hoomd.device.CPU()


# In[42]:


import mbuild as mb
from mbuild.compound import Compound
from mbuild.lattice import Lattice

from flowermd.base import Molecule, System


class Graphene(System):
    def __init__(
        self,
        x_repeat,
        y_repeat,
        n_layers,
        base_units=dict(),
        periodicity=(True, True, False),
    ):
        surface = mb.Compound(periodicity=periodicity)
        scale = 0.246 # 1.0 for testing, 0.246 for graphene bond lengths
        spacings = [s * scale for s in [1.0, 1.0, 1.0]]
        points = [[0, 0, 0], [2/3, 1/3, 0]]
        lattice = Lattice(
            lattice_spacing=spacings,
            angles=[90, 90, 120],
            lattice_points={"A": points},
        )
        carbon = Compound(name="C", element="C")
        layers = lattice.populate(
            compound_dict={"A": carbon}, x=x_repeat, y=y_repeat, z=n_layers
        )
        surface.add(layers)
        surface.freud_generate_bonds("C", "C", dmin=0.3*scale, dmax=0.6*scale)
        surface_mol = Molecule(num_mols=1, compound=surface)
        super(Graphene, self).__init__(
            molecules=[surface_mol],
            base_units=base_units,
        )

    def _build_system(self):
        return self.all_molecules[0]


# In[43]:


#OK, so we want to initialize a system with some chains and some flakes
kg_chain = LJChain(lengths=10,num_mols=500)
sheet = Graphene(x_repeat=5, y_repeat=5, n_layers=1, periodicity=(False, False, False))
system = Pack(molecules=[Molecule(compound=sheet.all_molecules[0], num_mols=25), kg_chain], density=0.2, packing_expand_factor = 5)


# In[44]:


system.visualize()


# In[45]:


s = 0.246 # 1.0 for scaling


# In[46]:


ff = BeadSpring(r_cut=2.5,
        beads={"A": dict(epsilon=.05, sigma=1.0*s),
               "C": dict(epsilon=0.05, sigma=1.0*s)},
        bonds={"C-C": dict(r0=0.57736*s, k=500), "A-A": dict(r0=1.1, k=0.01)},
        angles={"A-A-A": dict(t0=np.pi, k=0.1),
                "C-C-C": dict(t0=2*np.pi/3., k=500)},
        dihedrals={"A-A-A-A": dict(phi0=0.0, k=0.1, d=-1, n=1), "C-C-C-C": dict(phi0=0.1, k=0.1, d=-1, n=1)})


# In[47]:


# this FF is for research question
ff = BeadSpring(
    r_cut=2**(1/6) * s,  # WCA cutoff (purely repulsive)
    beads={
        "A": dict(epsilon=1.0, sigma=1.0 * s),  # chains
        "C": dict(epsilon=1.0, sigma=1.0 * s),  # flakes
    },
    bonds={
        "C-C": dict(r0=0.57736 * s, k=500),
        "A-A": dict(r0=1.0, k=500.0),  # increased k to avoid chain collapse
    },
    angles={
        "A-A-A": dict(t0=np.pi, k=10.0),   # moderate stiffness for chains
        "C-C-C": dict(t0=2 * np.pi / 3., k=500),
    },
    dihedrals={
        "A-A-A-A": dict(phi0=0.0, k=1.0, d=-1, n=1),
        "C-C-C-C": dict(phi0=0.1, k=0.1, d=-1, n=1),
    }
)


# In[49]:


sim = Simulation(initial_state=system.hoomd_snapshot, forcefield=ff.hoomd_forces, device=cpu, dt = 0.00005, gsd_write_freq=int(5000), log_file_name = "ej_flakes.txt")
sim.run_NVT(n_steps=5000, kT=0.001, tau_kt=1.0)
sim.run_NVT(n_steps=5e8,kT=7, tau_kt = 0.1) #short for initial testing
sim.flush_writers()


# In[50]:


sim_visualizer = FresnelGSD(gsd_file="trajectory.gsd", frame=2, view_axis=(1, 1, 1))
sim_visualizer.view()
# dt of .0001 is as big as we can get with the current parameterization
#Issues:
# Graphene flakes are initialized way too small, or the bond lengths are wrong
# Need to double-check all the bonds, angle, dihedral parameters
# then need to go back and get some chain lengths and flake sizes right.


# In[51]:


sim_visualizer = FresnelGSD(gsd_file="trajectory.gsd", frame=56, view_axis=(1, 1, 1))
sim_visualizer.view()


# In[ ]:




