#!/usr/bin/env python
# coding: utf-8

# In[1]:


#imports
import warnings
warnings.filterwarnings('ignore')
import hoomd
import gsd
import matplotlib.pyplot as plt
import numpy as np
import gsd.hoomd
from flowermd.base import Pack,Lattice, Simulation
from flowermd.library import EllipsoidForcefield, EllipsoidChain
from flowermd.utils import get_target_box_number_density
from flowermd.utils.constraints import create_rigid_ellipsoid_chain
import unyt as u
import hoomd


# In[2]:


#defining
#no angles if lengths=1 -- will error
chain = EllipsoidChain(lengths=1, num_mols=128, lpar=1, bead_mass=1)
system = Pack(molecules=chain, density=1*u.Unit("nm**-3"), packing_expand_factor=5)
system.hoomd_snapshot.bonds.types
rigid_frame, d_constraint = create_rigid_ellipsoid_chain(snapshot = system.hoomd_snapshot)
forces = EllipsoidForcefield(epsilon=1.0, lpar=1.0, lperp=0.5, r_cut=2.0, bond_k=100, bond_r0=0)


# In[3]:


#define the sim
#change file name
ellipsoid_sim=Simulation(initial_state=rigid_frame, forcefield=forces.hoomd_forces, 
                         constraint=d_constraint, dt=0.001, gsd_write_freq=int(1e4), 
                         gsd_file_name="elipschal.gsd", log_write_freq=int(1e4), log_file_name="log-elipschal")


# In[5]:


#CAN CHANGE= density, kT(keep consistent in shrink + run, period


target_box = get_target_box_number_density(density=1*u.Unit("nm**-3"),n_beads=128)
ellipsoid_sim.run_update_volume(final_box_lengths=target_box, kT=.01, n_steps=1e5,tau_kt=100*ellipsoid_sim.dt,period=10,thermalize_particles=True)
print("shrink finished")
ellipsoid_sim.run_NVT(n_steps=1e6, kT=.001, tau_kt=10*ellipsoid_sim.dt)
ellipsoid_sim.save_restart_gsd("restart.gsd")
ellipsoid_sim.flush_writers()
ellipsoid_sim.save_simulation("sim.pickle")


# In[ ]:


#sim parmas
def ellipsoid_gsd(gsd_file, new_file, ellipsoid_types, lpar, lperp):
    with gsd.hoomd.open(new_file, "w") as new_t:
        with gsd.hoomd.open(gsd_file) as old_t:
            for snap in old_t:
                shape_dicts_list = []
                for ptype in snap.particles.types:
                    if ptype == ellipsoid_types or ptype in ellipsoid_types:
                        shapes_dict = {
                            "type": "Ellipsoid",
                            "a": lpar,
                            "b": lperp,
                            "c": lperp,
                        }
                    else:
                        shapes_dict = {"type": "Sphere", "diameter": 0.001}
                    shape_dicts_list.append(shapes_dict)
                snap.particles.type_shapes = shape_dicts_list
                snap.validate()
                new_t.append(snap)


# In[ ]:


#name + dimensions (again)
ellipsoid_gsd(gsd_file="ellipsoids1.gsd",new_file="ellipsoid-packing-ovitosim.gsd",lpar=1.0,lperp=0.5, ellipsoid_types="R")

