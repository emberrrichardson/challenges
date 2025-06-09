#!/usr/bin/env python
# coding: utf-8

# In[32]:


import warnings
warnings.filterwarnings('ignore')
import hoomd
import gsd
import matplotlib.pyplot as plt
import numpy as np
import gsd.hoomd
from flowermd.utils import get_target_box_number_density
import unyt as u
import hoomd


# In[33]:


from flowermd.base import Pack,Lattice, Simulation
from flowermd.library import EllipsoidChain

ellipsoid_chain = EllipsoidChain(lengths=1,num_mols=128,lpar=1.0,bead_mass=1.0)
system = Pack(molecules=ellipsoid_chain, density=0.05*u.Unit("nm**-3"), 
              packing_expand_factor=6,edge=2,overlap=1,fix_orientation=True)


# In[34]:


from flowermd.library import EllipsoidForcefield

ff = EllipsoidForcefield(epsilon=1.0,lpar=1.0,lperp=0.5,
                         r_cut=2.5,bond_k=100,bond_r0=0)


# In[35]:


from flowermd.utils.constraints import create_rigid_ellipsoid_chain

rigid_frame, rigid = create_rigid_ellipsoid_chain(system.hoomd_snapshot)



# In[36]:


gsd_path=('ellipsoid-chain-2mer.gsd')
ellipsoid_sim = Simulation(
    initial_state=rigid_frame,
    forcefield=ff.hoomd_forces,
    constraint=rigid,
    dt=0.001,
    gsd_write_freq=int(1e4),
    gsd_file_name=gsd_path,
    log_write_freq=int(1e4),
    log_file_name='log.txt')

target_box = get_target_box_number_density(density=0.8*u.Unit("nm**-3"),n_beads=128)
ellipsoid_sim.run_update_volume(final_box_lengths=target_box, kT=1.0, n_steps=1e5,tau_kt=100*ellipsoid_sim.dt,period=10,
                                thermalize_particles=True)
print("shrink finished")
ellipsoid_sim.run_NVT(n_steps=1e6, kT=1.0, tau_kt=100*ellipsoid_sim.dt)
#ellipsoid_sim.save_restart_gsd("restart.gsd")
ellipsoid_sim.flush_writers()
#ellipsoid_sim.save_simulation("sim.pickle")


# In[37]:


ff = EllipsoidForcefield(epsilon=1.0,lpar=1.0,lperp=0.5,r_cut=2.0,bond_k=100,bond_r0=0)
rigid = restart_rigid_ellipsoid_sim()
gsd_path = "sim_restart.gsd"

sim = Simulation(
    initial_state="restart.gsd",
    forcefield=ff.hoomd_forces,
    constraint=rigid,
    gsd_write_freq=int(1e3),
    gsd_file_name=gsd_path,
    log_write_freq=int(1e4),
    log_file_name='log2.txt',
)
sim.run_NVT(n_steps=1e5, kT=1.0, tau_kt=10*sim.dt)
sim.flush_writers()


# In[ ]:


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


ellipsoid_gsd(gsd_file=gsd_path,new_file="ovito-ellipsoid.gsd",lpar=1.0,lperp=0.5, ellipsoid_types="R")


# In[ ]:




