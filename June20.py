#!/usr/bin/env python
# coding: utf-8

# In[1]:


import warnings
warnings.filterwarnings('ignore')
import hoomd
import gsd
import matplotlib.pyplot as plt
import numpy as np
import gsd.hoomd
from flowermd.utils import get_target_box_number_density
import unyt as u
#sim_visualizer = FresnelGSD(
    #gsd_file="trajectory.gsd", frame=10, view_axis=(1, 1, 1))
import hoomd
#import signac
from flowermd.base import Pack,Lattice, Simulation
from flowermd.library import EllipsoidChain
from flowermd.library import EllipsoidForcefield

from flowermd.utils.constraints import create_rigid_ellipsoid_chain


# In[2]:


#changed edge value
ellipsoid_chain = EllipsoidChain(lengths=1,num_mols=128,lpar=1.0,bead_mass=1.0)
system = Pack(molecules=ellipsoid_chain, density=.25*u.Unit("nm**-3"), 
              packing_expand_factor=4,edge=1,overlap=1,fix_orientation=True)


# In[ ]:


ff = EllipsoidForcefield(epsilon=1.0,lpar=1.0,lperp=0.5,
                         r_cut=3.0,bond_k=100,bond_r0=0.1)


# In[ ]:


rigid_frame, rigid_constraint = create_rigid_ellipsoid_chain(system.hoomd_snapshot)


# In[ ]:


#Thermalize particles=True
#kT=2
#density = 0.3 nm**-3
#5e6 after shrink
gsd_path=('june20-am.gsd')
ellipsoid_sim = Simulation(
    initial_state=rigid_frame,
    forcefield = ff.hoomd_forces,
    constraint=rigid_constraint,
    dt=0.0001,
    gsd_write_freq=int(1e4),
    gsd_file_name=gsd_path,
    log_write_freq=int(1e4),
    log_file_name='june20-am.txt')

target_box = get_target_box_number_density(density=.5*u.Unit("nm**-3"),n_beads=128)
ellipsoid_sim.run_update_volume(final_box_lengths=target_box, kT=1.0, n_steps=5e6,tau_kt=5*ellipsoid_sim.dt,
                                period=5,thermalize_particles=True)
print("shrink finished")
ellipsoid_sim.run_NVT(n_steps=5e6, kT=1.0, tau_kt=5*ellipsoid_sim.dt)
#ellipsoid_sim.save_restart_gsd("restart.gsd")
ellipsoid_sim.flush_writers()
#ellipsoid_sim.save_simulation("sim.pickle")


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


ellipsoid_gsd(gsd_file='/Users/ember/Downloads/june17-11am.gsd',new_file="ovitoj17-11am.gsd",lpar=1.0,lperp=0.5,
              ellipsoid_types="R")


# In[ ]:


sim_data = np.genfromtxt("june20-am.txt", names=True)
PotentialEnergy = sim_data["mdcomputeThermodynamicQuantitiespotential_energy"]
Time = sim_data["flowermdbasesimulationSimulationtimestep"]
x= Time
y= PotentialEnergy
plt.plot(x,y/128)
plt.xlabel('Time')
plt.ylabel('PotentialEnergy')
plt.savefig("June20-am-PE.png")
plt.close()

KineticEnergy = sim_data["mdcomputeThermodynamicQuantitieskinetic_energy"]
Time = sim_data["flowermdbasesimulationSimulationtimestep"]
x= Time
y= KineticEnergy
plt.plot(x,y/128)
plt.xlabel('Time')
plt.ylabel('PotentialEnergy')
plt.savefig("June20-am-KE.png")
plt.close()

Volume = sim_data["mdcomputeThermodynamicQuantitiesvolume"]
Time = sim_data["flowermdbasesimulationSimulationtimestep"]
x= Time
y= Volume
plt.plot(x,y/128)
plt.xlabel('Time')
plt.ylabel('Volume')
plt.savefig("June20-am-Vol.png")
plt.close()

Temp = sim_data["mdcomputeThermodynamicQuantitieskinetic_temperature"]
Time = sim_data["flowermdbasesimulationSimulationtimestep"]
x= Time
y= Temp
plt.plot(x,y/128)
plt.xlabel('Time')
plt.ylabel('Temp')
plt.savefig("June20-am-Temp.png")
plt.close()

