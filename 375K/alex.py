import mbuild as mb
import foyer
import hoomd
import numpy as np
from hoomd import deprecated
from hoomd import md
from hoomd import group

class CH2(mb.Compound):
    """Defines a united-atom CH2 particle """
    def __init__(self):
        super(CH2, self).__init__()
        self.add(mb.Particle(name='_CH2'))

        self.add(mb.Port(anchor=self[0], orientation=[0, 1, 0],
                         separation=0.075), 'up')
        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0],
                         separation=0.075), 'down')

class CH3(mb.Compound):
    """Defines a united-atom CH3 particle """
    def __init__(self):
        super(CH3, self).__init__()
        self.add(mb.Particle(name='_CH3'))

        self.add(mb.Port(anchor=self[0], separation=0.075), 'up')


class Hydroxyl(mb.Compound):
    """Defines a hydroxyl group """
    def __init__(self):
        super(Hydroxyl, self).__init__()
        #self.add(mb.Particle(name='_OH'))
        oxygen = mb.Particle(name='O', pos=[0.0, 0.0, 0.0])
        hydrogen = mb.Particle(name='H', pos=[0.0945, 0.0, 0.0])
        self.add((oxygen, hydrogen))
        self.add_bond((oxygen, hydrogen))
        self.add(mb.Port(anchor=oxygen, orientation=[-1, 0, 0], separation=0.075),
                 label='up')


class Alcohol(mb.Compound):
    """A united-atom alcohol chain of variable length """
    def __init__(self, chainlength):
        """Initialize a united-atom alcohol chain

        Parameters
        ----------
        chainlength : int
            Number of *carbons* in the chain backbone. This does not include
            the OH group

        """
        super(Alcohol, self).__init__()

        hydroxyl = Hydroxyl()
        ch3 = CH3()
        self.add(hydroxyl)
        self.add(ch3)

        if chainlength < 1:
            raise ValueError('Chain length must be greater than one!')
        elif chainlength == 1:
            mb.force_overlap(ch3, ch3['up'], hydroxyl['up'])
        else:
            chain_middle = mb.recipes.Polymer(CH2(), n=chainlength-1,
                                      port_labels=('up', 'down'))
            self.add(chain_middle)
            mb.force_overlap(chain_middle, chain_middle['up'], hydroxyl['up'])
            mb.force_overlap(ch3, ch3['up'], chain_middle['down'])




methanol = Alcohol(chainlength = 1)
meth_box = mb.packing.fill_box(compound= methanol, n_compounds=2000, box=mb.Box(lengths=[8,8,8]))



hoomd.context.initialize("")
box_arr = [8,8,8]
wrap_xyz = meth_box.xyz - 1*np.floor_divide(meth_box.xyz, box_arr) * box_arr
meth_box.xyz = wrap_xyz
struc =meth_box.to_parmed()
meth_box.save('meth_box.hoomdxml', forcefield_name='trappe-ua', overwrite=True)
ff = foyer.Forcefield(name='trappe-ua')
param = ff.apply(struc)
#for atom in param.atoms:
    #atom.charge = 0

from mbuild.formats.hoomd_simulation import create_hoomd_simulation
create_hoomd_simulation(param, ref_distance=10, ref_energy=1/4.184)
nl = hoomd.md.nlist.cell()
nl.reset_exclusions(['1-2', '1-3', '1-4'])

all = hoomd.group.all()

#fire = md.integrate.mode_minimize_fire(dt=0.002)
#nve = md.integrate.nve(group=all, limit = 0.01)
#while not fire.has_converged():
#    print([a.get_energy(all) for a in hoomd.context.current.forces])
#    hoomd.run(1000)
#nve.disable()

md.integrate.mode_standard(dt=0.002)
nvt = md.integrate.nvt(group=all, kT=3.118, tau=1)
hoomd.analyze.log(filename="analyze2.log",
                  quantities=['temperature', 'potential_energy','kinetic_energy', 'bond_harmonic_energy', 'pppm_energy', 'pair_lj_energy', 'angle_harmonic_energy','volume', 'pressure'],
                  period=100,
                  overwrite=True)
dump = hoomd.dump.dcd('nvt.dcd', period = 100, overwrite = True)
hoomd.run(1e5)
dump.disable()
print("SWITCHING TO NPT")
nvt.disable()
npt = md.integrate.npt(group=all, kT=3.118, tau=1, P=0.62, tauP=10)
hoomd.analyze.log(filename="analyze.log",
                  quantities=['temperature', 'potential_energy','kinetic_energy', 'bond_harmonic_energy', 'pppm_energy', 'pair_lj_energy', 'angle_harmonic_energy','volume','pressure'],
                  period=100,
                  overwrite=True)
hoomd.dump.dcd('npt.dcd', period = 100, overwrite = True)
hoomd.run(1e5)

#lattice = hoomd.init.read_gsd(filename="./oct_box.gsd")
#lattice = deprecated.init.read_xml(filename='./meth_box.hoomdxml')
#
##nl = hoomd.md.nlist.cell(r_buff=2, check_period=1)
#nl = hoomd.md.nlist.cell()
#
#b = md.bond.harmonic()
#b.bond_coeff.set('H-O', k=502416, r0=0.0945)
##Ch2sp3 = md.bond.harmonic(name="CH2sp3")
#b.bond_coeff.set('CH3_O-O', k=502416, r0=0.143)
#
#a1 = md.angle.harmonic()
#a1.angle_coeff.set('CH3_O-O-H', k=460.6212, t0=1.89368)
#
#
#lj = md.pair.lj(r_cut=2, nlist=nl)
#
#lj.pair_coeff.set('H', 'H', epsilon=0, sigma=1)
#lj.pair_coeff.set('O', 'O', epsilon=0.7732, sigma=0.3020)
#lj.pair_coeff.set('CH3_O', 'CH3_O', epsilon=0.8148, sigma=0.3750)
#
#lj.pair_coeff.set('O', 'H', epsilon=0, sigma=0.651)
#lj.pair_coeff.set('H', 'CH3_O', epsilon=0, sigma=0.6875)
#
#lj.pair_coeff.set('O', 'CH3_O', epsilon=0.7937275099, sigma=0.3385)
#
#
#charged = group.charged()
#pppm = md.charge.pppm(group=charged, nlist = nl)
#pppm.set_params(Nx=5, Ny=5, Nz=5, order=1, rcut=2)
'''

import random
random.seed(3.9241564195)
T_init = 0.01

for p in lattice.particles:
    p.velocity = (random.gauss(0, T_init), random.gauss(0, T_init), random.gauss(0, T_init))

md.integrate.mode_standard(dt=0.005)
fire = md.integrate.mode_minimize_fire(dt=0.05)
all = group.all()
nve = md.integrate.nve(group=all, limit = 0.01)

while not(fire.has_converged()):
    hoomd.run(100)

hoomd.analyze.log(filename="analyze.log",
                  quantities=['temperature', 'potential_energy','kinetic_energy', 'bond_harmonic_energy', 'pppm_energy', 'pair_lj_energy', 'angle_harmonic_energy'],
                  period=100,
                  overwrite=True)

#dist =md.constrain.distance()
hoomd.dump.dcd('nve.dcd', period = 100, overwrite = True)
hoomd.run(2e4)

nve.disable()
nve2 = md.integrate.nve(group=all, limit = 0.1)
hoomd.dump.dcd('nve2.dcd', period = 100, overwrite = True)
hoomd.analyze.log(filename="analyze1.log",
                  quantities=['temperature', 'potential_energy','kinetic_energy', 'bond_harmonic_energy', 'pppm_energy', 'pair_lj_energy', 'angle_harmonic_energy'],
                  period=100,
                  overwrite=True)
hoomd.run(1e5)
nve2.disable()

nve3 = md.integrate.nve(group=all, limit = 1)
hoomd.dump.dcd('nve3.dcd', period = 100, overwrite = True)
hoomd.run(2e4)
nve3.disable()

npt = md.integrate.nvt(group=all, kT=1.5, tau=1)
hoomd.analyze.log(filename="analyze2.log",
                  quantities=['temperature', 'potential_energy','kinetic_energy', 'bond_harmonic_energy', 'pppm_energy', 'pair_lj_energy', 'angle_harmonic_energy'],
                  period=100,
                  overwrite=True)
hoomd.dump.dcd('npt.dcd', period = 100, overwrite = True)
hoomd.run(1e5)

'''
