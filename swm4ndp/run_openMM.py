from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import sys, inspect, os

pdb = PDBFile('h2o_220.pdb')
strdir = 'test_sim/'
if not os.path.isdir(strdir): os.mkdir(strdir)

integ_md = DrudeSCFIntegrator(0.0005*picoseconds)
integ_md.setMinimizationErrorTolerance(1e-3)

modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField('swm4ndp.xml')
modeller.addExtraParticles(forcefield)

system = forcefield.createSystem(modeller.topology, nonbondedCutoff=1.0*nanometer, constraints=None, rigidWater=True)
nbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == NonbondedForce][0]

drudeForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == DrudeForce][0]
nbondedForce.setNonbondedMethod(NonbondedForce.PME)

hyper = CustomBondForce('step(r-rhyper)*((r-rhyper)*khyper)^powh')
hyper.addGlobalParameter('khyper', 100.0)
hyper.addGlobalParameter('rhyper', 0.02)
hyper.addGlobalParameter('powh', 6)
system.addForce(hyper)

for i in range(drudeForce.getNumParticles()):
    param = drudeForce.getParticleParameters(i)
    drude = param[0]
    parent = param[1]
    hyper.addBond(drude, parent)

for i in range(system.getNumForces()):
    f = system.getForce(i)
    type(f)
    f.setForceGroup(i)

#platform = Platform.getPlatformByName('CUDA')
#platform=Platform.getPlatformByName('Reference')
platform = Platform.getPlatformByName('OpenCL')
properties = {'OpenCLPrecision': 'double'}
#simmd = Simulation(modeller.topology, system, integ_md, platform)
simmd = Simulation(modeller.topology, system, integ_md, platform, properties)
simmd.context.setPositions(modeller.positions)

simmd.step(1)
# initial energies
state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True)
position = state.getPositions()
print(str(state.getKineticEnergy()))
print(str(state.getPotentialEnergy()))
for j in range(system.getNumForces()):
    f = system.getForce(j)
    print(type(f), str(simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))

# write initial pdb file with drude oscillators
PDBFile.writeFile(simmd.topology, position, open(strdir+'start_drudes.pdb', 'w'))
simmd.reporters = []
simmd.reporters.append(DCDReporter(strdir+'md_npt_equil.dcd', 1000))
simmd.reporters.append(CheckpointReporter(strdir+'md_npt_equil.chk', 1000))
simmd.reporters[1].report(simmd,state)
simmd.context.setVelocitiesToTemperature(300)

out = open("energy.txt", 'w')
out.write("Time  Total  Kinetic  Potential\n")
out.close()
for i in range(20000):
    out = open("energy.txt", 'a+')
    simmd.step(1)
    print(i,strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    print(i,datetime.now())
    state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True)
    print(str(state.getKineticEnergy()))
    print(str(state.getPotentialEnergy()))
    kin = state.getKineticEnergy()/kilojoule_per_mole
    pot = state.getPotentialEnergy()/kilojoule_per_mole
    tot = kin + pot
    time = 0.5*i/1000
    out.write(f"{time} {tot} {kin} {pot}\n")
    out.close()
    for j in range(system.getNumForces()):
        f = system.getForce(j)
        print(type(f), str(simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))

# print equilibrated pdb file
state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True)
position = state.getPositions()
simmd.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
PDBFile.writeFile(simmd.topology, position, open(strdir+'final.pdb', 'w'))

exit()
