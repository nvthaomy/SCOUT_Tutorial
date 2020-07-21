# -*- coding: utf-8 -*-
#!/usr/bin/env python2
"""
Created on Wed May 29 09:29:13 2019

@author: My Nguyen 

Simulating Gaussian polymer solution in NVT/NPT ensemble
With or without the external potential

Flow of the script

Provide inputs in:
    0) define energy, length, mass scales
    1) dimensionless parameters defining system topology and interactions
       d) only if apply an sinusoidal external potential on the polymers, set Uext to 0 to disable this
    2) integration options
    3) set up reports

4)Convert dimensionless parameters to real units based on the scales provided in 0) since openmm only handle real units
5)-9) Create system and interactions
10)-12) Simulation
"""
#import mdtraj as md
from simtk import openmm, unit
from simtk.unit import *
from simtk.openmm import app
import numpy as np
import simtk.openmm.app.topology as topology

#========================================
#0) DEFINE ENERGY, LENGTH & MASS SCALES###
#========================================
epsilon  = 1 * kilojoules_per_mole
sigma = 0.1 * nanometer
mass = 12 * amu
N_av = constants.AVOGADRO_CONSTANT_NA 
kb = constants.BOLTZMANN_CONSTANT_kB

#====================================
#1) SYSTEM DIMENSIONLESS PARAMETERS###
#====================================
#a)Molecules:
NA = 1000
NB = 1000
N = NA + NB

boxLs = [10.,10.,20.]

#c)Pair potential: u = B exp(-Kr^2)
K = 0.25
#interaction energy scale for each pair type
#polymer-polymer
B_AA =  1.0
#solvent-solvent
B_BB =  1.0
#polymer-solvent
B_AB =  1.0
reduced_Cutoff = 2.5/np.sqrt(K) 
nonbondedMethod = openmm.CustomNonbondedForce.CutoffPeriodic
#whether to add correction to the energy beyond nonbonded cutoff
UseLongRangeCorrection = False

#d) External potential: applying sinusoidal external potential on bead A, set Uext = 0 to disable 
Uext = 0.5 #amplitude of Uext in kBT
Nperiod = 1
axis  = 2 #0:x, 1:y, 2:z
reduced_planeLoc = 0 

#======================
#2) Integration Options
#======================
ensemble = 'NVT' # 'NVT' or 'NPT'
reduced_timestep = 0.05
reduced_temp = 1
reduced_Tdamp = 100*reduced_timestep #time units
reduced_pressure = 1
reduced_Pdamp = 1000*reduced_timestep #time units

steps = 10000/reduced_timestep
equilibrationSteps = 1000/reduced_timestep
stride = 100
#if platform is not set, openmm will try to select the fastest available Platform
platformName = None #'CPU','CUDA','OpenCL'
platformProperties = None

#platformProperties = {'Precision': 'mixed'}

#=====================
#3) Reports
#=====================
#set up data and trajectory reports:
#if use mdtraj's reporter:
#dcdReporter = mdtraj.reporters.DCDReporter(traj, 5000)

#if don't have mdtraj, can use openmm's reporter to write out trajectory in pdb or dcd format
dcdReporter = app.dcdreporter.DCDReporter('trajectory.dcd',stride,enforcePeriodicBox=True)
#pdbReporter = app.pdbreporter.PDBReporter( 'trajectory.dcd', 5000)
dataReporter = app.statedatareporter.StateDataReporter('log.txt', stride, totalSteps=steps,
    step=True, speed=True, potentialEnergy=True, kineticEnergy=True, 
    totalEnergy=True, temperature=True, volume=True, density=True, separator='\t')

dcdReporter_warmup = app.dcdreporter.DCDReporter('trajectory_warmup.dcd',stride,enforcePeriodicBox=True)
dataReporter_warmup= app.statedatareporter.StateDataReporter('log_warmup.txt', stride, totalSteps=equilibrationSteps,
    step=True, speed=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, volume=True, density=True, separator='\t')

#=============================
#4) Converting to real units
#=============================
boxLs = np.array(boxLs) * sigma
planeLoc =  reduced_planeLoc*sigma
nonbondedCutoff = reduced_Cutoff*sigma
dt = reduced_timestep* (mass*sigma*sigma/epsilon)**(1/2)
temperature = reduced_temp * epsilon/kb
friction = 1/(reduced_Tdamp) * (epsilon/(mass*sigma*sigma))**(1/2)
pressure = reduced_pressure * epsilon/(sigma**3) * N_av**-1
barostatInterval = int(reduced_Pdamp/reduced_timestep)
print ('\n== Parameters in real units ==')
print ('temperature:{}'.format(temperature))
print ('pressure:{}'.format(pressure))
print ('friction:{}'.format(friction))
print ('barostat interval:{}'.format(barostatInterval))
print ('nonbondedCutoff: {}'.format(nonbondedCutoff))

#========================================
#5) Create a system and add particles to it
#========================================
system = openmm.System()

# Particles are added one at a time
# Their indices in the System will correspond with their indices in the Force objects we will add later
print ("\nAdding {} A into system".format(NA))
for index in range(NA):
    system.addParticle(mass)
print ("Adding {} B into system".format(NB))
for index in range(NB):
    system.addParticle(mass)

# Set the periodic box vectors:
print ('Box dimensions {}'.format(boxLs))
box_vectors = np.diag([L/angstrom for L in boxLs]) * angstroms
system.setDefaultPeriodicBoxVectors(*box_vectors)

#==================
#6) Create topology
#==================
#Topology consists of a set of Chains 
#Each Chain contains a set of Residues, 
#and each Residue contains a set of Atoms.
nmols = [NA, NB]
residues = [["A"], ["B"]]
atomList = {"A": ['A'], "B":['B']} #list of atoms in each residue
elements = {"A":app.element.Element(200, 'A', 'AA', mass),
            "B":app.element.Element(201, 'B', 'BB', mass)}
def makeTop(NA,NB):
    nmols = [NA, NB]
    top = topology.Topology()
    for ispec in range(len(nmols)): #loop over each species
        for imol in range(nmols[ispec]): #loop over each molecule within species
            chain = top.addChain() #create chain for each molecule
            for ires in range( len(residues[ispec]) ): #polymer and solvent each has one residue
                resname = residues[ispec][ires]
                res = top.addResidue( resname, chain)
                atoms = atomList[resname]
                for atomInd,atomName in enumerate(atoms):
                    el = elements[atomName]
                    if atomInd > 0:
                        previousAtom = atom
                    atom = top.addAtom( atomName, el, res )
                    if atomInd > 0:
                        top.addBond(previousAtom,atom)
    return top
print ("\nCreating topology")
top = makeTop(NA,NB)
#the following line ensures the output trajectory has unit cell info
top.setPeriodicBoxVectors(system.getDefaultPeriodicBoxVectors())

#================================
#8) create custom nonbonded force
#================================
#gaussianFunc = 'B*exp(-K*r^2)'
#for each pair interaction type, need to add ALL atoms in the system to the force object, only paricles in the InteractionGroup will interact with each other
Ai = set()
Bi = set()
for atom in top.atoms():
	if atom.residue.name in ['A']:
		Ai.add(atom.index)
	elif atom.residue.name in ['B']:
		Bi.add(atom.index)
all_atoms = Ai.union(Bi)

#A-A:
AA_nonbondedForce = openmm.CustomNonbondedForce('BAA*exp(-K*r^2)')
AA_nonbondedForce.setNonbondedMethod(nonbondedMethod)
AA_nonbondedForce.addGlobalParameter('K',K/(sigma*sigma)) #length^-2
AA_nonbondedForce.addGlobalParameter('BAA',B_AA*epsilon) #energy/mol
AA_nonbondedForce.addInteractionGroup(Ai,Ai)
for i in range(system.getNumParticles()):
	AA_nonbondedForce.addParticle()
AA_nonbondedForce.setCutoffDistance(nonbondedCutoff)
AA_nonbondedForce.setUseLongRangeCorrection(UseLongRangeCorrection)
system.addForce(AA_nonbondedForce)
print ("\nNumber of particles in AA nonbonded force:{}".format(AA_nonbondedForce.getNumParticles()))
#B-B:
BB_nonbondedForce = openmm.CustomNonbondedForce('BBB*exp(-K*r^2)')
BB_nonbondedForce.setNonbondedMethod(nonbondedMethod)
BB_nonbondedForce.addGlobalParameter('K',K/(sigma*sigma)) #length^-2
BB_nonbondedForce.addGlobalParameter('BBB',B_BB*epsilon) #energy/mol
BB_nonbondedForce.addInteractionGroup(Bi,Bi)
for i in range(system.getNumParticles()):
        BB_nonbondedForce.addParticle()
BB_nonbondedForce.setCutoffDistance(nonbondedCutoff)
BB_nonbondedForce.setUseLongRangeCorrection(UseLongRangeCorrection)
system.addForce(BB_nonbondedForce)
print ("Number of particles in BB nonbonded force:{}".format(BB_nonbondedForce.getNumParticles()))

#A-B:
AB_nonbondedForce = openmm.CustomNonbondedForce('BAB*exp(-K*r^2)')
AB_nonbondedForce.setNonbondedMethod(nonbondedMethod)
AB_nonbondedForce.addGlobalParameter('K',K/(sigma*sigma)) #length^-2
AB_nonbondedForce.addGlobalParameter('BAB',B_AB*epsilon) #energy/mol
AB_nonbondedForce.addInteractionGroup(Ai,Bi)
for i in range(system.getNumParticles()):
        AB_nonbondedForce.addParticle()
AB_nonbondedForce.setCutoffDistance(nonbondedCutoff)
AB_nonbondedForce.setUseLongRangeCorrection(UseLongRangeCorrection)
system.addForce(AB_nonbondedForce)
print ("Number of particles in AB nonbonded force:{}".format(AB_nonbondedForce.getNumParticles()))

#=======================================
#9) create external potential on bead A
#=======================================
external={"U":Uext*epsilon,"NPeriod":Nperiod,"axis":axis ,"planeLoc":planeLoc}
direction=['x','y','z']
ax = external["axis"]
if external["U"] > 0.0 * kilojoules_per_mole:
    print('Creating sinusoidal external potential in the {} direction'.format(direction[axis]))
    energy_function = 'U*sin(2*pi*NPeriod*({axis}-r0)/L)'.format(axis=direction[ax])
    fExt = openmm.CustomExternalForce(energy_function)
    fExt.addGlobalParameter("U", external["U"])
    fExt.addGlobalParameter("NPeriod", external["NPeriod"])
    fExt.addGlobalParameter("pi",np.pi)
    fExt.addGlobalParameter("r0",external["planeLoc"])
    fExt.addGlobalParameter("L",boxLs[ax])
    for i in Ai:
        fExt.addParticle( i,[] )
    system.addForce(fExt)

#==============================
#10) Prepare the Simulation ##
#==============================
if ensemble == 'NPT':
    system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = openmm.LangevinIntegrator(temperature, friction, dt)
#use platform if specified, otherwise let omm decide
if platformName != None:
    platform = openmm.Platform.getPlatformByName(platformName)
    if platformProperties:
        simulation = app.Simulation(top,system, integrator, platform,platformProperties)
    else:
        simulation = app.Simulation(top,system, integrator, platform)
else:
    simulation = app.Simulation(top,system, integrator)

#initialize positions
positions = [boxLs[0]/angstrom,boxLs[1]/angstrom,boxLs[2]/angstrom]*angstrom * np.random.rand(N,3)
simulation.context.setPositions(positions)
#write initial positions to pdb
initialpdb = 'trajectory_init.pdb'
app.PDBFile.writeModel(simulation.topology, positions, open(initialpdb,'w'))

#Restart and Check point:
#to load a saved state: simulation.loadState('output.xml') or simulation.loadCheckpoint('checkpnt.chk')
simulation.reporters.append(app.checkpointreporter.CheckpointReporter('checkpnt.chk', 1000))

#============================
#11) Minimize and Equilibrate
#============================
print ("\nInitialize {} simulation".format(ensemble))
print ("Running energy minimization")
simulation.minimizeEnergy(maxIterations=1000)
simulation.context.setVelocitiesToTemperature(temperature*3)
print ("Running equilibration")
simulation.reporters.append(dataReporter_warmup)
simulation.reporters.append(dcdReporter_warmup)
simulation.step(equilibrationSteps)
simulation.saveState('output_warmup.xml')

#=============
#12) Simulate
#=============
print ("Running production")
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.currentStep = 0
simulation.step(steps)
simulation.saveState('output.xml')

