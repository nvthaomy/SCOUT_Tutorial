import sys, os
sys.path.insert(1,'/home/kshen/mylib/sim.20190727')
print(sys.path)
import sim, pickleTraj
import mappoly
#import system, optimizer
import numpy as np

units = sim.units.DimensionlessUnits

# === 0) Define topology ===
atom_names = ['A','B']
atom_charges = {'A':0., 'B':0.}

molecule_names = ['A','B']
molecule_defs = {'A':'A', 'B':'B'}

# note: more work is needed to add in bonding information

# === 1) Define system we want to simulate ===
sys_name = 'ABfluid'
num_mols = {'A':1000, 'B':1000}
num_mols = {'A':3375, 'B':3375}

interactions = { 
        ('A','A'):[ {'type':'LJG','Label':'AA','B':1.,'Kappa':0.25,'Dist0':0.,'Sigma':1.,'Epsilon':0.,'Cut':5., 'fixB':True} ],
        ('B','B'):[ {'type':'LJG','Label':'BB','B':1.,'Kappa':0.25,'Dist0':0.,'Sigma':1.,'Epsilon':0.,'Cut':5., 'fixB':True} ],
        ('A','B'):[ {'type':'LJG','Label':'AB','B':1.,'Kappa':0.25,'Dist0':0.,'Sigma':1.,'Epsilon':0.,'Cut':5., 'fixB':False} ],
        ('A'):[ {'type':'external_sine','Label':'Uext','UConst':0.5,'NPeriods':1,'PlaneAxis':0,'PlaneLoc':0.} ]
        }

default_fix_params = True


#integration parameters
dt = 0.1
temp = 1.
langevin_gamma = 1./(100.*dt)

steps_equil = 100./dt
steps_prod = 10000./dt
steps_stride = 10./dt

#srel parameters
use_openmm = False
use_lammps = False
sim.export.omm.platformName = 'OpenCL'
sim.export.omm.device = -1 #-1 is default, let openmm choose its own platform. otherwise is GPU device #
sim.export.omm.NPairPotentialKnots = 500 #number of points used to spline-interpolate the potential
sim.export.omm.InnerCutoff = 0.001 #0.001 is default. Note that a small value is not necessary, like in the lammps export, because the omm export interpolates to zero
sim.srel.optimizetrajomm.OpenMMStepsMin = 0 #number of steps to minimize structure, 0 is default
sim.srel.optimizetrajomm.OpenMMDelTempFiles = False #False is Default
sim.export.omm.UseTabulated = True

#other inputs
datadir='.'
traj_file = os.path.abspath('{}/L30-15-15_du0_Uext0.0.lammpstrj'.format(datadir))
traj = pickleTraj(traj_file)
box = traj.FrameData['BoxL']

# === 2) Create system ===
sim_atom_types = {}
sim_mol_types = []
# 2a) define atoms
for atom_name in atom_names:
    charge = atom_charges[atom_name]
    atom_type = sim.chem.AtomType(atom_name, Mass=1., Charge = charge)
    sim_atom_types[atom_name] = atom_type
print('atom types: {}'.format(sim_atom_types))
# 2b) define molecules
for mol_name in molecule_names:
    if num_mols[mol_name]>0:
        atoms_in_mol = []
        
        for atom_name in molecule_defs[mol_name]:
            atoms_in_mol.append( sim_atom_types[atom_name] )

        mol_type = sim.chem.MolType(mol_name, atoms_in_mol)
        sim_mol_types.append(mol_type)
print('molecule types: {}'.format(sim_mol_types))
# 2c) create 'World' and system
World = sim.chem.World(sim_mol_types, Dim=3, Units=units)
Sys = sim.system.System(World, Name = sys_name)
Sys.BoxL = box
print('Setting box: {}'.format(Sys.BoxL))

for i, mol_type in enumerate(sim_mol_types):
    num_mol = num_mols[mol_type.Name]

    print("Adding {} {} molecules to system".format(num_mol, mol_type.Name))
    for j in range(num_mol):
        Sys += mol_type.New()


# === 3) Create force field ===
filters={}
"""
for ii,name1 in enumerate(atom_names):
    filters[{name1}] = sim.atomselect.PolyFilter([sim_atom_types[name1]])
    for jj in range(ii,len(atom_names)):
        name2 = atom_names[jj]
        filters[{name1,name2}] = sim.atomselect.PolyFilter([sim_atom_types[name1],sim_atom_types[name2]])
"""

filters={}
force_field = []
for atom_tuple,interaction_list in interactions.items():
    for params in interaction_list:
        if params['type'].lower() in ['ljg']: #Lennard Jones Gaussian
            f = sim.atomselect.PolyFilter( [sim_atom_types[atom_tuple[0]],sim_atom_types[atom_tuple[1]]] )
            filters[atom_tuple] = f
            print('Adding LJG for {}'.format(atom_tuple))
            p = sim.potential.LJGaussian( Sys, Filter = f, Cut = params['Cut'], Fixed = default_fix_params,
                    B = params['B'], Kappa = params['Kappa'], Dist0 = params['Dist0'], Sigma = params['Sigma'], Epsilon = params['Epsilon'], Label = params['Label'] )
            p.Param.B.Fixed = params['fixB']
            force_field.append(p)
        if params['type'].lower() in ['extsin','ext_sin','external_sinusoid','external_sine']:
            f = sim.atomselect.PolyFilter( [sim_atom_types[atom_tuple[0]]] )
            filters[atom_tuple] = f
            print('Adding external sinusoid for {}'.format(atom_tuple))
            p = sim.potential.ExternalSinusoid(Sys, Filter = f, UConst=params['UConst'], NPeriods=params['NPeriods'], PlaneAxis=params['PlaneAxis'], PlaneLoc=params['PlaneLoc'], Label=params['Label'])
            force_field.append(p)

print('Force Field Summary: {}'.format(force_field))

# === 4) Set other settings ===
Sys.Load()
Int = Sys.Int
Int.Method = Int.Methods.VVIntegrate
Int.Method.Thermostat = Int.Method.ThermostatLangevin
Int.Method.Langevin = langevin_gamma
Int.Method.TimeStep = dt
Sys.TempSet = temp

# === 5) Create optimizer ===
print('Making optimizer')
# === create a mapping object (1:1, since actual mapped traj created earlier) ===
Map = sim.atommap.PosMap()
for (i, a) in enumerate(Sys.Atom):
    Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]

# read in the trajectory
# set box-len same as AA traj
Sys.BoxL = traj.FrameData['BoxL']
print("Box: {}".format(Sys.BoxL))

# create optimizer

ewaldPotentials = [ P for P in Sys.ForceField if P.Names[0] == 'ewald' ]
UseTarHists = False
if not ewaldPotentials: #if no ewald
    UseTarHists = True

print("Use target histograms: {}".format(UseTarHists))
Opt = sim.srel.OptimizeTrajClass(Sys, Map, Beta = 1./Sys.TempSet, Traj = traj, FilePrefix = Sys.Name, 
                                 SaveLoadArgData = True, TempFileDir = os.getcwd(), UseTarHists = UseTarHists)

# make opt use lammps for on-the-fly MD runs
if use_openmm:
    print("Using OpenMM")
    Opt = sim.srel.UseOpenMM(Opt)
elif use_lammps:
    print("Using Lammps")
    Opt = sim.srel.UseLammps(Opt)
else:
    print("Using internal Sim engine")

# prep the optimizer object
sim.srel.optimizetraj.PlotFmt = 'png'
Opt.MinReweightFrames = 10
#Opt.MinReweightFrac = 0.15

# md iterations
Opt.StepsEquil = int(steps_equil/Sys.Int.Method.TimeStep)
Opt.StepsProd  = int(steps_prod/Sys.Int.Method.TimeStep)
Opt.StepsStride = int(steps_stride/Sys.Int.Method.TimeStep)

print('Check fixed parameters status:')
print "i,Fixed"
for pair in Sys.ForceField:
    print("--- {} ---".format(pair.Label))
    for (i,Fixed) in enumerate(pair.Param.Fixed):
        print i,Fixed


# === 6) Run optimization ===
Opt.RunConjugateGradient()

print("\n\n\n")
print("=== Hessian/DDSrel ===")
print("{}".format(Opt.DDObj))
np.savetxt("Hessian.dat",Opt.DDObj)
import pickle
#pickle.dump(Opt, open("obj_Opt.p","wb"))
pickle.dump(Sys, open("obj_Sys.p","wb"))







