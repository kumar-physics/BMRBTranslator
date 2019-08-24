#
# Number of structures.
#
nstructures = 100

#
# Read PSF file(s).
#
import protocol
protocol.initStruct('protein.psf')

#
# Load paramaters.
# ---
# Read covalent and nonbonded parameters from parameter file(s).
# Note that only covalent parameters for bond lengths, bond angles, and
# impropers dihedrals are used in this script.  The torsion angle parameters,
# if any, are ommited because they are provided by a statistical potential
# below.
#
protocol.initParams('protein')


#
# Set random seed.
#
protocol.initRandomSeed(3421)   # by seed number

#
# Generate extended structure.
#
protocol.genExtendedStructure() 

#
# Create a potList.PotList() to contain the energy terms that will be active
# during structure calculations.
#
from potList import PotList
etotal = PotList()

#
# Lists highTempParams and rampedParams will hold simulationTools.StaticRamp and
# simulationTools.MultRamp objects to handle (i) changes in both energy scales
# and atomic repulsion setups (e.g., CA-only vs. all-atom interactions)
# between the high temp. and simulated annealing stages, and (ii) energy scale
# ramping withing the annealing stage itself.
#
from simulationTools import MultRamp, StaticRamp, InitialParams

highTempParams = []
rampedParams = []


# Next, the entire setup of each energy term (including their treatment during
# the high temperature and annealing stages) is performed in self-contained
# sections, so that removal of a term or addition of a new one can be done
# simply by commenting out or adding the corresponding section, respcetively.


# #
# # Set up RDC potential.
# #
#
# # Orientation tensor(s).
# # ---
# # Define one tensor per medium.
# # For each medium, specify an arbitrary name for the medium (a string), and
# # initial values of the tensor's magnitude (Da) and rhombicity (Rh) (e.g, as
# # estimated from the powder pattern approach).
# #
# from varTensorTools import create_VarTensor, calcTensorOrientation
# tensors = {}
# #                          medium     Da   Rh
# for (medium, Da, Rh) in [('tmv107',  -6.5, 0.62),
#                          ('bicelle', -9.9, 0.23)]:
#     tensor = create_VarTensor(medium)
#     tensor.setDa(Da)
#     tensor.setRh(Rh)
#     tensors[medium] = tensor
#
# highTempParams.append(StaticRamp("""
# for medium in tensors.values():
#     calcTensorOrientation(medium)
# """) )
#
# # List with RDC data.
# # ---
# # Each item in the list will be used to create an RDC energy term.  Each item
# # is a tuple that contains (in this order): the medium's name (as defined above
# # in the tensor set up), an arbitrary name for the experiment (a string, e.g.,
# # "NH" for N-H RDCs), the path of the corresponding restraint table (a string),
# # and a relative scale factor.  Here, the input RDCs have been normalized
# # relative to the N-H spin pair; the relative scale factor is the square of the
# # inverse of the normalization factor.
# #            medium    expt.   restraint file     scale
# rdcData = [('tmv107',  'NH' , 'tmv107_nh.tbl',    1),
#            ('tmv107',  'NCO', 'tmv107_nc.tbl',    0.05),
#            ('tmv107',  'HNC', 'tmv107_hnc.tbl',   0.108),
#            ('bicelle', 'NH' , 'bicelles_nh.tbl',  1),
#            ('bicelle', 'NCO', 'bicelles_nc.tbl',  0.05),
#            ('bicelle', 'HNC', 'bicelles_hnc.tbl', 0.108)]
# # RDC potential per se.
# from rdcPotTools import create_RDCPot, scale_toNH
# rdcs = PotList('rdc')
# for (medium, exp, table, scale) in rdcData:
#     name = '%s_%s' % (exp, medium)
#     rdc = create_RDCPot(name, table, tensors[medium])
#     rdc.setScale(scale)
#     # scale_toNH(rdc) # uncomment if unnormalized restraints
#     rdcs.append(rdc)
# etotal.append(rdcs)
# rampedParams.append(MultRamp(0.01, 1.0, "rdcs.setScale(VALUE)"))
#

#
# Set up distance restraint potential (e.g., from NOEs).
#
import noePotTools      
noe = noePotTools.create_NOEPot(name='noe', file='noe.tbl')
etotal.append(noe)
rampedParams.append(MultRamp(2, 30, "noe.setScale(VALUE)"))

#
# Set up torsion angle restraint potential (e.g., from J-couplings).
#
from xplorPot import XplorPot
dihedralTable = 'dihedral.tbl'
protocol.initDihedrals(dihedralTable)
etotal.append(XplorPot('CDIH'))
highTempParams.append(StaticRamp("etotal['CDIH'].setScale(10)"))
rampedParams.append(StaticRamp("etotal['CDIH'].setScale(200)"))


#
# Set up statistical backbone H-bond potential (HBDB).
#
protocol.initHBDB()
etotal.append(XplorPot('HBDB'))
#
# Set up statistical torsion angle potential (torsionDB).
#
import torsionDBPotTools
torsiondb = torsionDBPotTools.create_TorsionDBPot(name='torsiondb',
                                                  system='protein')
etotal.append(torsiondb)
rampedParams.append(MultRamp(0.002, 2, "torsiondb.setScale(VALUE)"))


#
# Setup interatomic repulsion potential (van der Waals-like term).
#
from repelPotTools import create_RepelPot, initRepel
repel = create_RepelPot('repel')
etotal.append(repel)
# CA-only interactions with large atomic radius to improve sampling.
highTempParams.append( StaticRamp("""initRepel(repel,
                                               use14=True,
                                               scale=0.004,
                                               repel=1.2,
                                               moveTol=45,
                                               interactingAtoms='name CA'
                                               )""") )
# All interactions active with more realistic atomic radii.
rampedParams.append(StaticRamp("initRepel(repel, use14=False)"))
rampedParams.append(MultRamp(0.004, 4, "repel.setScale(VALUE)"))

# Selected 1-4 interactions.
import torsionDBPotTools
repel14 = torsionDBPotTools.create_Terminal14Pot('repel14')
etotal.append(repel14)
highTempParams.append(StaticRamp("repel14.setScale(0)"))
rampedParams.append(MultRamp(0.004, 4, "repel14.setScale(VALUE)"))


# Bond Length Energy Term
etotal.append(XplorPot('BOND'))

# Bond Angle Energy Term
etotal.append(XplorPot('ANGL'))
rampedParams.append(MultRamp(0.4, 1.0, "etotal['ANGL'].setScale(VALUE)"))


# Improper Dihedral Angle Energy Term
etotal.append(XplorPot('IMPR'))
rampedParams.append(MultRamp(0.1, 1.0, "etotal['IMPR'].setScale(VALUE)"))

#
# Done with energy terms.


#
# Set up IVM object(s).
# ---
# The IVM module is used for performing dynamics and minimization in torsion-
# angle and in Cartesian space.
#

# IVM object for torsion-angle dynamics/minimization.
#
from ivm import IVM
dyn = IVM()

# Alignment tensor setup - fix tensor Rh and Da, vary orientation.
# for tensor in tensors.values():
#     tensor.setFreedom("fixDa, fixRh")

protocol.torsionTopology(dyn)



# IVM object for final Cartesian minimization.
#
minc = IVM()

# Alingment tensor setup - allow all tensor parameters to float.
# for tensor in tensors.values():
#     tensor.setFreedom("varyDa, varyRh")
    
protocol.cartesianTopology(minc)


#
# Give atoms uniform weights, except for pseudoatoms
#
protocol.massSetup()



#
# Temperature set up.
#
temp_ini = 3500.0   # initial temperature
temp_fin = 25.0     # final temperature


def calcOneStructure(loopInfo):
    """Calculate a single structure.

    """
    # Randomize torsion angles.
    from monteCarlo import randomizeTorsions
    randomizeTorsions(dyn)

    # Set torsion angles from restraints.
    # (They start satisfied, allowing the shortening of high temp dynamics.)
    import torsionTools 
    torsionTools.setTorsionsFromTable(dihedralTable)

    # The torsion restraints may include ring torsions and distort geometry.
    protocol.fixupCovalentGeom(maxIters=100, useVDW=True)


    #
    # High Temperature Dynamics Stage.
    #

    # Initialize parameters for high temperature dynamics.
    InitialParams(rampedParams)
    InitialParams(highTempParams) 
                            

    # Set up IVM object and run.
    protocol.initDynamics(dyn,
                          potList=etotal, 
                          bathTemp=temp_ini,
                          initVelocities=True,
                          finalTime=100,   
                          numSteps=1000,   
                          printInterval=100)

    dyn.setETolerance(temp_ini/100)# used to set step size (default: temp/1000) 

    dyn.run()

    #
    # Simulated Annealing Stage.
    #

    # Set up IVM object for annealing.
    protocol.initDynamics(dyn,
                          finalTime=0.2,  
                          numSteps=100,       
                          printInterval=100)



    # Set up cooling loop and run.
    from simulationTools import AnnealIVM
    AnnealIVM(initTemp=temp_ini,
              finalTemp=temp_fin,
              tempStep=12.5,
              ivm=dyn,
              rampedParams=rampedParams).run()
              
    #
    # Torsion angle minimization.
    #
    protocol.initMinimize(dyn,
                          printInterval=50)
    dyn.run()


    #
    # Cartesian minimization.
    #
    protocol.initMinimize(minc,
                          potList=etotal,
                          dEPred=10)
    minc.run()


#
# Calculate all structures and perform analysis.
#
from simulationTools import StructureLoop
StructureLoop(numStructures=nstructures,
              structLoopAction=calcOneStructure,
              doWriteStructures=True,
              # Arguments for generating structure statistics:
              genViolationStats=True,
              averageSortPots=[etotal['BOND'], # terms for structure sorting.
                               etotal['ANGL'], 
                               etotal['IMPR'], noe, etotal['CDIH']],
              #averageSortPots=[etotal['BOND'], # terms for structure sorting.
              #                 etotal['ANGL'],
              #                 etotal['IMPR'], noe, rdcs, etotal['CDIH']],
              averageTopFraction=0.1, # top fraction of structs. to report on.
              averagePotList=etotal, # terms analyzed.
              averageFitSel='not (name H* or PSEUDO)', # selection to fit...
              ).run()                               # to average structure...
                                                    # and report precision.

