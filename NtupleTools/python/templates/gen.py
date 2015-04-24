'''

Ntuple branch template sets for gen level candidates

Author: Devin N. Taylor

'''

from FinalStateAnalysis.Utilities.cfgtools import PSet

kinematics = PSet(
    objectPt = '{object}.pt',
    objectEta = '{object}.eta',
    objectAbsEta = 'abs({object}.eta)',
    objectPhi = '{object}.phi',
    objectCharge = '{object}.charge',
    objectMass = '{object}.mass',
    objectVZ = '{object}.vz',
    objectEnergy = '{object}.energy'
)

geninfo = PSet(
    objectPdgId = '{object}.pdgId',
)
