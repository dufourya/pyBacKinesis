"""
Library of classes defining behavioral attributes for cell agent for simulations.
The main CellBody class must contain the following methods:
'update' defines how the cell bodies behave at each time step.
'record' defines which paramters are recorded.
Subclasses of CellBody must be define to specify the parameters of different types of cell agents.
"""

"""External libaries"""
import numpy as np

rng = np.random.default_rng()

class CellBody:

    def __init__(self,CellNumber,SampleNumber,SampleInterval):
        self.CellNumber = CellNumber
        self.Records = BodyRecords(CellNumber,self.LigandNumber,SampleNumber,SampleInterval) # class to record positions and cell state

        self.Position = np.ones((CellNumber,3))
        self.RandomDiffusion = np.zeros((CellNumber,3))
        self.QuaternionOrientation = np.zeros((CellNumber,4)) #qr,qi,qj,qk
        self.QuaternionOrientation[:,0] = 1
        self.QuaternionRotation = np.zeros((CellNumber,4))
        self.Speed = np.zeros(CellNumber)
        self.OrientationXYZ = np.zeros((CellNumber,3)) #vector describing a single step of movement
        self.OrientationXYZ[:,0] = 1   
        self.Collision = np.zeros(CellNumber,dtype = np.bool)

        # Cell body state
        self.MotorG = np.zeros(CellNumber)
        self.MotorEnergyDelta = np.zeros(CellNumber)
        self.MotorState = np.zeros(CellNumber,dtype=np.bool) # 0 for ccw,1 for cw grow array when more Motors are implemented
        self.MotorHasSwitched = np.zeros(CellNumber,dtype=np.bool)
        self.InterruptedRun = np.zeros(CellNumber,dtype = np.bool)
        self.LigandConsumed = np.zeros((CellNumber,self.LigandNumber))
        self.LigandConsumedIntegrated = np.zeros((CellNumber,self.LigandNumber))

        # Cell Morphology
        self.CellSize = np.ones(CellNumber)
        self.CellPoleAge = np.zeros((CellNumber,2))
        self.MotorNumber = np.ones(CellNumber)

    def update(self,CellMind,SimulationTimeStep):
        self.update_cell_morphology(CellMind,SimulationTimeStep)
        self.update_motor_state(CellMind,SimulationTimeStep)
        self.update_swimming_state()
        self.update_position(SimulationTimeStep)
        self.update_body_orientation(SimulationTimeStep)
        
    def update_position(self,SimulationTimeStep):
        self.Position += self.Speed[:,None] * self.OrientationXYZ * SimulationTimeStep

    def randomize_orientation(self):
        u = rng.random((self.CellNumber,3))
        self.QuaternionRotation = np.vstack((np.sqrt(1-u[:,0])*np.sin(2*np.pi*u[:,1]),\
            np.sqrt(1-u[:,0])*np.cos(2*np.pi*u[:,1]),\
            np.sqrt(u[:,0])*np.sin(2*np.pi*u[:,2]),\
            np.sqrt(u[:,0])*np.cos(2*np.pi*u[:,2]))).T
        self.update_quaternion()
        self.update_orientation()

    def update_quaternion(self):
        q0 = self.QuaternionRotation[:,0]*self.QuaternionOrientation[:,0] - self.QuaternionRotation[:,1]*self.QuaternionOrientation[:,1] - self.QuaternionRotation[:,2]*self.QuaternionOrientation[:,2] - self.QuaternionRotation[:,3]*self.QuaternionOrientation[:,3]
        q1 = self.QuaternionRotation[:,1]*self.QuaternionOrientation[:,0] + self.QuaternionRotation[:,0]*self.QuaternionOrientation[:,1] + self.QuaternionRotation[:,3]*self.QuaternionOrientation[:,2] - self.QuaternionRotation[:,2]*self.QuaternionOrientation[:,3]
        q2 = self.QuaternionRotation[:,2]*self.QuaternionOrientation[:,0] - self.QuaternionRotation[:,3]*self.QuaternionOrientation[:,1] + self.QuaternionRotation[:,0]*self.QuaternionOrientation[:,2] + self.QuaternionRotation[:,1]*self.QuaternionOrientation[:,3]
        q3 = self.QuaternionRotation[:,3]*self.QuaternionOrientation[:,0] + self.QuaternionRotation[:,2]*self.QuaternionOrientation[:,1] - self.QuaternionRotation[:,1]*self.QuaternionOrientation[:,2] + self.QuaternionRotation[:,0]*self.QuaternionOrientation[:,3]
        self.QuaternionOrientation[:,0] = q0
        self.QuaternionOrientation[:,1] = q1
        self.QuaternionOrientation[:,2] = q2
        self.QuaternionOrientation[:,3] = q3

    def update_orientation(self):
        self.OrientationXYZ[:,0] = 1 - 2 * (self.QuaternionOrientation[:,2]*self.QuaternionOrientation[:,2] + self.QuaternionOrientation[:,3]*self.QuaternionOrientation[:,3])
        self.OrientationXYZ[:,1] =     2 * (self.QuaternionOrientation[:,1]*self.QuaternionOrientation[:,2] - self.QuaternionOrientation[:,0]*self.QuaternionOrientation[:,3])
        self.OrientationXYZ[:,2] =     2 * (self.QuaternionOrientation[:,0]*self.QuaternionOrientation[:,2] + self.QuaternionOrientation[:,1]*self.QuaternionOrientation[:,3])
    
    def update_cell_morphology(self,CellMind,SimulationTimeStep):
        self.LigandConsumed = SimulationTimeStep * (self.LigandConsumptionRate[None,:] * \
            (1 / (1 + self.LigandMichaelisConstant[None,:] / CellMind.LigandConcentration)) \
                - self.LigandSecretionRate[None,:])
        self.LigandConsumedIntegrated += self.LigandConsumed
        self.CellPoleAge += SimulationTimeStep
        self.CellSize += np.sum(self.LigandConsumed * self.ConversionLigandSize[None,:],axis=1)

    def record(self,Clock):
        self.Records.Time[self.Records.Iterator] = Clock
        self.Records.MotorG[:,self.Records.Iterator] = self.MotorG
        self.Records.Position[:,:,self.Records.Iterator] = self.Position
        self.Records.Speed[:,self.Records.Iterator] = self.Speed
        self.Records.OrientationXYZ[:,:,self.Records.Iterator] = self.OrientationXYZ
        self.Records.InterruptedRun[:,self.Records.Iterator] = self.InterruptedRun
        self.InterruptedRun = np.zeros(self.CellNumber,dtype=np.bool)
        self.Records.LigandConsumed[:,:,self.Records.Iterator] = self.LigandConsumed
        self.Records.LigandConsumedIntegrated[:,:,self.Records.Iterator] = self.LigandConsumedIntegrated
        self.Records.Iterator += 1

    def update_motor_state(self,CellMind,SimulationTimeStep):
        self.MotorG = (self.MotorEnergy / 4) - (((self.MotorEnergy - self.MotorEnergyDelta)/ 2) / (1 + self.MotorCheYBindingConstant / CellMind.CheYp))
        MotorExpG = np.exp(self.MotorG)
        MotorCheYBindingConstantp = self.MotorSwitchingRate * MotorExpG
        MotorCheYBindingConstantm = self.MotorSwitchingRate / MotorExpG
        MotorRandom = rng.random(self.CellNumber)
        NewMotorState = np.logical_or( \
            np.logical_and(self.MotorState,MotorRandom < np.exp(-SimulationTimeStep * MotorCheYBindingConstantp)),\
                np.logical_and(~self.MotorState,MotorRandom >= np.exp(-SimulationTimeStep * MotorCheYBindingConstantm)))
        self.MotorHasSwitched = np.logical_xor(self.MotorState,NewMotorState)
        self.MotorState = NewMotorState
        
    def update_swimming_state(self):
        # maybe add flagella bundle dynamics
        self.Speed = (1 * np.logical_not(self.MotorState) - 1 * self.MotorState) * self.SpeedMaximum * np.logical_not(self.MotorHasSwitched) + \
            self.MotorHasSwitched * self.SpeedMinimum
        self.InterruptedRun = np.logical_or(self.InterruptedRun,self.MotorHasSwitched)

    def update_body_orientation(self,SimulationTimeStep):
        self.RandomDiffusion = SimulationTimeStep / 2 * (self.RotationalDiffusion[None,:] + \
            (np.logical_and(self.MotorHasSwitched,self.MotorState)[:,None] * self.RotationMotorSwitchFR[None,:]) + \
                (np.logical_and(self.MotorHasSwitched,np.logical_not(self.MotorState))[:,None] * self.RotationMotorSwitchRF[None,:]) + \
                    (self.Collision[:,None] * self.RotationMotorSwitchRF[None,:]))
        self.QuaternionRotation[:,1:] = (rng.random((self.CellNumber,3)) - 0.5) * np.sqrt(self.RandomDiffusion * 12) # 12 is necessary to correct for the uniform distribution approximation to get variance of 1
        self.QuaternionRotation[:,0] = np.sqrt(1 - self.RandomDiffusion.sum(axis = 1))
        self.QuaternionRotation /= np.sqrt(np.sum(self.QuaternionRotation**2,axis = 1))[:,None]
        self.update_quaternion()
        self.update_orientation()

    def override_parameters(self, Override):
        for item in vars(Override):
            assert hasattr(self,item)
            setattr(self,item,getattr(Override,item))


class BodyRecords:

    def __init__(self,CellNumber,LigandNumber,SampleNumber,SampleInterval):
        self.SampleInterval = SampleInterval
        self.Time = np.zeros(SampleNumber)
        self.MotorG = np.zeros((CellNumber,SampleNumber))
        self.Position = np.zeros((CellNumber,3,SampleNumber))
        self.Speed = np.zeros((CellNumber,SampleNumber))
        self.OrientationXYZ = np.zeros((CellNumber,3,SampleNumber))
        self.InterruptedRun = np.zeros((CellNumber,SampleNumber),dtype=np.bool)
        self.LigandConsumed = np.zeros((CellNumber,LigandNumber,SampleNumber))
        self.LigandConsumedIntegrated = np.zeros((CellNumber,LigandNumber,SampleNumber))
        self.CellSize = np.zeros((CellNumber,SampleNumber))
        self.CellPoleAge = np.zeros((CellNumber,2,SampleNumber))
        self.MotorNumber = np.zeros((CellNumber,SampleNumber))
        self.Iterator = 0
    
    def reset(self,SampleNumber):
        CellNumber = self.Position.shape[0]
        LigandNumber = self.LigandConsumed.shape[1]
        self.Time = np.zeros(SampleNumber)
        self.MotorG = np.zeros((CellNumber,SampleNumber))
        self.Position = np.zeros((CellNumber,3,SampleNumber))
        self.Speed = np.zeros((CellNumber,SampleNumber))
        self.OrientationXYZ = np.zeros((CellNumber,3,SampleNumber))
        self.InterruptedRun = np.zeros((CellNumber,SampleNumber),dtype=np.bool)
        self.LigandConsumed = np.zeros((CellNumber,LigandNumber,SampleNumber))
        self.LigandConsumedIntegrated = np.zeros((CellNumber,LigandNumber,SampleNumber))
        self.CellSize = np.zeros((CellNumber,SampleNumber))
        self.CellPoleAge = np.zeros((CellNumber,2,SampleNumber))
        self.MotorNumber = np.zeros((CellNumber,SampleNumber))
        self.Iterator = 0


class EcoliBody(CellBody):

    def __init__(self,CellNumber,SampleNumber,SampleInterval,Override=None):

        ## cell parameters
        self.Name = 'Escherichia coli'
        self.RotationalDiffusion = np.ones(3) * 0.062
        self.RotationMotorSwitchFR = np.zeros(3)
        self.RotationMotorSwitchRF = np.ones(3) * 6.2

        ## Motor
        self.SpeedMaximum = 25 #cell speed um/s
        self.SpeedMinimum = 0.25 #cell speed um/s
        self.MotorSwitchingRate = 1.3
        self.MotorEnergy = 40
        # self.MotorAdaptationRate = 1
        self.MotorCheYBindingConstant = 3.06
        self.MotorAverageNumber = 3  #number of flagellar Motors
        self.MotorRelativeCost = 0.01

        # self.MotorBiasEquilibrium = 0.3
        # self.MotorGEquilibrium = np.log(1 / self.MotorBiasEquilibrium - 1) / 2

        ## Ligand consumption and secretion
        self.LigandNumber = 2
        self.LigandConsumptionRate = np.array([1,0])
        self.LigandMichaelisConstant = np.array([1e3,1e3])
        self.LigandSecretionRate = np.array([0,0.1])

        # cell growth and aging
        self.CellPoleMaximumAge = 60*60*1 #seconds
        self.CellMaximunSize = 2 #um3
        
        self.ConversionLigandSize = np.array([2e-3,2e-3])
        # self.MotorSynthesisRate = self.MotorAverageNumber / self.LigandConsumptionRate

        if not(Override is None):
            print('Overriding parameters in CellBody!')
            self.override_parameters(Override)

        super(EcoliBody,self).__init__(CellNumber,SampleNumber,SampleInterval)

        self.CellPoleAge = rng.random((CellNumber,2)) * self.CellPoleMaximumAge
        self.CellSize = 1 + rng.random(CellNumber)
        # self.MotorNumber = rng.poisson(self.MotorAverageNumber,CellNumber)
