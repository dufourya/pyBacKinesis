"""
Library of classes defining the signaling pathway of cell agent for simulations.
The main CellMind class must contain the following methods:
'update' defines how the dynamics of the signaling pathway changes at each time step.
'record' defines which paramters are recorded.
Subclasses of CellMind must be define to specify the parameters of different types of cell agents.
"""

"""External libaries"""
import numpy as np

rng = np.random.default_rng()

class CellMind:
  def __init__(self, CellNumber, SampleNumber, SampleInterval):
    self.CellNumber = CellNumber
    self.Records = MindRecords(CellNumber, self.ReceptorNumber, SampleNumber, SampleInterval)

    #initialize cell state
    self.LigandConcentration = np.zeros((CellNumber, self.ReceptorNumber))
    self.Methylation = np.zeros((CellNumber, self.ReceptorNumber))
    self.MethylationAdapted = np.zeros((CellNumber, self.ReceptorNumber))
    self.ResetReceptorMethylation = np.ones((CellNumber, self.ReceptorNumber), dtype = np.bool)
    self.ClusterFreeEnergy = np.zeros((CellNumber, self.ReceptorNumber))
    self.ClusterRelativeActivity = np.zeros((CellNumber, self.ReceptorNumber))
    self.ClusterTotalRelativeActivity = np.zeros(CellNumber)
    self.CheYp = np.zeros(CellNumber)
    
  def update(self, CellBody, World, SimulationTimeStep):
    self.update_world_signal(CellBody, World)
    self.update_signaling_activity()
    self.update_receptor_methylation(SimulationTimeStep)

  def record(self, Clock):
    self.Records.Time[self.Records.Iterator] = Clock
    self.Records.LigandConcentration[:,:,self.Records.Iterator] = self.LigandConcentration
    self.Records.MethylationState[:,:,self.Records.Iterator] = self.Methylation
    self.Records.MethylationAdapted[:,:,self.Records.Iterator] = self.MethylationAdapted
    self.Records.ClusterRelativeActivity[:,:,self.Records.Iterator] = self.ClusterRelativeActivity
    self.Records.ClusterTotalRelativeActivity[:,self.Records.Iterator] = self.ClusterTotalRelativeActivity
    self.Records.ReceptorMethylationRate[:,:,self.Records.Iterator] = self.ReceptorMethylationRate
    self.Records.CheY[:,self.Records.Iterator] = self.CheY
    self.Records.Iterator += 1

  def update_world_signal(self, CellBody, World):
    self.LigandConcentration = World.get_ligand_concentration(CellBody)

  def update_signaling_activity(self):
    self.MethylationAdapted = self.ReceptorMethylationEquilibrium[None,:] - (self.ReceptorCooperativity/self.ReceptorEnergyActive[None,:]) \
      * (self.ReceptorProportions[None,:] * np.sum(np.log((1 + self.LigandConcentration[:,:,None] / self.ReceptorInactiveBindingConst[None,:,:]) / \
        (1 + self.LigandConcentration[:,:,None] / self.ReceptorActiveBindingConst[None,:,:])), axis = 2))
    #reset methylation state if required
    self.Methylation[self.ResetReceptorMethylation] = self.MethylationAdapted[self.ResetReceptorMethylation]
    self.ResetReceptorMethylation = self.ResetReceptorMethylation * False
    self.ClusterFreeEnergy = self.ReceptorProportions[None,:] * \
      (self.ReceptorEnergyInactive[None,:] + self.ReceptorEnergyActive[None,:] * \
        (self.ReceptorMethylationEquilibrium[None,:] + self.Methylation - self.MethylationAdapted))
    self.ClusterRelativeActivity = 1 / (1 + np.exp(self.ClusterFreeEnergy))
    self.ClusterTotalRelativeActivity = 1 / (1 + np.exp(np.sum(self.ClusterFreeEnergy, axis = 1)))
    self.CheYp = self.CheY * self.ClusterTotalRelativeActivity
  
  def update_receptor_methylation(self, SimulationTimeStep):
    #update receptor Methylation adaptation
    self.Methylation += SimulationTimeStep * ( \
      self.ReceptorMethylationRate * (1 - self.ClusterRelativeActivity) - \
        self.ClusterRelativeActivity * self.ReceptorDemethylationRate * \
          (self.ClusterTotalRelativeActivity[:,None] * (self.DemethylationActivationFactor - 1) + 1))
    # add Methylation noise
    if not(self.MethylationNoiseSigma is None):
      self.Methylation += SimulationTimeStep * self.MethylationNoiseSigma * (rng.random((self.CellNumber,self.ReceptorNumber)) - 0.5)
    self.Methylation[self.Methylation < self.ReceptorMethylationMin] = self.ReceptorMethylationMin
    self.Methylation[self.Methylation > self.ReceptorMethylationMax] = self.ReceptorMethylationMax

  def override_parameters(self, Override):
    for item in vars(Override):
      assert hasattr(self,item)
      setattr(self,item,getattr(Override,item))
        


class MindRecords:

  def __init__(self, CellNumber, ReceptorNumber, SampleNumber, SampleInterval):
    self.SampleInterval = SampleInterval
    self.Time = np.zeros(SampleNumber)
    self.LigandConcentration = np.zeros((CellNumber, ReceptorNumber, SampleNumber))
    self.MethylationState = np.zeros((CellNumber, ReceptorNumber, SampleNumber))
    self.MethylationAdapted = np.zeros((CellNumber, ReceptorNumber, SampleNumber))
    self.ClusterRelativeActivity = np.zeros((CellNumber, ReceptorNumber, SampleNumber))
    self.ClusterTotalRelativeActivity = np.zeros((CellNumber, SampleNumber))
    self.ReceptorMethylationRate = np.zeros((CellNumber, ReceptorNumber, SampleNumber))
    self.CheY = np.zeros((CellNumber, SampleNumber))
    self.Iterator = 0
  
  def reset(self, SampleNumber):
    CellNumber = self.LigandConcentration.shape[0]
    ReceptorNumber = self.LigandConcentration.shape[1]
    self.Time = np.zeros(SampleNumber)
    self.LigandConcentration = np.zeros((CellNumber, ReceptorNumber, SampleNumber))
    self.MethylationState = np.zeros((CellNumber, ReceptorNumber, SampleNumber))
    self.MethylationAdapted = np.zeros((CellNumber, ReceptorNumber, SampleNumber))
    self.ClusterRelativeActivity = np.zeros((CellNumber, ReceptorNumber, SampleNumber))
    self.ClusterTotalRelativeActivity = np.zeros((CellNumber, SampleNumber))
    self.CheY = np.zeros((CellNumber, SampleNumber))
    self.Iterator = 0


class EcoliMind(CellMind):
  def __init__(self, CellNumber, SampleNumber, SampleInterval, Override=None):

    self.Name = 'Escherichia coli'

    # Biochemical constants
    self.ReceptorNumber = 2
    self.ReceptorProportions = np.array([1, 1])
    self.ReceptorProportions = self.ReceptorProportions / np.sum(self.ReceptorProportions)
    self.ReceptorEnergyInactive = np.ones(self.ReceptorNumber) * 6 #kbT
    self.ReceptorEnergyActive   = np.ones(self.ReceptorNumber) * -1 #kbT
    self.ReceptorInactiveBindingConst = np.array([[1e1, 1e3], [1e3, 1e1]]) #uM
    self.ReceptorActiveBindingConst   = np.array([[1e4, 1e6], [1e6, 1e4]]) #uM
    self.ClusterRelativeActivityEquilibrium = np.array([0.37, 0.37])

    self.ReceptorCooperativity = 6
    self.ReceptorMethylationMin = 0   #minimum value allowed in Methylation
    self.ReceptorMethylationMax = 48  #maximum value allowed in Methylation
    self.ReceptorMethylationRate = np.power(10, -2 + 2 * rng.random((CellNumber, self.ReceptorNumber)))
    # self.ReceptorMethylationRate = 0.1

    self.DemethylationActivationFactor = 70
    self.MethylationNoiseSigma = None
    
    self.CheY = rng.random(CellNumber) + 5.5 #uM

    if not(Override is None):
      print('Overriding parameters in CellMind!')
      self.override_parameters(Override)

    # calculated demethylation rate to set receptor resting activity
    self.ReceptorMethylationEquilibrium = (np.log(1 / self.ClusterRelativeActivityEquilibrium - 1) - self.ReceptorEnergyInactive) / self.ReceptorEnergyActive
    self.ReceptorDemethylationRate = self.ReceptorMethylationRate * (1 - self.ClusterRelativeActivityEquilibrium[None,:]) / \
      (self.ClusterRelativeActivityEquilibrium[None,:] * (self.ClusterRelativeActivityEquilibrium[None,:] * (self.DemethylationActivationFactor - 1) + 1))

    # intialize cells
    super(EcoliMind, self).__init__(CellNumber, SampleNumber, SampleInterval)