"""
Library of classes defining environments for agent-based simulations.
Each class must contain the following methods:
'update' defines how the world evolves at each time step.
'record' defines which paramters are recorded.
'get_ligand_concentration' communicates ligand concentrations with the CellBody class.
'initialize_cell_positions' defines how cells are seeded into the world.
'check_world_boundaries' determines the physical interactions between cells and ophysical obstacles or boundaries.
"""

"""External libaries"""
import numpy as np
from scipy import ndimage

rng = np.random.default_rng()

class WorldRecords:
    """Stores chemical concentrations over time"""
    def __init__(self,World,SampleNumber,SampleInterval):
        self.SampleInterval = SampleInterval
        self.Time = np.zeros(SampleNumber)
        self.LigandConcentration = np.zeros((World.DimensionsScaled[0],World.DimensionsScaled[1],World.LigandNumber,SampleNumber))
        self.Iterator = 0
    
    def reset(self,SampleNumber):
        self.Time = np.zeros(SampleNumber)
        self.LigandConcentration = np.zeros(self.LigandConcentration.shape[:-1] + (SampleNumber,))
        self.Iterator = 0


class BoxExponential:
    """Exponential gradients in a closed box"""
    def __init__(self,SampleNumber,SampleInterval):
        self.LigandNumber = 2
        self.LigandGradientLengthScale = np.array([4e2,4e2])
        self.LigandBackgroundConcentration = np.array([1,1])
        self.Dimensions = np.array([1e4,1e3,1e3],dtype = np.int16)
        self.WorldScaling = 1
        self.DimensionsScaled = np.uint16(self.Dimensions / self.WorldScaling)
        self.Records = WorldRecords(self,SampleNumber,SampleInterval)
        self.GradientDirection = np.array([1,1])
        self.GradientOffset = np.array([0,self.Dimensions[0]])

        
    def update(self,CellBody,SimulationTimeStep):
        self.check_world_boundaries(CellBody)

    def record(self,Clock):
        self.Records.Time[self.Records.Iterator] = Clock
        self.Records.LigandConcentration[:,:,:,self.Records.Iterator] = np.ones(self.DimensionsScaled[1])[None,:,None] * self.LigandBackgroundConcentration[None,None,:] * \
            np.exp(np.arange(self.DimensionsScaled[0])[:,None,None] / self.LigandGradientLengthScale[None,None,:])
        self.Records.Iterator += 1

    def get_ligand_concentration(self,CellBody):
        return self.LigandBackgroundConcentration[None,:] * np.exp(((CellBody.Position[:,0,None] * self.GradientDirection[None,:] + self.GradientOffset[None,:])) / self.LigandGradientLengthScale[None,:])

    def check_world_boundaries(self,CellBody):
        CellBody.Collision = np.logical_or(np.sum(CellBody.Position > self.Dimensions[None,:],axis=1),np.sum(CellBody.Position < np.array([0,0,0]),axis=1))
        CellBody.Position[CellBody.Position > self.Dimensions[None,:]] = (2 * self.Dimensions[None,:] - CellBody.Position)[CellBody.Position > self.Dimensions[None,:]]
        CellBody.Position[CellBody.Position < np.array([0,0,0])] = - CellBody.Position[CellBody.Position < np.array([0,0,0])]

    def initialize_cell_positions(self,CellBody):
        CellBody.Position = CellBody.Position * self.Dimensions * rng.random((CellBody.CellNumber,3))
        # CellBody.Position[:,0] = CellBody.Position[:,0] * self.Dimensions[0] / 4
        # CellBody.Position[:,1:] = CellBody.Position[:,1:] * self.Dimensions[1:] * rng.random((CellBody.CellNumber,2))


class BoxHillEquation:
    """Sigmoidal gradients in a closed box"""
    def __init__(self,SampleNumber,SampleInterval):
        self.LigandNumber = 2
        self.LigandMaxConcentration = np.array([1e4,1e4])
        self.LigandHalfGradient = np.array([5e3,5e3])
        self.LigandGradientSlope = np.array([3,3])
        self.Dimensions = np.array([1e4,1e3,1e3],dtype = np.int16)
        self.WorldScaling = 1
        self.DimensionsScaled = np.uint16(self.Dimensions / self.WorldScaling)
        self.Records = WorldRecords(self,SampleNumber,SampleInterval)
        self.GradientDirection = np.array([1,-1])
        self.GradientOffset = np.array([0,self.Dimensions[0]])
        
    def update(self,CellBody,SimulationTimeStep):
        self.check_world_boundaries(CellBody)

    def record(self,Clock):
        self.Records.Time[self.Records.Iterator] = Clock
        self.Records.LigandConcentration[:,:,:,self.Records.Iterator] = np.ones(self.DimensionsScaled[1])[None,:,None] * \
            (self.LigandMaxConcentration[None,:] / (1 + (self.LigandHalfGradient[None,None,:] / (0.01 + np.arange(self.DimensionsScaled[0])[:,None,None])) ** self.LigandGradientSlope[None,None,:])) # added 0.01 to avoid divide by 0
        self.Records.Iterator += 1

    def get_ligand_concentration(self,CellBody):
        return self.LigandMaxConcentration[None,:] / (1 + (self.LigandHalfGradient[None,:] / ((0.01 + CellBody.Position[:,0,None]) * self.GradientDirection[None,:] + self.GradientOffset[None,:])) ** self.LigandGradientSlope[None,:]) # added 0.01 to cell positions to avoid divide by 0

    def check_world_boundaries(self,CellBody):
        CellBody.Collision = np.logical_or(np.sum(CellBody.Position > self.Dimensions[None,:],axis=1),np.sum(CellBody.Position < np.array([0,0,0]),axis=1))
        CellBody.Position[CellBody.Position > self.Dimensions[None,:]] = (2 * self.Dimensions[None,:] - CellBody.Position)[CellBody.Position > self.Dimensions[None,:]]
        CellBody.Position[CellBody.Position < np.array([0,0,0])] = - CellBody.Position[CellBody.Position < np.array([0,0,0])]

    def initialize_cell_positions(self,CellBody):
        CellBody.Position = CellBody.Position * self.Dimensions * rng.random((CellBody.CellNumber,3))
        # CellBody.Position[:,0] = CellBody.Position[:,0] * self.Dimensions[0] / 4
        # CellBody.Position[:,1:] = CellBody.Position[:,1:] * self.Dimensions[1:] * rng.random((CellBody.CellNumber,2))


class BoxLinear:
    """Linear exponential gradients in a closed box"""
    def __init__(self,SampleNumber,SampleInterval):
        self.GradientSlope = np.array([0.1,-0.1])
        self.LigandBackgroundConcentration = np.array([1,1e3+1])
        self.Dimensions = np.array([10000,1000,100])
        self.LigandNumber = 2
        self.WorldScaling = 100
        self.DimensionsScaled = np.uint16(self.Dimensions / self.WorldScaling)
        self.Records = WorldRecords(self,SampleNumber,SampleInterval)
        
    def update(self,CellBody,SimulationTimeStep):
        self.check_world_boundaries(CellBody)

    def record(self,Clock):
        self.Records.Time[self.Records.Iterator] = Clock
        self.Records.LigandConcentration[:,:,:,self.Records.Iterator] = np.ones(self.DimensionsScaled[1])[None,:,None] * (self.LigandBackgroundConcentration[None,None,:] + \
            np.arange(self.DimensionsScaled[0])[:,None,None] * self.GradientSlope[None,None,:])
        self.Records.Iterator += 1

    def get_ligand_concentration(self,CellBody):
        return self.LigandBackgroundConcentration[None,:] + CellBody.Position[:,0,None] * self.GradientSlope[None,:]
    
    def check_world_boundaries(self,CellBody):
        CellBody.Collision = np.logical_or(np.sum(CellBody.Position > self.Dimensions[None,:],axis=1),np.sum(CellBody.Position < np.array([0,0,0]),axis=1))
        CellBody.Position[CellBody.Position > self.Dimensions[None,:]] = (2 * self.Dimensions[None,:] - CellBody.Position)[CellBody.Position > self.Dimensions[None,:]]
        CellBody.Position[CellBody.Position < np.array([0,0,0])] = - CellBody.Position[CellBody.Position < np.array([0,0,0])]
    
    def initialize_cell_positions(self,CellBody):
        # CellBody.Position[:,0] = CellBody.Position[:,0] * self.Dimensions[0] / 4
        CellBody.Position = CellBody.Position * self.Dimensions * rng.random((CellBody.CellNumber,3))
        # CellBody.Position[:,1:] = CellBody.Position[:,1:] * self.Dimensions[1:] * rng.random((CellBody.CellNumber,2))
 

class ThreadmillLinear:
    """Linear gradients in a closed box with bulk flow"""
    def __init__(self,SampleNumber,SampleInterval):
        self.GradientSlope = np.array([5e-3,5e-3])
        self.LigandBackgroundConcentration = np.array([1e-6,1e-6])
        self.Dimensions = np.array([50000,5000,5000])
        self.LigandNumber = 2
        self.WorldScaling = 10
        self.DimensionsScaled = np.uint16(self.Dimensions / self.WorldScaling)
        self.Records = WorldRecords(self,SampleNumber,SampleInterval)
        self.FlowRate = -20
        
    def update(self,CellBody,SimulationTimeStep):
        self.check_world_boundaries(CellBody,SimulationTimeStep)

    def record(self,Clock):
        self.Records.Time[self.Records.Iterator] = Clock
        self.Records.LigandConcentration[:,:,:,self.Records.Iterator] = np.ones(self.DimensionsScaled[1])[None,:,None] * (self.LigandBackgroundConcentration[None,None,:] + \
            np.arange(self.DimensionsScaled[0])[:,None,None] * self.GradientSlope[None,None,:])
        self.Records.Iterator += 1

    def get_ligand_concentration(self,CellBody):
        return self.LigandBackgroundConcentration[None,:] + CellBody.Position[:,0,None] * self.GradientSlope[None,:]
    
    def check_world_boundaries(self,CellBody,SimulationTimeStep):
        CellBody.Collision = np.logical_or(np.sum(CellBody.Position > self.Dimensions[None,:],axis=1),np.sum(CellBody.Position < np.array([0,0,0]),axis=1))
        CellBody.Position[:,0] = CellBody.Position[:,0] + self.FlowRate * SimulationTimeStep
        CellBody.Position[CellBody.Position > self.Dimensions[None,:]] = (2 * self.Dimensions[None,:] - CellBody.Position)[CellBody.Position > self.Dimensions[None,:]]
        CellBody.Position[CellBody.Position < np.array([0,0,0])] = - CellBody.Position[CellBody.Position < np.array([0,0,0])]
    
    def initialize_cell_positions(self,CellBody):
        CellBody.Position[:,0] = CellBody.Position[:,0] * self.Dimensions[0] * 0.95
        CellBody.Position[:,1:] = CellBody.Position[:,1:] * self.Dimensions[1:] * rng.random((CellBody.CellNumber,2))


class ThreadmillExponential:
    """Exponential gradients in a closed box with bulk flow"""
    def __init__(self,SampleNumber,SampleInterval):
        self.LigandNumber = 2
        self.LigandGradientLengthScale = np.array([1e3,1e3])
        self.LigandBackgroundConcentration = np.array([1,1])
        self.Dimensions = np.array([1.5e4,1.5e3,1.5e3],dtype = np.int16)
        self.WorldScaling = 10
        self.DimensionsScaled = np.uint16(self.Dimensions / self.WorldScaling)
        self.Records = WorldRecords(self,SampleNumber,SampleInterval)
        self.FlowRate = -1
        
    def update(self,CellBody,SimulationTimeStep):
        self.check_world_boundaries(CellBody,SimulationTimeStep)

    def record(self,Clock):
        self.Records.Time[self.Records.Iterator] = Clock
        self.Records.LigandConcentration[:,:,:,self.Records.Iterator] = np.ones(self.DimensionsScaled[1])[None,:,None] * self.LigandBackgroundConcentration[None,None,:] * \
            np.exp(np.arange(self.DimensionsScaled[0])[:,None,None] / self.LigandGradientLengthScale[None,None,:])
        self.Records.Iterator += 1

    def get_ligand_concentration(self,CellBody):
        return self.LigandBackgroundConcentration[None,:] * np.exp(CellBody.Position[:,0,None] / self.LigandGradientLengthScale[None,:])
    
    def check_world_boundaries(self,CellBody,SimulationTimeStep):
        CellBody.Collision = np.logical_or(np.sum(CellBody.Position > self.Dimensions[None,:],axis=1),np.sum(CellBody.Position < np.array([0,0,0]),axis=1))
        CellBody.Position[:,0] = CellBody.Position[:,0] + self.FlowRate * SimulationTimeStep
        CellBody.Position[CellBody.Position > self.Dimensions[None,:]] = (2 * self.Dimensions[None,:] - CellBody.Position)[CellBody.Position > self.Dimensions[None,:]]
        CellBody.Position[CellBody.Position < np.array([0,0,0])] = - CellBody.Position[CellBody.Position < np.array([0,0,0])]
    
    def initialize_cell_positions(self,CellBody):
        CellBody.Position[:,0] = CellBody.Position[:,0] * self.Dimensions[0] * 0.5
        CellBody.Position[:,1:] = CellBody.Position[:,1:] * self.Dimensions[1:] * rng.random((CellBody.CellNumber,2))


class StepLigand:
    """Uniform chemical concentrations with sudden changes at specified times"""
    def __init__(self,SampleNumber,SampleInterval):
        self.Dimensions = np.array([5000,5000,50])
        self.LigandConcentration = np.array([1,1])
        self.Time = 0
        self.LigandNumber = 2
        self.WorldScaling = 100
        self.DimensionsScaled = np.uint16(self.Dimensions / self.WorldScaling)
        self.Records = WorldRecords(self,SampleNumber,SampleInterval)

    def update(self,CellBody,SimulationTimeStep):
        self.Time += SimulationTimeStep
        if self.Time > 10:
            self.LigandConcentration = np.array([6,1])
        if self.Time > 30:
            self.LigandConcentration = np.array([6,4])
        if self.Time > 50:
            self.LigandConcentration = np.array([1,4])
        if self.Time > 70:
            self.LigandConcentration = np.array([1,1])
        self.check_world_boundaries(CellBody)

    def record(self,Clock):
        self.Records.Time[self.Records.Iterator] = Clock
        self.Records.LigandConcentration[:,:,:,self.Records.Iterator] = np.ones(self.DimensionsScaled[:2])[:,:,None] * self.LigandConcentration[None,None,:]
        self.Records.Iterator += 1

    def get_ligand_concentration(self,CellBody):
        return np.ones(CellBody.CellNumber)[:,None] * self.LigandConcentration[None,:]
    
    def check_world_boundaries(self,CellBody):
        CellBody.Collision = np.logical_or(np.sum(CellBody.Position > self.Dimensions[None,:],axis=1),np.sum(CellBody.Position < np.array([0,0,0]),axis=1))
        CellBody.Position[CellBody.Position > self.Dimensions[None,:]] = (2 * self.Dimensions[None,:] - CellBody.Position)[CellBody.Position > self.Dimensions[None,:]]
        CellBody.Position[CellBody.Position < np.array([0,0,0])] = - CellBody.Position[CellBody.Position < np.array([0,0,0])]
    
    def initialize_cell_positions(self,CellBody):
        CellBody.Position = CellBody.Position * self.Dimensions * rng.random((CellBody.CellNumber,3))


class Maze_and_Diffusion:
    """Reaction-diffusion gradients with physical obstacles"""
    def __init__(self,SampleNumber,SampleInterval):
        self.Dimensions = np.array([5000,500,100])
        self.WorldScaling = 10 #scale gradient matrix to speed up reaction-diffusion updates
        self.DimensionsScaled = np.uint16(self.Dimensions / self.WorldScaling)

        self.Obstacles = np.zeros(self.DimensionsScaled[:-1],dtype = np.bool)
        self.ObstaclesIndices = np.zeros(((np.ndim(self.Obstacles),) + self.Obstacles.shape),dtype=np.int32)
        self.ObstaclesDensity = 0.2
        self.ObstaclesSize = 100
        # self.initialize_square_positions()
        # self.initialize_circle_positions()
        self.initialize_u_shape_positions()
        # self.initialize_cross_shape_positions()

        self.LigandNumber = 2
        self.LigandSourceConcentration = np.array([1000,0])
        self.LigandInitialConcentration = np.array([1000,1e-6])
        self.LigandDiffusionCoefficient = np.array([5e2,4e2])
        self.LigandConcentration = np.ones(np.append(self.DimensionsScaled,self.LigandNumber)) * self.LigandInitialConcentration[None,None,None,:]

        self.Records = WorldRecords(self,SampleNumber,SampleInterval)
    
    def update(self,CellBody,SimulationTimeStep):
        self.check_world_boundaries(CellBody)
        self.update_consumption_diffusion(CellBody,SimulationTimeStep)

    def record(self,Clock):
        self.Records.Time[self.Records.Iterator] = Clock
        self.Records.LigandConcentration[:,:,:,self.Records.Iterator] = np.mean(self.LigandConcentration,axis = 2).squeeze()
        self.Records.Iterator += 1

    def get_ligand_concentration(self,CellBody):
        return self.LigandConcentration[CellBody.PositionScaledWorld[:,0],CellBody.PositionScaledWorld[:,1],CellBody.PositionScaledWorld[:,2],:]
    
    def check_world_boundaries(self,CellBody):
        CellBody.Collision = np.logical_or(np.sum(CellBody.Position > self.Dimensions[None,:],axis=1),np.sum(CellBody.Position < np.array([0,0,0]),axis=1))
        CellBody.Position[CellBody.Position > self.Dimensions[None,:]] = (2 * self.Dimensions[None,:] - CellBody.Position)[CellBody.Position > self.Dimensions[None,:]]
        CellBody.Position[CellBody.Position < np.array([0,0,0])] = - CellBody.Position[CellBody.Position < np.array([0,0,0])]
        CellBody.PositionScaledWorld = np.uint16(CellBody.Position / self.WorldScaling)
        Collisions = self.Obstacles[CellBody.PositionScaledWorld[:,0],CellBody.PositionScaledWorld[:,1]]
        CellBody.Position[Collisions,:2] = CellBody.Position[Collisions,:2] - \
            ((CellBody.PositionScaledWorld[Collisions,:2] - self.ObstaclesIndices[:,CellBody.PositionScaledWorld[Collisions,0],CellBody.PositionScaledWorld[Collisions,1]].T) * \
                (self.WorldScaling * 0.5 - np.absolute(CellBody.Position[Collisions,:2] - (0.5 + CellBody.PositionScaledWorld[Collisions,:2]) * self.WorldScaling)))
        CellBody.Collision = np.logical_or(CellBody.Collision,Collisions)
    
    def initialize_cell_positions(self,CellBody):
        CellBody.Position = CellBody.Position * self.Dimensions * rng.random((CellBody.CellNumber,3))
        CellBody.PositionScaledWorld = np.uint16(CellBody.Position / self.WorldScaling)
    
    def update_consumption_diffusion(self,CellBody,SimulationTimeStep):
        for k in range(self.LigandNumber):
            np.subtract.at(self.LigandConcentration[:,:,:,k],(CellBody.PositionScaledWorld[:,0],CellBody.PositionScaledWorld[:,1],CellBody.PositionScaledWorld[:,2]),CellBody.LigandConsumed[:,k] / self.WorldScaling**3)
            self.LigandConcentration[self.LigandConcentration < 1e-12] = 1e-12
            if self.LigandSourceConcentration[k] > 0:
                self.LigandConcentration[-1,:,:,k] = np.ones(self.DimensionsScaled[[1,2]]) * self.LigandSourceConcentration[k]
            self.LigandConcentration[self.ObstaclesIndices[0,:],self.ObstaclesIndices[1,:],:,k] = self.LigandConcentration[:,:,:,k]
            self.LigandConcentration[:,:,:,k] = ndimage.gaussian_filter(self.LigandConcentration[:,:,:,k],\
                sigma = np.sqrt(self.LigandDiffusionCoefficient[k] / self.WorldScaling * SimulationTimeStep))

    def initialize_circle_positions(self):
        Number = np.uint16(self.ObstaclesDensity * self.DimensionsScaled[0] * self.DimensionsScaled[1] / (np.pi * (self.ObstaclesSize/2/self.WorldScaling)**2))
        self.Obstacles[np.unravel_index(rng.integers(np.prod(self.Obstacles.shape),size = Number),self.Obstacles.shape)] = True
        Distances = ndimage.distance_transform_edt(np.logical_not(self.Obstacles))
        self.Obstacles = Distances < self.ObstaclesSize/self.WorldScaling/2
        ndimage.distance_transform_edt(self.Obstacles,return_indices=True,indices=self.ObstaclesIndices)

    def initialize_square_positions(self):
        Shape = np.ones(np.uint16(np.ones(2) * self.ObstaclesSize / self.WorldScaling))
        Number = np.uint16(self.ObstaclesDensity * self.DimensionsScaled[0] * self.DimensionsScaled[1] / np.sum(Shape))
        self.Obstacles[np.unravel_index(rng.integers(np.prod(self.Obstacles.shape),size = Number),self.Obstacles.shape)] = True
        self.Obstacles = ndimage.morphology.binary_dilation(self.Obstacles,Shape)
        self.Obstacles = ndimage.binary_fill_holes(self.Obstacles).astype(bool)
        ndimage.distance_transform_edt(self.Obstacles,return_indices=True,indices=self.ObstaclesIndices)

    def initialize_u_shape_positions(self):
        Shape = np.zeros(np.uint16(np.ones(2) * self.ObstaclesSize / self.WorldScaling),dtype = np.bool)
        Shape[:,:3] = True
        Shape[-3:,:] = True
        Shape[:,-3:] = True
        Number = np.uint16(self.ObstaclesDensity * self.DimensionsScaled[0] * self.DimensionsScaled[1] / np.sum(Shape))
        self.Obstacles[np.unravel_index(rng.integers(np.prod(self.Obstacles.shape),size = Number),self.Obstacles.shape)] = True
        self.Obstacles = ndimage.morphology.binary_dilation(self.Obstacles,Shape)
        self.Obstacles = ndimage.morphology.binary_fill_holes(self.Obstacles)
        ndimage.distance_transform_edt(self.Obstacles,return_indices=True,indices=self.ObstaclesIndices)

    def initialize_cross_shape_positions(self):
        Shape = np.zeros(np.uint16(np.ones(2) * self.ObstaclesSize / self.WorldScaling),dtype = np.bool)
        Shape[:,np.uint16(Shape.shape[1]/2 - 1.5):np.uint16(Shape.shape[1]/2 + 1.5)] = True
        Shape[np.uint16(Shape.shape[1]/2 - 1.5):np.uint16(Shape.shape[1]/2 + 1.5),:] = True
        Number = np.uint16(self.ObstaclesDensity * self.DimensionsScaled[0] * self.DimensionsScaled[1] / np.sum(Shape))
        self.Obstacles[np.unravel_index(rng.integers(np.prod(self.Obstacles.shape),size = Number),self.Obstacles.shape)] = True
        self.Obstacles = ndimage.morphology.binary_dilation(self.Obstacles,Shape)
        self.Obstacles = ndimage.morphology.binary_fill_holes(self.Obstacles)
        ndimage.distance_transform_edt(self.Obstacles,return_indices=True,indices=self.ObstaclesIndices)


class OpenCapillary:
    """Reaction-diffusion gradients with a constant source"""
    def __init__(self,SampleNumber,SampleInterval):
        self.Dimensions = np.array([10000,1000,100])
        self.WorldScaling = 25 #scale gradient matrix to speed up reaction-diffusion updates
        self.DimensionsScaled = np.uint16(self.Dimensions / self.WorldScaling)
        self.LigandNumber = 2
        self.LigandSourceConcentration = np.array([1e3,1e3])
        self.LigandInitialConcentration = np.array([1e3,1e3])
        self.LigandDiffusionCoefficient = np.array([1e3,1e3])
        self.LigandConcentration = np.ones(np.append(self.DimensionsScaled,self.LigandNumber)) * self.LigandInitialConcentration[None,None,None,:]
        self.Records = WorldRecords(self,SampleNumber,SampleInterval)
    
    def update(self,CellBody,SimulationTimeStep):
        self.check_world_boundaries(CellBody)
        self.update_consumption_diffusion(CellBody,SimulationTimeStep)

    def record(self,Clock):
        self.Records.Time[self.Records.Iterator] = Clock
        self.Records.LigandConcentration[:,:,:,self.Records.Iterator] = np.mean(self.LigandConcentration,axis = 2).squeeze()
        self.Records.Iterator += 1

    def get_ligand_concentration(self,CellBody):
        return self.LigandConcentration[CellBody.PositionScaledWorld[:,0],CellBody.PositionScaledWorld[:,1],CellBody.PositionScaledWorld[:,2],:]
    
    def check_world_boundaries(self,CellBody):
        CellBody.Collision = np.logical_or(np.sum(CellBody.Position > self.Dimensions[None,:],axis=1),np.sum(CellBody.Position < np.array([0,0,0]),axis=1))
        CellBody.Position[CellBody.Position > self.Dimensions[None,:]] = (2 * self.Dimensions[None,:] - CellBody.Position)[CellBody.Position > self.Dimensions[None,:]]
        CellBody.Position[CellBody.Position < np.array([0,0,0])] = - CellBody.Position[CellBody.Position < np.array([0,0,0])]
        CellBody.PositionScaledWorld = np.uint16(CellBody.Position / self.WorldScaling)
    
    def initialize_cell_positions(self,CellBody):
        CellBody.Position = CellBody.Position * self.Dimensions * rng.random((CellBody.CellNumber,3))
        CellBody.PositionScaledWorld = np.uint16(CellBody.Position / self.WorldScaling)
    
    def update_consumption_diffusion(self,CellBody,SimulationTimeStep):
        for k in range(self.LigandNumber):
            np.subtract.at(self.LigandConcentration[:,:,:,k],(CellBody.PositionScaledWorld[:,0],CellBody.PositionScaledWorld[:,1],CellBody.PositionScaledWorld[:,2]),CellBody.LigandConsumed[:,k] / self.WorldScaling**3)
            self.LigandConcentration[self.LigandConcentration < 1e-12] = 1e-12
            if self.LigandSourceConcentration[k] > 0:
                self.LigandConcentration[-1,:,:,k] = np.ones(self.DimensionsScaled[[1,2]]) * self.LigandSourceConcentration[k]
            self.LigandConcentration[:,:,:,k] = ndimage.gaussian_filter(self.LigandConcentration[:,:,:,k],\
                sigma = np.sqrt(self.LigandDiffusionCoefficient[k] / self.WorldScaling * SimulationTimeStep))


class Maze_and_Exponential:
    """Exponential gradients with physical obstacles"""
    def __init__(self,SampleNumber,SampleInterval):
        self.Dimensions = np.array([1000,1000,100])
        self.WorldScaling = 25 #scale gradient matrix to speed up reaction-diffusion updates
        self.DimensionsScaled = np.uint16(self.Dimensions / self.WorldScaling)
        self.Obstacles = np.zeros(self.DimensionsScaled[:-1],dtype = np.bool)
        self.ObstaclesIndices = np.zeros(((np.ndim(self.Obstacles),) + self.Obstacles.shape),dtype=np.int32)
        self.ObstaclesDensity = 0.2
        self.ObstaclesSize = 100
        # self.initialize_square_positions()
        # self.initialize_circle_positions()
        self.initialize_u_shape_positions()
        # self.initialize_cross_shape_positions()
        self.LigandNumber = 2
        self.LigandBackgroundConcentration = np.array([1e2,1e2])
        self.LigandGradientLengthScale = np.array([3e3,3e3])
        self.Records = WorldRecords(self,SampleNumber,SampleInterval)
    
    def update(self,CellBody,SimulationTimeStep):
        self.check_world_boundaries(CellBody)

    def record(self,Clock):
        self.Records.Time[self.Records.Iterator] = Clock
        self.Records.LigandConcentration[:,:,:,self.Records.Iterator] = np.ones(self.DimensionsScaled[1])[None,:,None] * self.LigandBackgroundConcentration[None,None,:] * \
            np.exp(np.arange(self.DimensionsScaled[0])[:,None,None] / self.LigandGradientLengthScale[None,None,:])
        self.Records.Iterator += 1

    def get_ligand_concentration(self,CellBody):
        return self.LigandBackgroundConcentration[None,:] * np.exp(CellBody.Position[:,0,None] / self.LigandGradientLengthScale[None,:])
    
    def initialize_cell_positions(self,CellBody):
        CellBody.Position = CellBody.Position * self.Dimensions * rng.random((CellBody.CellNumber,3))
        # CellBody.Position[:,0] = CellBody.Position[:,0] * self.Dimensions[0] / 4
        # CellBody.Position[:,1:] = CellBody.Position[:,1:] * self.Dimensions[1:] * rng.random((CellBody.CellNumber,2))
        CellBody.PositionScaledWorld = np.uint16(CellBody.Position / self.WorldScaling)
    
    def check_world_boundaries(self,CellBody):
        CellBody.Collision = np.logical_or(np.sum(CellBody.Position > self.Dimensions[None,:],axis=1),np.sum(CellBody.Position < np.array([0,0,0]),axis=1))
        CellBody.Position[CellBody.Position > self.Dimensions[None,:]] = (2 * self.Dimensions[None,:] - CellBody.Position)[CellBody.Position > self.Dimensions[None,:]]
        CellBody.Position[CellBody.Position < np.array([0,0,0])] = - CellBody.Position[CellBody.Position < np.array([0,0,0])]
        CellBody.PositionScaledWorld = np.uint16(CellBody.Position / self.WorldScaling)
        Collisions = self.Obstacles[CellBody.PositionScaledWorld[:,0],CellBody.PositionScaledWorld[:,1]]
        CellBody.Position[Collisions,:2] = CellBody.Position[Collisions,:2] - \
            ((CellBody.PositionScaledWorld[Collisions,:2] - self.ObstaclesIndices[:,CellBody.PositionScaledWorld[Collisions,0],CellBody.PositionScaledWorld[Collisions,1]].T) * \
                (self.WorldScaling * 0.5 - np.absolute(CellBody.Position[Collisions,:2] - (0.5 + CellBody.PositionScaledWorld[Collisions,:2]) * self.WorldScaling)))
        CellBody.Collision = np.logical_or(CellBody.Collision,Collisions)

    def initialize_circle_positions(self):
        Number = np.uint16(self.ObstaclesDensity * self.DimensionsScaled[0] * self.DimensionsScaled[1] / (np.pi * (self.ObstaclesSize/2/self.WorldScaling)**2))
        self.Obstacles[np.unravel_index(rng.integers(np.prod(self.Obstacles.shape),size = Number),self.Obstacles.shape)] = True
        Distances = ndimage.distance_transform_edt(np.logical_not(self.Obstacles))
        self.Obstacles = Distances < self.ObstaclesSize/self.WorldScaling/2
        ndimage.distance_transform_edt(self.Obstacles,return_indices=True,indices=self.ObstaclesIndices)

    def initialize_square_positions(self):
        Shape = np.ones(np.uint16(np.ones(2) * self.ObstaclesSize / self.WorldScaling))
        Number = np.uint16(self.ObstaclesDensity * self.DimensionsScaled[0] * self.DimensionsScaled[1] / np.sum(Shape))
        self.Obstacles[np.unravel_index(rng.integers(np.prod(self.Obstacles.shape),size = Number),self.Obstacles.shape)] = True
        self.Obstacles = ndimage.morphology.binary_dilation(self.Obstacles,Shape)
        self.Obstacles = ndimage.binary_fill_holes(self.Obstacles).astype(bool)
        ndimage.distance_transform_edt(self.Obstacles,return_indices=True,indices=self.ObstaclesIndices)

    def initialize_u_shape_positions(self):
        Shape = np.zeros(np.uint16(np.ones(2) * self.ObstaclesSize / self.WorldScaling),dtype = np.bool)
        Shape[:,:5]  = True
        Shape[-5:,:] = True
        Shape[:,-5:] = True
        Number = np.uint16(self.ObstaclesDensity * self.Dimensions[0] * self.Dimensions[1] / self.ObstaclesSize**2)
        self.Obstacles[np.unravel_index(rng.integers(np.prod(self.Obstacles.shape),size = Number),self.Obstacles.shape)] = True
        self.Obstacles = ndimage.morphology.binary_dilation(self.Obstacles,Shape)
        self.Obstacles = ndimage.morphology.binary_fill_holes(self.Obstacles)
        ndimage.distance_transform_edt(self.Obstacles,return_indices=True,indices=self.ObstaclesIndices)

    def initialize_cross_shape_positions(self):
        Shape = np.zeros(np.uint16(np.ones(2) * self.ObstaclesSize / self.WorldScaling),dtype = np.bool)
        Shape[:,np.uint16(Shape.shape[1]/2 - 1.5):np.uint16(Shape.shape[1]/2 + 1.5)] = True
        Shape[np.uint16(Shape.shape[1]/2 - 1.5):np.uint16(Shape.shape[1]/2 + 1.5),:] = True
        Number = np.uint16(self.ObstaclesDensity * self.DimensionsScaled[0] * self.DimensionsScaled[1] / np.sum(Shape))
        self.Obstacles[np.unravel_index(rng.integers(np.prod(self.Obstacles.shape),size = Number),self.Obstacles.shape)] = True
        self.Obstacles = ndimage.morphology.binary_dilation(self.Obstacles,Shape)
        self.Obstacles = ndimage.morphology.binary_fill_holes(self.Obstacles)
        ndimage.distance_transform_edt(self.Obstacles,return_indices=True,indices=self.ObstaclesIndices)