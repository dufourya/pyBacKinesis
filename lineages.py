"""
Library of classes defining how population of agent cells evolve during the simulation.
Each class must contain the following methods:
'update' defines how the population evolves at each time step.
'record' defines which paramters are recorded.
"""

"""External libaries"""
import numpy as np

rng = np.random.default_rng()

class LineageRecords:

    def __init__(self, CellNumber, SampleNumber, SampleInterval):
        self.SampleInterval = SampleInterval
        self.Time = np.zeros(SampleNumber)
        self.Iterator = 0
    
    def reset(self, SampleNumber):
        self.Time = np.zeros(SampleNumber)
        self.Iterator = 0


class SimpleEvolution:

    def __init__(self, CellNumber, SampleNumber, SampleInterval):
        self.CellNumber = CellNumber
        self.Records = LineageRecords(CellNumber, SampleNumber, SampleInterval)
        self.BirthIndices = None
        self.DeathIndices = None
        self.SdCheY = 0.5
        self.SdReceptorMethylationRate = 0.25
        # self.RespawnPosition = 0.25

    def update(self, CellBody, CellMind, World):
        self.update_birth_death(CellBody, World)
        if np.any(self.BirthIndices):
            self.evolve_population(CellMind, CellBody, World)

    def record(self, Clock):
        self.Records.Time[self.Records.Iterator] = Clock
        self.Records.Iterator += 1

    def update_birth_death(self, CellBody, World):
        # Winners = CellBody.Position[:,0] > (World.Dimensions[0] - 50)
        DividingCells = CellBody.CellSize > CellBody.CellMaximunSize
        # self.BirthIndices = np.nonzero(np.logical_or(Winners, DividingCells))[0]
        self.BirthIndices = np.nonzero(DividingCells)[0]
        # Losers = CellBody.Position[:,0] < 50
        DyingCells = np.any(CellBody.CellPoleAge > CellBody.CellPoleMaximumAge, axis=1)
        # self.DeathIndices = np.nonzero(np.logical_or(Losers, DyingCells))[0]
        self.DeathIndices = np.nonzero(DyingCells)[0]
        Difference = self.BirthIndices.__len__() - self.DeathIndices.__len__()
        if Difference > 0:
            self.DeathIndices = np.concatenate([self.DeathIndices, rng.integers(0, CellBody.CellNumber, Difference)])
        elif Difference < 0:
            self.BirthIndices = np.concatenate([self.BirthIndices, rng.integers(0, CellBody.CellNumber, -Difference)])

    def evolve_population(self, CellMind, CellBody, World):
        # CellBody.Position[self.BirthIndices, 0] = World.Dimensions[0] * self.RespawnPosition
        # CellBody.Position[self.DeathIndices, 0] = World.Dimensions[0] * self.RespawnPosition
        CellBody.Position[self.DeathIndices, :] = CellBody.Position[self.BirthIndices, :]
        CellMind.CheY[self.DeathIndices] = CellMind.CheY[self.BirthIndices] + \
            (rng.random(self.BirthIndices.__len__()) - 0.5) * self.SdCheY
        CellMind.CheY[CellMind.CheY < 1e-12] = 1e-12
        CellMind.ReceptorMethylationRate[self.DeathIndices,:] = CellMind.ReceptorMethylationRate[self.BirthIndices,:] * \
            np.exp((rng.random((self.BirthIndices.__len__(),CellMind.ReceptorNumber)) - 0.5) * self.SdReceptorMethylationRate)
        CellMind.ReceptorMethylationRate[CellMind.ReceptorMethylationRate < 0] = 0
        # CellMind.ResetReceptorMethylation[self.DeathIndices] = True
        # CellMind.ResetReceptorMethylation[self.BirthIndices] = True
        CellBody.CellPoleAge[self.DeathIndices,:] = CellBody.CellPoleAge[self.BirthIndices,:] * 1
        CellBody.CellPoleAge[self.DeathIndices,0] = 0
        CellBody.CellPoleAge[self.BirthIndices,1] = 0
        CellBody.LigandConsumedIntegrated[self.BirthIndices,:] = 0
        CellBody.LigandConsumedIntegrated[self.DeathIndices,:] = 0
        CellBody.CellSize[self.BirthIndices] = 1
        CellBody.CellSize[self.DeathIndices] = 1

class NoEvolution:

    def __init__(self, CellNumber, SampleNumber, SampleInterval):
        self.Records = LineageRecords(CellNumber, SampleNumber, SampleInterval)

    def update(self, CellBody, CellMind, World):
        pass

    def record(self, Clock):
        self.Records.Time[self.Records.Iterator] = Clock
        self.Records.Iterator += 1