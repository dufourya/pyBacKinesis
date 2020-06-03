import time
import pickle

import cell_minds
import cell_bodies
import worlds
import lineages
import simulate

#%%
CELL_NUMBER = 1000
TOTAL_SIMULATION_TIME = 2000
SIMULATION_TIME_STEP = 0.01
COLLECTION_TIME_STEP = 10

#%%
SampleInterval = int(COLLECTION_TIME_STEP/SIMULATION_TIME_STEP)
SampleNumber = int(TOTAL_SIMULATION_TIME/COLLECTION_TIME_STEP)
StepNumber = int(TOTAL_SIMULATION_TIME/SIMULATION_TIME_STEP)

CellMind = cell_minds.EcoliMind(CELL_NUMBER, SampleNumber, SampleInterval)
CellBody = cell_bodies.EcoliBody(CELL_NUMBER, SampleNumber, SampleInterval)

# World    = worlds.OpenCapillary(int(SampleNumber/100), int(SampleInterval*100))
# World    = worlds.StepLigand(int(SampleNumber/100), int(SampleInterval*100))
World    = worlds.BoxExponential(SampleNumber, SampleInterval)
# World    = worlds.BoxLinear(int(SampleNumber/100), int(SampleInterval*100))
# World    = worlds.Maze_and_Diffusion(int(SampleNumber/100), int(SampleInterval*100))
# World    = worlds.Maze_and_Exponential(int(SampleNumber/100), int(SampleInterval*100))

# Lineage = lineages.SimpleEvolution(CELL_NUMBER, SampleNumber, SampleInterval)

#%%
t0 = time.time()
Results_0 = simulate.start_simulation(CellMind, CellBody, World, Lineage, StepNumber, SIMULATION_TIME_STEP)
t1 = str(int(time.time() - t0))
print(t1 + " seconds")

pickle.dump(Results_0,open( "Results_0.p", "wb" ))

#%%
t0 = time.time()
Results_1 = simulate.continue_simulation(Results_0, StepNumber)
t1 = str(int(time.time() - t0))
print(t1 + " seconds")

pickle.dump(Results_1,open( "Results_1.p", "wb" ))

#%%
t0 = time.time()
Results_2 = simulate.continue_simulation(Results_1, StepNumber)
t1 = str(int(time.time() - t0))
print(t1 + " seconds")

pickle.dump(Results_2,open( "Results_2.p", "wb" ))

#%%
t0 = time.time()
Results_3 = simulate.continue_simulation(Results_2, StepNumber)
t1 = str(int(time.time() - t0))
print(t1 + " seconds")

pickle.dump(Results_3,open( "Results_3.p", "wb" ))

#%%
t0 = time.time()
Results_4 = simulate.continue_simulation(Results_3, StepNumber)
t1 = str(int(time.time() - t0))
print(t1 + " seconds")

pickle.dump(Results_4,open( "Results_4.p", "wb" ))