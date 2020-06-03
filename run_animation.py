import cell_minds
import cell_bodies
import worlds
import lineages
import animate

#%%
CELL_NUMBER = 100
SIMULATION_TIME_STEP = 0.01

#%%
SampleInterval = 0
SampleNumber = 1

CellMind = cell_minds.EcoliMind(CELL_NUMBER, SampleNumber, SampleInterval)
CellBody = cell_bodies.EcoliBody(CELL_NUMBER, SampleNumber, SampleInterval)

Lineage = lineages.SimpleEvolution(CELL_NUMBER, SampleNumber, SampleInterval)

# World    = worlds.OpenCapillary(SampleNumber, SampleInterval)
# World    = worlds.StepLigand(int(SampleNumber/100), int(SampleInterval*100))
# World    = worlds.BoxExponential(SampleNumber, SampleInterval)
World    = worlds.ThreadmillExponential(SampleNumber, SampleInterval)
# World    = worlds.BoxLinear(int(SampleNumber/100), int(SampleInterval*100))
# World    = worlds.Maze_and_Diffusion(int(SampleNumber/100), int(SampleInterval*100))
# World    = worlds.Maze_and_Exponential(int(SampleNumber/100), int(SampleInterval*100))


#%%
FRAME_PERIOD = 10
animate.start_animation(CellMind, CellBody, World, Lineage, SIMULATION_TIME_STEP, FRAME_PERIOD)
