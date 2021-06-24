#!/usr/bin/env python

"""
Setup and run a simulation of agent-based bacterial chemotaxis in physically and chemically defined environments.
'cell_minds' defines the structure and dynamics of the chemotaxis signaling pathway.
'cell_bodies' defines the physical parameters and behavioral response of bacterial cell agents.
'worlds' defines the phyical and chemical structure of the simulated environment.
'lineages' defines the rule for evolving the cell population based on each agent performance.
'simulate' runs the simulation loop to advance time steps and record the state of the simulation.
The output of the simulation, which contains records of many parameters, is saved in a data structure for further analysis.
"""

"""External libraries"""
import time
import pickle

"""Internal libraries"""
import cell_minds
import cell_bodies
import worlds
import lineages
import simulate


def main():

    """Parameters for simulation"""    
    CELL_NUMBER = 10
    TOTAL_SIMULATION_TIME = 500
    SIMULATION_TIME_STEP = 0.01
    COLLECTION_TIME_STEP = 1

    """Parameters for data recording"""    
    SampleInterval = int(COLLECTION_TIME_STEP/SIMULATION_TIME_STEP)
    SampleNumber = int(TOTAL_SIMULATION_TIME/COLLECTION_TIME_STEP)
    StepNumber = int(TOTAL_SIMULATION_TIME/SIMULATION_TIME_STEP)

    """Initialize classes to build the simulation"""
    CellMind = cell_minds.EcoliMind(CELL_NUMBER,SampleNumber,SampleInterval)        # chemotaxis signaling pathway based on Escherichia coli
    CellBody = cell_bodies.EcoliBody(CELL_NUMBER,SampleNumber,SampleInterval)       # agent cell behavior based on Escherichia coli
    # Lineage = lineages.SimpleEvolution(CELL_NUMBER,SampleNumber,SampleInterval)   # simple rules for evolution by natural selection based on cell performance
    Lineage = lineages.NoEvolution(CELL_NUMBER,SampleNumber,SampleInterval)         # no evolution of the original agent cell population
    # World    = worlds.OpenCapillary(SampleNumber,SampleInterval)                  # defines a capillary environment with 1 open side and reaction-diffusion dynamics
    # World    = worlds.StepLigand(SampleNumber,SampleInterval)                     # defines simple environment with a sudden change of chemical signals after some time
    # World    = worlds.BoxExponential(SampleNumber,SampleInterval)                 # defines a closed box with exponentail chemical gradients
    World    = worlds.ThreadmillExponential(SampleNumber,SampleInterval)            # defines a closed box with exponentail chemical gradients and bulk flow
    # World    = worlds.BoxLinear(SampleNumber,SampleInterval)                      # defines a closed box with linear chemical gradients
    # World    = worlds.BoxHillEquation(SampleNumber,SampleInterval)                # defines a closed box with sigmoidal chemical gradients
    # World    = worlds.Maze_and_Diffusion(SampleNumber,SampleInterval)             # defines a closed box with reaction-diffusion dynamics and physical obstacles
    # World    = worlds.Maze_and_Exponential(SampleNumber,SampleInterval)           # defines a closed box with exponentail chemical gradients and physical obstacles

    """Run the simulation and save records"""
    t0 = time.time()
    Results_0 = simulate.start_simulation(CellMind, CellBody, World, Lineage, StepNumber, SIMULATION_TIME_STEP)
    t1 = str(int(time.time() - t0))
    print(t1 + " seconds")
    pickle.dump(Results_0,open( "Results_0.p", "wb" ))

    """Continue the simulation from the previous state and save records"""
    # t0 = time.time()
    # Results_1 = simulate.continue_simulation(Results_0, StepNumber)
    # t1 = str(int(time.time() - t0))
    # print(t1 + " seconds")
    # pickle.dump(Results_1,open( "Results_1.p", "wb" ))


if __name__ == "__main__":
    main()