#!/usr/bin/env python

"""
Setup and run a simulation of agent-based bacterial chemotaxis in physically and chemically defined environments.
'cell_minds' defines the structure and dynamics of the chemotaxis signaling pathway.
'cell_bodies' defines the physical parameters and behavioral response of bacterial cell agents.
'worlds' defines the phyical and chemical structure of the simulated environment.
'lineages' defines the rule for evolving the cell population based on each agent performance.
'animate' runs the simulation loop to advance time steps and display the state of the simulation in live plots.
The output of the simulation is displayed but not recorded and saved.
"""

"""Internal libraries"""
import cell_minds
import cell_bodies
import worlds
import lineages
import animate
import parse_override_parameters

"""External libraries"""
import sys

def main():
    
    """Parameters for simulation"""  
    CELL_NUMBER = 1000
    SIMULATION_TIME_STEP = 0.01

    """Parameters for data recording"""    
    SampleInterval = 0
    SampleNumber = 1

    """Import additional parameters"""
    override_parameters = parse_override_parameters.parse_json(sys.argv[1:])
    
    """Initialize classes to build the simulation"""
    CellMind = cell_minds.EcoliMind(CELL_NUMBER,SampleNumber,SampleInterval,override_parameters.CellMind)        # chemotaxis signaling pathway based on Escherichia coli
    CellBody = cell_bodies.EcoliBody(CELL_NUMBER,SampleNumber,SampleInterval,override_parameters.CellBody)       # agent cell behavior based on Escherichia coli
    # Lineage = lineages.SimpleEvolution(CELL_NUMBER,SampleNumber,SampleInterval)   # simple rules for evolution by natural selection based on cell performance
    Lineage = lineages.NoEvolution(CELL_NUMBER,SampleNumber,SampleInterval,override_parameters.Lineage)         # no evolution of the original agent cell population
    # World    = worlds.OpenCapillary(SampleNumber,SampleInterval)                  # defines a capillary environment with 1 open side and reaction-diffusion dynamics
    # World    = worlds.StepLigand(SampleNumber,SampleInterval)                     # defines simple environment with a sudden change of chemical signals after some time
    # World    = worlds.BoxExponential(SampleNumber,SampleInterval)                 # defines a closed box with exponentail chemical gradients
    # World    = worlds.ThreadmillExponential(SampleNumber,SampleInterval)            # defines a closed box with exponentail chemical gradients and bulk flow
    # World    = worlds.BoxLinear(SampleNumber,SampleInterval)                      # defines a closed box with linear chemical gradients
    World    = worlds.BoxHillEquation(SampleNumber,SampleInterval)                # defines a closed box with sigmoidal chemical gradients
    # World    = worlds.Maze_and_Diffusion(SampleNumber,SampleInterval)             # defines a closed box with reaction-diffusion dynamics and physical obstacles
    # World    = worlds.Maze_and_Exponential(SampleNumber,SampleInterval)           # defines a closed box with exponentail chemical gradients and physical obstacles

    """Run live animation of the simulation"""
    FRAME_PERIOD = round(100/10)
    animate.start_animation(CellMind,CellBody,World,Lineage,SIMULATION_TIME_STEP,FRAME_PERIOD)


if __name__ == "__main__":
    main()
