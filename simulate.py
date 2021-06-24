"""
This library defines the functions necessary to initialize, start, run, or continue the simulation.
'run_sim_loop' runs the main simulation loop and returns all the objects defining the simulation with the records.
"""


class results:
    """Stores objects defining the simulation"""
    def __init__(self,CellMind,CellBody,World,Lineage,StepNumber,SimulationTimeStep,Clock):
        self.CellMind = CellMind
        self.CellBody = CellBody
        self.World = World
        self.Lineage = Lineage
        self.StepNumber = StepNumber
        self.SimulationTimeStep = SimulationTimeStep
        self.Clock = Clock
        

def start_simulation(CellMind,CellBody,World,Lineage,StepNumber,SimulationTimeStep):
    """Initialize cell positions and the state of the signaling pathway"""
    World.initialize_cell_positions(CellBody)
    CellBody.randomize_orientation()
    CellMind.update(CellBody,World,SimulationTimeStep)    
    return run_sim_loop(CellMind,CellBody,World,Lineage,StepNumber,SimulationTimeStep)


def continue_simulation(PreviousResults,StepNumber):
    """Takes the state of a finished simulation loop to continue"""
    SimulationTimeStep = PreviousResults.SimulationTimeStep
    CellMind = PreviousResults.CellMind
    CellBody = PreviousResults.CellBody
    World    = PreviousResults.World
    Lineage  = PreviousResults.Lineage
    CellMind.Records.reset(int(StepNumber/CellMind.Records.SampleInterval))
    CellBody.Records.reset(int(StepNumber/CellBody.Records.SampleInterval))
    World.Records.reset(int(StepNumber/World.Records.SampleInterval))
    Lineage.Records.reset(int(StepNumber/Lineage.Records.SampleInterval))
    return run_sim_loop(CellMind,CellBody,World,Lineage,StepNumber,SimulationTimeStep)


def run_sim_loop(CellMind,CellBody,World,Lineage,StepNumber,SimulationTimeStep):
    printProgressBar(0,StepNumber)
    """Main simulation loop"""
    for Step in range(StepNumber):        
        if (Step+1) % int(StepNumber/100) == 0:
            printProgressBar(Step+1,StepNumber)
        """Commit state to record if needed"""
        if (Step+1) % CellMind.Records.SampleInterval == 0:
            Clock = Step * SimulationTimeStep
            CellMind.record(Clock)
        if (Step+1) % CellBody.Records.SampleInterval == 0:
            Clock = Step * SimulationTimeStep    
            CellBody.record(Clock)
        if (Step+1) % World.Records.SampleInterval == 0:
            Clock = Step * SimulationTimeStep
            World.record(Clock)        
        if (Step+1) % Lineage.Records.SampleInterval == 0:
            Clock = Step * SimulationTimeStep
            Lineage.record(Clock)
        """Update all the classes by 1 step"""
        CellMind.update(CellBody,World,SimulationTimeStep)
        CellBody.update(CellMind,SimulationTimeStep)
        World.update(CellBody,SimulationTimeStep)
        Lineage.update(CellBody,CellMind,World) 
    """Return the final state of the simulation and all the records"""  
    return results(CellMind,CellBody,World,Lineage,StepNumber,SimulationTimeStep,Clock)


def printProgressBar (iteration,total,prefix = 'Progress:',suffix = 'Complete',decimals = 0,length = 50,fill = '█',printEnd = "\r"):
    """Print the progress of the simulation loop"""
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix,bar,percent,suffix),end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()