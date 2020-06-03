# Module for simulation protocols
class results:

    def __init__(self, CellMind, CellBody, World, Lineage, StepNumber, SimulationTimeStep, Clock):
        self.CellMind = CellMind
        self.CellBody = CellBody
        self.World = World
        self.Lineage = Lineage
        self.StepNumber = StepNumber
        self.SimulationTimeStep = SimulationTimeStep
        self.Clock = Clock
        

def start_simulation(CellMind, CellBody, World, Lineage, StepNumber, SimulationTimeStep):

    # Initialize cell position
    World.initialize_cell_positions(CellBody)
    CellBody.randomize_orientation()
    CellMind.update(CellBody, World, SimulationTimeStep)
    
    printProgressBar(0, StepNumber)

    # Simulation loop
    for Step in range(StepNumber):
        
        if (Step+1) % int(StepNumber/100) == 0:
            printProgressBar(Step+1, StepNumber)

        # Data collection
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

        # Update classes
        CellMind.update(CellBody, World, SimulationTimeStep)
        CellBody.update(CellMind, SimulationTimeStep)
        World.update(CellBody, SimulationTimeStep)
        Lineage.update(CellBody, CellMind, World)
    
    Clock = StepNumber * SimulationTimeStep

    # Return a result class with simulation samples
    return results(CellMind, CellBody, World, Lineage, StepNumber, SimulationTimeStep, Clock)


def continue_simulation(PreviousResults, StepNumber):

    SimulationTimeStep = PreviousResults.SimulationTimeStep

    CellMind = PreviousResults.CellMind
    CellBody = PreviousResults.CellBody
    World    = PreviousResults.World
    Lineage  = PreviousResults.Lineage
    
    CellMind.Records.reset(int(StepNumber/CellMind.Records.SampleInterval))
    CellBody.Records.reset(int(StepNumber/CellBody.Records.SampleInterval))
    World.Records.reset(int(StepNumber/World.Records.SampleInterval))
    Lineage.Records.reset(int(StepNumber/Lineage.Records.SampleInterval))

    printProgressBar(0, StepNumber)

    # Simulation loop
    for Step in range(StepNumber):
        
        if (Step+1) % int(StepNumber/100) == 0:
            printProgressBar(Step+1, StepNumber)

        # Data collection
        if (Step+1) % CellMind.Records.SampleInterval == 0:
            Clock = Step * SimulationTimeStep + PreviousResults.Clock
            CellMind.record(Clock)

        if (Step+1) % CellBody.Records.SampleInterval == 0:
            Clock = Step * SimulationTimeStep + PreviousResults.Clock
            CellBody.record(Clock)

        if (Step+1) % World.Records.SampleInterval == 0:
            Clock = Step * SimulationTimeStep + PreviousResults.Clock
            World.record(Clock)

        if (Step+1) % Lineage.Records.SampleInterval == 0:
            Clock = Step * SimulationTimeStep + PreviousResults.Clock
            Lineage.record(Clock)
        
        # Update classes
        CellMind.update(CellBody, World, SimulationTimeStep)
        CellBody.update(CellMind, World, SimulationTimeStep)
        World.update(CellBody, SimulationTimeStep)
        Lineage.update(CellBody, CellMind, World)
    
    Clock = StepNumber * SimulationTimeStep + PreviousResults.Clock

    # Return a result class with simulation samples
    return results(CellMind, CellBody, World, Lineage, StepNumber, SimulationTimeStep, Clock)


# Print iterations progress
def printProgressBar (iteration, total, prefix = 'Progress:', suffix = 'Complete', decimals = 0, length = 50, fill = '█', printEnd = "\r"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()