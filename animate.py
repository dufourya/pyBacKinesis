import sys 
import numpy as np
import pyqtgraph as pg


class MainWindow(pg.GraphicsWindow):

    def __init__(self):
        super(MainWindow, self).__init__(title = 'Chemotaxis')

        self.ScreenGeometry = pg.QtGui.QDesktopWidget().screenGeometry()
        self.resize(self.ScreenGeometry.width()*0.75, self.ScreenGeometry.height()*0.75)
        self.move((self.ScreenGeometry.width() / 2) - (self.frameSize().width() / 2),
                  (self.ScreenGeometry.height() / 2) - (self.frameSize().height() / 2))
    
    def update(self, CellMind, CellBody, World, Lineage):
        self.scatter_Cell_Position.setData(CellBody.Position[:,0], CellBody.Position[:,1])
        self.scatter_CheY_MethylationRate1.setData(CellMind.CheY,np.log10(CellMind.ReceptorMethylationRate[:,0]))
        self.scatter_CheY_MethylationRate2.setData(CellMind.CheY,np.log10(CellMind.ReceptorMethylationRate[:,1]))
        self.scatter_MethylationRate_MethylationRate.setData(np.log10(CellMind.ReceptorMethylationRate[:,0]),np.log10(CellMind.ReceptorMethylationRate[:,1]))
        self.scatter_XPosition_Methylation1.setData(CellBody.Position[:,0],CellMind.Methylation[:,0])
        self.scatter_XPosition_Methylation2.setData(CellBody.Position[:,0],CellMind.Methylation[:,1])
        self.scatter_CheY_CheYp.setData(CellMind.CheY,CellMind.CheYp)
        # self.scatter_Age_LigandConsumedIntegrated1.setData(CellBody.CellPoleAge.max(axis=1), CellBody.LigandConsumedIntegrated[:,0])
        # self.scatter_Age_LigandConsumedIntegrated2.setData(CellBody.CellPoleAge.max(axis=1), CellBody.LigandConsumedIntegrated[:,1])
        self.scatter_PoleAge_CellSize.setData(CellBody.CellPoleAge.max(axis=1), CellBody.CellSize)
        pg.QtGui.QApplication.processEvents()
    
    def add_Cell_Position(self, CellBody, World):
        self.plot_Cell_Position = self.addPlot(row=2, col=0, colspan=5, title='Cells')
        self.plot_Cell_Position.setRange(xRange = [0, World.Dimensions[0]], yRange = [0, World.Dimensions[1]])
        self.scatter_Cell_Position = pg.ScatterPlotItem(x = CellBody.Position[:,0], y = CellBody.Position[:,1])
        self.scatter_Cell_Position.setSize(3)
        self.scatter_Cell_Position.setPen(None)
        self.scatter_Cell_Position.setBrush(255, 255, 255, 168)
        self.plot_Cell_Position.addItem(self.scatter_Cell_Position)
        self.plot_Cell_Position.setAspectLocked()
        self.plot_Cell_Position.hideAxis('left')
        self.plot_Cell_Position.hideAxis('bottom')
        # self.plot_Cell_Position.setMenuEnabled(enableMenu=False, enableViewBoxMenu='same')
        # self.plot_Cell_Position.hideButtons()

    def add_gradients(self, World):
        self.plot_gradient1 = self.addPlot(row=0, col=0, colspan=5, title='Gradient 1')
        self.plot_gradient1.hideAxis('left')
        self.plot_gradient1.hideAxis('bottom')
        self.image_gradient1 = pg.ImageItem(image = np.log(World.Records.LigandConcentration[:,:,0,0]), autoLevels = True)
        self.plot_gradient1.addItem(self.image_gradient1)
        self.plot_gradient1.setAspectLocked()

        self.plot_gradient2 = self.addPlot(row=1, col=0, colspan=5, title='Gradient 2')
        self.plot_gradient2.hideAxis('left')
        self.plot_gradient2.hideAxis('bottom')
        self.image_gradient2 = pg.ImageItem(image = np.log(World.Records.LigandConcentration[:,:,1,0]), autoLevels = True)
        self.plot_gradient2.addItem(self.image_gradient2)
        self.plot_gradient2.setAspectLocked()

    def add_CheY_MethylationRate(self, CellMind):
        self.plot_CheY_MethylationRate = self.addPlot(row=3,col=0,colspan=1, rowspan=1)
        self.plot_CheY_MethylationRate.setLabel('bottom', text='CheY', units = 'uM')
        self.plot_CheY_MethylationRate.setLabel('left', text='Methylation rate', units = '1/s')

        self.scatter_CheY_MethylationRate1 = pg.ScatterPlotItem(x = CellMind.CheY, y = np.log10(CellMind.ReceptorMethylationRate[:,0]))
        self.scatter_CheY_MethylationRate1.setSize(3)
        self.scatter_CheY_MethylationRate1.setPen(None)
        self.scatter_CheY_MethylationRate1.setBrush(255, 102, 0, 168)

        self.scatter_CheY_MethylationRate2 = pg.ScatterPlotItem(x = CellMind.CheY, y = np.log10(CellMind.ReceptorMethylationRate[:,1]))
        self.scatter_CheY_MethylationRate2.setSize(3)
        self.scatter_CheY_MethylationRate2.setPen(None)
        self.scatter_CheY_MethylationRate2.setBrush(0, 98, 255, 168)

        self.plot_CheY_MethylationRate.addItem(self.scatter_CheY_MethylationRate1)
        self.plot_CheY_MethylationRate.addItem(self.scatter_CheY_MethylationRate2)

    def add_MethylationRate_MethylationRate(self, CellMind):
        self.plot_MethylationRate_MethylationRate = self.addPlot(row=3,col=1,colspan=1, rowspan=1)
        self.plot_MethylationRate_MethylationRate.setLabel('bottom', text='Methylation rate 1', units = '1/s')
        self.plot_MethylationRate_MethylationRate.setLabel('left', text='Methylation rate 2', units = '1/s')

        self.scatter_MethylationRate_MethylationRate = pg.ScatterPlotItem(x = np.log10(CellMind.ReceptorMethylationRate[:,0]), y = np.log10(CellMind.ReceptorMethylationRate[:,1]))
        self.scatter_MethylationRate_MethylationRate.setSize(3)
        self.scatter_MethylationRate_MethylationRate.setPen(None)
        self.scatter_MethylationRate_MethylationRate.setBrush(255, 255, 255, 168)

        self.plot_MethylationRate_MethylationRate.addItem(self.scatter_MethylationRate_MethylationRate)

    def add_XPosition_Methylation(self, CellBody, CellMind, World):
        self.plot_XPosition_Methylation = self.addPlot(row=3, col=2, colspan=1, rowspan=1)
        # self.plot_CheY_MethylationRate.setMenuEnabled(enableMenu=False, enableViewBoxMenu='same')
        # self.plot_CheY_MethylationRate.hideButtons()
        self.plot_XPosition_Methylation.setLabel('bottom', text='X Position', units = 'um')
        self.plot_XPosition_Methylation.setLabel('left', text='Receptor Methylation', units = '#')
        self.plot_XPosition_Methylation.setRange(xRange = [0, World.Dimensions[0]], yRange = [CellMind.ReceptorMethylationMin, CellMind.ReceptorMethylationMax])

        self.scatter_XPosition_Methylation1 = pg.ScatterPlotItem(x = CellBody.Position[:,0], y = CellMind.Methylation[:,0])
        self.scatter_XPosition_Methylation1.setSize(3)
        self.scatter_XPosition_Methylation1.setPen(None)
        self.scatter_XPosition_Methylation1.setBrush(255, 102, 0, 168)

        self.scatter_XPosition_Methylation2 = pg.ScatterPlotItem(x = CellBody.Position[:,0], y = CellMind.Methylation[:,1])
        self.scatter_XPosition_Methylation2.setSize(3)
        self.scatter_XPosition_Methylation2.setPen(None)
        self.scatter_XPosition_Methylation2.setBrush(0, 98, 255, 168)

        self.plot_XPosition_Methylation.addItem(self.scatter_XPosition_Methylation1)
        self.plot_XPosition_Methylation.addItem(self.scatter_XPosition_Methylation2)

    def add_CheY_CheYp(self, CellMind):
        self.plot_CheY_CheYp = self.addPlot(row=3, col=3, colspan=1, rowspan=1)
        self.plot_CheY_CheYp.setLabel('bottom', text='CheY', units = 'uM')
        self.plot_CheY_CheYp.setLabel('left', text='CheYp', units = 'uM')

        self.scatter_CheY_CheYp = pg.ScatterPlotItem(x = CellMind.CheY, y = CellMind.CheYp)
        self.scatter_CheY_CheYp.setSize(3)
        self.scatter_CheY_CheYp.setPen(None)
        self.scatter_CheY_CheYp.setBrush(255, 255, 255, 168)

        self.plot_CheY_CheYp.addItem(self.scatter_CheY_CheYp)

    def add_Age_LigandConsumedIntegrated(self, CellBody):
        self.plot_Age_LigandConsumedIntegrated = self.addPlot(row=3, col=4, colspan=1, rowspan=1)
        self.plot_Age_LigandConsumedIntegrated.setLabel('bottom', text='Age', units = 's')
        self.plot_Age_LigandConsumedIntegrated.setLabel('left', text='Ligand consumed', units = 'mol')

        self.scatter_Age_LigandConsumedIntegrated1 = pg.ScatterPlotItem(x = CellBody.CellPoleAge.max(axis=1), y = CellBody.LigandConsumedIntegrated[:,0])
        self.scatter_Age_LigandConsumedIntegrated1.setSize(3)
        self.scatter_Age_LigandConsumedIntegrated1.setPen(None)
        self.scatter_Age_LigandConsumedIntegrated1.setBrush(255, 102, 0, 168)

        self.scatter_Age_LigandConsumedIntegrated2 = pg.ScatterPlotItem(x = CellBody.CellPoleAge.max(axis=1), y = CellBody.LigandConsumedIntegrated[:,1])
        self.scatter_Age_LigandConsumedIntegrated2.setSize(3)
        self.scatter_Age_LigandConsumedIntegrated2.setPen(None)
        self.scatter_Age_LigandConsumedIntegrated2.setBrush(0, 98, 255, 168)

        self.plot_Age_LigandConsumedIntegrated.addItem(self.scatter_Age_LigandConsumedIntegrated1)
        self.plot_Age_LigandConsumedIntegrated.addItem(self.scatter_Age_LigandConsumedIntegrated2)

    def add_PoleAge_CellSize(self, CellBody):
        self.plot_PoleAge_CellSize = self.addPlot(row=3, col=4, colspan=1, rowspan=1)
        self.plot_PoleAge_CellSize.setLabel('bottom', text='Pole Age', units = 's')
        self.plot_PoleAge_CellSize.setLabel('left', text='Cell Size', units = 'um3')

        self.scatter_PoleAge_CellSize = pg.ScatterPlotItem(x = CellBody.CellPoleAge.max(axis=1), y = CellBody.CellSize)
        self.scatter_PoleAge_CellSize.setSize(3)
        self.scatter_PoleAge_CellSize.setPen(None)
        self.scatter_CheY_CheYp.setBrush(255, 255, 255, 168)

        self.plot_PoleAge_CellSize.addItem(self.scatter_PoleAge_CellSize)

def start_animation(CellMind, CellBody, World, Lineage, SimulationTimeStep, FramePeriod):

    FramePeriod = np.around(FramePeriod/SimulationTimeStep).astype(int)
    # Initialize cell position
    World.initialize_cell_positions(CellBody)
    World.record(0)
    CellBody.randomize_orientation()
    CellMind.update(CellBody, World, SimulationTimeStep)

    # Create plot window
    App = pg.Qt.QtWidgets.QApplication(sys.argv)
    Main = MainWindow()
    Main.add_gradients(World)
    Main.add_Cell_Position(CellBody, World)
    Main.add_CheY_MethylationRate(CellMind)
    Main.add_MethylationRate_MethylationRate(CellMind)
    Main.add_XPosition_Methylation(CellBody, CellMind, World)
    Main.add_CheY_CheYp(CellMind)
    # Main.add_Age_LigandConsumedIntegrated(CellBody)
    Main.add_PoleAge_CellSize(CellBody)
    Main.show()

    # Simulation loop
    Clock = 0

    while True:

        # Update classes
        CellMind.update(CellBody, World, SimulationTimeStep)
        CellBody.update(CellMind, SimulationTimeStep)
        World.update(CellBody, SimulationTimeStep)
        Lineage.update(CellBody, CellMind, World)

        # Update plots
        if Clock % FramePeriod == 0:
            Main.update(CellMind, CellBody, World, Lineage)

        Clock += 1

    sys.exit(App.exec_())
