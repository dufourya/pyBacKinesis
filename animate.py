"""
This library defines the graphical window and different plots to animate the simulation.
"""

"""External libraries"""
import sys 
import numpy as np
import pyqtgraph as pg


class MainWindow(pg.GraphicsWindow):

    def __init__(self,World,CellBody,CellMind,Lineage,SimulationTimeStep,FramePeriod):
        super(MainWindow,self).__init__(title = 'Chemotaxis')
        self.World = World
        self.CellBody = CellBody
        self.CellMind = CellMind
        self.Lineage = Lineage
        self.SimulationTimeStep = SimulationTimeStep
        self.FramePeriod = FramePeriod
        self.has_ligand = hasattr(self.World,'LigandConcentration')
        self.has_obstacles = hasattr(self.World,'Obstacles')
        self.clock = 0
        self.ScreenGeometry = pg.QtGui.QDesktopWidget().screenGeometry()
        self.resize(self.ScreenGeometry.width()*0.75,self.ScreenGeometry.height()*0.75)
        self.move((self.ScreenGeometry.width() / 2) - (self.frameSize().width() / 2),
                  (self.ScreenGeometry.height() / 2) - (self.frameSize().height() / 2))
        self.initialize_simulation()
        self.initialize_plots()
        self.timer = pg.QtCore.QTimer()
        self.timer.timeout.connect(self.update_simulation)
        self.timer.start()

    def initialize_simulation(self):
        self.World.initialize_cell_positions(self.CellBody)
        self.World.record(0)
        self.CellBody.randomize_orientation()
        self.CellMind.update(self.CellBody,self.World,self.SimulationTimeStep)
    
    def initialize_plots(self):
        self.add_gradients()
        self.add_Cell_Position()
        self.add_CheY_MethylationRate()
        self.add_MethylationRate_MethylationRate()
        self.add_XPosition_Methylation()
        self.add_CheY_CheYp()
        # self.add_Age_LigandConsumedIntegrated()
        self.add_PoleAge_CellSize()

    def update_simulation(self):
        self.CellMind.update(self.CellBody,self.World,self.SimulationTimeStep)
        self.CellBody.update(self.CellMind,self.SimulationTimeStep)
        self.World.update(self.CellBody,self.SimulationTimeStep)
        self.Lineage.update(self.CellBody,self.CellMind,self.World)
        self.clock += 1
        if (self.clock % self.FramePeriod == 0):
            self.update_plots()

    def update_plots(self):
        self.scatter_Cell_Position.setData(self.CellBody.Position[:,0],self.CellBody.Position[:,1])
        self.scatter_CheY_MethylationRate1.setData(self.CellMind.CheY,np.log10(self.CellMind.ReceptorMethylationRate[:,0]))
        self.scatter_CheY_MethylationRate2.setData(self.CellMind.CheY,np.log10(self.CellMind.ReceptorMethylationRate[:,1]))
        self.scatter_MethylationRate_MethylationRate.setData(np.log10(self.CellMind.ReceptorMethylationRate[:,0]),np.log10(self.CellMind.ReceptorMethylationRate[:,1]))
        self.scatter_XPosition_Methylation1.setData(self.CellBody.Position[:,0],self.CellMind.Methylation[:,0])
        self.scatter_XPosition_Methylation2.setData(self.CellBody.Position[:,0],self.CellMind.Methylation[:,1])
        self.scatter_CheY_CheYp.setData(self.CellMind.CheY,self.CellMind.CheYp)
        # self.scatter_Age_LigandConsumedIntegrated1.setData(self.CellBody.CellPoleAge.max(axis=1),self.CellBody.LigandConsumedIntegrated[:,0])
        # self.scatter_Age_LigandConsumedIntegrated2.setData(self.CellBody.CellPoleAge.max(axis=1),self.CellBody.LigandConsumedIntegrated[:,1])
        self.scatter_PoleAge_CellSize.setData(self.CellBody.CellPoleAge.max(axis=1),self.CellBody.CellSize)
        if self.has_ligand:
            self.image_gradient1.setImage(image = np.log(np.mean(self.World.LigandConcentration[:,:,:,0],axis=2)))
            self.image_gradient2.setImage(image = np.log(np.mean(self.World.LigandConcentration[:,:,:,1],axis=2)))
    
    def add_Cell_Position(self):
        self.plot_Cell_Position = self.addPlot(row=2,col=0,colspan=5,title='Cells')
        self.plot_Cell_Position.setRange(xRange = [0,self.World.Dimensions[0]],yRange = [0,self.World.Dimensions[1]])
        if self.has_obstacles:
            self.image_obstacles = pg.ImageItem(image = self.World.Obstacles,levels=(0,5))      
            self.image_obstacles.setRect(pg.QtCore.QRect(0,0,self.World.Dimensions[0],self.World.Dimensions[1]))
            self.plot_Cell_Position.addItem(self.image_obstacles)
        self.scatter_Cell_Position = pg.ScatterPlotItem(x = self.CellBody.Position[:,0],y = self.CellBody.Position[:,1])
        self.scatter_Cell_Position.setSize(3)
        self.scatter_Cell_Position.setPen(None)
        self.scatter_Cell_Position.setBrush(100,200,150,200)
        self.plot_Cell_Position.addItem(self.scatter_Cell_Position)
        self.plot_Cell_Position.setAspectLocked()
        self.plot_Cell_Position.hideAxis('left')
        self.plot_Cell_Position.hideAxis('bottom')

    def add_gradients(self):
        self.plot_gradient1 = self.addPlot(row=0,col=0,colspan=5,title='Gradient 1')
        self.plot_gradient1.hideAxis('left')
        self.plot_gradient1.hideAxis('bottom')
        self.image_gradient1 = pg.ImageItem(image = np.log(self.World.Records. LigandConcentration[:,:,0,0]),autoLevels = True)
        self.plot_gradient1.addItem(self.image_gradient1)
        self.plot_gradient1.setAspectLocked()
        self.plot_gradient2 = self.addPlot(row=1,col=0,colspan=5,title='Gradient 2')
        self.plot_gradient2.hideAxis('left')
        self.plot_gradient2.hideAxis('bottom')
        self.image_gradient2 = pg.ImageItem(image = np.log(self.World.Records.LigandConcentration[:,:,1,0]),autoLevels = True)
        self.plot_gradient2.addItem(self.image_gradient2)
        self.plot_gradient2.setAspectLocked()

    def add_CheY_MethylationRate(self):
        self.plot_CheY_MethylationRate = self.addPlot(row=3,col=0,colspan=1,rowspan=1)
        self.plot_CheY_MethylationRate.setLabel('bottom',text='CheY',units = 'uM')
        self.plot_CheY_MethylationRate.setLabel('left',text='Methylation rate',units = '1/s')
        self.scatter_CheY_MethylationRate1 = pg.ScatterPlotItem(x = self.CellMind.CheY,y = np.log10(self.CellMind.ReceptorMethylationRate[:,0]))
        self.scatter_CheY_MethylationRate1.setSize(3)
        self.scatter_CheY_MethylationRate1.setPen(None)
        self.scatter_CheY_MethylationRate1.setBrush(255,102,0,200)
        self.scatter_CheY_MethylationRate2 = pg.ScatterPlotItem(x = self.CellMind.CheY,y = np.log10(self.CellMind.ReceptorMethylationRate[:,1]))
        self.scatter_CheY_MethylationRate2.setSize(3)
        self.scatter_CheY_MethylationRate2.setPen(None)
        self.scatter_CheY_MethylationRate2.setBrush(0,98,255,200)
        self.plot_CheY_MethylationRate.addItem(self.scatter_CheY_MethylationRate1)
        self.plot_CheY_MethylationRate.addItem(self.scatter_CheY_MethylationRate2)

    def add_MethylationRate_MethylationRate(self):
        self.plot_MethylationRate_MethylationRate = self.addPlot(row=3,col=1,colspan=1,rowspan=1)
        self.plot_MethylationRate_MethylationRate.setLabel('bottom',text='Methylation rate 1',units = '1/s')
        self.plot_MethylationRate_MethylationRate.setLabel('left',text='Methylation rate 2',units = '1/s')
        self.scatter_MethylationRate_MethylationRate = pg.ScatterPlotItem(x = np.log10(self.CellMind.ReceptorMethylationRate[:,0]),y = np.log10(self.CellMind.ReceptorMethylationRate[:,1]))
        self.scatter_MethylationRate_MethylationRate.setSize(3)
        self.scatter_MethylationRate_MethylationRate.setPen(None)
        self.scatter_MethylationRate_MethylationRate.setBrush(200,100,100,200)
        self.plot_MethylationRate_MethylationRate.addItem(self.scatter_MethylationRate_MethylationRate)

    def add_XPosition_Methylation(self):
        self.plot_XPosition_Methylation = self.addPlot(row=3,col=2,colspan=1,rowspan=1)
        self.plot_XPosition_Methylation.setLabel('bottom',text='X Position',units = 'um')
        self.plot_XPosition_Methylation.setLabel('left',text='Receptor Methylation',units = '#')
        self.plot_XPosition_Methylation.setRange(xRange = [0,self.World.Dimensions[0]],yRange = [self.CellMind.ReceptorMethylationMin,self.CellMind.ReceptorMethylationMax])
        self.scatter_XPosition_Methylation1 = pg.ScatterPlotItem(x = self.CellBody.Position[:,0],y = self.CellMind.Methylation[:,0])
        self.scatter_XPosition_Methylation1.setSize(3)
        self.scatter_XPosition_Methylation1.setPen(None)
        self.scatter_XPosition_Methylation1.setBrush(255,102,0,200)
        self.scatter_XPosition_Methylation2 = pg.ScatterPlotItem(x = self.CellBody.Position[:,0],y = self.CellMind.Methylation[:,1])
        self.scatter_XPosition_Methylation2.setSize(3)
        self.scatter_XPosition_Methylation2.setPen(None)
        self.scatter_XPosition_Methylation2.setBrush(0,98,255,200)
        self.plot_XPosition_Methylation.addItem(self.scatter_XPosition_Methylation1)
        self.plot_XPosition_Methylation.addItem(self.scatter_XPosition_Methylation2)

    def add_CheY_CheYp(self):
        self.plot_CheY_CheYp = self.addPlot(row=3,col=3,colspan=1,rowspan=1)
        self.plot_CheY_CheYp.setLabel('bottom',text='CheY',units = 'uM')
        self.plot_CheY_CheYp.setLabel('left',text='CheYp',units = 'uM')
        self.scatter_CheY_CheYp = pg.ScatterPlotItem(x = self.CellMind.CheY,y = self.CellMind.CheYp)
        self.scatter_CheY_CheYp.setSize(3)
        self.scatter_CheY_CheYp.setPen(None)
        self.scatter_CheY_CheYp.setBrush(200,100,100,200)
        self.plot_CheY_CheYp.addItem(self.scatter_CheY_CheYp)

    def add_Age_LigandConsumedIntegrated(self):
        self.plot_Age_LigandConsumedIntegrated = self.addPlot(row=3,col=4,colspan=1,rowspan=1)
        self.plot_Age_LigandConsumedIntegrated.setLabel('bottom',text='Age',units = 's')
        self.plot_Age_LigandConsumedIntegrated.setLabel('left',text='Ligand consumed',units = 'mol')
        self.scatter_Age_LigandConsumedIntegrated1 = pg.ScatterPlotItem(x = self.CellBody.CellPoleAge.max(axis=1),y = self.CellBody.LigandConsumedIntegrated[:,0])
        self.scatter_Age_LigandConsumedIntegrated1.setSize(3)
        self.scatter_Age_LigandConsumedIntegrated1.setPen(None)
        self.scatter_Age_LigandConsumedIntegrated1.setBrush(255,102,0,200)
        self.scatter_Age_LigandConsumedIntegrated2 = pg.ScatterPlotItem(x = self.CellBody.CellPoleAge.max(axis=1),y = self.CellBody.LigandConsumedIntegrated[:,1])
        self.scatter_Age_LigandConsumedIntegrated2.setSize(3)
        self.scatter_Age_LigandConsumedIntegrated2.setPen(None)
        self.scatter_Age_LigandConsumedIntegrated2.setBrush(0,98,255,200)
        self.plot_Age_LigandConsumedIntegrated.addItem(self.scatter_Age_LigandConsumedIntegrated1)
        self.plot_Age_LigandConsumedIntegrated.addItem(self.scatter_Age_LigandConsumedIntegrated2)

    def add_PoleAge_CellSize(self):
        self.plot_PoleAge_CellSize = self.addPlot(row=3,col=4,colspan=1,rowspan=1)
        self.plot_PoleAge_CellSize.setLabel('bottom',text='Pole Age',units = 's')
        self.plot_PoleAge_CellSize.setLabel('left',text='Cell Size',units = 'um3')
        self.scatter_PoleAge_CellSize = pg.ScatterPlotItem(x = self.CellBody.CellPoleAge.max(axis=1),y = self.CellBody.CellSize)
        self.scatter_PoleAge_CellSize.setSize(3)
        self.scatter_PoleAge_CellSize.setPen(None)
        self.scatter_CheY_CheYp.setBrush(200,100,100,200)
        self.plot_PoleAge_CellSize.addItem(self.scatter_PoleAge_CellSize)


def start_animation(CellMind,CellBody,World,Lineage,SimulationTimeStep,FramePeriod):
    """Create plots and start simulation"""
    App = pg.Qt.QtWidgets.QApplication(sys.argv)
    Main = MainWindow(World,CellBody,CellMind,Lineage,SimulationTimeStep,FramePeriod)
    Main.show()
    sys.exit(App.exec_())
