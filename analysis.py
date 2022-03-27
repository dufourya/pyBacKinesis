"""
This library contains various functions to analyze the results of simulations
"""

"""External libary"""
import numpy as np
import matplotlib.pyplot as plt


def plot_results(Results):

    plt.style.use('fast')
    #%%
    fig = plt.figure()
    ax = fig.add_subplot()
    for k in range(0,Results.CellMind.CellNumber,10):
        x = Results.CellMind.Records.Time
        for kk in range(0,Results.CellMind.ReceptorNumber):
            y = Results.CellMind.Records.LigandConcentration[k,kk,:]
            ax.plot(x, y, alpha=0.5, linewidth = 0.5)
    #ax.set_aspect(1.0)
    ax.set_xlim([0, Results.CellMind.Records.Time[-1]])
    ax.set_yscale('log')
    #ax.set_ylim([0, 1])
    ax.axis('on')
    fig.savefig('plot_LigandConcentration.png', dpi=300, bbox_inches='tight')

    fig = plt.figure()
    ax = fig.add_subplot()
    for k in range(0,Results.CellMind.CellNumber,10):
        x = Results.CellMind.Records.Time
        for kk in range(0,Results.CellMind.ReceptorNumber):
            y = Results.CellMind.Records.MethylationAdapted[k,kk]
            ax.plot(x, y, alpha=0.5, linewidth = 0.5)
    #ax.set_aspect(1.0)
    #ax.set_xlim([0, Results.CellMind.Records.Time[-1]])
    #ax.set_ylim([0, 1])
    ax.axis('on')
    fig.savefig('plot_MethylationAdapted.png', dpi=300, bbox_inches='tight')


    fig = plt.figure()
    ax = fig.add_subplot()
    for k in range(0,Results.CellMind.CellNumber,10):
        x = Results.CellMind.Records.Time
        for kk in range(0,Results.CellMind.ReceptorNumber):
            y = Results.CellMind.Records.MethylationState[k,kk]
            ax.plot(x, y, alpha=0.5, linewidth = 0.5)
    #ax.set_aspect(1.0)
    #ax.set_xlim([0, Results.CellMind.Records.Time[-1]])
    #ax.set_ylim([0, 1])
    ax.axis('on')
    fig.savefig('plot_MethylationState.png', dpi=300, bbox_inches='tight')


    fig = plt.figure()
    ax = fig.add_subplot()
    for k in range(0,Results.CellMind.CellNumber,10):
        x = Results.CellMind.Records.Time
        for kk in range(0,Results.CellMind.ReceptorNumber):
            y = Results.CellMind.Records.ClusterRelativeActivity[k,kk]
            ax.plot(x, y, alpha=0.5, linewidth = 0.5)
    #ax.set_aspect(1.0)
    #ax.set_xlim([0, Results.CellMind.Records.Time[-1]])
    #ax.set_ylim([0, 1])
    ax.axis('on')
    fig.savefig('plot_ClusterRelativeActivity.png', dpi=300, bbox_inches='tight')

    fig = plt.figure()
    ax = fig.add_subplot()
    for k in range(0,Results.CellMind.CellNumber,10):
        x = Results.CellMind.Records.Time
        y = Results.CellMind.Records.ClusterTotalRelativeActivity[k,:]
        ax.plot(x, y, alpha=0.5, linewidth = 0.5)
    #ax.set_aspect(1.0)
    #ax.set_xlim([0, Results.CellMind.Records.Time[-1]])
    #ax.set_ylim([0, 1])
    ax.axis('on')
    fig.savefig('plot_ClusterTotalActivity.png', dpi=300, bbox_inches='tight')

    fig = plt.figure()
    ax = fig.add_subplot()
    for k in range(0,Results.CellBody.CellNumber,10):
        x = Results.CellBody.Records.Time
        y = Results.CellBody.Records.MotorG[k,:]
        ax.plot(x, y, alpha=0.5, linewidth = 0.5)
    #ax.set_aspect(1.0)
    #ax.set_xlim([0, Results.CellMind.Records.Time[-1]])
    #ax.set_ylim([0, 1])
    ax.axis('on')
    fig.savefig('plot_MortorG.png', dpi=300, bbox_inches='tight')

    fig = plt.figure()
    ax = fig.add_subplot()
    for k in range(0,Results.CellBody.CellNumber):
        x = Results.CellBody.Records.Position[k,0,:]
        y = Results.CellBody.Records.Position[k,1,:]
        ax.plot(x, y, alpha=0.5, linewidth = 0.5)
    ax.set_aspect(1.0)
    ax.set_xlim([0, Results.World.Dimensions[0]])
    ax.set_ylim([0, Results.World.Dimensions[1]])
    ax.axis('off')
    fig.savefig('plot_CellTrajectory.png', dpi=300, bbox_inches='tight')

    fig = plt.figure()
    ax = fig.add_subplot()
    x = Results.CellBody.Records.Position[:,0,-1]
    y = Results.CellBody.Records.Position[:,1,-1] 
    ax.scatter(x, y, s=2, alpha=0.5)
    x = Results.CellBody.Records.Position[:,0,1]
    y = Results.CellBody.Records.Position[:,1,1]
    ax.scatter(x, y, s=1, alpha=0.5)
    ax.set_aspect(1.0)
    ax.set_xlim([0, Results.World.Dimensions[0]])
    ax.set_ylim([0, Results.World.Dimensions[1]])
    ax.axis('off')
    fig.savefig('plot_FinalCellPosition.png', dpi=300, bbox_inches='tight')

    if hasattr(Results.World.Records, 'LigandConcentration'):
        fig = plt.figure()
        ax = fig.add_subplot(121)
        ax.imshow(np.log(Results.World.Records.LigandConcentration[:,:,0,-1].T), origin='lower')
        ax.axis('off')
        ax = fig.add_subplot(122)
        ax.imshow(np.log(Results.World.Records.LigandConcentration[:,:,1,-1].T), origin='lower')
        ax.axis('off')
        fig.savefig('plot_WorldLigandConcentration.png', dpi=300, bbox_inches='tight')

    if hasattr(Results.World, 'Obstacles'):
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.imshow(Results.World.Obstacles[:,:].T, origin='lower')
        ax.axis('off')
        fig.savefig('plot_WorldObstacles.png', dpi=300, bbox_inches='tight')

    fig = plt.figure()
    ax = fig.add_subplot()
    x = Results.CellMind.Records.CheY[:,0]
    ax.hist(x,np.linspace(0,20,40),histtype='step')
    # x = Results.CellMind.Records.CheY[:,50]
    # ax.hist(x,np.linspace(0,20,50))
    x = Results.CellMind.Records.CheY[:,-1]
    ax.hist(x,np.linspace(0,20,40),histtype='step')
    fig.savefig('plot_Distribution_CheY.png', dpi=300, bbox_inches='tight')

    fig = plt.figure()
    ax = fig.add_subplot()
    x = Results.CellMind.Records.ReceptorMethylationRate[:,:,0].flatten()
    ax.hist(x,np.linspace(0,1,40),histtype='step')
    x = Results.CellMind.Records.ReceptorMethylationRate[:,:,-1].flatten()
    ax.hist(x,np.linspace(0,1,40),histtype='step')
    fig.savefig('plot_Distribution_ReceptorMethylationRate.png', dpi=300, bbox_inches='tight')

    fig = plt.figure()
    ax = fig.add_subplot()
    x = Results.CellMind.Records.CheY[:,0]
    y = Results.CellMind.Records.ReceptorMethylationRate[:,0,0]
    ax.scatter(x,y)
    x = Results.CellMind.Records.CheY[:,-1]
    y = Results.CellMind.Records.ReceptorMethylationRate[:,0,-1]
    ax.scatter(x,y)
    fig.savefig('plot_Scatter_CheY_ReceptorMethylationRate.png', dpi=300, bbox_inches='tight')


def concatenate_results(*ResultsList):
    CatResults = ResultsList[-1]
    for Results in reversed(ResultsList[:-1]):
        for key in Results.CellMind.Records.__dict__.keys():
            AttrVal = getattr(Results.CellMind.Records, key)
            if type(AttrVal) is np.ndarray:
                CatAttrVal = np.append(AttrVal, getattr(CatResults.CellMind.Records, key), axis = AttrVal.ndim-1)
                setattr(CatResults.CellMind.Records, key, CatAttrVal)
        for key in Results.CellBody.Records.__dict__.keys():
            AttrVal = getattr(Results.CellBody.Records, key)
            if type(AttrVal) is np.ndarray:
                CatAttrVal = np.append(AttrVal, getattr(CatResults.CellBody.Records, key), axis = AttrVal.ndim-1)
                setattr(CatResults.CellBody.Records, key, CatAttrVal)
        for key in Results.World.Records.__dict__.keys():
            AttrVal = getattr(Results.World.Records, key)
            if type(AttrVal) is np.ndarray:
                CatAttrVal = np.append(AttrVal, getattr(CatResults.World.Records, key), axis = AttrVal.ndim-1)
                setattr(CatResults.World.Records, key, CatAttrVal)    
    return CatResults


def calculate_mean_drift_velocity(Results):
    return np.mean((Results.CellBody.Records.Position[:,:,-1] - Results.CellBody.Records.Position[:,:,0]) / (Results.Clock), axis = 0)