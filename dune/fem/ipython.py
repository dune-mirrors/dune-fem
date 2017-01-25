def plotPointData(grid, solution, level=0, xlim=None, ylim=None, cmap=None):
    import matplotlib
    from matplotlib import pyplot
    from numpy import amin, amax, linspace
    from IPython.core.display import display

    triangulation = grid.triangulation(level)
    data = solution.pointData(level)

    levels = linspace(amin(data[:,0]), amax(data[:,0]), 256)

    fig = pyplot.figure()
    fig.gca().set_aspect('equal')
    if xlim:
        fig.gca().set_xlim(xlim)
    if ylim:
        fig.gca().set_ylim(ylim)
    pyplot.triplot(grid.triangulation(), antialiased=True, linewidth=0.2, color='black')
    pyplot.tricontourf(triangulation, data[:,0], cmap=cmap, levels=levels)
    display(fig)
