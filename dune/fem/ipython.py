import matplotlib
from matplotlib import pyplot
from numpy import amin, amax, linspace, linalg
from IPython.core.display import display

def plotPointData(grid, solution, level=0, component=0, showGrid=True, xlim=None, ylim=None, cmap=None):
    triangulation = grid.triangulation(level)
    data = solution.pointData(level)

    if component == -1:
        levels = linspace(amin(linalg.norm(data,axis=1)),
                          amax(linalg.norm(data,axis=1)), 256)
    else:
        levels = linspace(amin(data[:,component]), amax(data[:,component]), 256)

    fig = pyplot.figure()
    fig.gca().set_aspect('equal')
    if xlim:
        fig.gca().set_xlim(xlim)
    if ylim:
        fig.gca().set_ylim(ylim)
    if showGrid:
        pyplot.triplot(grid.triangulation(), antialiased=True, linewidth=0.2, color='black')

    if component == -1:
        pyplot.tricontourf(triangulation, linalg.norm(data,axis=1), cmap=cmap, levels=levels)
    elif component == -2:
        print(component)
        pyplot.quiver(triangulation.x, triangulation.y, data[:,0], data[:,1],
          units='xy', scale=10., zorder=3, color='blue',
          width=0.007, headwidth=3., headlength=4.)
    else:
        pyplot.tricontourf(triangulation, data[:,component], cmap=cmap, levels=levels)
    display(fig)


def plotComponents(grid, solution, show=None, showGrid=True):
    if not show:
        show = range(solution.dimRange)

    fig = pyplot.figure()
    triangulation = grid.triangulation()
    if showGrid:
        offset = 1
    else:
        offset = 0
    subfig = 101+(len(show)+offset)*10

    if showGrid:
        # plot the grid
        pyplot.subplot(subfig)
        pyplot.gca().set_aspect('equal')
        pyplot.gca().locator_params(tight=True, nbins=3)
        pyplot.triplot(triangulation, antialiased=True, linewidth=0.2, color='black')

    # add the data
    for p in show:
        pyplot.subplot(subfig+offset+p)
        pyplot.gca().set_aspect('equal')
        pyplot.gca().locator_params(tight=True, nbins=3)
        data = solution.pointData()
        levels = linspace(amin(data[:,p]), amax(data[:,p]), 256)
        pyplot.tricontourf(triangulation, data[:,p], cmap=pyplot.cm.rainbow, levels=levels)

    display(fig)

def

def mayaviPointData(grid, solution, level=0, component=0):
    from mayavi import mlab
    triangulation = grid.triangulation(level)
    z = uh.pointData(level)[:,component]
    s = mlab.triangular_mesh(triangulation.x, triangulation.y, z,
                                triangulation.triangles)
    mlab.show()
