import matplotlib
from matplotlib import pyplot
from numpy import amin, amax, linspace, linalg
from matplotlib.collections import PolyCollection

from dune.plotting import block, disable
globalBlock = block

def _plotPointData(fig, grid, solution, level=0, gridLines="black", vectors=None,
        xlim=None, ylim=None, clim=None, cmap=None, colorbar="vertical", triplot=False):

    if colorbar == True:
        colorbar = "vertical"
    elif colorbar == False:
        colorbar = None

    if (gridLines is not None) and (gridLines != ""):
        polys = grid.polygons()
        for p in polys:
            coll = PolyCollection(p,facecolor='none',edgecolor=gridLines,linewidth=0.5,zorder=2)
            pyplot.gca().add_collection(coll)

    if not solution == None:
        triangulation = grid.triangulation(level)
        data = solution.pointData(level)
        try:
            x1 = vectors[0]
            x2 = vectors[1]
            if x1 >= solution.dimRange or x2 >= solution.dimRange:
                vectors = None
        except:
            vectors = None

        if not vectors == None:
            pyplot.quiver(triangulation.x, triangulation.y, data[:,x1], data[:,x2],
                      units='xy', scale=10., zorder=3, color='blue',
                      width=0.007, headwidth=3., headlength=4.)
        else:
            if solution.dimRange > 1:
                data = linalg.norm(data,axis=1)
            else:
                data = data[:,0]
            minData = amin(data)
            maxData = amax(data)
            if clim == None:
                clim = [minData, maxData]
            levels = linspace(clim[0], clim[1], 256, endpoint=True)
            if triplot == True:
                pyplot.triplot(triangulation, antialiased=True, linewidth=0.2, color='black')
            else:
                try:
                    pyplot.tricontourf(triangulation, data, cmap=cmap, levels=levels, extend="both")
                except:
                    pyplot.tricontourf(triangulation, data, cmap=cmap, extend="both")
            if colorbar is not None:
                # having extend not 'both' does not seem to work (needs fixing)...
                if clim[0] > minData and clim[1] < maxData:
                    extend = 'both'
                elif clim[0] > minData:
                    extend = 'min'
                elif clim[1] < maxData:
                    extend = 'max'
                else:
                    extend = 'neither'
                v = linspace(clim[0], clim[1], 10, endpoint=True)
                norm = matplotlib.colors.Normalize(vmin=clim[0], vmax=clim[1])
                if not isinstance(colorbar,dict):
                    cbar = {}
                    cbar.setdefault("orientation",colorbar)
                else:
                    cbar = colorbar
                    cbar.setdefault("orientation","vertical")
                cbar.setdefault("shrink",1.0)
                cbar.setdefault("extend",extend)
                cbar.setdefault("norm",norm)
                cbar.setdefault("ticks",v)
                cbar = pyplot.colorbar(**cbar)
                # cbar = pyplot.colorbar(orientation="vertical",shrink=1.0, extend=extend, norm=norm, ticks=v)
                cbar.ax.tick_params(labelsize=10)

    fig.gca().set_aspect('equal')
    fig.gca().autoscale()
    if xlim:
        fig.gca().set_xlim(xlim)
    if ylim:
        fig.gca().set_ylim(ylim)

from ufl.core.expr import Expr
from dune.ufl import expression2GF
def plotPointData(solution, figure=None,
        level=0, gridLines="black", vectors=False,
        xlim=None, ylim=None, clim=None, cmap=None,
        colorbar="vertical", grid=None, triplot=False,
        block=globalBlock):
    if disable: return
    try:
        grid = solution.grid
    except AttributeError:
        if isinstance(solution, Expr) and grid is not None:
            assert grid, "need to provide a named grid argument to plot a ufl expression directly"
            solution = expression2GF(grid, solution, 1)
        else:
            grid = solution
            solution = None
    if not grid.dimension == 2:
        print("inline plotting so far only available for 2d grids")
        return

    if figure is None:
        figure = pyplot.figure()
        newFig = True
    else:
        try:
            subPlot = figure[1]
            figure = figure[0]
            pyplot.subplot(subPlot)
        except:
            pass
        newFig = False
    _plotPointData(figure, grid, solution, level, gridLines,
                    vectors, xlim, ylim, clim, cmap, colorbar, triplot)

    if newFig and block:
        pyplot.show(block=block)
    # return figure

def plotComponents(solution, level=0, show=None, gridLines="black",
        xlim=None, ylim=None, clim=None, cmap=None,
        block=globalBlock, **kwargs):
    if disable: return
    try:
        grid = solution.grid
    except AttributeError:
        if isinstance(solution, Expr):
            grid = kwargs.get("grid",None)
            assert grid, "need to provide a named grid argument to plot a ufl expression directly"
            solution = expression2GF(grid,solution,1)
        else:
            grid = solution
            solution = None
    if not grid.dimension == 2:
        print("inline plotting so far only available for 2d grids")
        return

    if not show:
        show = range(solution.dimRange)

    fig = pyplot.figure()
    if (gridLines is not None) and (gridLines != ""):
        offset = 1
    else:
        offset = 0
    subfig = 101+(len(show)+offset)*10

    # first the grid if required
    if (gridLines is not None) and (gridLines != ""):
        pyplot.subplot(subfig)
        _plotPointData(fig,grid,None,level,gridLines,False,xlim,ylim,clim,cmap)

    # add the data
    for p in show:
        pyplot.subplot(subfig+offset+p)
        _plotPointData(fig,grid,solution[p],level,"",False,xlim,ylim,clim,cmap,False)

    pyplot.show(block=globalBlock)
    # return fig

def mayaviPointData(grid, solution, level=0, component=0, block=globalBlock):
    if disable: return
    from mayavi import mlab
    triangulation = grid.triangulation(level)
    z = uh.pointData(level)[:,component]
    s = mlab.triangular_mesh(triangulation.x, triangulation.y, z,
                                triangulation.triangles)
    mlab.show(block=globalBlock)
