import matplotlib
from matplotlib import pyplot
from numpy import amin, amax, linspace, linalg, log
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm

from dune.plotting import block, disable
from ufl import as_vector
globalBlock = block

def _plotPointData(fig, grid, solution, level=0, gridLines="black", vectors=None,
        onlyContours=False, contours=None, contourWidth=2, contourColor="black",
        xlim=None, ylim=None, clim=None, cmap=None, colorbar="vertical",
        triplot=False, logscale=False):

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

        if not vectors is None:
            pyplot.quiver(triangulation.x, triangulation.y, data[:,x1], data[:,x2],
                      units='xy', scale=10., zorder=3, color='blue',
                      width=0.007, headwidth=3., headlength=4.)
            return
        else:
            if solution.dimRange > 1:
                data = linalg.norm(data,axis=1)
            else:
                data = data[:,0]
        if not onlyContours:
            if logscale:
                data = abs(data)
                logNorm = {"norm":LogNorm()}
            else:
                logNorm = {}
            minData = amin(data)
            maxData = amax(data)
            if clim == None:
                clim = [minData, maxData]
            # having extend not 'both' does not seem to work (needs fixing)...
            if clim[0] > minData and clim[1] < maxData:
                extend = 'both'
            elif clim[0] > minData:
                extend = 'min'
            elif clim[1] < maxData:
                extend = 'max'
            else:
                extend = 'neither'
            norm = matplotlib.colors.Normalize(vmin=clim[0], vmax=clim[1])
            levels = linspace(clim[0], clim[1], 256, endpoint=True)
            if triplot == True:
                pyplot.triplot(triangulation, antialiased=True,
                        linewidth=0.2, color='black')
            else:
                try:
                    pyplot.tricontourf(triangulation, data, cmap=cmap, levels=levels,
                            **logNorm,extend=extend)
                except:
                    pyplot.tricontourf(triangulation, data, cmap=cmap,
                            **logNorm,extend="both",norm=norm)
                if colorbar is not None:
                    v = linspace(clim[0], clim[1], 10, endpoint=True)
                    if not isinstance(colorbar,dict):
                        cbar = {}
                        cbar.setdefault("orientation",colorbar)
                    else:
                        cbar = colorbar
                        cbar.setdefault("orientation","vertical")
                    cbar.setdefault("shrink",1.0)
                    cbar.setdefault("ticks",v)
                    cbar = pyplot.colorbar(**cbar)
                    cbar.ax.tick_params(labelsize=10)
        if contours is not None:
            pyplot.tricontour(triangulation, data, levels=contours,
                    colors=contourColor, linewidths=contourWidth)

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
        onlyContours=False, contours=None, contourWidth=2, contourColor="black",
        xlim=None, ylim=None, clim=None, cmap=None,
        colorbar="vertical", grid=None, triplot=False,
        block=globalBlock,
        logscale=False):
    if disable: return
    try:
        grid = solution.grid
    except AttributeError:
        if isinstance(solution, list) or isinstance(solution,tuple):
            solution = as_vector(solution)
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
                    vectors, onlyContours, contours, contourWidth, contourColor,
                    xlim, ylim, clim, cmap, colorbar, triplot,
                    logscale)

    if newFig and block:
        pyplot.show(block=block)
    # return figure

def plotComponents(solution, figure=None, level=0, show=None, gridLines="black",
        onlyContours=False, contours=None, contourWidth=2, contourColor="black",
        xlim=None, ylim=None, clim=None, cmap=None,
        block=globalBlock, grid=None, colorbar=None):
    if disable: return
    try:
        grid = solution.grid
    except AttributeError:
        if isinstance(solution, Expr):
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

    if figure is None:
        figure = pyplot.figure()
        newFig = True
    else:
        newFig = False
    if (gridLines is not None) and (gridLines != ""):
        offset = 1
    else:
        offset = 0
    subfig = 101+(len(show)+offset)*10

    # first the grid if required
    if (gridLines is not None) and (gridLines != ""):
        pyplot.subplot(subfig)
        _plotPointData(figure,grid,None,level,gridLines,False,
                onlyContours, contours, contourWidth, contourColor,
                xlim,ylim,clim,cmap,colorbar)

    # add the data
    for p in show:
        pyplot.subplot(subfig+offset+p)
        _plotPointData(figure,grid,solution[p],level,"",False,
                onlyContours, contours, contourWidth, contourColor,
                xlim,ylim,clim,cmap,colorbar)

    if newFig and block:
        pyplot.show(block=block)
    # return fig

def mayaviPointData(grid, solution, level=0, component=0, block=globalBlock):
    if disable: return
    from mayavi import mlab
    triangulation = grid.triangulation(level)
    z = uh.pointData(level)[:,component]
    s = mlab.triangular_mesh(triangulation.x, triangulation.y, z,
                                triangulation.triangles)
    mlab.show(block=globalBlock)
