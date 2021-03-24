import matplotlib
from matplotlib import pyplot
from numpy import amin, amax, linspace, linalg, log, append, zeros, int32
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm

from dune.plotting import block, disable
from ufl import as_vector
globalBlock = block

def triangulationOfNetwork(grid, level=0, linewidth=0.01):
    from matplotlib.tri import Triangulation
    x, t = grid.tesselate(level)
    n, m = len(x), len(t)
    t = append(t.T, [zeros(m, dtype=int32)], axis=0).T
    t = append(t, zeros((m,3), dtype=int32), axis=0)
    x = append(x, zeros((2*m,2)), axis=0)
    for i in range(m):
        i0 = n+2*i
        t[i][2] = i0
        t[m+i][0] = t[i][1]
        t[m+i][1] = i0
        t[m+i][2] = i0+1
        x0 = x[t[i][0]]
        x1 = x[t[i][1]]
        d = x1 - x0
        d *= 0.5 * linewidth / (d[0]**2+d[1]**2)**(1/2)
        d = [d[1],-d[0]]
        x[i0] = x0 + d
        x[i0+1] = x1 + d
        x[t[i][0]] -= d
        x[t[i][1]] -= d
    return Triangulation(x[:,0], x[:,1], t)

def _plotPointData(fig, grid, solution, level=0, gridLines="black", linewidth=0.2, vectors=None,
        onlyContours=False, contours=None, contourWidth=2, contourColor="black",
        xlim=None, ylim=None, clim=None, cmap=None, colorbar="vertical",
        triplot=False, logscale=False, ticks=11):

    if colorbar == True:
        colorbar = "vertical"
    elif colorbar == False:
        colorbar = None

    if grid.dimWorld==2 and (gridLines is not None) and (gridLines != ""):
        polys = grid.polygons()
        for p in polys:
            coll = PolyCollection(p,facecolor='none',edgecolor=gridLines,linewidth=linewidth,zorder=2)
            pyplot.gca().add_collection(coll)

    if not solution == None:
        data = solution.pointData(level)

        if grid.dimGrid == 1 and grid.dimension == 1:
            if solution.dimRange > 1:
                data = linalg.norm(solution.pointData(level),axis=1)
            else:
                data = solution.pointData(level)[:,0]
            pyplot.plot(grid.tesselate(level)[0], data, '-p')
            if xlim:
                fig.gca().set_xlim(xlim)
            fig.tight_layout()
            return

        if grid.dimGrid == 1:
            triangulation = triangulationOfNetwork(grid, level, linewidth)
            n = len(data)
            data = append(data.T, [zeros(len(triangulation.triangles))], axis=1).T
            for i in range(n//2):
                for j in range(2):
                    data[n+2*i+j] = data[2*i+j]
        else:
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
                      units='xy', scale=None, zorder=3, color='black')
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
                        linewidth=linewidth, color='black')
            else:
                try:
                    pyplot.tricontourf(triangulation, data, cmap=cmap, levels=levels,
                            **logNorm,extend=extend)
                except:
                    pyplot.tricontourf(triangulation, data, cmap=cmap,
                            **logNorm,extend=extend,norm=norm)
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
                    cbar = pyplot.colorbar(**cbar, fraction=0.046, pad=0.04)
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
    fig.tight_layout()

from ufl.core.expr import Expr
from dune.ufl import expression2GF
def plotPointData(solution, figure=None, linewidth=0.1,
        level=0, gridLines="black", vectors=False,
        onlyContours=False, contours=None, contourWidth=2, contourColor="black",
        xlim=None, ylim=None, clim=None, cmap=None,
        colorbar="vertical", grid=None, triplot=False,
        block=globalBlock,
        logscale=False, ticks=11):
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
    if not grid.dimWorld == 2 and solution is None:
        print("inline plotting of grid only available for 2d grids")
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
    _plotPointData(figure, grid, solution, level, gridLines, linewidth,
                    vectors, onlyContours, contours, contourWidth, contourColor,
                    xlim, ylim, clim, cmap, colorbar, triplot,
                    logscale, ticks)

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
    if not grid.dimWorld == 2 and solution is None:
        print("inline plotting of grid only available for 2d grids")
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
        _plotPointData(figure,grid,None,level,gridLines,0.2,False,
                onlyContours, contours, contourWidth, contourColor,
                xlim,ylim,clim,cmap,colorbar)

    # add the data
    for p in show:
        pyplot.subplot(subfig+offset+p)
        _plotPointData(figure,grid,solution[p],level,"",0.2,False,
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
