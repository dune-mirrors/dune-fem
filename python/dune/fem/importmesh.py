import sys
import numpy as np
from dune.grid import reader
from dune.fem.function import gridFunction

def meshDim(mesh):
    cell_names = {c.type.lower() for c in mesh.cells}
    three_d = {"tetra", "hexahedron", "pyramid", "wedge", "tetra10", "hexa20"}
    two_d = {"triangle", "quad", "triangle6", "quad9"}
    elementType="general"
    if cell_names & three_d:
        if len(cell_names & three_d) == 1:
            if cell_names & {"tetra"}: elementType="simplex"
            elif cell_names & {"hexahedron"}: elementType="cube"
        return 3, elementType
    if cell_names and (cell_names <= two_d or (cell_names & two_d and not (cell_names & three_d))):
        if len(cell_names & two_d) == 1:
            if cell_names & {"triangle"}: elementType="simplex"
            elif cell_names & {"quad"}: elementType="cube"
        return 2, elementType
    return 1, "simplex"

def mesh2DGF(points, cells, bndDomain = None, bndSegments = None, periodic = None, dim = None,
             computeFirstIndex = True, defaultBndId = None):
    """
    Parameter:
       points      array(list) of points of length dim
       cells       dict containing element vertex numbers
       bndDomain   dict id -> list[lower,upper] (or id -> str) where lower and upper describe the bounding box of the boundary section
       bndSegments dict id -> list of lists containing vertex numbers of boundary segments
       periodic    string containing periodic boundary transformation
       dim         dimension of grid
       computeFirstIndex if true, first vertex index is computed bases on all cells vertices (default is on)

    Returns
        String containing the mesh description in DGF format.
    """
    import numpy as np

    if defaultBndId:
        assert not bndDomain
        bndDomain = {defaultBndId:"default"}

    if dim is None:
        if "tetra" in cells or "hexahedron" in cells:
            dim = 3
        elif "quad" in cells or "triangle" in cells:
            dim = 2
        else:
            dim = len(points[0])

    simplex = "triangle" if "triangle" in cells else None
    if dim == 3 and "tetra" in cells:
        simplex = "tetra"

    # default numbering is 0 based
    firstVertexIndex = 0
    if computeFirstIndex:
        for key, cellVertices in cells.items():
            firstVertexIndex = min(firstVertexIndex, np.min(cellVertices))

    dgf="DGF\nVertex\n"
    # if first vertex index is not 0 we need to add this here
    if firstVertexIndex > 0:
        dgf += f"firstindex {firstVertexIndex}\n"

    for p in points:
        for i in range(dim):
            dgf += str(p[i]) + " "
        dgf += "\n"
    dgf += "#\n\n"

    if simplex is not None:
        dgf += "Simplex\n"
        for t in cells[simplex]:
            for v in t:
                dgf += str(v) + " "
            dgf += "\n"
        dgf += "#\n\n"

    if "quad" in cells and dim==2:
        dgf += "Cube\n"
        # gmsh has a different reference quadrilateral
        vxmap = [0,1,3,2] # flip vertex 2 and 3
        for t in cells["quad"]:
            for i in range(4):
                dgf += str(t[vxmap[i]]) + " "
            dgf += "\n"
        dgf += "#\n\n"
    if "hexahedron" in cells and dim==3:
        dgf += "Cube\n"
        # gmsh has a different reference hexahedron
        vxmap = [0,1,3,2,4,5,7,6] # flip vertex 2,3 and 6,7
        for t in cells["hexahedron"]:
            for i in range(8):
                dgf += str(t[vxmap[i]]) + " "
            dgf += "\n"
        dgf += "#\n\n"

    # boundary segments
    if bndSegments is not None:
        assert isinstance(bndSegments, dict), "Expecting a dictionary for boundary domain"
        dgf += "BoundarySegments\n"
        for bndid,bndsegs in bndSegments.items():
            for segment in bndsegs:
                dgf += str(bndid)
                for vx in segment:
                    dgf += " " + str(vx)
                dgf += "\n"
        dgf += "#\n\n"

    # boundary domain section
    if bndDomain is not None:
        assert isinstance(bndDomain, dict), "Expecting a dictionary for boundary domain"
        dgf += "BoundaryDomain\n"
        for bndid,bnd in bndDomain.items():
            if isinstance(bnd,str):
                if bnd == 'default':
                    dgf += bnd + " " + str(bndid) +"\n"
                else:
                    dgf += str(bndid) + " " + bnd +"\n"
            else: # tuple or list
                # has to be either list or tuple here
                assert isinstance(bnd,(tuple,list))
                assert len(bnd) == 2 # a lower left and upper right corner
                dgf += str(bndid)
                for coord in bnd:
                    assert len(coord) == dim # should a coordinate in the domain
                    for c in coord:
                        dgf += " " + str(c)
                dgf += "\n"

        dgf += "#\n"

    # periodic boundaries
    if periodic is not None:
        dgf += periodic

    return dgf

# remark: currently does not support surface grids
# Nice to have: would be good to have a way of setting a 'default' boundary id on missing boundaries
#               Also setting other properties, i.e., longest edge, type of refinement etc.
#               Possibly add a 'parameter' entry to the constructor dict with this information?
def importMesh(msh, ignoreInternalId=0, defaultBndId=None, bndDomain=None, generateDGF=False):
    try:
        import meshio
    except ImportError:
        raise ImportError("Function `importMesh` uses the `meshio` package - run `pip install meshio`")
    mesh = meshio.read(msh)
    dim, elementType = meshDim(mesh)

    if dim == 2:
        zCoord = mesh.points[:,2]
        assert np.isclose(min(zCoord),max(zCoord)) # check it's not a surface grid
        points = np.delete( mesh.points, 2, 1) # remove the z component from the points
        bndCells = ["line"]
    else:
        points = mesh.points
        bndCells = {"triangle", "quad", "triangle6", "quad9"}

    cells = mesh.cells_dict
    points = points.astype("float")

    segments = {}
    segs = []
    for cell in bndCells:
        # we prefer 'physical' tagging
        try:
            bnd = list(cells[cell])
        except KeyError:
            continue
        try:
            for i,(line,id) in enumerate( zip(bnd,mesh.cell_data_dict["gmsh:physical"][cell]) ):
                if id == ignoreInternalId: continue # inside skeleton tag
                if id <= 0: continue
                if id in segments:
                    segments[id] += [line]
                else:
                    segments[id] = [line]
                segs += [[id,*line]]
                bnd[i] = None
        except KeyError:
           pass
        # we prefer 'physical' but will try 'geometrical' tagging as well
        try:
            for line,id in zip(bnd,mesh.cell_data_dict["gmsh:geometrical"][cell]):
                if line is None: continue
                if id == ignoreInternalId: continue # inside skeleton tag
                if id <= 0: continue
                if id in segments:
                    segments[id] += [line]
                else:
                    segments[id] = [line]
                segs += [[id,*line]]
        except KeyError:
           pass
    segs = np.array(segs)

    if not generateDGF:
        if elementType == "simplex":
            cells = cells["triangle"] if dim==2 else cells["tetra"]
            minIndex = np.inf
            for c in cells:
                minIndex = min(minIndex,min(c))
            for c in cells:
                c -= minIndex
            for s in segs:
                s -= minIndex
            domain = {"vertices":points, "simplices":cells, "boundaries":segs}
        elif elementType == "cube":
            cells = cells["quad"] if dim==2 else cells["hexahedron"]
            if dim==2:
                for c in cells: # gmsh reader has different cube ordering (this could be done in alugrid's gridfactory)
                    c[2],c[3] = c[3],c[2]
            elif dim==3:
                for c in cells: # gmsh reader has different cube ordering (this could be done in alugrid's gridfactory)
                    c[2],c[3] = c[3],c[2]
                    c[6],c[7] = c[7],c[6]
            minIndex = np.inf
            for c in cells:
                minIndex = min(minIndex,min(c))
            for c in cells:
                c -= minIndex
            for s in segs:
                s -= minIndex
            domain = {"vertices":points, "cubes":cells, "boundaries":segs}
            if defaultBndId:
                domain["defaultBndId"] = defaultBndId;
        return {"constructor":domain, "dimgrid":dim, "elementType":elementType}
    else:
        domain = mesh2DGF(points, cells, bndSegments=segments, bndDomain=bndDomain, defaultBndId=defaultBndId)
        try:
            with open (generateDGF, "w") as file:
                file.write(domain)
        except:
            pass
        return {"constructor":(reader.dgfString, domain),
                "dimgrid":dim, "elementType":elementType}

# a simple projection for boundary ids
def projectBoundaryIds( gridView ):
    """
    Parameter:
       gridView     a grid view to project boundary ids for

    Returns:
       a piecewise discrete function containing values corresponding to
       adjacent boundaries. `0` refers to interior elements.
    """
    import io, sys
    from dune.fem.space import finiteVolume
    from dune.generator import algorithm

    code = """
    #include <dune/fem/misc/boundaryidprovider.hh>

    template <class GridView, class Intersection>
    int boundaryId( const GridView& gv, const Intersection& i )
    {
      return Dune::Fem::boundaryId( gv, i );
    }
    """
    bndId = None

    @gridFunction(gridView, order=1, name="bndId")
    def bndFunc(e,x):
        bndId = None
        maxId = 0
        if e.hasBoundaryIntersections:
            for i in gridView.intersections( e ):
                dist = np.linalg.norm( np.array(i.geometryInInside.toGlobal(i.geometryInInside.toLocal(x)))
                                      - np.array(x) )
                if bndId is None:
                    bndId = algorithm.load("boundaryId", io.StringIO(code), gridView, i )
                if i.boundary:
                    id = bndId( gridView, i )
                else:
                    id = 0
                if dist <= 1e-12:
                    maxId = max(maxId,id)
        return maxId

    return bndFunc

def main() -> int:
    from dune.alugrid import aluGrid as leafGridView
    import matplotlib.pyplot as plt
    msh = sys.argv[1]

    """
    # test the gmsh reader - boundary ids are missing and 2d cubes fail (3d cubes work...)
    # gridView = leafGridView((reader.gmsh,msh)) # , defaultBndId=7, dimgrid=2, elementType="simplex")
    gridView = leafGridView((reader.gmsh,msh), dimgrid=3, elementType="cube")
    bndIds = projectBoundaryIds(gridView)
    if gridView.dimension == 2:
        bndIds.plot(level=0, gridLines="white", cmap=plt.cm.jet, block=False)
    else:
        gridView.writeVTK(msh, pointdata=[bndIds], subsampling=0)
    """

    gridView = leafGridView((reader.meshio,msh), defaultBndId=5)
    bndIds = projectBoundaryIds(gridView)
    if gridView.dimension == 2:
        bndIds.plot(level=0, gridLines="white", cmap=plt.cm.jet, block=False)
    else:
        gridView.writeVTK(msh, pointdata=[bndIds], subsampling=0)

    domain = importMesh(msh, defaultBndId=5)
    if domain["elementType"] == "general":
        print("can't read in grid with general element type")
        return
    gridView = leafGridView(domain)
    bndIds = projectBoundaryIds(gridView)
    if domain["dimgrid"] == 2:
        bndIds.plot(level=0, gridLines="white", cmap=plt.cm.jet, block=False)
    else:
        gridView.writeVTK(msh, pointdata=[bndIds], subsampling=0)

    domain = importMesh(msh, generateDGF="test.dgf", defaultBndId=5)
    if domain["elementType"] == "general":
        print("can't read in grid with general element type")
        return
    gridView = leafGridView(domain)
    bndIds = projectBoundaryIds(gridView)
    if domain["dimgrid"] == 2:
        bndIds.plot(level=0, gridLines="white", cmap=plt.cm.jet, block=False)
    else:
        gridView.writeVTK(msh, pointdata=[bndIds], subsampling=0)

    # without a default - this will use '1'
    gridView = leafGridView((reader.meshio,msh))
    bndIds = projectBoundaryIds(gridView)
    if domain["dimgrid"] == 2:
        bndIds.plot(level=0, gridLines="white", cmap=plt.cm.jet, block=False)
    else:
        gridView.writeVTK(msh, pointdata=[bndIds], subsampling=0)

    if False: # without a default using DGF - this should segfault
        domain = importMesh(msh, generateDGF="test.dgf")
        if domain["elementType"] == "general":
            print("can't read in grid with general element type")
            return
        gridView = leafGridView(domain)
        bndIds = projectBoundaryIds(gridView)
        if domain["dimgrid"] == 2:
            bndIds.plot(level=0, gridLines="white", cmap=plt.cm.jet, block=False)
        else:
            gridView.writeVTK(msh, pointdata=[bndIds], subsampling=0)

    plt.show()

if __name__ == "__main__":
    raise SystemExit(main())
