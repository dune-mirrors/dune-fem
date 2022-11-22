# Sources:
# https://kitware.github.io/paraview-docs/latest/python/paraview.util.vtkAlgorithm.html
# https://www.litianyi.me/2019/10/27/paraview-meshio-plugin/
# https://github.com/tianyikillua/paraview-meshio/blob/master/meshioPlugin.py
import numpy as np
import os,sys
import vtk
from paraview.util.vtkAlgorithm import VTKPythonAlgorithmBase
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

def setDuneModulePaths():
    # in older paraview versions there is no way to set the
    # virtual environment to use - use a environment variable
    # to set it before starting paraview:
    try:
        envdir = os.path.realpath(os.environ['VIRTUAL_ENV'])
        dunePaths = find_egglinks(os.path.join(envdir,"lib"))
        sys.path += dunePaths
        if not "DUNE_PY_DIR" in os.environ:
            os.environ["DUNE_PY_DIR"] = os.path.join(envdir,".cache")
        # print(os.environ["DUNE_PY_DIR"], dunePaths)
    except KeyError:
        # print("no virtual env path found!")
        pass


def find_egglinks(directory_name):
    dune_found = []
    for path, subdirs, files in os.walk(directory_name):
        if not path.endswith("site-packages"):
            continue
        dune_found.append(path)
        for name in files:
            if not "dune" in name:
                continue
            ext = os.path.splitext(name)[1]
            if ext == ".egg-link":
                file_path = os.path.join(path,name)
                with open(file_path,"r") as f:
                    dune_found.append(f.read().split()[0])
    return dune_found

dune_extensions = ["dbf"] # ,"dgf"]
@smproxy.reader(
    label="Dune Reader",
    extensions=dune_extensions,
    file_description="dune-supported files",
)
class DuneReader(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType="vtkUnstructuredGrid"
        )
        self._filename = None
        self._level = 0
        self._gridView = None
        setDuneModulePaths()
        # do not yet know how to fix the problem in paraview that
        # dune.generated in not found - the following hack works for some reason
        import dune.common
        dune.common.FieldVector([1])
        import dune.common.pickle
        self.load = dune.common.pickle.load

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(
        extensions=dune_extensions, file_description="dune supported files"
    )
    def SetFileName(self, filename):
        if self._filename != filename:
            self._filename = filename
            self.Modified()

    @smproperty.intvector(name="Level", default_values="0")
    @smdomain.intrange(min=0, max=5)
    def SetLevel(self, level):
        self._level = level
        self.Modified()

    def loadData(self):
        ext = os.path.splitext(self._filename)[1]
        if ext == ".dgf":
            print("Still need to implement dgf reading")
            print("Which grid to use with which dimensions?")
        else:
            with open(self._filename,"rb") as f:
                df = self.load(f)
            self._df = [d for d in df if hasattr(d,"gridView")]
            self._gridView = self._df[0].gridView

    def RequestData(self, request, inInfo, outInfo):
        if self._gridView is None:
            self.loadData()

        points, cells = self._gridView.tessellate(self._level)
        output = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(outInfo))

        # points need to be 3d:
        if self._gridView.dimWorld == 2:
            vtk_type = vtk.VTK_TRIANGLE
            points = np.hstack([points, np.zeros((len(points), 1))])
        elif self._gridView.dimWorld == 3:
            if self._gridView.dimGrid == 2:
                vtk_type = vtk.VTK_TRIANGLE
            else:
                vtk_type = vtk.VTK_TETRAHEDRON
        output.SetPoints(points)

        cell_types = np.array([], dtype=np.ubyte)
        cell_offsets = np.array([], dtype=int)
        cell_conn = np.array([], dtype=int)
        ncells, npoints = cells.shape
        cell_types = np.hstack(
                       [cell_types, np.full(ncells, vtk_type, dtype=np.ubyte)]
                     )
        offsets = len(cell_conn) + (1 + npoints) * np.arange(ncells, dtype=int)
        cell_offsets = np.hstack([cell_offsets, offsets])
        conn = np.hstack(
                   [npoints * np.ones((ncells, 1), dtype=int), cells]
               ).flatten()
        cell_conn = np.hstack([cell_conn, conn])
        output.SetCells(cell_types, cell_offsets, cell_conn)  # cell connectivities

        # data
        for df in self._df:
            array = df.pointData(self._level)
            output.PointData.append(array, df.name)  # point data
            # output.CellData.append(array, name)  # cell data
            # output.FieldData.append(array, name)  # field data

        return 1
