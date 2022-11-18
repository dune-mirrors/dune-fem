# Sources:
# https://kitware.github.io/paraview-docs/latest/python/paraview.util.vtkAlgorithm.html
# https://www.litianyi.me/2019/10/27/paraview-meshio-plugin/
# https://github.com/tianyikillua/paraview-meshio/blob/master/meshioPlugin.py

import numpy as np
import vtk
import pickle
from paraview.util.vtkAlgorithm import VTKPythonAlgorithmBase
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

import dune.common
from dune.alugrid import aluConformGrid as view
from dune.grid import cartesianDomain

dune_extensions = ["dbf","dgf"]
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
        # need to load at least one module - then dune.generated is
        # available. Actually not sure why this is needed or works...
        dune.common.FieldVector([1])
        view( cartesianDomain([-2,-2],[2,2],[10,10]) )

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

    def load(self):
        with open(self._filename,"rb") as f:
            self._df = pickle.load(f)
            """
            self._gridView = pickle.load(f).leafView
            space = lagrange(self._gridView, order=4)
            self._df = []
            self._df += [space.interpolate(0,name="test")]
            self._df[-1].read( pickle.load(f) )
            self._df += [space.interpolate(0,name="test2")]
            self._df[-1].read( pickle.load(f) )
            """
        self._gridView = self._df[0].gridView

    def RequestData(self, request, inInfo, outInfo):
        if self._gridView is None:
            self.load()

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
