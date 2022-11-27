# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

###############################################################################
# This paraview reader adds support for 'dune binary format' # files (dbf).
# The file is assumed to be written using 'dune.common.pickle.dump'. It
# therefore consists of two parts (the required jit module source code and
# a pickled list of objects). This list is searched for objects containing
# a 'gridView' attribute - these are all assumed to be grid functions
# over the same grid view and with a 'pointData' attribute.
# If no entry in the list with a 'gridView' attribute is found the first
# entry is assumed to be a grid view and only the grid is plotted.
###############################################################################

import numpy as np
import os,sys,vtk,importlib
from paraview.util.vtkAlgorithm import VTKPythonAlgorithmBase
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

# In older paraview versions there is no way to set the
# virtual environment to use - use a environment variable
# to set it before starting paraview.
# This should be improved to take other usages into account.

# This finds all 'egg-link' files in a given folder structure.
# These correspond to packages installed 'editable' and need to
# be added by hand to the Python search path:
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
# We find an active virtual env by checking if the environment variable
# 'VIRTUAL_ENV' is set - this work at least with activated environments
# setup with 'venv' on Linux:
def setDuneModulePaths():
    try:
        envdir = os.path.realpath(os.environ['VIRTUAL_ENV'])
        dunePaths = find_egglinks(os.path.join(envdir,"lib"))
        sys.path += dunePaths
        if not "DUNE_PY_DIR" in os.environ:
            os.environ["DUNE_PY_DIR"] = os.path.join(envdir,".cache")
        sys.path += os.path.join(os.environ["DUNE_PY_DIR"],"python","dune","generated")
        # print(os.environ["DUNE_PY_DIR"], dunePaths)
    except KeyError:
        # print("no virtual env path found!")
        pass

############################################



# Actual reader
# -------------
# Some documentation
# https://kitware.github.io/paraview-docs/latest/python/paraview.util.vtkAlgorithm.html
# https://github.com/Kitware/ParaView/blob/master/Examples/Plugins/PythonAlgorithm/PythonAlgorithmExamples.py

dune_extensions = ["dbf"] # ,"dgf"]
@smproxy.reader(
    label="Dune Reader",
    extensions=dune_extensions,
    file_description="dune-supported files",
)
class DuneReader(VTKPythonAlgorithmBase):
    def __init__(self, parent=None):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType="vtkUnstructuredGrid"
        )
        self._filename = None
        self._level = 0
        self._transform = None
        self._transformFcts = []
        self._transformFct = ""
        self._gridView = None
        setDuneModulePaths()
        import dune.common.pickle
        self.load = dune.common.pickle.load
        # self.AddObserver('ModifiedEvent', MyObserver())
        # self.AddObserver(vtkCommand.KeyPressEvent,MyObserver())

    def keyPress(self, obj, event):
        key = obj.GetKeySym()
        print( self.GetClassName(), event, key)

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(
        extensions=dune_extensions, file_description="dune supported files"
    )
    def SetFileName(self, filename):
        if self._filename != filename:
            self._filename = filename
            self.Modified()

    @smproperty.stringvector(name="Transform")
    def SetTransform(self, transform):
        if transform == "None":
            return
        try:
            mod = importlib.import_module(transform)
            transformFcts = [m.__name__ for m in mod.register]
            transformFcts[:0] = ["None"]
            transformFct  = transformFcts[0]
            self._transform = mod
            self._transformFcts = transformFcts
            self._transformFct  = transformFct
            self.Modified()
        except:
            print("Failed to import script",transform)

    @smproperty.stringvector(name="TransformFct", information_only="1")
    def getTransformFcts(self):
        return self._transformFcts
        """
        if self._transform is None:
            return []
        return [ d for d in dir(self._transform)
                 if d.startswith("trans_") and callable(getattr(self._transform,d)) ]
        """
    @smproperty.stringvector(name="transfnc", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="TransformFct" function="TransformFct"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def setTransformFcts(self, value):
        self._transformFct = value
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
            if len(self._df) > 0:
                self._gridView = self._df[0].gridView
            else:
                self._gridView = df[0]
            # make some checks:
            assert hasattr(self._gridView,"dimension"), "file read contains no valid grid view"
            assert all( [hasattr(d,"pointData") for d in self._df] ), "found a non valid grid function (no 'pointData' attribute"
            assert all( [self._gridView == d.gridView for d in self._df] ), "all grid function must be over the same gridView"

    def RequestData(self, request, inInfo, outInfo):
        if self._gridView is None:
            self.loadData()
        # data
        if self._transform is not None and self._transformFct != "None":
            gfs = getattr(self._transform,self._transformFct)(self._gridView, self._df)
        else:
            gfs = self._df

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
                vtk_type = vtk.VTK_TETRA
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

        for df in gfs:
            array = df.pointData(self._level)
            output.PointData.append(array, df.name)  # point data
            # output.CellData.append(array, name)  # cell data
            # output.FieldData.append(array, name)  # field data

        return 1
