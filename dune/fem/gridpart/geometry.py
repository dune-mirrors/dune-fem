from __future__ import absolute_import, division, print_function, unicode_literals

import hashlib

from ...generator.generator import SimpleGenerator
from . import addAttr

generator = SimpleGenerator("GridPart", "dune/fempy/py", "Dune::FemPy")

def create(coordFunction):
    """create a GeometryGridPart.

    Args:
        coordFunction: coordinate function of the constructed GeometryGridPart

    Returns:
        GridPart: the constructed GridPart
    """
    includes = coordFunction._module._includes + ["dune/fem/gridpart/geometrygridpart.hh"]
    coordFunctionType = coordFunction._module._typeName
    typeName = "Dune::Fem::GeometryGridPart< " + coordFunctionType + " >"

    constructor = ["[] ( " + typeName + " &self" + ", " + coordFunctionType + " &coordFunction ) {",
                   "    new (&self) " + typeName + "( coordFunction );",
                   "  }, pybind11::keep_alive< 1, 2 >()"]

    typeHash = "gridpart_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, typeHash, [constructor])
    addAttr(module, module.GridPart)
    return module.GridPart(coordFunction)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
