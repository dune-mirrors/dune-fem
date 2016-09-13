from __future__ import absolute_import, division, print_function, unicode_literals

import hashlib

from ...generator.generator import SimpleGenerator
from . import addAttr

generator = SimpleGenerator("GridPart", "dune/fempy/py", "Dune::FemPy")

def cppBool(value):
    return "true" if value else "false"

def create(hostGridPart, contains, useFilteredIndexSet=False):
    """create a FilteredGridPart.

    Args:
        hostGridPart:        grid part to filter
        contains:            function (Element -> bool) returning whether an element is contained in the resulting grid part
        useFilteredIndexSet: build index set containing only filtered entites? (defaults to false)

    Returns:
        GridPart: the constructed GridPart
    """
    includes = hostGridPart._module._includes + ["dune/fem/gridpart/filteredgridpart.hh", "dune/fem/gridpart/filter/simple.hh"]
    hostGridPartType = hostGridPart._module._typeName
    filterType = "Dune::Fem::SimpleFilter< " + hostGridPartType + " >"
    typeName = "Dune::Fem::FilteredGridPart< " + hostGridPartType + ", " + filterType + ", " + cppBool(useFilteredIndexSet) + " >"
    typeHash = "gridpart_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    constructors = [["[] ( " + typeName + " &self" + ", " + hostGridPartType + " &hostGridPart, pybind11::function contains ) {",
                     "    auto containsCpp = [ contains ] ( const " + hostGridPartType + "::Codim< 0 >::EntityType &e ) {",
                     "        return contains( e ).template cast< bool >();",
                     "      };",
                     "    new (&self) " + typeName + "( hostGridPart, " + filterType + "( hostGridPart, containsCpp ) );",
                     "  }, pybind11::keep_alive< 1, 2 >()"]]
    module = generator.load(includes, typeName, typeHash, constructors)
    addAttr(module, module.GridPart)
    return module.GridPart(hostGridPart, contains)

#############################################
if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
