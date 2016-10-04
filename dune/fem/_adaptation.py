from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib

from dune.generator.generator import SimpleGenerator

generator = SimpleGenerator("GridAdaptation", "Dune::FemPy")

modules = {}

def module(grid):
    try:
        return modules[grid._typeName]
    except KeyError:
        pass

    typeName = "Dune::FemPy::GridAdaptation< " + grid._typeName + " >"
    includes = grid._includes + ["dune/fempy/py/grid/adaptation.hh"]
    moduleName = "adapt_" + hashlib.md5(typeName.encode('utf8')).hexdigest()

    module = generator.load(includes, typeName, moduleName)
    modules[grid._typeName] = module
    return module


def adapt(grid, *args, **kwargs):
    module(grid).gridAdaptation(grid).adapt(*args, **kwargs)


def loadBalance(grid, *args, **kwargs):
    module(grid).gridAdaptation(grid).loadBalance(*args, **kwargs)


if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
