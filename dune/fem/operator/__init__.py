"""Functions for creating python modules and C++ classes for operators.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import hashlib

import sys
import logging
logger = logging.getLogger(__name__)

from ufl.equation import Equation
from ufl import Form

from dune.generator import Constructor, Method
from dune.generator.generator import SimpleGenerator

generator = SimpleGenerator("Operator", "Dune::FemPy")

def load(includes, typeName, *args):
    includes = includes + ["dune/fempy/py/operator.hh"]
    moduleName = "femoperator" + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, *args)
    return module


def galerkin(integrands, domainSpace, rangeSpace=None):
    if rangeSpace is None:
        rangeSpace = domainSpace

    modelParam = None
    if isinstance(integrands, (list, tuple)):
        modelParam = integrands[1:]
        integrands = integrands[0]
    if isinstance(integrands,Form):
        integrands = integrands == 0
    if isinstance(integrands,Equation):
        from dune.fem.model._models import integrands as makeIntegrands
        if modelParam:
            integrands = makeIntegrands(domainSpace.grid,integrands,*modelParam)
        else:
            integrands = makeIntegrands(domainSpace.grid,integrands)

    domainSpaceType = domainSpace._typeName
    rangeSpaceType = rangeSpace._typeName

    _, domainFunctionIncludes, domainFunctionType, _, _, _ = domainSpace.storage
    _, rangeFunctionIncludes, rangeFunctionType, _, _, _ = rangeSpace.storage

    includes = ["dune/fem/schemes/galerkin.hh", "dune/fempy/py/grid/gridpart.hh"]
    includes += domainSpace._includes + domainFunctionIncludes
    includes += rangeSpace._includes + rangeFunctionIncludes

    integrandsType = 'Dune::Fem::VirtualizedIntegrands< typename ' + rangeSpaceType + '::GridPartType, ' + integrands._domainValueType + ', ' + integrands._rangeValueType + ' >'

    typeName = 'Dune::Fem::GalerkinOperator< ' + integrandsType + ', ' + domainFunctionType + ', ' + rangeFunctionType + ' >'

    constructor = Constructor(['pybind11::object gridView', integrandsType + ' &integrands'],
                              ['return new ' + typeName + '( Dune::FemPy::gridPart< typename ' + rangeSpaceType + '::GridPartType::GridViewType >( gridView ), integrands );'],
                              ['"grid"_a', '"integrands"_a', 'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()'])

    return load(includes, typeName, constructor).Operator(rangeSpace.grid, integrands)
