"""Functions for creating python modules and C++ classes for operators.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import hashlib

import sys
import logging
logger = logging.getLogger(__name__)

from dune.generator.generator import SimpleGenerator

generator = SimpleGenerator("Operator", "Dune::FemPy")

def load(includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/fempy/py/operator.hh"]
    moduleName = "femoperator" + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, constructors, methods)
    return module


def galerkin(integrands, domainSpace, rangeSpace=None):
    if rangeSpace is None:
        rangeSpace = domainSpace

    domainSpaceType = domainSpace._typeName
    rangeSpaceType = rangeSpace._typeName

    _, domainFunctionIncludes, domainFunctionType, _, _ = domainSpace.storage
    _, rangeFunctionIncludes, rangeFunctionType, _, _ = rangeSpace.storage

    includes = ["dune/fem/schemes/galerkin.hh", "dune/fempy/py/grid/gridpart.hh"]
    includes += domainSpace._includes + domainFunctionIncludes
    includes += rangeSpace._includes + rangeFunctionIncludes

    domainValueType = 'std::tuple< typename ' + domainSpaceType + '::RangeType, typename ' + domainSpaceType + '::JacobianRangeType >'
    rangeValueType = 'std::tuple< typename ' + rangeSpaceType + '::RangeType, typename ' + rangeSpaceType + '::JacobianRangeType >'
    integrandsType = 'Dune::Fem::VirtualizedIntegrands< typename ' + rangeSpaceType + '::GridPartType, ' + domainValueType + ', ' + rangeValueType + ' >'

    typeName = 'Dune::Fem::GalerkinOperator< ' + integrandsType + ', ' + domainFunctionType + ', ' + rangeFunctionType + ' >'

    ctors = []
    ctors.append(['[] ( ' + typeName + ' &self, pybind11::object gridView, ' + integrandsType + ' &integrands ) {',
                  '    new (&self) ' + typeName + '( Dune::FemPy::gridPart< typename ' + rangeSpaceType + '::GridPartType::GridViewType >( gridView ), std::ref( integrands ) );',
                  '  }, "grid"_a, "integrands"_a, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >()'])

    return load(includes, typeName, ctors).Operator(rangeSpace.grid, integrands)
