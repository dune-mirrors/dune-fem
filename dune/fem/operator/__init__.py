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
    module = generator.load(includes, typeName, moduleName, *args, dynamicAttr=True)
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

    storage, domainFunctionIncludes, domainFunctionType, _, _, _ = domainSpace.storage
    _, rangeFunctionIncludes, rangeFunctionType, _, _, _ = rangeSpace.storage

    includes = ["dune/fem/schemes/galerkin.hh", "dune/fempy/py/grid/gridpart.hh"]
    includes += domainSpace._includes + domainFunctionIncludes
    includes += rangeSpace._includes + rangeFunctionIncludes

    integrandsType = 'Dune::Fem::VirtualizedIntegrands< typename ' + rangeSpaceType + '::GridPartType, ' + integrands._domainValueType + ', ' + integrands._rangeValueType + ' >'

    import dune.create as create
    linearOperator = create.discretefunction(storage)(domainSpace,rangeSpace)[3]

    # typeName = 'Dune::Fem::GalerkinOperator< ' + integrandsType + ', ' + domainFunctionType + ', ' + rangeFunctionType + ' >'
    typeName = 'Dune::Fem::DifferentiableGalerkinOperator< ' + integrandsType + ', ' + linearOperator + ' >'

    constructor = Constructor(['const '+domainSpaceType+'& dSpace','const '+rangeSpaceType+' &rSpace', integrandsType + ' &integrands'],
                              ['return new ' + typeName + '( dSpace, rSpace, integrands );'],
                              ['pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])

    op = load(includes, typeName, constructor).Operator(domainSpace,rangeSpace, integrands)
    op.model = integrands
    return op

def h1(model, domainSpace, rangeSpace=None):
    if rangeSpace is None:
        rangeSpace = domainSpace

    modelParam = None
    if isinstance(model, (list, tuple)):
        modelParam = model[1:]
        model = model[0]
    if isinstance(model,Form):
        model = model == 0
    if isinstance(model,Equation):
        from dune.fem.model._models import elliptic as makeElliptic
        if modelParam:
            model = makeElliptic(domainSpace.grid,model,*modelParam)
        else:
            model = makeElliptic(domainSpace.grid,model)

    domainSpaceType = domainSpace._typeName
    rangeSpaceType = rangeSpace._typeName

    storage, domainFunctionIncludes, domainFunctionType, _, _, _ = domainSpace.storage
    _, rangeFunctionIncludes, rangeFunctionType, _, _, _ = rangeSpace.storage

    includes = ["dune/fem/schemes/elliptic.hh", "dune/fempy/py/grid/gridpart.hh"]
    includes += domainSpace._includes + domainFunctionIncludes
    includes += rangeSpace._includes + rangeFunctionIncludes
    includes += ["dune/fem/schemes/diffusionmodel.hh", "dune/fempy/parameter.hh"]

    import dune.create as create
    linearOperator = create.discretefunction(storage)(domainSpace,rangeSpace)[3]

    modelType = "DiffusionModel< " +\
          "typename " + domainSpaceType + "::GridPartType, " +\
          domainSpaceType + "::dimRange, " +\
          rangeSpaceType + "::dimRange, " +\
          "typename " + domainSpaceType + "::RangeFieldType >"
    typeName = "DifferentiableEllipticOperator< " + linearOperator + ", " + modelType + ">"

    constructor = Constructor(['const '+domainSpaceType+'& dSpace','const '+rangeSpaceType+' &rSpace', modelType + ' &model'],
                              ['return new ' + typeName + '( dSpace, rSpace, model );'],
                              ['pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])

    scheme = load(includes, typeName, constructor).Operator(domainSpace,rangeSpace, model)
    scheme.model = model
    return scheme
