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

def load(includes, typeName, *args, baseClasses=None, preamble=None,
         codegen=None):
    moduleName = hashlib.md5(typeName.encode('utf-8')).hexdigest()
    if baseClasses is None:
        baseClasses = []
    if codegen is not None:
        if codegen[0].codegenStorage:
            includesExt, moduleNameExt = codegen[0].codegen(
                "op"+ "_" + moduleName,
                interiorQuadratureOrders=codegen[1],
                skeletonQuadratureOrders=codegen[2] )
            includes = includesExt + includes
            moduleName = moduleNameExt + "_" + moduleName
    includes = includes + ["dune/fempy/py/operator.hh"]
    moduleName = "femoperator" + "_" + moduleName
    module = generator.load(includes, typeName, moduleName, *args,
                  preamble=preamble, dynamicAttr=True, baseClasses=baseClasses)
    module.Operator.codegen = codegen
    return module

linearGenerator = SimpleGenerator("LinearOperator", "Dune::FemPy")

def loadLinear(includes, typeName, *args, backend=None, preamble=None):
    from dune.fem.space import addBackend
    includes = includes + ["dune/fempy/py/operator.hh"]
    moduleName = "femoperator" + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = linearGenerator.load(includes, typeName, moduleName, *args, preamble=preamble, dynamicAttr=True)
    LinearOperator = module.LinearOperator
    try:
        backend = backend[0] if backend[0] == backend[1] else None
    except:
        pass
    if hasattr(LinearOperator,"_backend") and backend is not None:
        addBackend(LinearOperator,backend)
    return module


def _galerkin(integrands, domainSpace=None, rangeSpace=None,
              virtualize=None, communicate=True, operatorPrefix = '' ):
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
        if rangeSpace == None:
            try:
                rangeSpace = integrands.lhs.arguments()[0].ufl_function_space()
            except AttributeError:
                raise ValueError("no range space provided and could not deduce from form provided")
        if domainSpace == None:
            try:
                domainSpace = integrands.lhs.arguments()[1].ufl_function_space()
            except AttributeError:
                raise ValueError("no domain space provided and could not deduce from form provided")
        if modelParam:
            integrands = makeIntegrands(domainSpace.grid,integrands,*modelParam)
        else:
            integrands = makeIntegrands(domainSpace.grid,integrands)

    if virtualize is None:
        virtualize = integrands.virtualized

    if not hasattr(rangeSpace,"interpolate"):
        raise ValueError("given range space has to be a discrete space")
    if not hasattr(domainSpace,"interpolate"):
        raise ValueError("given domain space has to be a discrete space")

    domainSpaceType = domainSpace._typeName
    rangeSpaceType  = rangeSpace._typeName

    storage, domainFunctionIncludes, domainFunctionType, _, _, dbackend = domainSpace.storage
    rstorage, rangeFunctionIncludes,  rangeFunctionType,  _, _, rbackend = rangeSpace.storage

    includes = ["dune/fem/schemes/galerkin.hh", "dune/fempy/py/grid/gridpart.hh"]
    if operatorPrefix == 'MOL':
        includes += ["dune/fem/schemes/molgalerkin.hh"]

    includes += domainSpace._includes + domainFunctionIncludes
    includes += rangeSpace._includes + rangeFunctionIncludes

    if virtualize:
        integrandsType = 'Dune::Fem::VirtualizedIntegrands< typename ' + rangeSpaceType + '::GridPartType, ' + integrands._domainValueType + ", " + integrands._rangeValueType+ ' >'
    else:
        includes += integrands._includes
        integrandsType = integrands._typeName

    if not rstorage == storage:
        typeName = 'Dune::Fem::' + operatorPrefix + 'GalerkinOperator< ' + integrandsType + ', ' + domainFunctionType + ', ' + rangeFunctionType + ' >'
        constructor = Constructor(['pybind11::object gridView', integrandsType + ' &integrands'],
                                  ['return new DuneType( Dune::FemPy::gridPart< typename ' + rangeSpaceType + '::GridPartType::GridViewType >( gridView ), integrands );'],
                                  ['"grid"_a', '"integrands"_a', 'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()'])
        constructor = Constructor(['const '+domainSpaceType+'& dSpace','const '+rangeSpaceType+' &rSpace', integrandsType + ' &integrands'],
                                  ['return new DuneType( dSpace.gridPart(), integrands );'],
                                  ['pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()'])
    else:
        import dune.create as create
        linearOperator = create.discretefunction(storage)(domainSpace,rangeSpace)[3]
        typeName = 'Dune::Fem::' + operatorPrefix + 'DifferentiableGalerkinOperator< ' + integrandsType + ', ' + linearOperator + ' >'
        constructor = Constructor(['const '+domainSpaceType+'& dSpace','const '+rangeSpaceType+' &rSpace', integrandsType + ' &integrands'],
                                  ['return new DuneType( dSpace, rSpace, integrands );'],
                                  ['pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])
    if integrands.hasDirichletBoundary:
        typeName = 'DirichletWrapperOperator< ' + typeName + ' >'

    setCommunicate = Method('setCommunicate', '''[]( DuneType &self, const bool communicate ) { self.setCommunicate( communicate ); }''' )
    op = load(includes, typeName, setCommunicate, constructor).Operator(domainSpace,rangeSpace,integrands)
    op.model = integrands
    # apply communicate flag
    op.setCommunicate( communicate )
    return op

# galerkin operator
def galerkin(integrands, domainSpace=None, rangeSpace=None,
             virtualize=None, communicate=True):
    return _galerkin(integrands, domainSpace=domainSpace, rangeSpace=rangeSpace,
                     virtualize=virtualize, communicate=communicate)

# method of lines galerkin operator (applies inverse mass matrix)
def molGalerkin(integrands, domainSpace=None, rangeSpace=None,
                virtualize=None, communicate=True):
    return _galerkin(integrands, domainSpace=domainSpace, rangeSpace=rangeSpace,
                     virtualize=virtualize, communicate=communicate,
                     operatorPrefix='MOL')

def h1(model, domainSpace=None, rangeSpace=None):
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
        if rangeSpace == None:
            try:
                rangeSpace = model.lhs.arguments()[0].ufl_function_space()
            except AttributeError:
                raise ValueError("no range space provided and could not deduce from form provided")
        if domainSpace == None:
            try:
                domainSpace = model.lhs.arguments()[1].ufl_function_space()
            except AttributeError:
                raise ValueError("no domain space provided and could not deduce from form provided")
        if modelParam:
            model = makeElliptic(domainSpace.grid,model,*modelParam)
        else:
            model = makeElliptic(domainSpace.grid,model)

    if not hasattr(rangeSpace,"interpolate"):
        raise ValueError("wrong range space")
    if not hasattr(domainSpace,"interpolate"):
        raise ValueError("wrong domain space")

    domainSpaceType = domainSpace._typeName
    rangeSpaceType = rangeSpace._typeName

    storage,  domainFunctionIncludes, domainFunctionType, _, _, dbackend = domainSpace.storage
    rstorage, rangeFunctionIncludes,  rangeFunctionType,  _, _, rbackend = rangeSpace.storage
    if not rstorage == storage:
        raise ValueError("storage for both spaces must be identical to construct operator")

    includes = ["dune/fem/schemes/elliptic.hh", "dune/fem/schemes/dirichletwrapper.hh", "dune/fempy/py/grid/gridpart.hh"]
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
    if model.hasDirichletBoundary:
        typeName = 'DirichletWrapperOperator< ' + typeName + ' >'

    constructor = Constructor(['const '+domainSpaceType+'& dSpace','const '+rangeSpaceType+' &rSpace', modelType + ' &model'],
                              ['return new ' + typeName + '( dSpace, rSpace, model );'],
                              ['pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])

    scheme = load(includes, typeName, constructor).Operator(domainSpace,rangeSpace, model)
    scheme.model = model
    return scheme

# TODO: linear could also return the 'affine shift' and the point of
#       linearization (0 now) could also be passed in
#       The returning operator could have a 'update' method and then we
#       remove the 'jacobian' from the official interface
def linear(operator, ubar=None,parameters={}):
    from dune.fem.function import discreteFunction
    assert hasattr(operator,"jacobian"), "operator does not allow assembly"
    rangeSpace  = operator.rangeSpace
    domainSpace = operator.domainSpace

    domainSpaceType = domainSpace._typeName
    rangeSpaceType = rangeSpace._typeName

    storage,  domainFunctionIncludes, domainFunctionType, _, _, dbackend = domainSpace.storage
    rstorage, rangeFunctionIncludes,  rangeFunctionType,  _, _, rbackend = rangeSpace.storage
    if not rstorage == storage:
        raise ValueError("storage for both spaces must be identical to construct operator")

    includes = ["dune/fempy/py/grid/gridpart.hh"]
    includes += domainSpace._includes + domainFunctionIncludes
    includes += rangeSpace._includes + rangeFunctionIncludes

    import dune.create as create
    typeName = create.discretefunction(storage)(domainSpace,rangeSpace)[3]

    constructor = Constructor(['const '+domainSpaceType+'& dSpace','const '+rangeSpaceType+' &rSpace',
                               'const pybind11::dict &parameters'],
                              ['return new ' + typeName + '( "tmp", dSpace, rSpace, '+
                                       'Dune::FemPy::pyParameter( "fem.solver.", parameters, std::make_shared< std::string >() ));'],
                              ['pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])

    lin = loadLinear(includes, typeName, constructor, backend=(dbackend,rbackend)).LinearOperator(domainSpace,rangeSpace,parameters)
    if ubar is None:
        # operator.jacobian(domainSpace.interpolate([0,]*domainSpace.dimRange,"tmp"), lin)
        operator.jacobian(discreteFunction(domainSpace,"tmp"), lin)
    else:
        operator.jacobian(domainSpace.interpolate(ubar,"tmp"), lin)
    return lin
