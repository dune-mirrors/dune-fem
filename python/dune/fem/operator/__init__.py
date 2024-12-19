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
from dune.common.utility import isString
from dune.fem.deprecated import deprecated
from dune.fem.discretefunction import _storage

_defaultGenerator = SimpleGenerator("Operator", "Dune::FemPy")

def load(includes, typeName, *args, baseClasses=None, preamble=None,
         codegen=None, generator = _defaultGenerator ):
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

_linearGenerator = SimpleGenerator("LinearOperator", "Dune::FemPy")

def loadLinear(includes, typeName, *args, backend=None, preamble=None,
               generator=_linearGenerator):
    from dune.fem.space import addBackend
    includes = includes + ["dune/fempy/py/operator.hh"]
    moduleName = "femoperator" + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, *args, preamble=preamble, dynamicAttr=True)
    LinearOperator = module.LinearOperator
    try:
        backend = backend[0] if backend[0] == backend[1] else None
    except:
        pass
    if hasattr(LinearOperator,"_backend") and backend is not None:
        # the as_numpy backend needs scipy.sparse.csr_matrix
        # in the C++ code this will simply fail if scipy is not available
        if backend == 'as_numpy':
            from scipy.sparse import csr_matrix
        addBackend(LinearOperator,backend)
    return module

def _opLinear(self, assemble=True, parameters=None):
    A = _linear([self.domainSpace,self.rangeSpace], parameters)
    if assemble:
        self.jacobian(self.domainSpace.zero, A)
    return A
def _opDirichletIndices(self, id=None):
    if id is None:
        return [i*len(block)+j for i,block in enumerate(self.dirichletBlocks)
                               for j,b in enumerate(block) if b>0]
    else:
        return [i*len(block)+j for i,block in enumerate(self.dirichletBlocks)
                               for j,b in enumerate(block) if b==id]

def _galerkin(integrands, domainSpace=None, rangeSpace=None,
              virtualize=None, communicate=True,
              jacobianStorage=None,
              operatorPrefix = '' ):
    """
    Parameter:
        integrands       Model representing the operator B(u,v): V -> W, i.e. weak form of a PDE
        domainSpace      V, discrete function space representing the domain
        rangeSpace       W, discrete function space representing the range
        virtualize       virtualize the integrands model passed
        communicate      If true, a synchronization will be done during the
                         application of the operator
        jacobianStorage  Storage for the Jacobian linear operator (default: same as spaces)

    Return:
        A Galerkin operator implementing the integration of the weak form of the
        given PDE.
    """
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
            integrands = makeIntegrands(domainSpace.gridView,integrands,*modelParam)
        else:
            integrands = makeIntegrands(domainSpace.gridView,integrands)

    if virtualize is None:
        virtualize = integrands.virtualized

    if not hasattr(rangeSpace,"interpolate"):
        raise ValueError("given range space has to be a discrete space")
    if not hasattr(domainSpace,"interpolate"):
        raise ValueError("given domain space has to be a discrete space")

    domainSpaceType = domainSpace.cppTypeName
    rangeSpaceType  = rangeSpace.cppTypeName

    storage                = domainSpace.storage.name
    domainFunctionIncludes = domainSpace.storage.includes
    domainFunctionType     = domainSpace.storage.type
    dbackend               = domainSpace.storage.backend

    rstorage               = rangeSpace.storage.name
    rangeFunctionIncludes  = rangeSpace.storage.includes
    rangeFunctionType      = rangeSpace.storage.type
    rbackend               = rangeSpace.storage.backend

    # use storage of discrete function if not specified
    if jacobianStorage is None:
        jacobianStorage = storage
    else:
        assert isString(jacobianStorage)

    includes = ["dune/fem/schemes/galerkin.hh", "dune/fempy/py/grid/gridpart.hh"]
    if operatorPrefix == 'MOL':
        includes += ["dune/fem/schemes/molgalerkin.hh"]

    includes += domainSpace.cppIncludes + domainFunctionIncludes
    includes += rangeSpace.cppIncludes + rangeFunctionIncludes

    if virtualize:
        integrandsType = 'Dune::Fem::VirtualizedIntegrands< typename ' + rangeSpaceType + '::GridPartType, ' + integrands._domainValueType + ", " + integrands._rangeValueType+ ' >'
    else:
        includes += integrands.cppIncludes
        integrandsType = integrands.cppTypeName

    # get storage depending on choices
    solverstorage  = _storage( dfStorage=storage, solverStorage=jacobianStorage)(domainSpace,rangeSpace)
    linearOperator = solverstorage.linopType

    # add extra includes for linear operator
    includes += solverstorage.linopIncludes

    if not rstorage == storage:
        typeName = 'Dune::Fem::' + operatorPrefix + 'GalerkinOperator< ' + integrandsType + ', ' + domainFunctionType + ', ' + rangeFunctionType + ' >'
        constructor = Constructor(['pybind11::object gridView', integrandsType + ' &integrands'],
                                  ['return new DuneType( Dune::FemPy::gridPart< typename ' + rangeSpaceType + '::GridPartType::GridViewType >( gridView ), integrands );'],
                                  ['"grid"_a', '"integrands"_a', 'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()'])
        constructor = Constructor(['const '+domainSpaceType+'& dSpace','const '+rangeSpaceType+' &rSpace', integrandsType + ' &integrands'],
                                  ['return new DuneType( dSpace.gridPart(), integrands );'],
                                  ['pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()'])
    else:
        # get storage depending on choices
        typeName = 'Dune::Fem::' + operatorPrefix + 'DifferentiableGalerkinOperator< ' + integrandsType + ', ' + linearOperator + ' >'
        constructor = Constructor(['const '+domainSpaceType+'& dSpace','const '+rangeSpaceType+' &rSpace', integrandsType + ' &integrands'],
                                  ['return new DuneType( dSpace, rSpace, integrands );'],
                                  ['pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])
    if integrands.hasDirichletBoundary:
        typeName = 'DirichletWrapperOperator< ' + typeName + ' >'

    setCommunicate = Method('setCommunicate', '''[]( DuneType &self, const bool communicate ) { self.setCommunicate( communicate ); }''' )
    gridSizeInterior = Method('gridSizeInterior', '''[]( DuneType &self ) { return self.gridSizeInterior(); }''' )
    op = load(includes, typeName, setCommunicate, gridSizeInterior, constructor).Operator(domainSpace,rangeSpace,integrands)
    op.model = integrands
    op.__class__.linear = _opLinear
    op.__class__.dirichletIndices = _opDirichletIndices
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
            model = makeElliptic(domainSpace.gridView,model,*modelParam)
        else:
            model = makeElliptic(domainSpace.gridView,model)

    if not hasattr(rangeSpace,"interpolate"):
        raise ValueError("wrong range space")
    if not hasattr(domainSpace,"interpolate"):
        raise ValueError("wrong domain space")

    domainSpaceType = domainSpace.cppTypeName
    rangeSpaceType = rangeSpace.cppTypeName

    storage                = domainSpace.storage.name
    domainFunctionIncludes = domainSpace.storage.includes
    domainFunctionType     = domainSpace.storage.type
    dbackend               = domainSpace.storage.backend

    rstorage               = rangeSpace.storage.name
    rangeFunctionIncludes  = rangeSpace.storage.includes
    rangeFunctionType      = rangeSpace.storage.type
    rbackend               = rangeSpace.storage.backend

    if not rstorage == storage:
        raise ValueError("storage for both spaces must be identical to construct operator")

    includes = ["dune/fem/schemes/elliptic.hh", "dune/fem/schemes/dirichletwrapper.hh", "dune/fempy/py/grid/gridpart.hh"]
    includes += domainSpace.cppIncludes + domainFunctionIncludes
    includes += rangeSpace.cppIncludes + rangeFunctionIncludes
    includes += ["dune/fem/schemes/conservationlawmodel.hh", "dune/fempy/parameter.hh"]

    import dune.create as create
    linearOperator = create.discretefunction(storage)(domainSpace,rangeSpace).linopType

    modelType = "ConservationLawModel< " +\
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
    op.__class__.linear = lambda self, parameters=None: (
        _linear([self.domainSpace,self.rangeSpace], parameters)
      )
    return scheme

def _linear(operator, ubar=None,parameters=None,affineShift=False):
    """ operator can be either a dune.fem.operator/scheme or a tuple/list of spaces.
        In the second case ubar needs to be None and the operator is unassembled.
        In the first case the operator will be assembled around 'ubar' or around zero if ubar is None.

        If affineShift is true the function returns both the linearization
        A and the affine shift, i.e., both L[0] and DL[0]. Note that the
        `right hand side` returned by `affine` is `b=-L[0]`.
    """
    if parameters is None:
        parameters = {}
    try:
        domainSpace = operator[0]
        rangeSpace  = operator[1]
        assemble = False
    except TypeError:
        assert hasattr(operator,"jacobian"), "operator does not allow assembly"
        rangeSpace  = operator.rangeSpace
        domainSpace = operator.domainSpace
        assemble = True

    domainSpaceType = domainSpace.cppTypeName
    rangeSpaceType = rangeSpace.cppTypeName

    storage                = domainSpace.storage.name
    domainFunctionIncludes = domainSpace.storage.includes
    domainFunctionType     = domainSpace.storage.type
    dbackend               = domainSpace.storage.backend

    rstorage               = rangeSpace.storage.name
    rangeFunctionIncludes  = rangeSpace.storage.includes
    rangeFunctionType      = rangeSpace.storage.type
    rbackend               = rangeSpace.storage.backend

    if not rstorage == storage:
        raise ValueError("storage for both spaces must be identical to construct operator")

    includes = ["dune/fempy/py/grid/gridpart.hh"]
    includes += domainSpace.cppIncludes + domainFunctionIncludes
    includes += rangeSpace.cppIncludes + rangeFunctionIncludes
    includes += domainSpace.storage.linopIncludes

    import dune.create as create
    typeName = create.discretefunction(storage)(domainSpace,rangeSpace).linopType

    constructor = Constructor(['const '+domainSpaceType+'& dSpace','const '+rangeSpaceType+' &rSpace',
                               'const pybind11::dict &parameters'],
                              ['return new ' + typeName + '( "tmp", dSpace, rSpace, '+
                                       'Dune::FemPy::pyParameter( "fem.solver.", parameters, std::make_shared< std::string >() ));'],
                              ['pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])

    lin = loadLinear(includes, typeName, constructor, backend=(dbackend,rbackend)).LinearOperator(domainSpace,rangeSpace,parameters)

    if assemble:
        try:
            lin.dirichletBlocks = operator.dirichletBlocks
        except:
            pass
        if ubar is None:
            operator.jacobian(domainSpace.zero, lin)
        else:
            operator.jacobian(ubar, lin)
    if affineShift:
        b = rangeSpace.zero.copy()
        operator(rangeSpace.zero,b)
        return lin,b
    else:
        return lin
def linear(operator, ubar=None,parameters=None):
    deprecated("dune.fem.operator is deprecated use the ``linear`` method on the scheme/operator  instead.")
    return _linear(operator,ubar,parameters)
