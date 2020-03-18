from __future__ import division, print_function, unicode_literals

from ufl import FiniteElementBase, FunctionSpace, dx
from ufl import Coefficient, FacetNormal, Form, SpatialCoordinate
from ufl import CellVolume, MinCellEdgeLength, MaxCellEdgeLength
from ufl import FacetArea, MinFacetEdgeLength, MaxFacetEdgeLength
from ufl import action, derivative, as_vector
from ufl.classes import Indexed
from ufl.core.multiindex import FixedIndex, MultiIndex
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.algorithms.replace import Replacer
from ufl.constantvalue import IntValue, Zero
from ufl.corealg.map_dag import map_expr_dags
from ufl.differentiation import Grad
from ufl.equation import Equation
from ufl import UFLException
from dune.ufl.tensors import ExprTensor
from dune.source.cplusplus import UnformattedExpression, SwitchStatement, Declaration, UnformattedBlock, assign, Block, ConstructExpression

from dune.source.builtin import make_pair
from dune.source.cplusplus import InitializerList, Variable
from dune.source.cplusplus import construct, lambda_, makeExpression, maxEdgeLength, minEdgeLength, return_
from dune.source.cplusplus import SourceWriter
from dune.source.algorithm.extractvariables import extractVariablesFromExpressions, extractVariablesFromStatements

from dune.common.hashit import hashIt

from dune.ufl import codegen, DirichletBC
from dune.ufl.gatherderivatives import gatherDerivatives
from dune.ufl.linear import splitForm
import dune.ufl.tensors as tensors

from .model import Integrands

def generateCode(predefined, testFunctions, tensorMap, tempVars=True):
    # build list of all expressions to compile
    expressions = []
    for phi in testFunctions:
        if (phi,) not in tensorMap:
            continue
        tensor = tensorMap[(phi,)]
        keys = tensor.keys()
        expressions += [tensor[i] for i in keys]

    # compile all expressions at once
    preamble, results = codegen.generateCode(predefined, expressions, tempVars=tempVars)

    # extract generated code for expressions and build values
    values = []
    for phi in testFunctions:
        value = tensors.fill(phi.ufl_shape, makeExpression(0))
        if (phi,) in tensorMap:
            tensor = tensorMap[(phi,)]
            keys = tensor.keys()
            for i, r in zip(keys, results[:len(keys)]):
                value = tensors.setItem(value, i, r)
            results = results[len(keys):]
        values += [tensors.reformat(lambda row: InitializerList(*row), phi.ufl_shape, value)]

    return preamble, values

# used for Dirichlet conditions - taken from elliptic
def generateDirichletCode(predefined, tensor, tempVars=True):
    keys = tensor.keys()
    expressions = [tensor[i] for i in keys]
    preamble, results = codegen.generateCode(predefined, expressions, tempVars=tempVars)
    result = Variable('auto', 'result')
    return preamble + [assign(result[i], r) for i, r in zip(keys, results)]

def generateDirichletDomainCode(predefined, tensor, tempVars=True):
    # predefined={}
    keys = tensor.keys()
    expressions = [tensor[i] for i in keys]
    preamble, results = codegen.generateCode(predefined, expressions, tempVars=tempVars)
    result = Variable('int', 'domainId')
    return [preamble,assign(result, results[0])]

def generateLinearizedCode(predefined, testFunctions, trialFunctionMap, tensorMap, tempVars=True):
    """generate code for a bilinear form

    Args:
        predefined:       list of predefined arguments or coefficients
        testFunctions:    list of arguments to interpret as test functions
        trialFunctionMap: map of variable to list of arguments to interpret as trial functions
        tensorMap:        map of expression tensors of shape (testFunction x trialFunction)
        tempVars:         introduce temporary variables during code generation
    """

    # build list of all expressions to compile
    expressions = []
    for var, trialFunctions in trialFunctionMap.items():
        for phi in testFunctions:
            for psi in trialFunctions:
                if (phi, psi) not in tensorMap:
                    continue
                tensor = tensorMap[(phi, psi)]
                keys = tensor.keys()
                expressions += [tensor[i] for i in keys]

    # compile all expressions at once
    preamble, results = codegen.generateCode(predefined, expressions, tempVars=tempVars)

    # extract generated code for expressions and build values
    values = {}
    for var, trialFunctions in trialFunctionMap.items():
        values[var] = []
        for phi in testFunctions:
            value = tensors.fill(phi.ufl_shape, None)
            for idx in range(len(trialFunctions)):
                psi = trialFunctions[idx]
                if (phi, psi) in tensorMap:
                    tensor = tensorMap[(phi, psi)]
                    keys = tensor.keys()
                    for ij, r in zip(keys, results[:len(keys)]):
                        if isinstance(tensor[ij], Zero):
                            continue
                        i = ij[:len(phi.ufl_shape)]
                        j = ij[len(phi.ufl_shape):]
                        if isinstance(tensor[ij], IntValue) and int(tensor[ij]) == 1:
                            r = var[idx][j]
                        else:
                            r = r * var[idx][j]
                        s = tensors.getItem(value, i)
                        s = r if s is None else s + r
                        value = tensors.setItem(value, i, s)
                    results = results[len(keys):]
            value = tensors.apply(lambda v : makeExpression(0) if v is None else v, phi.ufl_shape, value)
            values[var] += [tensors.reformat(lambda row: InitializerList(*row), phi.ufl_shape, value)]

    return preamble, values


def generateUnaryCode(predefined, testFunctions, tensorMap, tempVars=True):
    preamble, values = generateCode(predefined, testFunctions, tensorMap, tempVars=tempVars)
    return preamble + [return_(construct('RangeValueType', *values, brace=True))]


def generateUnaryLinearizedCode(predefined, testFunctions, trialFunctions, tensorMap, tempVars=True):
    var = Variable('std::tuple< RangeType, JacobianRangeType >', 'phi')
    if tensorMap is None:
        values = []
        for phi in testFunctions:
            value = tensors.fill(phi.ufl_shape, None)
            value = tensors.apply(lambda v : makeExpression(0), phi.ufl_shape, value)
            values += [tensors.reformat(lambda row: InitializerList(*row), phi.ufl_shape, value)]
        return [return_(lambda_(args=['const DomainValueType &phi'],\
                  code=return_(construct('RangeValueType',*values,
                    brace=True)) ))]

    preamble, values = generateLinearizedCode(predefined, testFunctions, {var: trialFunctions}, tensorMap, tempVars=tempVars)
    capture = extractVariablesFromExpressions(values[var]) - {var}
    return preamble + [return_(lambda_(capture=capture, args=['const DomainValueType &phi'], code=return_(construct('RangeValueType', *values[var], brace=True))))]


def generateBinaryCode(predefined, testFunctions, tensorMap, tempVars=True):
    restrictedTestFunctions = [phi('+') for phi in testFunctions] + [phi('-') for phi in testFunctions]
    preamble, values = generateCode(predefined, restrictedTestFunctions, tensorMap, tempVars=tempVars)
    return preamble + [return_(make_pair(construct('RangeValueType', *values[:len(testFunctions)], brace=True), construct('RangeValueType', *values[len(testFunctions):], brace=True)))]


def generateBinaryLinearizedCode(predefined, testFunctions, trialFunctions, tensorMap, tempVars=True):
    restrictedTestFunctions = [phi('+') for phi in testFunctions] + [phi('-') for phi in testFunctions]

    trialFunctionsIn = [psi('+') for psi in trialFunctions]
    trialFunctionsOut = [psi('-') for psi in trialFunctions]

    if tensorMap is None:
        value = construct('RangeValueType', *[0 for i in range(len(testFunctions))], brace=True)
        tensorIn = lambda_(args=['const DomainValueType &phiIn'], code=return_(make_pair(value, value)))
        tensorOut = lambda_(args=['const DomainValueType &phiOut'], code=return_(make_pair(value, value)))
        return [return_(make_pair(tensorIn, tensorOut))]

    varIn = Variable('std::tuple< RangeType, JacobianRangeType >', 'phiIn')
    varOut = Variable('std::tuple< RangeType, JacobianRangeType >', 'phiOut')
    preamble, values = generateLinearizedCode(predefined, restrictedTestFunctions, {varIn: trialFunctionsIn, varOut: trialFunctionsOut}, tensorMap, tempVars=tempVars)

    captureIn = extractVariablesFromExpressions(values[varIn]) - {varIn}
    captureOut = extractVariablesFromExpressions(values[varOut]) - {varOut}

    tensorIn = lambda_(capture=captureIn, args=['const DomainValueType &phiIn'], code=return_(make_pair(construct('RangeValueType', *values[varIn][:len(testFunctions)], brace=True), construct('RangeValueType', *values[varIn][len(testFunctions):], brace=True))))
    tensorOut = lambda_(capture=captureOut, args=['const DomainValueType &phiOut'], code=return_(make_pair(construct('RangeValueType', *values[varOut][:len(testFunctions)], brace=True), construct('RangeValueType', *values[varOut][len(testFunctions):], brace=True))))

    return preamble + [return_(make_pair(tensorIn, tensorOut))]


def toFileName(value):
    import unicodedata
    value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore')
    value = unicode(re.sub('[^\w\s-]', '', value).strip().lower())
    value = unicode(re.sub('[-\s]+', '-', value))
    return value


def _compileUFL(integrands, form, *args, tempVars=True):
    if isinstance(form, Equation):
        form = form.lhs - form.rhs
    if not isinstance(form, Form):
        raise ValueError("ufl.Form or ufl.Equation expected.")

    # added for dirichlet treatment same as elliptic model
    dirichletBCs = [arg for arg in args if isinstance(arg, DirichletBC)]

    uflExpr = [form] + [bc.ufl_value for bc in dirichletBCs]
    if len(form.arguments()) < 2:
        raise ValueError("Integrands model requires form with at least two arguments.")

    x = SpatialCoordinate(form.ufl_cell())
    n = FacetNormal(form.ufl_cell())

    cellVolume = CellVolume(form.ufl_cell())
    maxCellEdgeLength = MaxCellEdgeLength(form.ufl_cell())
    minCellEdgeLength = MinCellEdgeLength(form.ufl_cell())

    facetArea = FacetArea(form.ufl_cell())
    maxFacetEdgeLength = MaxFacetEdgeLength(form.ufl_cell())
    minFacetEdgeLength = MinFacetEdgeLength(form.ufl_cell())

    phi, u = form.arguments()
    ubar = Coefficient(u.ufl_function_space())

    derivatives = gatherDerivatives(form, [phi, u])

    derivatives_phi = derivatives[0]
    derivatives_u = derivatives[1]
    derivatives_ubar = map_expr_dags(Replacer({u: ubar}), derivatives_u)

    try:
        integrands.field = u.ufl_function_space().field
    except AttributeError:
        pass

    integrals = splitForm(form, [phi])

    dform = apply_derivatives(derivative(action(form, ubar), ubar, u))
    linearizedIntegrals = splitForm(dform, [phi, u])

    if not set(integrals.keys()) <= {'cell', 'exterior_facet', 'interior_facet'}:
        raise Exception('unknown integral encountered in ' + str(set(integrals.keys())) + '.')

    if 'cell' in integrals.keys():
        arg = Variable(integrands.domainValueTuple, 'u')

        predefined = {derivatives_u[i]: arg[i] for i in range(len(derivatives_u))}
        predefined[x] = integrands.spatialCoordinate('x')
        predefined[cellVolume] = integrands.cellVolume()
        predefined[maxCellEdgeLength] = maxEdgeLength(integrands.cellGeometry())
        predefined[minCellEdgeLength] = minEdgeLength(integrands.cellGeometry())
        integrands.predefineCoefficients(predefined, False)
        integrands.interior = generateUnaryCode(predefined, derivatives_phi, integrals['cell'], tempVars=tempVars)

        predefined = {derivatives_ubar[i]: arg[i] for i in range(len(derivatives_u))}
        predefined[x] = integrands.spatialCoordinate('x')
        predefined[cellVolume] = integrands.cellVolume()
        predefined[maxCellEdgeLength] = maxEdgeLength(integrands.cellGeometry())
        predefined[minCellEdgeLength] = minEdgeLength(integrands.cellGeometry())
        integrands.predefineCoefficients(predefined, False)
        integrands.linearizedInterior = generateUnaryLinearizedCode(predefined, derivatives_phi, derivatives_u, linearizedIntegrals.get('cell'), tempVars=tempVars)

    if 'exterior_facet' in integrals.keys():
        arg = Variable(integrands.domainValueTuple, 'u')

        predefined = {derivatives_u[i]: arg[i] for i in range(len(derivatives_u))}
        predefined[x] = integrands.spatialCoordinate('x')
        predefined[n] = integrands.facetNormal('x')
        predefined[cellVolume] = integrands.cellVolume()
        predefined[maxCellEdgeLength] = maxEdgeLength(integrands.cellGeometry())
        predefined[minCellEdgeLength] = minEdgeLength(integrands.cellGeometry())
        predefined[facetArea] = integrands.facetArea()
        predefined[maxFacetEdgeLength] = maxEdgeLength(integrands.facetGeometry())
        predefined[minFacetEdgeLength] = minEdgeLength(integrands.facetGeometry())
        integrands.predefineCoefficients(predefined, False)
        integrands.boundary = generateUnaryCode(predefined, derivatives_phi, integrals['exterior_facet'], tempVars=tempVars);

        predefined = {derivatives_ubar[i]: arg[i] for i in range(len(derivatives_u))}
        predefined[x] = integrands.spatialCoordinate('x')
        predefined[n] = integrands.facetNormal('x')
        predefined[cellVolume] = integrands.cellVolume()
        predefined[maxCellEdgeLength] = maxEdgeLength(integrands.cellGeometry())
        predefined[minCellEdgeLength] = minEdgeLength(integrands.cellGeometry())
        predefined[facetArea] = integrands.facetArea()
        predefined[maxFacetEdgeLength] = maxEdgeLength(integrands.facetGeometry())
        predefined[minFacetEdgeLength] = minEdgeLength(integrands.facetGeometry())
        integrands.predefineCoefficients(predefined, False)
        integrands.linearizedBoundary = generateUnaryLinearizedCode(predefined, derivatives_phi, derivatives_u, linearizedIntegrals.get('exterior_facet'), tempVars=tempVars)

    if 'interior_facet' in integrals.keys():
        argIn = Variable(integrands.domainValueTuple, 'uIn')
        argOut = Variable(integrands.domainValueTuple, 'uOut')

        predefined = {derivatives_u[i](s): arg[i] for i in range(len(derivatives_u)) for s, arg in (('+', argIn), ('-', argOut))}
        predefined[x] = integrands.spatialCoordinate('xIn')
        predefined[n('+')] = integrands.facetNormal('xIn')
        predefined[cellVolume('+')] = integrands.cellVolume('Side::in')
        predefined[cellVolume('-')] = integrands.cellVolume('Side::out')
        predefined[maxCellEdgeLength('+')] = maxEdgeLength(integrands.cellGeometry('Side::in'))
        predefined[maxCellEdgeLength('-')] = maxEdgeLength(integrands.cellGeometry('Side::out'))
        predefined[minCellEdgeLength('+')] = minEdgeLength(integrands.cellGeometry('Side::in'))
        predefined[minCellEdgeLength('-')] = minEdgeLength(integrands.cellGeometry('Side::out'))
        predefined[facetArea] = integrands.facetArea()
        predefined[maxFacetEdgeLength] = maxEdgeLength(integrands.facetGeometry())
        predefined[minFacetEdgeLength] = minEdgeLength(integrands.facetGeometry())
        integrands.predefineCoefficients(predefined, True)
        integrands.skeleton = generateBinaryCode(predefined, derivatives_phi, integrals['interior_facet'], tempVars=tempVars)

        predefined = {derivatives_ubar[i](s): arg[i] for i in range(len(derivatives_u)) for s, arg in (('+', argIn), ('-', argOut))}
        predefined[x] = integrands.spatialCoordinate('xIn')
        predefined[n('+')] = integrands.facetNormal('xIn')
        predefined[cellVolume('+')] = integrands.cellVolume('Side::in')
        predefined[cellVolume('-')] = integrands.cellVolume('Side::out')
        predefined[maxCellEdgeLength('+')] = maxEdgeLength(integrands.cellGeometry('Side::in'))
        predefined[maxCellEdgeLength('-')] = maxEdgeLength(integrands.cellGeometry('Side::out'))
        predefined[minCellEdgeLength('+')] = minEdgeLength(integrands.cellGeometry('Side::in'))
        predefined[minCellEdgeLength('-')] = minEdgeLength(integrands.cellGeometry('Side::out'))
        predefined[facetArea] = integrands.facetArea()
        predefined[maxFacetEdgeLength] = maxEdgeLength(integrands.facetGeometry())
        predefined[minFacetEdgeLength] = minEdgeLength(integrands.facetGeometry())
        integrands.predefineCoefficients(predefined, True)
        integrands.linearizedSkeleton = generateBinaryLinearizedCode(predefined, derivatives_phi, derivatives_u, linearizedIntegrals.get('interior_facet'), tempVars=tempVars)

    if dirichletBCs:
        integrands.hasDirichletBoundary = True

        predefined = {}
        # predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + integrands.arg_x.name + ' ) )')
        predefined[x] = UnformattedExpression('auto', 'intersection.geometry().center( )')
        integrands.predefineCoefficients(predefined, False)

        maxId = 0
        codeDomains = []
        bySubDomain = dict()
        neuman = []
        wholeDomain = None
        for bc in dirichletBCs:
            if bc.subDomain in bySubDomain:
                raise Exception('Multiply defined Dirichlet boundary for subdomain ' + str(bc.subDomain))
            if not isinstance(bc.functionSpace, (FunctionSpace, FiniteElementBase)):
                raise Exception('Function space must either be a ufl.FunctionSpace or a ufl.FiniteElement')
            if isinstance(bc.functionSpace, FunctionSpace) and (bc.functionSpace != u.ufl_function_space()):
                raise Exception('Space of trial function and dirichlet boundary function must be the same - note that boundary conditions on subspaces are not available, yet')
            if isinstance(bc.functionSpace, FiniteElementBase) and (bc.functionSpace != u.ufl_element()):
                raise Exception('Cannot handle boundary conditions on subspaces, yet')

            if isinstance(bc.value, list):
                neuman = [i for i, x in enumerate(bc.value) if x == None]
            else:
                neuman = []

            value = ExprTensor(u.ufl_shape)
            for key in value.keys():
                value[key] = Indexed(bc.ufl_value, MultiIndex(tuple(FixedIndex(k) for k in key)))
            if bc.subDomain is None:
                wholeDomain = value,neuman
            elif isinstance(bc.subDomain,int):
                bySubDomain[bc.subDomain] = value,neuman
                maxId = max(maxId, bc.subDomain)
            else:
                domain = ExprTensor(())
                for key in domain.keys():
                    domain[key] = Indexed(bc.subDomain, MultiIndex(tuple(FixedIndex(k) for k in key)))
                codeDomains.append( (value,neuman,domain) )
        defaultCode = []
        if len(codeDomains)>0:
            defaultCode.append(Declaration(Variable('int', 'domainId')))
        # defaultCode.append(Declaration(Variable('auto', 'x'),
        #     initializer=UnformattedExpression('auto','intersection.geometry().center()')))
        for i,v in enumerate(codeDomains):
            block = Block()
            block.append(generateDirichletDomainCode(predefined, v[2], tempVars=tempVars))
            block.append('if (domainId)')
            ifBlock = UnformattedBlock()
            ifBlock.append('std::fill( dirichletComponent.begin(), dirichletComponent.end(), ' + str(maxId+i+2) + ' );')
            if len(v[1])>0:
                [ifBlock.append('dirichletComponent[' + str(c) + '] = 0;') for c in v[1]]
            ifBlock.append('return true;')
            block.append(ifBlock)
            defaultCode.append(block)
        if wholeDomain is not None:
            block = UnformattedBlock()
            block.append('std::fill( dirichletComponent.begin(), dirichletComponent.end(), ' + str(maxId+1) + ' );')
            if len(wholeDomain[1])>0:
                [block.append('dirichletComponent[' + str(c) + '] = 0;') for c in wholeDomain[1]]
            block.append('return true;')
            defaultCode.append(block)
        defaultCode.append(return_(False))

        bndId = Variable('const int', 'bndId')
        getBndId = UnformattedExpression('int', 'BoundaryIdProviderType::boundaryId( ' + integrands.arg_i.name + ' )')
        # getBndId = UnformattedExpression('int', 'boundaryIdGetter_.boundaryId( ' + integrands.arg_i.name + ' )')
        switch = SwitchStatement(bndId, default=defaultCode)
        for i,v in bySubDomain.items():
            code = []
            if len(v[1])>0:
                [code.append('dirichletComponent[' + str(c) + '] = 0;') for c in v[1]]
            code.append(return_(True))
            switch.append(i, code)
        integrands.isDirichletIntersection = [Declaration(bndId, initializer=getBndId),
                                         UnformattedBlock('std::fill( dirichletComponent.begin(), dirichletComponent.end(), ' + bndId.name + ' );'),
                                         switch
                                        ]

        predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + integrands.arg_x.name + ' ) )')
        if wholeDomain is None:
            defaultCode = assign(integrands.arg_r, construct("RRangeType", 0))
        else:
            defaultCode = generateDirichletCode(predefined, wholeDomain[0], tempVars=tempVars)
        switch = SwitchStatement(integrands.arg_bndId, default=defaultCode)
        for i, v in bySubDomain.items():
            switch.append(i, generateDirichletCode(predefined, v[0], tempVars=tempVars))
        for i,v in enumerate(codeDomains):
            switch.append(i+maxId+2, generateDirichletCode(predefined, v[0], tempVars=tempVars))
        integrands.dirichlet = [switch]

    return integrands

def compileUFL(form, *args, tempVars=True, virtualize=True):
    if isinstance(form, Equation):
        form = form.lhs - form.rhs

    if not isinstance(form, Form):
        raise ValueError("ufl.Form or ufl.Equation expected.")

    # added for dirichlet treatment same as elliptic model
    dirichletBCs = [arg for arg in args if isinstance(arg, DirichletBC)]

    uflExpr = [form] + [bc.ufl_value for bc in dirichletBCs]
    if len(form.arguments()) < 2:
        raise ValueError("Integrands model requires form with at least two arguments.")

    phi, u = form.arguments()

    derivatives = gatherDerivatives(form, [phi, u])

    derivatives_phi = derivatives[0]
    derivatives_u = derivatives[1]

    integrands = Integrands((d.ufl_shape for d in derivatives_u), (d.ufl_shape for d in derivatives_phi),
                            uflExpr,virtualize)

    return _compileUFL(integrands,form,*args,tempVars)
