from __future__ import division, print_function, unicode_literals

from ufl import Coefficient, FacetNormal, Form, SpatialCoordinate
from ufl import CellVolume, MinCellEdgeLength, MaxCellEdgeLength
from ufl import FacetArea, MinFacetEdgeLength, MaxFacetEdgeLength
from ufl import action, derivative
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.algorithms.replace import Replacer
from ufl.constantvalue import IntValue, Zero
from ufl.corealg.map_dag import map_expr_dags
from ufl.differentiation import Grad
from ufl.equation import Equation

from dune.source.builtin import make_pair
from dune.source.cplusplus import InitializerList, Variable
from dune.source.cplusplus import construct, lambda_, makeExpression, maxEdgeLength, minEdgeLength, return_
from dune.source.algorithm.extractvariables import extractVariablesFromExpressions, extractVariablesFromStatements

from dune.ufl import codegen
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
    if tensorMap is None:
        return [return_(lambda_(args=['const DomainValueType &phi'], code=return_(construct('RangeValueType', *[0 for i in range(len(testFunctions))], brace=True))))]

    var = Variable('std::tuple< RangeType, JacobianRangeType >', 'phi')
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
        # value = construct('RangeValueType', *[0 for i in range(len(testFunctions)), brace=True])
        value = construct('RangeValueType', *[0 for i in range(len(testFunctions))])
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


def compileUFL(form, tempVars=True):
    if isinstance(form, Equation):
        form = form.lhs - form.rhs
    if not isinstance(form, Form):
        raise Exception("ufl.Form expected.")
    if len(form.arguments()) < 2:
        raise Exception("Elliptic model requires form with at least two arguments.")

    x = SpatialCoordinate(form.ufl_cell())
    n = FacetNormal(form.ufl_cell())

    cellVolume = CellVolume(form.ufl_cell())
    maxCellEdgeLength = MaxCellEdgeLength(form.ufl_cell())
    minCellEdgeLength = MinCellEdgeLength(form.ufl_cell())

    facetArea = FacetArea(form.ufl_cell())
    maxFacetEdgeLength = MaxFacetEdgeLength(form.ufl_cell())
    minFacetEdgeLength = MinFacetEdgeLength(form.ufl_cell())

    phi = form.arguments()[0]

    u = form.arguments()[1]
    ubar = Coefficient(u.ufl_function_space())

    derivatives = gatherDerivatives(form, [phi, u])

    derivatives_phi = derivatives[0]
    derivatives_u = derivatives[1]
    derivatives_ubar = map_expr_dags(Replacer({u: ubar}), derivatives_u)

    integrands = Integrands(form.signature(), (d.ufl_shape for d in derivatives_u), (d.ufl_shape for d in derivatives_phi))
    try:
        integrands.field = u.ufl_function_space().field()
    except AttributeError:
        pass

    constants = dict()
    coefficients = dict()
    for coefficient in set(form.coefficients()):
        if coefficient.is_cellwise_constant():
            dimRange = (1 if len(coefficient.ufl_shape) == 0 else coefficient.ufl_shape[0])
            constants[coefficient] = integrands.addConstant('Dune::FieldVector< double, ' + str(dimRange) + ' >')
        else:
            try:
                coefficients[coefficient] = integrands.addCoefficient(coefficient.ufl_shape[0], coefficient.ufl_function_space().field())
            except AttributeError:
                coefficients[coefficient] = integrands.addCoefficient(coefficient.ufl_shape[0])

    integrals = splitForm(form, [phi])

    dform = apply_derivatives(derivative(action(form, ubar), ubar, u))
    linearizedIntegrals = splitForm(dform, [phi, u])

    if not set(integrals.keys()) <= {'cell', 'exterior_facet', 'interior_facet'}:
        raise Exception('unknown integral encountered in ' + str(set(integrals.keys())) + '.')

    def predefineCoefficients(predefined, x, side=None):
        for coefficient, idx in coefficients.items():
            for derivative in integrands.coefficient(idx, x, side=side):
                if side is None:
                    predefined[coefficient] = derivative
                elif side == 'Side::in':
                    predefined[coefficient('+')] = derivative
                elif side == 'Side::out':
                    predefined[coefficient('-')] = derivative
                coefficient = Grad(coefficient)

    if 'cell' in integrals.keys():
        arg = Variable(integrands.domainValueTuple(), 'u')

        predefined = {derivatives_u[i]: arg[i] for i in range(len(derivatives_u))}
        predefined[x] = integrands.spatialCoordinate('x')
        predefined[cellVolume] = integrands.cellVolume()
        predefined[maxCellEdgeLength] = maxEdgeLength(integrands.cellGeometry())
        predefined[minCellEdgeLength] = minEdgeLength(integrands.cellGeometry())
        predefined.update({c: integrands.constant(i) for c, i in constants.items()})
        predefineCoefficients(predefined, 'x')
        integrands.interior = generateUnaryCode(predefined, derivatives_phi, integrals['cell'], tempVars=tempVars)

        predefined = {derivatives_ubar[i]: arg[i] for i in range(len(derivatives_u))}
        predefined[x] = integrands.spatialCoordinate('x')
        predefined[cellVolume] = integrands.cellVolume()
        predefined[maxCellEdgeLength] = maxEdgeLength(integrands.cellGeometry())
        predefined[minCellEdgeLength] = minEdgeLength(integrands.cellGeometry())
        predefined.update({c: integrands.constant(i) for c, i in constants.items()})
        predefineCoefficients(predefined, 'x')
        integrands.linearizedInterior = generateUnaryLinearizedCode(predefined, derivatives_phi, derivatives_u, linearizedIntegrals.get('cell'), tempVars=tempVars)

    if 'exterior_facet' in integrals.keys():
        arg = Variable(integrands.domainValueTuple(), 'u')

        predefined = {derivatives_u[i]: arg[i] for i in range(len(derivatives_u))}
        predefined[x] = integrands.spatialCoordinate('x')
        predefined[n] = integrands.facetNormal('x')
        predefined[cellVolume] = integrands.cellVolume()
        predefined[maxCellEdgeLength] = maxEdgeLength(integrands.cellGeometry())
        predefined[minCellEdgeLength] = minEdgeLength(integrands.cellGeometry())
        predefined[facetArea] = integrands.facetArea()
        predefined[maxFacetEdgeLength] = maxEdgeLength(integrands.facetGeometry())
        predefined[minFacetEdgeLength] = minEdgeLength(integrands.facetGeometry())
        predefined.update({c: integrands.constant(i) for c, i in constants.items()})
        predefineCoefficients(predefined, 'x')
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
        predefined.update({c: integrands.constant(i) for c, i in constants.items()})
        predefineCoefficients(predefined, 'x')
        integrands.linearizedBoundary = generateUnaryLinearizedCode(predefined, derivatives_phi, derivatives_u, linearizedIntegrals.get('exterior_facet'), tempVars=tempVars)

    if 'interior_facet' in integrals.keys():
        argIn = Variable(integrands.domainValueTuple(), 'uIn')
        argOut = Variable(integrands.domainValueTuple(), 'uOut')

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
        predefined.update({c: integrands.constant(i) for c, i in constants.items()})
        predefineCoefficients(predefined, 'xIn', 'Side::in')
        predefineCoefficients(predefined, 'xOut', 'Side::out')
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
        predefined.update({c: integrands.constant(i) for c, i in constants.items()})
        predefineCoefficients(predefined, 'xIn', 'Side::in')
        predefineCoefficients(predefined, 'xOut', 'Side::out')
        integrands.linearizedSkeleton = generateBinaryLinearizedCode(predefined, derivatives_phi, derivatives_u, linearizedIntegrals.get('interior_facet'), tempVars=tempVars)

    coefficients.update(constants)
    return integrands, coefficients
