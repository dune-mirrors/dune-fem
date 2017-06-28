from __future__ import absolute_import, division, print_function, unicode_literals

from ufl import Coefficient, Form, FiniteElementBase, FunctionSpace, SpatialCoordinate
from ufl import action, adjoint, as_vector, derivative, div, dx, inner, replace
from ufl.algorithms import expand_compounds, expand_derivatives, expand_indices
from ufl.algorithms.analysis import extract_arguments_and_coefficients
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.classes import Indexed
from ufl.differentiation import Grad
from ufl.equation import Equation
from ufl.core.multiindex import FixedIndex, MultiIndex

from dune.ufl import DirichletBC, GridFunction
from dune.ufl import codegen
from dune.ufl.tensors import ExprTensor
from dune.ufl.linear import splitMultiLinearExpr

from dune.source.cplusplus import UnformattedExpression
from dune.source.cplusplus import Declaration, NameSpace, SwitchStatement, TypeAlias, UnformattedBlock, Variable
from dune.source.cplusplus import assign, construct, return_

from .model import EllipticModel


def splitUFLForm(form):
    phi = form.arguments()[0]
    dphi = Grad(phi)

    source = ExprTensor(phi.ufl_shape)
    diffusiveFlux = ExprTensor(dphi.ufl_shape)
    boundarySource = ExprTensor(phi.ufl_shape)

    form = expand_indices(expand_derivatives(expand_compounds(form)))
    for integral in form.integrals():
        if integral.integral_type() == 'cell':
            fluxExprs = splitMultiLinearExpr(integral.integrand(), [phi])
            for op in fluxExprs:
                if op[0] == phi:
                    source = source + fluxExprs[op]
                elif op[0] == dphi:
                    diffusiveFlux = diffusiveFlux + fluxExprs[op]
                else:
                    raise Exception('Invalid derivative encountered in bulk integral: ' + str(op[0]))
        elif integral.integral_type() == 'exterior_facet':
            fluxExprs = splitMultiLinearExpr(integral.integrand(), [phi])
            for op in fluxExprs:
                if op[0] == phi:
                    boundarySource = boundarySource + fluxExprs[op]
                else:
                    raise Exception('Invalid derivative encountered in boundary integral: ' + str(op[0]))
        else:
            raise NotImplementedError('Integrals of type ' + integral.integral_type() + ' are not supported.')

    return source, diffusiveFlux, boundarySource


#def splitUFL2(u,du,d2u,tree):
#    tree0 = ExprTensor(u.ufl_shape)
#    tree1 = ExprTensor(u.ufl_shape)
#    tree2 = ExprTensor(u.ufl_shape)
#
#    for index in tree.keys():
#        q = splitMultiLinearExpr(tree[index], [u])
#        if q==0: continue
#        for op in q:
#            if not isinstance(op, tuple) or (len(op) != 1):
#                raise Exception('Missing trial function in bulk integral')
#            if op[0] == u:
#                tree0[index] += sum(i[0]*i[1] for i in zip(q[op].as_ufl(),u))
#            elif op[0] == du:
#                for r in range(du.ufl_shape[0]):
#                    for d in range(du.ufl_shape[1]):
#                        tree1[index] += q[op].as_ufl()[r,d]*du[r,d]
#            elif op[0] == d2u:
#                for r in range(d2u.ufl_shape[0]):
#                    for d1 in range(d2u.ufl_shape[1]):
#                        for d2 in range(d2u.ufl_shape[2]):
#                            tree2[index] += q[op].as_ufl()[r,d1,d2]*d2u[r,d1,d2]
#            else:
#                raise Exception('Invalid trial function derivative encountered in bulk integral: ' + str(op[0]))
#    return tree0, tree1, tree2


def generateCode(predefined, tensor, tempVars=True):
    keys = tensor.keys()
    expressions = [tensor[i] for i in keys]
    preamble, results = codegen.generateCode(predefined, expressions, tempVars=tempVars)
    result = Variable('auto', 'result')
    return preamble + [assign(result[i], r) for i, r in zip(keys, results)]


def compileUFL(form, *args, **kwargs):
    if isinstance(form, Equation):
        form = form.lhs - form.rhs
    if not isinstance(form, Form):
        raise Exception("ufl.Form expected.")
    if len(form.arguments()) < 2:
        raise Exception("Elliptic model requires form with at least two arguments.")

    phi = form.arguments()[0]
    dimRange = phi.ufl_shape[0]

    u = form.arguments()[1]
    du = Grad(u)
    d2u = Grad(du)
    ubar = Coefficient(u.ufl_function_space())
    dubar = Grad(ubar)
    d2ubar = Grad(dubar)

    x = SpatialCoordinate(form.ufl_cell())

    try:
        field = u.ufl_function_space().field()
    except AttributeError:
        field = "double"

    # if exact solution is passed in subtract a(u,.) from the form
    if "exact" in kwargs:
        b = replace(form, {u: as_vector(kwargs["exact"])} )
        form = form - b

    dform = apply_derivatives(derivative(action(form, ubar), ubar, u))

    source, diffusiveFlux, boundarySource = splitUFLForm(form)
    linSource, linDiffusiveFlux, linBoundarySource = splitUFLForm(dform)
    fluxDivergence, _, _ = splitUFLForm(inner(source.as_ufl() - div(diffusiveFlux.as_ufl()), phi) * dx(0))

    # split linNVSource off linSource
    # linSources = splitUFL2(u, du, d2u, linSource)
    # linNVSource = linSources[2]
    # linSource = linSources[0] + linSources[1]

    model = EllipticModel(dimRange, form.signature())

    model.hasNeumanBoundary = not boundarySource.is_zero()

    #expandform = expand_indices(expand_derivatives(expand_compounds(equation.lhs)))
    #if expandform == adjoint(expandform):
    #    model.symmetric = 'true'
    model.field = field

    dirichletBCs = [arg for arg in args if isinstance(arg, DirichletBC)]
    if "dirichlet" in kwargs:
        dirichletBCs += [DirichletBC(u.ufl_function_space(), as_vector(value), bndId) for bndId, value in kwargs["dirichlet"].items()]

    uflCoefficients = set(form.coefficients())
    for bc in dirichletBCs:
        _, c = extract_arguments_and_coefficients(bc.value)
        uflCoefficients |= set(c)

    constants = dict()
    coefficients = dict()
    for coefficient in uflCoefficients:
        try:
            name = getattr(coefficient, "name")
        except AttributeError:
            name = str(coefficient)
        if coefficient.is_cellwise_constant():
            if len(coefficient.ufl_shape) == 0:
                constants[coefficient] = model.addConstant('double', name)
            elif len(coefficient.ufl_shape) == 1:
                constants[coefficient] = model.addConstant('Dune::FieldVector< double, ' + str(coefficient.ufl_shape[0]) + ' >', name)
            else:
                Exception('Currently, only scalars and vectors are supported as constants')
        else:
            try:
                coefficients[coefficient] = model.addCoefficient(coefficient.ufl_shape[0], name, coefficient.ufl_function_space().field())
            except AttributeError:
                coefficients[coefficient] = model.addCoefficient(coefficient.ufl_shape[0], name)

    tempVars = kwargs.get("tempVars", True)

    def predefineCoefficients(predefined, x):
        for coefficient, idx in coefficients.items():
            for derivative in model.coefficient(idx, x):
                predefined[coefficient] = derivative
                coefficient = Grad(coefficient)
        predefined.update({c: model.constant(i) for c, i in constants.items()})

    predefined = {u: model.arg_u, du: model.arg_du}
    predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
    predefineCoefficients(predefined, model.arg_x)
    model.source = generateCode(predefined, source, tempVars=tempVars)
    model.diffusiveFlux = generateCode(predefined, diffusiveFlux, tempVars=tempVars)
    predefined.update({ubar: model.arg_ubar, dubar: model.arg_dubar})
    model.linSource = generateCode(predefined, linSource, tempVars=tempVars)
    model.linDiffusiveFlux = generateCode(predefined, linDiffusiveFlux, tempVars=tempVars)

    # model.linNVSource = generateCode({u: arg, du: darg, d2u: d2arg, ubar: argbar, dubar: dargbar, d2ubar: d2argbar}, linNVSource, model.coefficients, tempVars)

    predefined = {u: model.arg_u}
    predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
    predefineCoefficients(predefined, model.arg_x)
    model.alpha = generateCode(predefined, boundarySource, tempVars=tempVars)
    predefined.update({ubar: model.arg_ubar})
    model.linAlpha = generateCode(predefined, linBoundarySource, tempVars=tempVars)

    predefined = {u: model.arg_u, du: model.arg_du, d2u: model.arg_d2u}
    predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
    predefineCoefficients(predefined, model.arg_x)
    model.fluxDivergence = generateCode(predefined, fluxDivergence, tempVars=tempVars)

    if dirichletBCs:
        model.hasDirichletBoundary = True

        bySubDomain = dict()
        for bc in dirichletBCs:
            if bc.subDomain in bySubDomain:
                raise Exception('Multiply defined Dirichlet boundary for subdomain ' + str(bc.subDomain))

            if not isinstance(bc.functionSpace, (FunctionSpace, FiniteElementBase)):
                raise Exception('Function space must either be a ufl.FunctionSpace or a ufl.FiniteElement')
            if isinstance(bc.functionSpace, FunctionSpace) and (bc.functionSpace != u.ufl_function_space()):
                raise Exception('Cannot handle boundary conditions on subspaces, yet')
            if isinstance(bc.functionSpace, FiniteElementBase) and (bc.functionSpace != u.ufl_element()):
                raise Exception('Cannot handle boundary conditions on subspaces, yet')

            value = ExprTensor(u.ufl_shape)
            for key in value.keys():
                value[key] = Indexed(bc.value, MultiIndex(tuple(FixedIndex(k) for k in key)))
            bySubDomain[bc.subDomain] = value

        bndId = Variable('const int', 'bndId')
        getBndId = UnformattedExpression('int', 'Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType >::boundaryId( ' + model.arg_i.name + ' )')

        switch = SwitchStatement(bndId, default=return_(False))
        for i in bySubDomain:
            switch.append(i, return_(True))
        model.isDirichletIntersection = [Declaration(bndId, initializer=UnformattedExpression('int', 'BoundaryIdProviderType::boundaryId( ' + model.arg_i.name + ' )')),
                                         UnformattedBlock('std::fill( dirichletComponent.begin(), dirichletComponent.end(), ' + bndId.name + ' );'),
                                         switch
                                        ]

        switch = SwitchStatement(model.arg_bndId, default=assign(model.arg_r, construct("RangeType", 0)))
        predefined = {}
        predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
        predefineCoefficients(predefined, model.arg_x)
        for i, value in bySubDomain.items():
            switch.append(i, generateCode(predefined, value, tempVars=tempVars))
        model.dirichlet = [switch]

    coefficients.update(constants)
    return model, coefficients
