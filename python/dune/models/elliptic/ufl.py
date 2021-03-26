from __future__ import absolute_import, division, print_function, unicode_literals

from ufl import Coefficient, Form, FiniteElementBase, FunctionSpace, SpatialCoordinate
from ufl.core.expr import Expr
from ufl import action, adjoint, as_vector, derivative, div, dx, inner, replace
from ufl import replace, TestFunction, TrialFunction
from ufl.algorithms import expand_compounds, expand_derivatives, expand_indices
from ufl.algorithms.analysis import extract_arguments_and_coefficients
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.classes import Indexed
from ufl.differentiation import Grad
from ufl.equation import Equation
from ufl.core.multiindex import FixedIndex, MultiIndex
from ufl import UFLException

from dune.ufl import DirichletBC, GridFunction
from dune.ufl import codegen
from dune.ufl.tensors import ExprTensor
from dune.ufl.linear import splitMultiLinearExpr

from dune.common.hashit import hashIt

from dune.source.cplusplus import UnformattedExpression, Block
from dune.source.cplusplus import Declaration, NameSpace, SwitchStatement, TypeAlias, UnformattedBlock, Variable
from dune.source.cplusplus import assign, construct, return_

from .model import EllipticModel

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
    return [assign(result, results[0])]

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

def toFileName(value):
    import unicodedata
    value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore')
    value = unicode(re.sub('[^\w\s-]', '', value).strip().lower())
    value = unicode(re.sub('[-\s]+', '-', value))
    return value

def modelSignature(form,*args):
    dirichletBCs = [str(arg.ufl_value) for arg in args if isinstance(arg, DirichletBC)]
    expr = [str(arg) for arg in args if isinstance(arg, Expr)]
    sig = form.signature()
    if len(dirichletBCs) > 0:
        dirichletBCs.append(sig)
        sig = hashIt( dirichletBCs )
    if len(expr) > 0:
        expr.append(sig)
        sig = hashIt( expr )
    return sig

def compileUFL(form, patch, *args, **kwargs):
    if isinstance(form, Equation):
        form = form.lhs - form.rhs
    if not isinstance(form, Form):
        raise Exception("ufl.Form expected.")
    if len(form.arguments()) < 2:
        raise Exception("Elliptic model requires form with at least two arguments.")

    phi_, u_ = form.arguments()

    if phi_.ufl_function_space().scalar:
        phi = TestFunction(phi_.ufl_function_space().toVectorSpace())
        form = replace(form,{phi_:phi[0]})
    else:
        phi = phi_
    if u_.ufl_function_space().scalar:
        u = TrialFunction(u_.ufl_function_space().toVectorSpace())
        form = replace(form,{u_:u[0]})
    else:
        u = u_
    _, coeff_ = extract_arguments_and_coefficients(form)
    coeff_ = set(coeff_)

    # added for dirichlet treatment same as elliptic model
    dirichletBCs = [arg for arg in args if isinstance(arg, DirichletBC)]
    # remove the dirichletBCs
    arg = [arg for arg in args if not isinstance(arg, DirichletBC)]
    for dBC in dirichletBCs:
        _, coeff__ = extract_arguments_and_coefficients(dBC.ufl_value)
        coeff_ |= set(coeff__)
    if patch is not None:
        for a in patch:
            try:
                _, coeff__ = extract_arguments_and_coefficients(a)
                coeff_ |= set(coeff__)
            except:
                pass # a might be a float/int and not a ufl expression

    coeff = {c : c.toVectorCoefficient()[0] for c in coeff_ if len(c.ufl_shape) == 0 and not c.is_cellwise_constant()}
    form = replace(form,coeff)
    for bc in dirichletBCs:
        bc.ufl_value = replace(bc.ufl_value, coeff)
    if patch is not None:
        patch = [a if not isinstance(a, Expr) else replace(a,coeff) for a in patch]

    phi = form.arguments()[0]
    dimRange = phi.ufl_shape[0]

    u = form.arguments()[1]
    du = Grad(u)
    d2u = Grad(du)
    ubar = Coefficient(u.ufl_function_space())
    dubar = Grad(ubar)
    d2ubar = Grad(dubar)
    dimDomain = u.ufl_shape[0]

    x = SpatialCoordinate(form.ufl_cell())

    try:
        field = u.ufl_function_space().field
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

    if patch is not None:
        model = EllipticModel(dimDomain, dimRange, u, modelSignature(form,*patch,*args))
    else:
        model = EllipticModel(dimDomain, dimRange, u, modelSignature(form,None,*args))
    model._replaceCoeff = coeff

    model.hasNeumanBoundary = not boundarySource.is_zero()

    #expandform = expand_indices(expand_derivatives(expand_compounds(equation.lhs)))
    #if expandform == adjoint(expandform):
    #    model.symmetric = 'true'
    model.field = field

    dirichletBCs = [arg for arg in args if isinstance(arg, DirichletBC)]
    # deprecated
    # if "dirichlet" in kwargs:
    #     dirichletBCs += [DirichletBC(u.ufl_function_space(), as_vector(value), bndId) for bndId, value in kwargs["dirichlet"].items()]

    uflCoefficients = set(form.coefficients())
    for bc in dirichletBCs:
        _, c = extract_arguments_and_coefficients(bc.ufl_value)
        uflCoefficients |= set(c)
    if patch is not None:
        for a in patch:
            if isinstance(a, Expr):
                _, c = extract_arguments_and_coefficients(a)
                uflCoefficients |= set(c)

    constants = dict()
    coefficients = dict()

    for coefficient in uflCoefficients:
        try:
            name = getattr(coefficient, "name")
        except AttributeError:
            name = str(coefficient)
        if coefficient.is_cellwise_constant():
            try:
                parameter = getattr(coefficient, "parameter")
            except AttributeError:
                parameter = None
            if len(coefficient.ufl_shape) == 0:
                constants[coefficient] = model.addConstant('double', name=name, parameter=parameter)
            elif len(coefficient.ufl_shape) == 1:
                constants[coefficient] = model.addConstant('Dune::FieldVector< double, ' + str(coefficient.ufl_shape[0]) + ' >', name=name, parameter=parameter)
            else:
                Exception('Currently, only scalars and vectors are supported as constants')
        else:
            shape = coefficient.ufl_shape[0]
            try:
                coefficients[coefficient] = model.addCoefficient(
                        shape,
                        coefficient._typeName,
                        name=name,
                        field=coefficient.ufl_function_space().field)
            except AttributeError:
                coefficients[coefficient] = model.addCoefficient(
                        shape,
                        coefficient._typeName,
                        name=name)

    model.coefficients = coefficients
    model.constants = constants

    tempVars = kwargs.get("tempVars", True)

    predefined = {u: model.arg_u, du: model.arg_du, d2u: model.arg_d2u}
    predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
    model.predefineCoefficients(predefined,'x')
    model.source = generateCode(predefined, source, tempVars=tempVars)
    model.diffusiveFlux = generateCode(predefined, diffusiveFlux, tempVars=tempVars)
    predefined.update({ubar: model.arg_ubar, dubar: model.arg_dubar, d2ubar: model.arg_d2ubar})
    model.linSource = generateCode(predefined, linSource, tempVars=tempVars)
    model.linDiffusiveFlux = generateCode(predefined, linDiffusiveFlux, tempVars=tempVars)

    # model.linNVSource = generateCode({u: arg, du: darg, d2u: d2arg, ubar: argbar, dubar: dargbar, d2ubar: d2argbar}, linNVSource, model.coefficients, tempVars)

    predefined = {u: model.arg_u}
    predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
    model.predefineCoefficients(predefined,'x')
    model.alpha = generateCode(predefined, boundarySource, tempVars=tempVars)
    predefined.update({ubar: model.arg_ubar})
    model.linAlpha = generateCode(predefined, linBoundarySource, tempVars=tempVars)

    predefined = {u: model.arg_u, du: model.arg_du, d2u: model.arg_d2u}
    predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
    model.predefineCoefficients(predefined,'x')
    model.fluxDivergence = generateCode(predefined, fluxDivergence, tempVars=tempVars)

    if dirichletBCs:
        model.hasDirichletBoundary = True

        predefined = {}
        predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
        model.predefineCoefficients(predefined,'x')

        maxId = 0
        codeDomains = []
        bySubDomain = dict()
        neuman = []
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
            if isinstance(bc.subDomain,int):
                bySubDomain[bc.subDomain] = value,neuman
                maxId = max(maxId, bc.subDomain)
            else:
                domain = ExprTensor(())
                for key in domain.keys():
                    domain[key] = Indexed(bc.subDomain, MultiIndex(tuple(FixedIndex(k) for k in key)))
                codeDomains.append( (value,neuman,domain) )
        defaultCode = []
        defaultCode.append(Declaration(Variable('int', 'domainId')))
        defaultCode.append(Declaration(Variable('auto', 'tmp0'),
            initializer=UnformattedExpression('auto','intersection.geometry().center()')))
        for i,v in enumerate(codeDomains):
            block = Block()
            defaultCode.append(
                    generateDirichletDomainCode(predefined, v[2], tempVars=tempVars))
            defaultCode.append('if (domainId)')
            block = UnformattedBlock()
            block.append('std::fill( dirichletComponent.begin(), dirichletComponent.end(), ' + str(maxId+i+1) + ' );')
            if len(v[1])>0:
                [block.append('dirichletComponent[' + str(c) + '] = 0;') for c in v[1]]
            block.append('return true;')
            defaultCode.append(block)
        defaultCode.append(return_(False))

        bndId = Variable('const int', 'bndId')
        getBndId = UnformattedExpression('int', 'BoundaryIdProviderType::boundaryId( ' + model.arg_i.name + ' )')
        switch = SwitchStatement(bndId, default=defaultCode)
        for i,v in bySubDomain.items():
            code = []
            if len(v[1])>0:
                [code.append('dirichletComponent[' + str(c) + '] = 0;') for c in v[1]]
            code.append(return_(True))
            switch.append(i, code)
        model.isDirichletIntersection = [Declaration(bndId, initializer=getBndId),
                                         UnformattedBlock('std::fill( dirichletComponent.begin(), dirichletComponent.end(), ' + bndId.name + ' );'),
                                         switch
                                        ]

        switch = SwitchStatement(model.arg_bndId, default=assign(model.arg_r, construct("RRangeType", 0)))
        for i, v in bySubDomain.items():
            switch.append(i, generateDirichletCode(predefined, v[0], tempVars=tempVars))
        for i,v in enumerate(codeDomains):
            switch.append(i+maxId+1, generateDirichletCode(predefined, v[0], tempVars=tempVars))
        model.dirichlet = [switch]

    return model
