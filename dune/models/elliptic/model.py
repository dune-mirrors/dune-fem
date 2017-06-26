from __future__ import division, print_function, unicode_literals

from dune.common.compatibility import isInteger

from dune.source.builtin import get, make_shared
from dune.source.cplusplus import UnformattedExpression
from dune.source.cplusplus import AccessModifier, Constructor, Declaration, Function, Method, NameSpace, Struct, TypeAlias, UnformattedBlock, Variable
from dune.source.cplusplus import assign, construct, dereference, lambda_, nullptr, return_, this
from dune.source.cplusplus import SourceWriter
from dune.source.fem import declareFunctionSpace


class EllipticModel:
    def __init__(self, dimRange, signature):
        assert isInteger(dimRange)
        self.dimRange = dimRange
        self.init = []
        self.vars = None
        self.signature = signature
        self.field = "double"

        self._constants = []
        self._constantNames = {}
        self._coefficients = []

        self.arg_r = Variable("RangeType &", "result")
        self.arg_dr = Variable("JacobianRangeType &", "result")

        self.arg_x = Variable("const Point &", "x")
        self.arg_u = Variable("const RangeType &", "u")
        self.arg_du = Variable("const JacobianRangeType &", "du")
        self.arg_d2u = Variable("const HessianRangeType &", "d2u")
        self.arg_ubar = Variable("const RangeType &", "ubar")
        self.arg_dubar = Variable("const JacobianRangeType &", "dubar")
        self.arg_d2ubar = Variable("const HessianRangeType &", "d2ubar")

        self.arg_i = Variable("const IntersectionType &", "intersection")
        self.arg_bndId = Variable("int", "bndId")

        self.source = [assign(self.arg_r, construct("RangeType", 0))]
        self.linSource = [assign(self.arg_r, construct("RangeType", 0))]
        # self.linNVSource = [assign(self.arg_r, construct("RangeType", 0))]
        self.diffusiveFlux = [assign(self.arg_dr, construct("JacobianRangeType", 0))]
        self.linDiffusiveFlux = [assign(self.arg_dr, construct("JacobianRangeType", 0))]
        self.fluxDivergence = [assign(self.arg_r, construct("RangeType", 0))]
        self.alpha = [assign(self.arg_r, construct("RangeType", 0))]
        self.linAlpha = [assign(self.arg_r, construct("RangeType", 0))]

        self.hasDirichletBoundary = False
        self.hasNeumanBoundary = False
        self.isDirichletIntersection = [return_(False)]
        self.dirichlet = [assign(self.arg_r, construct("RangeType", 0))]
        self.symmetric = False

    def addCoefficient(self, dimRange, name=None, field="double"):
        idx = len(self._coefficients)
        self._coefficients.append({'dimRange': dimRange, 'name': name, 'field': field})
        return idx

    def addConstant(self, cppType, name=None):
        idx = len(self._constants)
        self._constants.append(cppType)
        if name is not None:
            self._constantNames[name] = idx
        return idx

    def constant(self, idx):
        return UnformattedExpression(self._constants[idx], 'constant< ' + str(idx) + ' >()')

    def constant(self, idx):
        return UnformattedExpression(self._constants[idx], 'constant< ' + str(idx) + ' >()')

    def coefficient(self, idx, x):
        coefficient = []
        for t, n in (('RangeType', 'evaluate'), ('JacobianRangeType', 'jacobian'), ('HessianRangeType', 'hessian')):
            result = Variable('typename std::tuple_element_t< ' + str(idx) + ', CoefficientFunctionSpaceTupleType >::' + t, 'result')
            code = [Declaration(result),
                    UnformattedExpression('void', 'std::get< ' + str(idx) + ' >( coefficients_ ).' + n + '( x, ' + result.name + ' )'),
                    return_(result)]
            coefficient += [lambda_(capture=[this], args=['auto x'], code=code)(x)]
        return coefficient

    @property
    def hasCoefficients(self):
        return bool(self._coefficients)

    @property
    def hasConstants(self):
        return bool(self._constants)

    def code(self, name='Model', targs=[]):
        constants_ = Variable('std::tuple< ' + ', '.join('std::shared_ptr< ' + c  + ' >' for c in self._constants) + ' >', 'constants_')
        coefficients_ = Variable('std::tuple< ' + ', '.join('Coefficient' + str(i) for i, c in enumerate(self._coefficients)) + ' >', 'coefficients_')
        entity_ = Variable('const EntityType *', 'entity_')

        code = Struct(name, targs=(['class GridPart'] + ['class Coefficient' + str(i) for i, c in enumerate(self._coefficients)] + targs))

        code.append(TypeAlias("GridPartType", "GridPart"))
        code.append(TypeAlias("EntityType", "typename GridPart::template Codim< 0 >::EntityType"))
        code.append(TypeAlias("IntersectionType", "typename GridPart::IntersectionType"))

        code.append(declareFunctionSpace("typename GridPartType::ctype", SourceWriter.cpp_fields(self.field), UnformattedExpression("int", "GridPartType::dimensionworld"), self.dimRange))
        code.append(Declaration(Variable("const int", "dimLocal"), initializer=UnformattedExpression("int", "GridPartType::dimension"), static=True))

        if self.hasConstants:
            code.append(TypeAlias("ConstantType", "typename std::tuple_element_t< i, " + constants_.cppType + " >::element_type", targs=["std::size_t i"]))
            code.append(Declaration(Variable("const std::size_t", "numConstants"), initializer=len(self._constants), static=True))

        if self.hasCoefficients:
            coefficientSpaces = ["Dune::Fem::FunctionSpace< DomainFieldType, " + SourceWriter.cpp_fields(c['field']) + ", dimDomain, " + str(c['dimRange']) + " >" for c in self._coefficients]
            code.append(TypeAlias("CoefficientFunctionSpaceTupleType", "std::tuple< " + ", ".join(coefficientSpaces) + " >"))
            code.append(TypeAlias('CoefficientType', 'std::tuple_element_t< i, ' + coefficients_.cppType + ' >', targs=['std::size_t i']))

        if self.hasCoefficients:
            args = [Variable("const Coefficient" + str(i) + " &", "coefficient" + str(i)) for i, c in enumerate(self._coefficients)]
            init = ["coefficients_( " + ", ".join("coefficient" + str(i) for i, c in enumerate(self._coefficients)) + " )"]
            constructor = Constructor(args=args, init=init)
        else:
            constructor = Constructor()
        if self.hasConstants:
            constructor.append([assign(get(str(i))(constants_), make_shared(c)()) for i, c in enumerate(self._constants)])
        code.append(constructor)

        init = ['entity_ = &entity;']
        init += ['std::get< ' + str(i) + ' >( ' + coefficients_.name + ' ).init( entity );' for i, c in enumerate(self._coefficients)]
        init = [UnformattedBlock(init)] + self.init + [return_(True)]
        code.append(Method('bool', 'init', args=['const EntityType &entity'], code=init, const=True))

        code.append(Method('const EntityType &', 'entity', code=return_(dereference(entity_)), const=True))
        code.append(Method('std::string', 'name', const=True, code=return_(UnformattedExpression('const char *', '"' + name + '"'))))

        code.append(TypeAlias("BoundaryIdProviderType", "Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType >"))
        code.append(Declaration(Variable("const bool", "symmetric"), initializer=self.symmetric, static=True))

        code.append(Method('void', 'source', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_r], code=self.source, const=True))
        code.append(Method('void', 'linSource', targs=['class Point'], args=[self.arg_ubar, self.arg_dubar, self.arg_x, self.arg_u, self.arg_du, self.arg_r], code=self.linSource, const=True))

        code.append(Method('void', 'diffusiveFlux', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_dr], code=self.diffusiveFlux, const=True))
        code.append(Method('void', 'linDiffusiveFlux', targs=['class Point'], args=[self.arg_ubar, self.arg_dubar, self.arg_x, self.arg_u, self.arg_du, self.arg_dr], code=self.linDiffusiveFlux, const=True))

        code.append(Method('void', 'fluxDivergence', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_d2u, self.arg_r], code=self.fluxDivergence, const=True))

        code.append(Method('void', 'alpha', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_r], code=self.alpha, const=True))
        code.append(Method('void', 'linAlpha', targs=['class Point'], args=[self.arg_ubar, self.arg_x, self.arg_u, self.arg_r], code=self.linAlpha, const=True))

        code.append(Method('bool', 'hasDirichletBoundary', const=True, code=return_(self.hasDirichletBoundary)))
        code.append(Method('bool', 'hasNeumanBoundary', const=True, code=return_(self.hasNeumanBoundary)))

        code.append(Method('bool', 'isDirichletIntersection', args=[self.arg_i, 'Dune::FieldVector< int, dimRange > &dirichletComponent'], code=self.isDirichletIntersection, const=True))
        code.append(Method('void', 'dirichlet', targs=['class Point'], args=[self.arg_bndId, self.arg_x, self.arg_r], code=self.dirichlet, const=True))

        if self.hasConstants:
            code.append(Method("const ConstantType< i > &", "constant", targs=["std::size_t i"], code=return_(dereference(get("i")(constants_))), const=True))
            code.append(Method("ConstantType< i > &", "constant", targs=["std::size_t i"], code=return_(dereference(get("i")(constants_)))))

        if self.hasCoefficients:
            code.append(Method("const CoefficientType< i > &", "coefficient", targs=["std::size_t i"], code=return_(get("i")(coefficients_)), const=True))
            code.append(Method("CoefficientType< i > &", "coefficient", targs=["std::size_t i"], code=return_(get("i")(coefficients_))))

        code.append(AccessModifier("private"))
        code.append(Declaration(entity_, nullptr, mutable=True))
        if self.hasConstants:
            code.append(Declaration(constants_, mutable=True))
        if self.hasCoefficients:
            code.append(Declaration(coefficients_, mutable=True))
        return code

    #def write(self, sourceWriter, name='Model', targs=[]):
    #    sourceWriter.emit(self.code(name=name, targs=targs))

    def exportSetConstant(self, sourceWriter, modelClass='Model', wrapperClass='ModelWrapper'):
        sourceWriter.openFunction('std::size_t renumberConstants', args=['pybind11::handle &obj'])
        sourceWriter.emit('std::string id = pybind11::str( obj );')
        sourceWriter.emit('if( obj.attr("name") ) id = pybind11::str(obj.attr("name"));')
        for name, number in self._constantNames.items():
            sourceWriter.emit('if (id == "' + name + '") return ' + str(number) + ';')
        sourceWriter.emit('throw pybind11::value_error("coefficient \'" + id + "\' has not been registered");')
        sourceWriter.closeFunction()

        sourceWriter.openFunction('void setConstant', targs=['std::size_t i'], args=[modelClass + ' &model', 'pybind11::handle value'])
        sourceWriter.emit('model.template constant< i >() = value.template cast< typename ' + modelClass + '::ConstantType< i > >();')
        sourceWriter.closeFunction()

        sourceWriter.openFunction('auto defSetConstant', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
        sourceWriter.emit(TypeAlias('Dispatch', 'std::function< void( ' + modelClass + ' &model, pybind11::handle ) >'))
        sourceWriter.emit('std::array< Dispatch, sizeof...( i ) > dispatch = {{ Dispatch( setConstant< i > )... }};')
        sourceWriter.emit('')
        sourceWriter.emit('return [ dispatch ] ( ' + wrapperClass + ' &model, pybind11::handle coeff, pybind11::handle value ) {')
        sourceWriter.emit('    std::size_t k = renumberConstants( coeff );')
        sourceWriter.emit('    if( k >= dispatch.size() )')
        sourceWriter.emit('      throw std::range_error( "No such coefficient: "+std::to_string(k)+" >= "+std::to_string(dispatch.size()) );' )
        sourceWriter.emit('    dispatch[ k ]( model.impl(), value );')
        sourceWriter.emit('    return k;')
        sourceWriter.emit('  };')
        sourceWriter.closeFunction()

    def export(self, sourceWriter, modelClass='Model', wrapperClass='ModelWrapper'):
        if self.hasConstants:
            sourceWriter.emit('cls.def( "setConstant", defSetConstant( std::make_index_sequence< ' + modelClass + '::numConstants >() ) );')
        coefficients = [('Dune::FemPy::VirtualizedGridFunction< GridPart, Dune::FieldVector< ' + SourceWriter.cpp_fields(c['field']) + ', ' + str(c['dimRange']) + ' > >') for c in self._coefficients]
        sourceWriter.emit('')
        sourceWriter.emit('cls.def( "__init__", [] ( ' + ', '.join([wrapperClass + ' &self'] + ['const ' + c + ' &coefficient' + str(i) for i, c in enumerate(coefficients)]) + ' ) {')
        if self.hasCoefficients:
            sourceWriter.emit('  new (&self) ' + wrapperClass + '( ' + ', '.join('coefficient' + str(i) + '.localFunction()' for i, c in enumerate(coefficients)) + ' );')
        else:
            sourceWriter.emit('  new (&self) ' + wrapperClass + '();')
        #if self.coefficients:
        #    sourceWriter.emit('  const int size = std::tuple_size<Coefficients>::value;')
        #    sourceWriter.emit('  auto dispatch = defSetCoefficient( std::make_index_sequence<size>() );' )
        #    sourceWriter.emit('  std::vector<bool> coeffSet(size,false);')
        #    sourceWriter.emit('  for (auto item : coeff) {')
        #    sourceWriter.emit('    int k = dispatch(instance, item.first, item.second); ')
        #    sourceWriter.emit('    coeffSet[k] = true;')
        #    sourceWriter.emit('  }')
        #    sourceWriter.emit('  if ( !std::all_of(coeffSet.begin(),coeffSet.end(),[](bool v){return v;}) )')
        #    sourceWriter.emit('    throw pybind11::key_error("need to set all coefficients during construction");')
        if self.hasCoefficients:
            sourceWriter.emit('  }, ' + ', '.join('pybind11::keep_alive< 1, ' + str(i) + ' >()' for i, c in enumerate(coefficients, start=2)) + ' );')
        else:
            sourceWriter.emit('  } );')

    def codeCoefficient(self, code, coefficients, constants):
        """find coefficients/constants in code string and do replacements
        """
        for name, value in coefficients.items():
            if not any(name == c['name'] for c in self._coefficients):
                self.addCoefficient(value.dimRange, name)
        for name, dimRange in constants.items():
            if name not in self._constantNames:
                self.addConstant('Dune::FieldVector< double, ' + str(dimRange) + ' >', name)

        if '@const:' in code:
            codeCst = code.split('@const:')
            import itertools
            for name in set([''.join(itertools.takewhile(str.isalpha, str(c))) for c in codeCst[1:]]):
                if name not in self._constantNames:
                    cname = '@const:' + name
                    afterName = code.split(cname)[1:]
                    if afterName[0][0] == '[':
                        beforeText = [an.split(']')[0].split('[')[1] for an in afterName]
                        dimRange = max( [int(bt) for bt in beforeText] ) + 1
                    else:
                        dimRange = 1
                    self.addConstant('Dune::FieldVector< double, ' + str(dimRange) + ' >', name)

        for i, c in enumerate(self._coefficients):
            jacname = '@jac:' + c['name']
            if jacname in code:
                varname = 'dc' + str(i)
                code = code.replace(jacname, varname)
                decl = 'CoefficientJacobianRangeType< ' + str(i) + ' > ' + varname + ';'
                if not decl in code:
                    code = decl + '\ncoefficient< ' + str(i) + ' >().jacobian( x, ' + varname + ' );' + code
            gfname = '@gf:' + c['name']
            if gfname in code:
                varname = 'c' + str(i)
                code = code.replace(gfname, varname)
                decl = 'CoefficientRangeType< ' + str(i) + ' > c' + str(i) + ';'
                if not decl in code:
                    code = decl + '\ncoefficient< ' + str(i) + ' >().evaluate( x, ' + varname + ' );' + code

        for name, i in self._constantNames.items():
            cname = '@const:' + name
            if cname in code:
                varname = 'cc' + str(i)
                code = code.replace(cname, varname)
                init = 'const ' + self._constants[i] + ' &' + varname + ' = constant< ' + str(i) + ' >();'
                if not init in code:
                    code = init + code

    def appendCode(self, key, code, **kwargs):
        function = getattr(self, key)
        coef = kwargs.pop("coefficients", {})
        const = kwargs.pop("constants", {})
        function.append(UnformattedBlock(self.codeCoefficient(code, coef, const)))
        setattr(self, key, function)


# splitUFLForm
# ------------

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



# splitUFL2
# ---------

def splitUFL2(u,du,d2u,tree):
    tree0 = ExprTensor(u.ufl_shape)
    tree1 = ExprTensor(u.ufl_shape)
    tree2 = ExprTensor(u.ufl_shape)

    for index in tree.keys():
        q = splitMultiLinearExpr(tree[index], [u])
        if q==0: continue
        for op in q:
            if not isinstance(op, tuple) or (len(op) != 1):
                raise Exception('Missing trial function in bulk integral')
            if op[0] == u:
                tree0[index] += sum(i[0]*i[1] for i in zip(q[op].as_ufl(),u))
            elif op[0] == du:
                for r in range(du.ufl_shape[0]):
                    for d in range(du.ufl_shape[1]):
                        tree1[index] += q[op].as_ufl()[r,d]*du[r,d]
            elif op[0] == d2u:
                for r in range(d2u.ufl_shape[0]):
                    for d1 in range(d2u.ufl_shape[1]):
                        for d2 in range(d2u.ufl_shape[2]):
                            tree2[index] += q[op].as_ufl()[r,d1,d2]*d2u[r,d1,d2]
            else:
                raise Exception('Invalid trial function derivative encountered in bulk integral: ' + str(op[0]))
    return tree0, tree1, tree2



# generateCode
# ------------

def generateCode(predefined, tensor, coefficients, tempVars=True):
    from dune.ufl import codegen
    keys = tensor.keys()
    expressions = [tensor[i] for i in keys]
    preamble, results = codegen.generateCode(predefined, expressions, coefficients, tempVars=tempVars)
    result = Variable('auto', 'result')
    return preamble + [assign(result[i], r) for i, r in zip(keys, results)]


#def compileUFL(equation, dirichlet = {}, exact = None, tempVars = True):
def compileUFL(equation, *args, **kwargs):
    form = equation.lhs - equation.rhs
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
        b = ufl.replace(form, {u: ufl.as_vector(kwargs["exact"])} )
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

    expandform = expand_indices(expand_derivatives(expand_compounds(equation.lhs)))
    if expandform == adjoint(expandform):
        model.symmetric = 'true'
    model.field = field

    dirichletBCs = [arg for arg in args if isinstance(arg, DirichletBC)]
    if "dirichlet" in kwargs:
        dirichletBCs += [DirichletBC(u.ufl_space(), ufl.as_vector(value), bndId) for bndId, value in kwargs["dirichlet"].items()]

    coefficients = set(form.coefficients())
    for bc in dirichletBCs:
        _, c = extract_arguments_and_coefficients(bc.value)
        coefficients |= set(c)

    idxConst = 0
    idxCoeff = 0
    for coefficient in coefficients:
        if coefficient.is_cellwise_constant():
            field = None  # must be improved for 'complex'
            idx = idxConst
            dimRange = 1 if coefficient.ufl_shape==() else coefficient.ufl_shape[0]
            idxConst += 1
        else:
            field = coefficient.ufl_function_space().field()
            dimRange = coefficient.ufl_shape[0]
            idx = idxCoeff
            idxCoeff += 1
        try:
            name = coefficient.str()
        except:
            name = str(coefficient)
        model.coefficients.append({ \
            'name' : name, \
            'number' : idx, \
            'counter' : coefficient.count(), \
            'dimRange' : dimRange,\
            'constant' : coefficient.is_cellwise_constant(),\
            'field': field } )

    tempVars = kwargs.get("tempVars", True)

    predefined = {u: model.arg_u, du: model.arg_du}
    predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
    model.source = generateCode(predefined, source, model.coefficients, tempVars)
    model.diffusiveFlux = generateCode(predefined, diffusiveFlux, model.coefficients, tempVars=tempVars)
    predefined.update({ubar: model.arg_ubar, dubar: model.arg_dubar})
    model.linSource = generateCode(predefined, linSource, model.coefficients, tempVars=tempVars)
    model.linDiffusiveFlux = generateCode(predefined, linDiffusiveFlux, model.coefficients, tempVars=tempVars)

    # model.linNVSource = generateCode({u: arg, du: darg, d2u: d2arg, ubar: argbar, dubar: dargbar, d2ubar: d2argbar}, linNVSource, model.coefficients, tempVars)

    predefined = {u: model.arg_u}
    predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
    model.alpha = generateCode(predefined, boundarySource, model.coefficients, tempVars)
    predefined.update({ubar: model.arg_ubar})
    model.linAlpha = generateCode(predefined, linBoundarySource, model.coefficients, tempVars)

    predefined = {u: model.arg_u, du: model.arg_du, d2u: model.arg_d2u}
    predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
    model.fluxDivergence = generateCode(predefined, fluxDivergence, model.coefficients, tempVars=tempVars)

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
        for i, value in bySubDomain.items():
            switch.append(i, generateCode(predefined, value, model.coefficients, tempVars=tempVars))
        model.dirichlet = [switch]

    return model


#def generateModel(grid, model, dirichlet = {}, exact = None, tempVars = True, header = False):
def generateModel(grid, model, *args, **kwargs):
    assert True
    from dune.common.hashit import hashIt
    start_time = timeit.default_timer()

    if isinstance(model, Equation):
        model = compileUFL(model, *args, **kwargs)

    # if not isinstance(grid, types.ModuleType):
    #     grid = grid._module
    name = 'ellipticmodel_' + model.signature + "_" + hashIt(grid._typeName)

    code = [Include(config.h)]

    code += [Include(i) for i in grid._includes]
    code.append(Include("dune/fem/misc/boundaryidprovider.hh>"))

    code.append(Include("dune/corepy/pybind11/pybind11.h"))
    code.append(Include("dune/corepy/pybind11/extensions.h"))
    code.append(Include("dune/fempy/py/grid/gridpart.hh"))
    if model.coefficients:
        code.append(Include("dune/fempy/function/virtualizedgridfunction.hh"))
    code.append(Include("dune/fem/schemes/diffusionmodel.hh"))

    nameSpace = NameSpace("ModelImpl_" + model.signature)
    nameSpace.append(model.code())
    code.append(nameSpace)

    code += [TypeAlias("GridPart", "typename Dune::FemPy::GridPart< " + grid._typeName + " >")]

    rangeTypes = ["Dune::FieldVector< " + SourceWriter.cpp_fields(c['field']) + ", " + str(c['dimRange']) + " >" for c in model.coefficients if not c['constant']]
    coefficients = ["Dune::FemPy::VirtualizedLocalFunction< GridPart, " + r + " >" for r in rangeTypes]
    code += [TypeAlias("Model", nameSpace.name + "::Model< " + ", ".join(["GridPart"] + coefficients) + " >")]

    code += [TypeAlias("ModelWrapper", "DiffusionModelWrapper< Model >"),
             TypeAlias("ModelBase", "typename ModelWrapper::Base")]

    writer = SourceWriter()
    writer.emit(code)

    if model.constants:
        model.exportSetConstant(writer)

    writer.openPythonModule(name)
    writer.emit('// export abstract base class')
    writer.emit('if( !pybind11::already_registered< ModelBase >() )')
    writer.emit('  pybind11::class_< ModelBase >( module, "ModelBase" );')
    writer.emit('')
    writer.emit('// actual wrapper class for model derived from abstract base')
    writer.emit('pybind11::class_< ModelWrapper > cls( module, "Model", pybind11::base< ModelBase >() );')
    writer.emit('cls.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return ' + str(model.dimRange) + '; } );')

    writer.emit('')
    model.export(writer, 'Model', 'ModelWrapper')
    writer.closePythonModule(name)

    if "header" in kwargs:
        with open(kwargs["header"], 'w') as modelFile:
            modelFile.write(writer.writer.getvalue())
    return writer, name


#def importModel(grid, model, dirichlet = {}, exact = None, tempVars = True, header = False):
def importModel(grid, model, *args, **kwargs):
    if isinstance(model, str):
        with open(model, 'r') as modelFile:
            data = modelFile.read()
        name = data.split('PYBIND11_PLUGIN( ')[1].split(' )')[0]
        builder.load(name, data, "ellipticModel")
        return importlib.import_module("dune.generated." + name)
    writer, name = generateModel(grid, model, *args, **kwargs)
    module = builder.load(name, writer.writer.getvalue(), "ellipticModel")
    writer.close()
    return module
