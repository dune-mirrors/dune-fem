from __future__ import division, print_function, unicode_literals

from ufl import Form, Coefficient, TrialFunction, TestFunction, replace
from ufl.core.expr import Expr
from ufl.algorithms.analysis import extract_arguments_and_coefficients
from ufl.equation import Equation

from dune.common.utility import isString

from dune.source.cplusplus import Include, NameSpace, TypeAlias
from dune.source.cplusplus import SourceWriter
from dune.source.fem import fieldTensorType
from dune.ufl import GridFunction, DirichletBC
from dune.ufl.gatherderivatives import gatherDerivatives
from dune.ufl.codegen import uflSignature

from .ufl import _compileUFL
from .model import Integrands


def init(integrands, source, *args, **kwargs):
    coefficients = kwargs.pop('coefficients', dict())
    coefficientNames = integrands._coefficientNames
    if len(args) == 1 and isinstance(args[0], dict):
        coefficients.update(args[0])
        args = []
    else:
        args = list(a for a in args if not isinstance(a,DirichletBC) and not a is None)
        dirichletBCs = list(a for a in args if isinstance(a,DirichletBC))
        if len(args) > len(coefficientNames) + len(dirichletBCs):
            raise ValueError('Too many coefficients passed.')

    args += [None] * (len(coefficientNames) - len(args))

    for name, value in kwargs.items():
        try:
            i = coefficientNames[name]
        except KeyError:
            raise ValueError('No such coefficent: ' + name + '.')

        if args[i] is not None:
            raise ValueError('Coefficient already given as positional argument: ' + name + '.')
        args[i] = value

    for key, value in coefficients.items():
        if isinstance(key, Coefficient):
            try:
                i = integrands._renumbering[key]
            except AttributeError:
                raise ValueError('Cannot map UFL coefficients, because model was not generated from UFL form.')
            except KeyError:
                raise ValueError('No such coefficient: ' + str(key) + '.')
        elif isString(key):
            try:
                i = coefficientNames[key]
            except KeyError:
                raise ValueError('No such coefficent: ' + key + '.')
        else:
            raise ValueError('Expecting keys of coefficient map to be strings or  intances of ufl.Coefficient.')
        if args[i] is not None:
            raise ValueError('Coefficient already given as positional or keyword argument: ' + str(key) + '.')
        args[i] = value

    if hasattr(integrands, '_renumbering'):
        for c, i in integrands._renumbering.items():
            if isinstance(c, GridFunction):
                if args[i] is None:
                    args[i] = c.gf

    if any(arg is None for arg in args):
        missing = [name for name, i in coefficientNames.items() if args[i] is None]
        raise ValueError('Missing coefficients: ' + ', '.join(missing) + '.')

    integrands.base.__init__(integrands, *args, **kwargs)

    for c in source.constantList:
        if hasattr(c,"name") and hasattr(c,"value"):
            assert hasattr(integrands,c.name)
            getattr(type(integrands), c.name).fset(integrands, c.value)


def setConstant(integrands, index, value):
    try:
        index = integrands._renumbering[index]
    except KeyError:
        pass
    integrands._setConstant(index, value)


class Source(object):
    version = "v1_1"
    def __init__(self, integrands, gridType, gridIncludes, modelIncludes, form, *args,
            tempVars=True, virtualize=True):
        self.gridType = gridType
        self.gridIncludes = gridIncludes
        if modelIncludes is not None:
            self.modelIncludes = modelIncludes
        else:
            self.modelIncludes = ["dune/fempy/py/integrands.hh"]
        self.integrands = integrands
        self.tempVars = tempVars
        self.virtualize = virtualize
        self.args = args
        self.form = form

    def signature(self):
        return uflSignature(self.form,
                *self.integrands._coefficients,
                *self.integrands._constantNames,
                *[a for a in self.args if isinstance(a,DirichletBC)],
                *self.integrands.baseSignature
                )+Source.version

    def name(self):
        from dune.common.hashit import hashIt
        if self.virtualize:
            return self.integrands.baseName + '_' + self.signature() + '_' + hashIt(self.gridType)
        else:
            return self.integrands.baseName + '_nv_' + self.signature() + '_' + hashIt(self.gridType)

    def valueTuples(self):
        if isinstance(self.form, Form):
            derivatives = gatherDerivatives(self.form)
            return ['std::tuple< ' + ', '.join(fieldTensorType(v.ufl_shape) for v in d) + ' >' for d in derivatives]
        else:
            return [self.form.rangeValueTuple, self.form.domainValueTuple]

    # actual code generation (the generator converts this object to a string)
    def __str__(self):
        if isinstance(self.form, Form):
            # actual code generation
            integrands = _compileUFL(self.integrands, self.form, *self.args, tempVars=self.tempVars)
        else:
            integrands = self.integrands

        code  = [Include('config.h')]
        code += [Include(i) for i in self.gridIncludes]

        code += integrands.includes()
        code.append(Include("dune/python/pybind11/pybind11.h"))
        code.append(Include("dune/python/pybind11/extensions.h"))
        code.append(Include("dune/fempy/py/grid/gridpart.hh"))

        if integrands._coefficients:
            if self.virtualize:
                code.append(Include("dune/fempy/function/virtualizedgridfunction.hh"))
                code.append(Include('dune/fempy/function/simplegridfunction.hh'))
                code.append(Include('dune/fem/misc/gridfunctionview.hh'))
            else:
                for c in integrands._coefficients:
                    for i in c._includes:
                        code.append(Include(i))
        for i in self.modelIncludes:
            code.append(Include(i))

        nameSpace = NameSpace('Integrands_' + self.signature())
        # add integrands class
        nameSpace.append(integrands.code(self.name(),integrands.targs))
        code.append(nameSpace)

        writer = SourceWriter()
        writer.emit("#ifndef GuardIntegrands_" + self.signature())
        writer.emit("#define GuardIntegrands_" + self.signature())
        writer.emit(code)

        name = self.name()
        writer.openPythonModule(name)

        coefficients = integrands.coefficientCppTypes
        integrandsName = nameSpace.name + '::Integrands< ' + ', '.join(['GridPart'] + coefficients) + ' >'
        writer.emit(TypeAlias('GridPart', 'typename Dune::FemPy::GridPart< ' + self.gridType + ' >'))
        writer.emit(TypeAlias('Integrands', integrandsName))
        writer.emit('auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("'+integrandsName+'"), Dune::Python::IncludeFiles({"python/dune/generated/'+name+'.cc"})).first;')
        writer.emit('Dune::FemPy::registerIntegrands< Integrands >( module, cls );')
        if coefficients:
            coefficientNames = integrands.coefficientNames
            initArgs = ', '.join('const ' + t + ' &' + n for t, n in zip(coefficients, coefficientNames))
            keepAlive = ', '.join('pybind11::keep_alive< 1, ' + str(i+2) + ' >()' for i in range(len(coefficientNames)))
            writer.emit('cls.def( pybind11::init( [] ( ' + initArgs + ' ) { return new Integrands( ' + ', '.join(coefficientNames) +  ' ); } ), ' + keepAlive + ' );')
        else:
            writer.emit('cls.def( pybind11::init( [] () { return new Integrands(); } ) );')

        for t, n, ns in zip(integrands.constantTypes, integrands.constantNames, integrands.constantShortNames):
            te = "Integrands::" + t
            writer.emit('cls.def_property( "' + ns + '", [] ( Integrands &self ) -> ' + te + ' { return self.' + n + '(); }, [] ( Integrands &self, const ' + te + ' &v ) { self.' + n + '() = v; } );')
        writer.emit('cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return '+str(self.virtualize).lower()+';});')
        hasDirichletBC = 'true' if integrands.hasDirichletBoundary else 'false'
        writer.emit('cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return '+hasDirichletBC+';});')

        writer.closePythonModule(name)
        writer.emit("#endif // GuardIntegrands_" + self.signature())

        source = writer.writer.getvalue()
        writer.close()
        return source


# Load the actual module - the code generation from the ufl form is done
# when the 'Source' class is converted to a string i.e. in Source.__str__
def load(grid, form, *args, renumbering=None, tempVars=True,
        virtualize=True, modelPatch=[None,None],
        includes=None):

    if not isinstance(modelPatch,list) and not isinstance(modelPatch,tuple):
        modelPatch = [modelPatch,None]

    if isinstance(form, Equation):
        form = form.lhs - form.rhs

    if isinstance(form, Integrands):
        integrands = form
    else:
        if len(form.arguments()) < 2:
            raise ValueError("Integrands model requires form with at least two arguments.")

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

        if not isinstance(form, Form):
            raise ValueError("ufl.Form or ufl.Equation expected.")

        _, coeff_ = extract_arguments_and_coefficients(form)
        coeff_ = set(coeff_)

        # added for dirichlet treatment same as elliptic model
        dirichletBCs = [arg for arg in args if isinstance(arg, DirichletBC)]
        # remove the dirichletBCs
        arg = [arg for arg in args if not isinstance(arg, DirichletBC)]
        for dBC in dirichletBCs:
            _, coeff__ = extract_arguments_and_coefficients(dBC.ufl_value)
            coeff_ |= set(coeff__)
        coeff = {c : c.toVectorCoefficient()[0] for c in coeff_ if len(c.ufl_shape) == 0 and not c.is_cellwise_constant()}

        form = replace(form,coeff)
        uflExpr = [form]
        for dBC in dirichletBCs:
            arg.append(dBC.replace(coeff))
            uflExpr += [dBC.ufl_value] # arg[-1].ufl_value]

        if modelPatch[1] is not None:
            uflExpr += modelPatch[1]

        derivatives = gatherDerivatives(form, [phi, u])

        derivatives_phi = derivatives[0]
        derivatives_u = derivatives[1]

        integrands = Integrands(u,
                                (d.ufl_shape for d in derivatives_u), (d.ufl_shape for d in derivatives_phi),
                                uflExpr,virtualize)

    if modelPatch[0] is not None:
        modelPatch[0](integrands)

    # set up the source class
    source = Source(integrands, grid._typeName, grid._includes, includes, form, *args,
             tempVars=tempVars,virtualize=virtualize)

    # ufl coefficient and constants only have numbers which depend on the
    # order in whch they were generated - we need to keep track of how
    # these numbers are translated into the tuple numbering in the
    # generated C++ code
    if isinstance(form, Form):
        coefficients = set(integrands.coefficientList+integrands.constantList)
        numCoefficients = len(coefficients)
        if renumbering is None:
            renumbering = dict()
            renumbering.update((c, i) for i, c in enumerate(sorted((c for c in coefficients if not c.is_cellwise_constant()), key=lambda c: c.count())))
            renumbering.update((c, i) for i, c in enumerate(c for c in coefficients if c.is_cellwise_constant()))
        coefficientNames = integrands._coefficientNames # ['coefficient' + str(i) if n is None else n for i, n in enumerate(getattr(c, 'name', None) for c in coefficients if not c.is_cellwise_constant())]
    else:
        coefficientNames = form.coefficientNames

    # call code generator
    from dune.generator import builder
    module = builder.load(source.name(), source, "Integrands")

    rangeValueTuple, domainValueTuple = source.valueTuples()
    setattr(module.Integrands, "_domainValueType", domainValueTuple)
    setattr(module.Integrands, "_rangeValueType", rangeValueTuple)
    # redirect the __init__ method to take care of setting coefficient and renumbering
    class Model(module.Integrands):
        def __init__(self, *args, **kwargs):
            self.base = module.Integrands
            init(self,integrands,*args,**kwargs)
            for c in integrands.constantList:
                c.registerModel(self)

    setattr(Model, '_coefficientNames', {n: i for i, n in enumerate(coefficientNames)})
    if renumbering is not None:
        setattr(Model, '_renumbering', renumbering)
        Model._setConstant = module.Integrands.__dict__['setConstant']
        setattr(Model, 'setConstant', setConstant)

    return Model
