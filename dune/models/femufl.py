r"""@package femufl
Provides the translation between a UFL interface file and a DUNE model file.

generate and generateFromExact are the main functions for creating the model. They essentially do the following:
1. Iterate over the expression trees for a, L, and the derivative F (computed using UFL).
2. 'Translate' the individual terms to strings compatible with DUNE.
3. Put the forms back together and sort them into variables (e.g. 'source' and 'flux'). This involves simplifying
the expressions using sympy.
4. Put these strings into the corresponding functions in the model.hh.in file.

In the case of the exact version, L is defined using the inputted exact expression.

The following shows an example usage for constructing a model class for
@f[ \int uv + \nabla u\cdot\nabla v @f]
Other such examples are contained within the src directory.

Examples:
    >>> from dune.fem.models.femufl import *
    >>> model    = DuneUFLModel(2,2,'LaplaceMass')

    >>> c0 = sympy.cos(2*math.pi*model.x0)
    >>> c1 = sympy.cos(2*math.pi*model.x1)
    >>> exact = [c0*c1,c0*c1]

    >>> u = model.trialFunction()
    >>> v = model.testFunction()
    >>> x = model.spatialCoordinate()
    >>> dx0 = dx(0)
    >>> a = (inner(u,v) + inner(grad(u),grad(v)))*dx0
    >>> model.generateFromExact(a,exact)
    model.hh has been changed
"""

########################

# Following code only compatible with UFL 1.6.0 (or later)
# http://fenicsproject.org/documentation/ufl/1.0-beta/user/form_language.html

from __future__ import print_function
import importlib
import fileinput
import math
import os
import subprocess
import sys
import timeit
from types import ModuleType
import sympy
from sympy.printing.ccode import CCodePrinter as SympyCCodePrinter

from mpi4py import MPI
from .. import femmpi

# Remark: should get rid of the * is apparently bad practice
from ufl import *
from ufl.algorithms import expand_compounds, expand_derivatives, expand_indices, ufl2dot
from ufl.algorithms.apply_derivatives import apply_derivatives

compilePath = os.path.join(os.path.dirname(__file__), "../generated")

# EXTENDING sympy.ccode
# make the sympy code generator use [i][j] notation instead of [i+n*j] as is used by default
class DuneUFLModel:
    """Model class that contains all of the following code, and is typically defined at the start of a UFL file.
    """
    def __init__(self, dimDomain, dimRange, modelName):
        """Initialise the model class.

        Args:
            dimDomain (int) : Dimension of the domain (i.e. the dimension of the spatial coordinate).
            dimRange (int) : Dimension of the range (i.e. the dimension of functions).
            modelName (str) : Name of the model (e.g. LaplaceMass, Transport).
        """
        self.dimD = dimDomain
        self.dimR = dimRange
        self.modelName = modelName
        # variables to later store the parsed ufl expression
        self.massCheck = 1
        self.boundCheck = 0
        self.diricCheck = 0
        self.massCoef = ''
        self.fluxCoef = ''
        self.alphaCoef = ''
        self.diricCoef = ''
        self.linSource = ''
        self.linFlux = ''
        self.linAlpha = ''
        self.source = ''
        self.flux = ''
        self.rhs = ''
        self.alpha = ''
        self.diric = ''
        self.rhsBound = ''
        self.modelTemplate = ''
        self.setCoef = ''
        self.initCoef = ''
        self.pyTemplate = ''
        self.pySetCoef = ''
        self.pyRangeType = ''
        # ... and the sympy equivalents
        self.matCodeSrc = sympy.zeros(self.dimR, 1)
        self.matCodeFlux = sympy.zeros(self.dimR, self.dimD)

        # ufl variables needed to define the forms
        self.cell = triangle # later add here something that depends on dimDomain
        self.element = VectorElement("Lagrange", self.cell, 1, self.dimR)
        self.vector = VectorElement("Lagrange", self.cell, 1, self.dimD)
        self.u_ = Coefficient(self.element)

        # sympy symbols needed for translating to sympy expressions
        self.u0 = sympy.IndexedBase('uBar')
        self.du0 = sympy.IndexedBase('gradientBar')
        self.v0 = sympy.IndexedBase('value')
        self.dv0 = sympy.IndexedBase('gradient')
        self.v1 = sympy.IndexedBase('v1')
        self.dv1 = sympy.IndexedBase('dv1')
        self.flux = sympy.IndexedBase('flux')
        self.x0, self.x1, self.x2 = sympy.symbols('x0 x1 x2')
        self.U = [sympy.Function('U'+str(i))(self.x0, self.x1, self.x2) for i in range(0, self.dimR)]
        self.residual = sympy.zeros(self.dimR, 1)
        self.coefficients = {}
        self.unsetCoefficients = {}
        self.path = os.path.join(os.path.dirname(__file__), "../generated")

    ########################
    # sympy helpers
    ########################
    class CCodePrinter(SympyCCodePrinter):
        """Modification of the sympy c-code printer to allow for [i][j] instead of [i,j].
        """
        def _print_Indexed(self, expr):
            output = self._print(expr.base.label) + ''.join(['[' + self._print(x) + ']' for x in expr.indices])
            return output
    def ccode(self, expr, assignTo=None, **settings):
        return self.CCodePrinter(settings).doprint(expr, assignTo)

    def sympyToString(self, sympyExpr):
        """Use this for function output of sympy expressions (should use # IndexedBase instead of x1,x2).
        """
        return str(sympyExpr).replace("x0", "x[0]").replace("x1", "x[1]")

    def stringTocFlux(self, expr):
        """Use sympy to generate code for a given flux expression (needs some common subexpression evaluation).
        """
        codeStr = ''
        self.matCodeFlux = sympy.zeros(self.dimR, self.dimD)
        for i in range(0, self.dimR):
            for d in range(0, self.dimD):
                exprSubs = str(expr)
                for name in self.coefficients:
                    coeff = self.coefficients[name]
                    if self.coefficients[name] != "":
                        for ii in range(0, len(coeff)):
                            exprSubs = exprSubs.replace(
                                name+"["+str(ii)+"]", "("+str(coeff[ii])+")")
                exprSubs = exprSubs.replace("value", "self.v0").replace("uBar", "self.u0").\
                    replace("gradient", "self.dv0").replace("gradBar", "self.du0").\
                    replace("gradphi", "self.dv1").replace("phi", "self.v1").\
                    replace("x0", "self.x0").replace("x1", "self.x1").replace("x2", "self.x2").\
                    replace("x[0]", "self.x0").replace("x[1]", "self.x1").replace("x[2]", "self.x2").\
                    replace("cos", "sympy.cos").replace("sin", "sympy.sin"). \
                    replace("sqrt", "sympy.sqrt")  # this is stupid and has to be improved!
                for ii in range(0, self.dimR):
                    if i == ii:
                        exprSubs = str(exprSubs).replace("self.v1["+str(ii)+"]", "1")
                    else:
                        exprSubs = str(exprSubs).replace("self.v1["+str(ii)+"]", "0")
                    for dd in range(0, self.dimD):
                        if d == dd and i == ii:
                            exprSubs = str(exprSubs).replace("self.dv1["+str(ii)+", "+str(dd)+"]", "1")
                        else:
                            exprSubs = str(exprSubs).replace("self.dv1["+str(ii)+", "+str(dd)+"]", "0")
                sympyExpr = eval(exprSubs)
                if codeStr != '':
                    codeStr = codeStr + '\n'
                codeStr = codeStr + '      flux[' + str(i) + '][' + str(d) + '] = ' + self.ccode(sympyExpr) + ';'
                self.matCodeFlux[i, d] = sympyExpr
        return self.sympyToString(codeStr)

    def stringToCSource(self, expr):
        """Use sympy to generate code for a given source expression.
        """
        codeStr = ''
        self.matCodeSrc = sympy.zeros(self.dimR, 1)
        for i in range(0, self.dimR):
            exprSubs = str(expr)
            for name in self.coefficients:
                coeff = self.coefficients[name]
                if self.coefficients[name] != "":
                    for ii in range(0, len(coeff)):
                        exprSubs = exprSubs.replace(name+"["+str(ii)+"]", "("+str(coeff[ii])+")")
            exprSubs = exprSubs.replace("value", "self.v0").replace("uBar", "self.u0").\
                    replace("gradient", "self.dv0").replace("gradBar", "self.du0").\
                    replace("gradphi", "self.dv1").replace("phi", "self.v1").\
                    replace("x0", "self.x0").replace("x1", "self.x1").replace("x2", "self.x2").\
                    replace("x[0]", "self.x0").replace("x[1]", "self.x1").replace("x[2]", "self.x2").\
                    replace("cos", "sympy.cos").replace("sin", "sympy.sin"). \
                    replace("sqrt", "sympy.sqrt")  # this is stupid and has to be improved!
            for ii in range(0, self.dimR):
                if i == ii:
                    exprSubs = str(exprSubs).replace("self.v1["+str(ii)+"]", "1")
                else:
                    exprSubs = str(exprSubs).replace("self.v1["+str(ii)+"]", "0")
            sympyExpr = eval(exprSubs)
            if codeStr != '':
                codeStr = codeStr + '\n'
            codeStr = codeStr + '      flux[' + str(i) + '] = ' + self.ccode(sympyExpr) + ';'
            self.matCodeSrc[i] = sympyExpr
        return self.sympyToString(codeStr)

    ################
    # UFL Helpers
    ################
    class NamedCoefficient(Coefficient):
        """Extend UFL "Coefficient" class to store a name.
        """
        def __init__(self, function_space, dimR, name="Coef"):
            Coefficient.__init__(self, function_space)
            self.name_ = name
            self.dimR = dimR
            #self.name_ = sympy.IndexedBase(name)
            command1 = 'global ' + name + '\n' + name + " = sympy.IndexedBase('" + name + "')"
            exec(command1)
        def name(self):
            return self.name_
        def dimRange(self):
            return self.dimR

    def expandExpr(self, a):
        """Expand expression in terms of in terms of indices.
        """
        ac = expand_compounds(a)
        acd = expand_derivatives(ac)
        acdi = expand_indices(acd)
        return acdi

    class Stack():
        """Define stack class.
        """
        def __init__(self):
            self.items = []
        def isEmpty(self):
            return self.items == []
        def push(self, item):
            return self.items.append(item)
        def pop(self):
            return self.items.pop()
        def peek(self):
            return self.items[len(self.items)-1]
        def size(self):
            return len(self.items)

    def tree(self, expr):
        """Take UFL expression tree and sort terms into source, flux etc.
        """
        if str(expr._ufl_class_.__name__) == "Sum":
            for o in expr.ufl_operands:
                self.tree(o)
        else:
            self.massCheck = 1
            stack = self.Stack()
            self.check(expr)
            self.stack_maker(expr, stack)
            if self.diricCheck == 1:
                if self.diricCoef:
                    self.diricCoef += " + "
                while stack.isEmpty() == 0:
                    self.diricCoef += stack.peek()
                    stack.pop()
            elif self.boundCheck == 1:
                if self.alphaCoef:
                    self.alphaCoef += " + "
                while stack.isEmpty() == 0:
                    self.alphaCoef += stack.peek()
                    stack.pop()
            else:
                if self.massCheck == 0:
                    if self.fluxCoef:
                        self.fluxCoef += " + "
                    while stack.isEmpty() == 0:
                        self.fluxCoef += stack.peek()
                        stack.pop()
                else:
                    if self.massCoef:
                        self.massCoef += " + "
                    while stack.isEmpty() == 0:
                        self.massCoef += stack.peek()
                        stack.pop()

    def stack_maker(self, expr, stack):
        """Iterate over expression tree and translate UFL variables to strings.
        """
        if str(expr._ufl_class_.__name__) == "Sin":
            arguments = ""
            sinStack = self.Stack()
            for op in expr.ufl_operands:
                self.stack_maker(op, sinStack)
            while sinStack.isEmpty() == False:
                arguments += sinStack.peek()
                sinStack.pop()
            string = "sin(" + arguments + ")"
            stack.push(string)
        elif str(expr._ufl_class_.__name__) == "Cos":
            arguments = ""
            cosStack = self.Stack()
            for op in expr.ufl_operands:
                self.stack_maker(op, cosStack)
            while cosStack.isEmpty() == False:
                arguments += cosStack.peek()
                cosStack.pop()
            string = "cos(" + arguments + ")"
            stack.push(string)
        elif str(expr._ufl_class_.__name__) == "Sqrt":
            arguments = ""
            sqrtStack = self.Stack()
            for op in expr.ufl_operands:
                self.stack_maker(op, sqrtStack)
            while sqrtStack.isEmpty() == False:
                arguments += sqrtStack.peek()
                sqrtStack.pop()
            string = "sqrt(" + arguments + ")"
            stack.push(string)
        elif str(expr) == "v_0":
            stack.push("phi")
        elif str(expr) == "w_0":
            stack.push("uBar")
        elif expr._ufl_class_.__name__ == "Coefficient":
            if not expr.name() in self.coefficients:
                stack.push(expr.name())
                if not expr.name() in self.unsetCoefficients:
                    self.unsetCoefficients[expr.name()] = expr.dimR
            else:
                stack.push(str(self.coefficients[expr.name()]))
        elif str(expr) == "v_1":
            stack.push("value")
        elif str(expr) == "x":
            stack.push('x')
        elif str(expr) == "grad(v_0)":
            stack.push("gradphi")
        elif str(expr) == "grad(v_1)":
            stack.push("gradient")
        elif str(expr) == "grad(w_0)":
            stack.push("gradBar")
        elif str(expr._ufl_class_.__name__) == "IntValue":
            stack.push(str(expr))
        elif str(expr._ufl_class_.__name__) == "FloatValue":
            stack.push(str(expr))
        elif str(expr._ufl_class_.__name__) == "Product":
            stack.push(")")
            self.stack_maker(expr.ufl_operands[0], stack)
            stack.push("*")
            self.stack_maker(expr.ufl_operands[1], stack)
            stack.push("(")
        elif str(expr._ufl_class_.__name__) == "Power":
            stack.push(")")
            self.stack_maker(expr.ufl_operands[1], stack)
            stack.push(",")
            self.stack_maker(expr.ufl_operands[0], stack)
            stack.push("pow(")
        elif str(expr._ufl_class_.__name__) == "Sum":
            stack.push(")")
            self.stack_maker(expr.ufl_operands[0], stack)
            stack.push("+")
            self.stack_maker(expr.ufl_operands[1], stack)
            stack.push("(")
        elif str(expr._ufl_class_.__name__) == "Division":
            stack.push(")")
            self.stack_maker(expr.ufl_operands[1], stack)
            stack.push("/")
            self.stack_maker(expr.ufl_operands[0], stack)
            stack.push("(")
        elif str(expr._ufl_class_.__name__) == "Dot":
            stack.push(")")
            self.stack_maker(expr.ufl_operands[1], stack)
            stack.push("*")
            self.stack_maker(expr.ufl_operands[0], stack)
            stack.push("(")
        elif  str(expr._ufl_class_.__name__) == "Indexed":
            stack.push("]")
            self.stack_maker(expr.ufl_operands[1], stack)
            stack.push("[")
            self.stack_maker(expr.ufl_operands[0], stack)
        elif  str(expr._ufl_class_.__name__) == "MultiIndex":
            stack.push(str(expr))
        else:
            for op in expr.ufl_operands:
                self.stack_maker(op, stack)

    def check(self, expr):
        """Check whether expression has a grad(v) term.
        """
        for op in expr.ufl_operands:
            self.check(op)
        if str(expr) == "grad(v_0)":
            self.massCheck = 0

    def formOutput(self, a):
        """Take a UFL form, puts it into tree() as a UFL expression.
        """
        acdi = self.expandExpr(a)
        self.massCoef = ""
        self.fluxCoef = ""
        self.alphaCoef = ""
        if isinstance(acdi, Form):
            integrals = acdi.integrals()
            for intg in integrals:
                if intg.integral_type() == 'exterior_facet':
                    self.boundCheck = 1
                    #if intg.subdomain_id() == x:
                else:
                    self.boundCheck = 0
                self.tree(intg.integrand())
        else: self.tree(acdi)

    def diricOutput(self, args):
        """Stores Dirichlet boundary data (from generate) in diric.
        """
        self.diricCheck = 1
        for arg in args:
            if hasattr(args, '_ufl_class_'):
                self.tree(arg)
            else:
                if self.diricCoef:
                    self.diricCoef += " + "
                self.diricCoef += str(arg)
        self.diricCheck = 0
        self.storeDiric()

    def modelPrint(self, exact):
        """Write expression strings into model file.
        """
        inputfile = self.path + '/model.hh.in'
        outputfile = self.path + '/' + self.modelName + "Model.hh"
        with open(outputfile, "wt") as fout:

            with open(inputfile, "rt") as fin:
                for line in fin:
                    if '#SOURCE' in line:
                        fout.write(self.source)
                    elif '#DIFFUSIVEFLUX' in line:
                        fout.write(self.flux)
                    elif '#LINSOURCE' in line:
                        fout.write(self.linSource)
                    elif '#LINDIFFUSIVEFLUX' in line:
                        fout.write(self.linFlux)
                    elif '#ALPHA' in line:
                        fout.write(self.alpha)
                    elif '#LINALPHA' in line:
                        fout.write(self.linAlpha)
                    elif '#HASDIRICHLET' in line:
                        if self.diric:
                            fout.write('      return true;\n')
                        else:
                            fout.write('      return false;\n')
                    elif '#HASNEUMAN' in line:
                        fout.write('      return true;\n')
                    elif '#ISDIRICHLETINTERSECTION' in line:
                        fout.write('      return false;\n')
                    elif '#FORCING' in line:
                        fout.write(self.rhs)
                    elif '#DIRICHLETDATA' in line:
                        if self.diric:
                            fout.write(self.diric)
                        else:
                            for i in range(0, self.dimR):
                                fout.write('      value['+str(i)+'] = ' \
                                           + self.sympyToString(self.ccode(exact[i])) + ';\n')
                    elif '#NEUMANDATA' in line:
                        fout.write(self.rhsBound)
                    elif "#JACOBIANEXACT" in line:
                        for i in range(0, self.dimR):
                            fout.write('      value['+str(i)+'][0] = ' \
                                       + self.sympyToString(self.ccode(sympy.diff(exact[i], self.x0))) + ';\n')
                            if self.dimD > 1:
                                fout.write('      value['+str(i)+'][1] = ' \
                                           + self.sympyToString(self.ccode(sympy.diff(exact[i], self.x1))) + ';\n')
                            if self.dimD > 2:
                                fout.write('      value['+str(i)+'][2] = ' \
                                           + self.sympyToString(self.ccode(sympy.diff(exact[i], self.x2))) + ';\n')
                    elif '#RESIDUAL' in line:
                        value = sympy.IndexedBase('value')
                        jacobian = sympy.IndexedBase('jacobian')
                        hessian = sympy.IndexedBase('hessian')
                        fout.write('      DomainType x = entity_->geometry()'
                                   + '.global( Dune::Fem::coordinate(point) );\n')
                        for i in range(0, self.dimR):
                            sympyResidual = self.residual[i]
                            for r in range(0, self.dimR):
                                sympyResidual = sympyResidual.subs(sympy.diff(self.U[r], self.x0, self.x0), hessian[r, 0, 0])
                                sympyResidual = sympyResidual.subs(sympy.diff(self.U[r], self.x0, self.x1), hessian[r, 0, 1])
                                sympyResidual = sympyResidual.subs(sympy.diff(self.U[r], self.x1, self.x1), hessian[r, 1, 1])
                                sympyResidual = sympyResidual.subs(sympy.diff(self.U[r], self.x1, self.x2), hessian[r, 1, 2])
                                sympyResidual = sympyResidual.subs(sympy.diff(self.U[r], self.x2, self.x2), hessian[r, 2, 2])
                                sympyResidual = sympyResidual.subs(sympy.diff(self.U[r], self.x0), jacobian[r, 0])
                                sympyResidual = sympyResidual.subs(sympy.diff(self.U[r], self.x1), jacobian[r, 1])
                                sympyResidual = sympyResidual.subs(sympy.diff(self.U[r], self.x2), jacobian[r, 2])
                                sympyResidual = sympyResidual.subs(self.U[r], value[r])
                            result = self.ccode(sympyResidual.simplify())
                            fout.write('      result['+str(i)+'] = ' + self.sympyToString(result)+';\n')
                    elif '#MODELTEMPLATE' in line:
                        line = line.replace('#MODELTEMPLATE', self.modelTemplate)
                        fout.write(line)
                    elif '#INIT' in line:
                        if self.unsetCoefficients:
                            line = line.replace('#INIT', self.initCoef)
                        else:
                            line = ''
                        fout.write(line)
                    elif '#SETCOEFFICIENT' in line:
                        if self.unsetCoefficients:
                            line = line.replace('#SETCOEFFICIENT', self.setCoef)
                        else:
                            line = ''
                        fout.write(line)
                    else:
                        line = line.replace('ModelTmp', self.modelName+'Model')
                        line = line.replace('DIMRANGE', str(self.dimR))
                        fout.write(line)
        print("model.hh has been changed")

    def storeSrc(self):
        """Store expression in source, flux and alpha.
        """
        # for source function
        if not self.massCoef:
            self.source = '      flux = RangeType(0);\n'
        else:
            cMass = self.stringToCSource(self.massCoef)
            if self.xTest(cMass) == 1:
                self.source += '      DomainType x = entity_->geometry().global( Dune::Fem::coordinate(point) );\n'
            for coef in self.unsetCoefficients:
                if coef in self.massCoef:
                    self.source += '      typename ' + coef + 'DiscreteFunction::RangeType ' + coef + ';\n      ' \
                                   + coef + 'Local_->evaluate(point,' + coef + ');\n'
            self.source += cMass + '\n'
        # for diffusiveFlux function
        if not self.fluxCoef:
            self.flux = '      flux = JacobianRangeType(0);\n'
        else:
            cFlux = self.stringTocFlux(self.fluxCoef)
            if self.xTest(cFlux) == 1:
                for coef in self.unsetCoefficients:
                    if coef in self.fluxCoef:
                        self.flux = '      DomainType x = entity_->geometry().global( Dune::Fem::coordinate(point) );' \
                                    + '\n      typename ' + coef + 'DiscreteFunction::RangeType ' + coef + ';\n      ' \
                                    + coef + 'Local_->evaluate(point,' + coef + ');\n' + cFlux + '\n'
                if str(self.flux) == 'flux':
                    self.flux = '      DomainType x = entity_->geometry().global( Dune::Fem::coordinate(point) );\n' \
                                + cFlux + '\n'
            else:
                for coef in self.unsetCoefficients:
                    if coef in self.fluxCoef:
                        self.flux = '      typename ' + coef + 'DiscreteFunction::RangeType ' + coef + ';\n      ' \
                                    + coef + 'Local_->evaluate(point,' + coef + ');\n' + cFlux + '\n'
                if str(self.flux) == 'flux':
                    self.flux = cFlux + '\n'
        # for alpha function
        if not self.alphaCoef:
            self.alpha = '      DomainType x = entity_->geometry().global( Dune::Fem::coordinate(point) );\n' \
                         '      linAlpha(value,x,value,val);\n'
        else:
            cAlpha = self.stringToCSource(self.alphaCoef)
            if self.xTest(cAlpha) == 1:
                self.alpha += '      DomainType x = entity_->geometry().global( Dune::Fem::coordinate(point) );\n'
            for coef in self.unsetCoefficients:
                if coef in self.alphaCoef:
                    self.alpha += '      typename ' + coef + 'DiscreteFunction::RangeType ' + coef + ';\n      ' \
                                  + coef + 'Local_->evaluate(point,' + coef + ');\n'
            self.alpha += cAlpha + '\n'
            self.alpha = self.alpha.replace('flux', 'val')

    def storeDiric(self):
        """Store Dirichlet boundary data in diric.
        """
        if self.diricCoef:
            cDiric = self.stringToCSource(self.diricCoef)
            if self.xTest(cDiric) == 1:
                self.diric += '      DomainType x = entity_->geometry().global( Dune::Fem::coordinate(point) );\n'
            for coef in self.unsetCoefficients:
                if coef in self.diricCoef:
                    self.diric += '      typename ' + coef + 'DiscreteFunction::RangeType ' + coef + ';\n      ' \
                                  + coef + 'Local_->evaluate(point,' + coef + ');\n'
            self.diric += cDiric + '\n'
            self.diric = self.diric.replace('flux', 'value')


    def storeLin(self):
        """Store linearised expression in linSource, linFlux and linAlpha.
        """
        # for linSource function
        if not self.massCoef:
            self.linSource = '      flux = RangeType(0);\n'
        else:
            cMass = self.stringToCSource(self.massCoef)
            if self.xTest(cMass) == 1:
                self.linSource += '      DomainType x = entity_->geometry().global( Dune::Fem::coordinate(point) );\n'
            for coef in self.unsetCoefficients:
                if coef in self.massCoef:
                    self.linSource += '      typename ' + coef + 'DiscreteFunction::RangeType ' + coef + ';\n      ' \
                                      + coef + 'Local_->evaluate(point,' + coef + ');\n'
            self.linSource += cMass + '\n'
        # for linDiffusiveFlux function
        if not self.fluxCoef:
            self.linFlux = '      flux = JacobianRangeType(0);\n'
        else:
            cFlux = self.stringTocFlux(self.fluxCoef)
            if self.xTest(cFlux) == 1:
                self.linFlux += '      DomainType x = entity_->geometry().global( Dune::Fem::coordinate(point) );\n'
            for coef in self.unsetCoefficients:
                if coef in self.fluxCoef:
                    self.linFlux += '      typename ' + coef + 'DiscreteFunction::RangeType ' + coef + ';\n      ' \
                                    + coef + 'Local_->evaluate(point,' + coef + ');\n'
            self.linFlux += cFlux + '\n'
        # for linAlpha function
        if not self.alphaCoef:
            self.linAlpha = '      RangeType alpha(0);\n      for (unsigned int i=0;i<val.size();++i)\n' \
                            '       val[i] = alpha[i]*value[i];\n'
        else:
            cAlpha = self.stringToCSource(self.alphaCoef)
            if self.xTest(cAlpha) == 1:
                self.linAlpha += '      DomainType x = entity_->geometry().global( Dune::Fem::coordinate(point) );\n'
            for coef in self.unsetCoefficients:
                if coef in self.alphaCoef:
                    self.linAlpha += '      typename ' + coef + 'DiscreteFunction::RangeType ' + coef + ';\n      ' \
                                     + coef + 'Local_->evaluate(point,' + coef + ');\n'
            self.linAlpha += cAlpha + '\n'
            self.linAlpha = self.linAlpha.replace('flux', 'val')

    def storeRhs(self):
        """Store rhs expression in rhs and rhsBound.
        """
        # for f function (forcing)
        if not self.massCoef:
            self.rhs = '      phi = RangeType(0);\n'
        else:
            cRhs = self.stringToCSource(self.massCoef)
            for coef in self.unsetCoefficients:
                if coef in self.massCoef:
                    self.rhs += '      typename ' + coef + 'DiscreteFunction::RangeType ' + coef + ';\n      ' \
                                + coef + 'Local_->evaluate(point,' + coef + ');\n'
            self.rhs += cRhs + '\n'
            self.rhs = self.rhs.replace('flux', 'phi')
        # for n function (neumann data)
        if not self.alphaCoef:
            self.rhsBound = '      value = RangeType(0);\n'
        else:
            cNeuman = self.stringToCSource(self.alphaCoef)
            if self.xTest(cNeuman) == 1:
                self.rhsBound += '      DomainType x = entity_->geometry().global( Dune::Fem::coordinate(point) );\n'
            for coef in self.unsetCoefficients:
                if coef in self.alphaCoef:
                    self.rhsBound += '      typename ' + coef + 'DiscreteFunction::RangeType ' + coef + ';\n      ' \
                                     + coef + 'Local_->evaluate(point,' + coef + ');\n'
            self.rhsBound += cNeuman + '\n'
            self.rhsBound = self.rhsBound.replace('flux', 'value')

    def storeCoef(self):
        """Store named coefficients (e.g. "velocity").
        """
        for i, (coef, dim) in enumerate(self.unsetCoefficients.items()):
            self.modelTemplate += ', class ' + coef + 'DiscreteFunction'
            self.setCoef += 'void set' + coef +'(const ' + coef + 'DiscreteFunction &' + coef + 'DF)\n    {\n' \
                            '      std::cout << "set ' + coef + '" << std::endl;\n      ' + coef + 'Local_ = ' \
                            'make_shared<typename ' + coef + 'DiscreteFunction::LocalFunctionType>( ' + coef + 'DF );' \
                            '\n    }\n    std::shared_ptr<typename ' + coef + 'DiscreteFunction::LocalFunctionType> ' \
                            + coef + 'Local_;'
            self.initCoef += 'if (!' + coef + 'Local_)\n      {\n        std::cout << "' + coef + ' not initialized -' \
                             ' call set' + coef + ' first" << std::endl;\n        abort();\n      }\n      ' \
                             + coef + 'Local_->init(entity);'
            self.pyTemplate += ', Dune::FemPy::VirtualizedGridFunction< GridPart, ' + coef + 'RangeType >'
            self.pySetCoef += '.def( "set' + coef + '", &PyModel::set' + coef + ' ) \\'
            self.pyRangeType += 'static const int ' + coef + 'dimRange = ' + str(dim) + ';\nstatic const int ' + coef +'dimDomain = ' \
                           'GridPart::dimensionworld;\ntypedef Dune::Fem::FunctionSpace< double, double, ' + coef + 'dimDomain, ' \
                           + coef + 'dimRange > ' + coef + 'FunctionSpaceType;\ntypedef typename ' + coef + 'FunctionSpaceType::RangeType ' \
                           + coef + 'RangeType;'
            if not (i + 1) == len(self.unsetCoefficients):
                self.setCoef += '\n    '
                self.initCoef += '\n      '
                self.pySetCoef += '\n  '
                self.pyRangeType += '\n'

    def uflToStrong(self):
        """Calculate strong form of an expression from the weak form.
        """
        # need to call formOutput(a), storeSrc() first for this to work
        self.residual = sympy.zeros(self.dimR, 1)
        subs = []
        for i in range(0, self.dimR):
            dexactdx = sympy.diff(self.U[i], self.x0)
            dexactdy = sympy.diff(self.U[i], self.x1)
            dexactdz = sympy.diff(self.U[i], self.x2)
            subs.extend([(self.v0[i], self.U[i]), \
                           (self.dv0[i, 0], dexactdx), \
                           (self.dv0[i, 1], dexactdy), \
                           (self.dv0[i, 2], dexactdz)])
        for i in range(0, self.dimR):
            if self.matCodeSrc[i]:
                self.residual[i] = self.residual[i] + self.matCodeSrc[i].subs(subs)
            if self.matCodeFlux[i, 0]:
                self.residual[i] = self.residual[i] - sympy.diff(self.matCodeFlux[i, 0].subs(subs), self.x0)
            if self.matCodeFlux[i, 1]:
                self.residual[i] = self.residual[i] - sympy.diff(self.matCodeFlux[i, 1].subs(subs), self.x1)
            # if self.matCodeFlux[i, 2]:
            #     self.residual[i] = self.residual[i] - sympy.diff(self.matCodeFlux[i, 2].subs(subs), self.x2)

    def computeForcing(self, exact):
        subs = []
        for i in range(0, self.dimR):
            subs.extend([(self.U[i], exact[i])])
        forcingStr = [self.sympyToString(self.residual[i].subs(subs).doit()) for i in range(0, self.dimR)]
        return forcingStr

    def xTest(self, string):
        """Test whether a C++ string has "x" as a variable.
        """
        for i in range(len(string)):
            if string[i] == 'x':
                if string[i-1].isalpha() == 0:
                    return True

    ########################
    # Main methods
    ########################
    def setCoefficient(self, name, expr):
        """For setting a named coefficient.
        """
        self.coefficients[name] = expr

    def generateFromExact(self, a, exact, *args):
        """Generate a DUNE model file using a UFL expression (this version uses 'exact' to calculate RHS).

        Args:
            a : UFL expression for the bilinear form.
            exact : Sympy expression for the exact solution (used for error analysis).
            *g : (optional) UFL expression for any Dirichlet conditions.

        Returns:
            Generates a DUNE file called Model.hh where "Model" is the name used to initialise DuneUFLModel.
        """
        # dirichlet conditions
        self.diricOutput(args)
        # define variables
        # store form a
        self.formOutput(a)
        self.storeSrc()
        # calculate strong form
        self.uflToStrong()
        # define rhs forcing using exact solution and print
        x = self.spatialCoordinate()
        forcing = self.computeForcing(exact)
        #uflForcing = as_vector( [eval(f) for f in forcing] ) # some problem with # python 3.4?
        # uflForcing = [eval(f) for f in forcing]
        uflForcing = []
        for f in forcing:
            uflForcing.extend([eval(f)])
        uflForcing = as_vector(uflForcing)
        L = inner(self.testFunction(), uflForcing)*dx(0)
        self.formOutput(L)
        self.storeRhs()
        # define linearization
        F = apply_derivatives(derivative(action(a, self.u_), self.u_, self.trialFunction()))
        self.formOutput(F)
        self.storeLin()
        if self.unsetCoefficients:
            self.storeCoef()
        # output model
        self.modelPrint(exact)

    def generate(self, a, L, exact, *args):
        """Generate a DUNE model file using a UFL expression.

        Args:
            a : UFL expression for the bilinear form.
            L : UFL expression for the RHS
            exact : Sympy expression for the exact solution (used for error analysis).
            *g : (optional) UFL expression for any Dirichlet conditions.

        Returns:
            Generates a DUNE file called Model.hh where "Model" is the name used to initialise DuneUFLModel.
        """
        # dirichlet conditions
        self.diricOutput(args)
        # calculate strong form
        self.uflToStrong()
        # store form a
        self.formOutput(a)
        self.storeSrc()
        # print (forcing)
        self.formOutput(L)
        self.storeRhs()
        # define linearization
        F = apply_derivatives(derivative(action(a, self.u_), self.u_, self.trialFunction()))
        self.formOutput(F)
        self.storeLin()
        if self.unsetCoefficients:
            self.storeCoef()
        # output model
        self.modelPrint(exact)

    def make(self, grid):
        """Create wrapper file.
        """
        startTime = timeit.default_timer()

        inputfile = self.path + '/modelimpl.hh.in'
        outputfile = self.path + '/modelimpl.hh'
        if isinstance(grid, ModuleType):
            module = grid
        else:
            module = grid._module
        name = self.modelName + "Model" + module._typeHash
        print(name)

        if femmpi.comm.rank == 0:
            with open(outputfile, "wt") as fout:
                with open(inputfile, "rt") as fin:
                    fout.write(module._includes)
                    for line in fin:
                        line = line.replace('ModelTmp', self.modelName + "Model")
                        if '#MODELNAME' in line:
                            line = line.replace('#MODELNAME', name)
                        if '#GRIDPARTCHOICE' in line:
                            line = line.replace('#GRIDPARTCHOICE', module._typeName)
                        if '#PYTEMPLATE' in line:
                            line = line.replace('#PYTEMPLATE', self.pyTemplate)
                        if '#DIMRANGE' in line:
                            line = line.replace('#DIMRANGE', str(self.dimR))
                        if '#PYSETCOEFFICIENT' in line:
                            if self.unsetCoefficients:
                                line = line.replace('#PYSETCOEFFICIENT', self.pySetCoef)
                            else:
                                line = ''
                        elif '#PYRANGETYPE' in line:
                            line = line.replace('#PYRANGETYPE', self.pyRangeType)
                        fout.write(line)
            print("modelimpl.hh has been changed")

            # the new model is constructed in the file modeltmp.cc for which make targets exist:
            cmake = subprocess.Popen(["cmake", "--build", "../../..", "--target", "modelimpl"], cwd=compilePath)
            cmake.wait()
            os.rename(os.path.join(compilePath, "modelimpl.so"), os.path.join(compilePath, name + ".so"))


            print(self.modelName + 'Model.cc and ' + self.modelName + 'Model.hh have been created')
        femmpi.comm.barrier()

        print("Building model took ", timeit.default_timer() - startTime, "s")

        return name

    def makeAndImport(self, grid):
        """Does make and imports the module.
        """
        name = self.make(grid)
        return importlib.import_module("dune.generated."+name)

    def exportTodot(self, a):
        """Draw Graph of a UFL expression into PDF format.
        """
        acdi = self.expandExpr(a)
        if isinstance(acdi, Form):
            integrals = acdi.integrals()
            for op in integrals:
                self.exportTodot(op.integrand())
        else:
            textFile = open("test.dot", "w")
            textFile.write("{}".format(ufl2dot(a, labeling="compact")[0]))
            textFile.write("{}".format(ufl2dot(a, labeling="repr")[0]))
            textFile.write("{}".format(ufl2dot(acdi, labeling="compact", object_names={})[0]))
            textFile.write("{}".format(ufl2dot(acdi, labeling="repr")[0]))
            textFile.close()
            os.system("dot -Tps test.dot -o test.pdf")

    def trialFunction(self):
        """Initialise trial function (typically "u").
        """
        return TrialFunction(self.element)

    def testFunction(self):
        """Initialise test function (typically "v").
        """
        return TestFunction(self.element)

    def spatialCoordinate(self):
        """Initialise spatial coordinate (typically "x").
        """
        return SpatialCoordinate(triangle)

    def coefficient(self, name):
        """Initialise named coefficient (e.g. "velocity" or "diffusion).
        """
        return self.NamedCoefficient(self.element, self.dimR, name)

#############################################
if __name__ == "__main__":
    import doctest
    doctest.testmod()
