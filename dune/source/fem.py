from __future__ import division, print_function, unicode_literals

from .cplusplus import Declaration, TypeAlias, Variable

def declareFunctionSpace(domainField, rangeField, dimDomain, dimRange, name="FunctionSpaceType"):
    code = [TypeAlias(name, "Dune::Fem::FunctionSpace< " + ", ".join([domainField, rangeField, str(dimDomain), str(dimRange)]) + " >")]
    code += [TypeAlias(t, "typename " + name + "::" + t) for t in ["DomainFieldType", "RangeFieldType", "DomainType", "RangeType", "JacobianRangeType", "HessianRangeType"]]
    code += [Declaration(Variable("const int", "dimDomain"), initializer=dimDomain, static=True)]
    code += [Declaration(Variable("const int", "dimRange"), initializer=dimRange, static=True)]
    return code
