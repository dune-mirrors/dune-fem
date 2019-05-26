from __future__ import division, print_function, unicode_literals

from .cplusplus import Declaration, TypeAlias, Variable, UnformattedExpression
from .expression import makeExpression
from .formatter.expression import formatExpression

def declareFunctionSpace(domainField, rangeField, dimDomain, dimRange,
        name="FunctionSpaceType", prefix="", **kwargs):
    dimDomain, dimRange = makeExpression(dimDomain), makeExpression(dimRange)
    code = [TypeAlias(name, "Dune::Fem::FunctionSpace< " + ", ".join([domainField, rangeField, " ".join(formatExpression(dimDomain)), " ".join(formatExpression(dimRange))]) + " >")]
    code += [TypeAlias(prefix+t, "typename " + name + "::" + t)\
            for t in ["DomainFieldType", "RangeFieldType", "DomainType", "RangeType", "JacobianRangeType", "HessianRangeType"]]
    dimDomainName = kwargs.get("dimDomainName",prefix.lower()+"dimDomain")
    if dimDomainName is not None:
        code += [Declaration(Variable("const int", dimDomainName), initializer=dimDomain, static=True)]
    dimRangeName = kwargs.get("dimRangeName",prefix.lower()+"dimRange")
    if dimRangeName is not None:
        code += [Declaration(Variable("const int", dimRangeName), initializer=dimRange, static=True)]
    return code


def fieldTensorType(shape, field = 'double'):
    field = 'std::complex< double >' if field == 'complex' else field
    if len(shape) == 0:
        return field
    elif len(shape) == 1:
        return 'Dune::FieldVector< ' + field + ', ' + str(shape[0]) + ' >'
    elif len(shape) == 2:
        return 'Dune::FieldMatrix< ' + field + ', ' + str(shape[0]) + ', ' + str(shape[1]) + ' >'
    elif len(shape) == 3:
        return 'Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< ' + field + ', ' + str(shape[1]) + ', ' + str(shape[2]) + ' >, ' + str(shape[0]) + ' >'
    else:
        raise ValueError('No C++ type defined for tensors of shape ' + str(shape) + '.')
