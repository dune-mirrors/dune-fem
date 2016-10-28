from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

import dune.common.checkconfiguration as checkconfiguration

def adaptive(space, name="tmp", instanciate=True, **unused):
    """create a discrete function - using the fem adaptive storage as linear algebra backend

    Args:
        space: discrete space

    Returns:
        DiscreteFunction: the constructed discrete function
    """

    from . import module
    includes = [ "dune/fem/function/adaptivefunction.hh" ] + space._module._includes
    spaceType = space._module._typeName
    typeName = "Dune::Fem::AdaptiveDiscreteFunction< " + spaceType + " >"

    if instanciate:
        return module("fem", includes, typeName).DiscreteFunction(space,name)
    else:
        linearOperatorType = "Dune::Fem::SparseRowLinearOperator< " + typeName + "," + typeName + ">"
        includes += ["dune/fem/operator/linear/spoperator.hh"]
        return includes, typeName, linearOperatorType

def eigen(space, name="tmp", instanciate=True, **unused):
    """create a discrete function - using the eigen library as linear algebra backend

    Args:
        space: discrete space

    Returns:
        DiscreteFunction: the constructed discrete function
    """

    try:
        checkconfiguration.preprocessorTest([ ("#if HAVE_EIGEN","Eigen package is not available") ])
    except checkconfiguration.ConfigurationError as err:
        if logger.getEffectiveLevel() == logging.DEBUG:
            raise
        else:
            print("configuration error while creating a discrete function with storage=eigen exiting...")
            sys.exit(-1)

    from . import module
    includes = [ "dune/fem/function/vectorfunction/managedvectorfunction.hh", "dune/fem/storage/eigenvector.hh" ] + space._module._includes
    spaceType = space._module._typeName
    field = space.field
    typeName = "Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< " +\
            spaceType + ", Dune::Fem::EigenVector< " + field + " > > >"

    if instanciate:
        return module("eigen", includes, typeName).DiscreteFunction(space,name)
    else:
        linearOperatorType = "Dune::Fem::EigenLinearOperator< " + typeName + "," + typeName + ">"
        includes += ["dune/fem/operator/linear/eigenoperator.hh"]
        return includes, typeName, linearOperatorType

def istl(space, name="tmp", instanciate=True,  **unused):
    """create a discrete function - using the fem istl storage as linear algebra backend

    Args:
        space: discrete space

    Returns:
        DiscreteFunction: the constructed discrete function
    """

    from . import module
    includes = [ "dune/fem/function/blockvectorfunction.hh" ] + space._module._includes
    spaceType = space._module._typeName
    typeName = "Dune::Fem::ISTLBlockVectorDiscreteFunction< " + spaceType + " >"

    if instanciate:
        return module("istl", includes, typeName).DiscreteFunction(space,name)
    else:
        linearOperatorType = "Dune::Fem::ISTLLinearOperator< " + typeName + "," + typeName + ">"
        includes += ["dune/fem/operator/linear/istloperator.hh"]
        return includes, typeName, linearOperatorType
