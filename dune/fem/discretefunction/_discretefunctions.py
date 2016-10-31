from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

import dune.common.checkconfiguration as checkconfiguration

from . import _solvers as solvers

def adaptive():
    dfType = lambda space: "Dune::Fem::AdaptiveDiscreteFunction< " + space._typeName + " >"
    return lambda space:[\
        "fem",\
        ["dune/fem/function/adaptivefunction.hh","dune/fem/operator/linear/spoperator.hh"] +\
              space._includes,\
        dfType(space),\
        "Dune::Fem::SparseRowLinearOperator< " + dfType(space) + "," + dfType(space) + ">",\
        solvers.pardgsolver
    ]

def eigen():
    try:
        checkconfiguration.preprocessorTest([ ("#if HAVE_EIGEN","Eigen package is not available") ])
    except checkconfiguration.ConfigurationError as err:
        if logger.getEffectiveLevel() == logging.DEBUG:
            raise
        else:
            print("configuration error while creating a discrete function with storage=eigen exiting...")
            sys.exit(-1)

    dfType = lambda space: "Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< " +\
                           space._typeName + ", Dune::Fem::EigenVector< " + space.field + " > > >"
    return lambda space:[\
        "eigen",\
        ["dune/fem/function/vectorfunction/managedvectorfunction.hh",\
                "dune/fem/storage/eigenvector.hh",\
                "dune/fem/operator/linear/eigenoperator.hh"] +\
              space._includes,\
        dfType(space),\
        "Dune::Fem::EigenLinearOperator< " + dfType(space) + "," + dfType(space) + ">",\
        solvers.eigensolver
    ]

def istl():
    dfType = lambda space: "Dune::Fem::ISTLBlockVectorDiscreteFunction< " + space._typeName + " >"
    return lambda space:[\
        "istl",\
        ["dune/fem/function/blockvectorfunction.hh", "dune/fem/operator/linear/istloperator.hh"] +\
              space._includes,\
        dfType(space),\
        "Dune::Fem::ISTLLinearOperator< " + dfType(space) + "," + dfType(space) + ">",
        solvers.istlsolver
    ]
