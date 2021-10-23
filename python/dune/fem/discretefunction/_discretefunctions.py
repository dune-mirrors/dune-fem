from __future__ import absolute_import, division, print_function, unicode_literals

import sys,os
import logging
logger = logging.getLogger(__name__)

from dune.common.checkconfiguration import assertCMakeHave, ConfigurationError
import dune.common.checkconfiguration as checkconfiguration
import dune.generator
from dune.deprecate import deprecated

from . import _solvers as solvers

@deprecated(name="storage=fem",msg="'storage=\"fem\" for the space is deprecated, use `storage=\"numpy\" instead")
def fem():
    dfType = lambda space: "Dune::Fem::AdaptiveDiscreteFunction< " + space.cppTypeName + " >"
    rdfType = lambda space,rspace: dfType(space if rspace is None else rspace)
    return lambda space, rspace=None:[\
        "fem",\
        ["dune/fem/function/adaptivefunction.hh","dune/fem/operator/linear/spoperator.hh"] +\
              space.cppIncludes,\
        dfType(space),\
        "Dune::Fem::SparseRowLinearOperator< " + dfType(space) + "," +\
        rdfType(space,rspace) + "," + "Dune::Fem::SparseRowMatrix<" + space.field + ",int>" ">",\
        solvers.femsolver,\
        "as_numpy"
    ]
@deprecated(name="storage=fem",msg="'storage=\"adaptive\" for the space is deprecated, use `storage=\"numpy\" instead")
def adaptive():
    dfType = lambda space: "Dune::Fem::AdaptiveDiscreteFunction< " + space.cppTypeName + " >"
    rdfType = lambda space,rspace: dfType(space if rspace is None else rspace)
    return lambda space, rspace=None:[\
        "fem",\
        ["dune/fem/function/adaptivefunction.hh","dune/fem/operator/linear/spoperator.hh"] +\
              space.cppIncludes,\
        dfType(space),\
        "Dune::Fem::SparseRowLinearOperator< " + dfType(space) + "," +\
        rdfType(space,rspace) + "," + "Dune::Fem::SparseRowMatrix<" + space.field + ",int>" ">",\
        solvers.femsolver,\
        "as_numpy"
    ]

def numpy():
    dfType = lambda space: "Dune::Fem::AdaptiveDiscreteFunction< " + space.cppTypeName + " >"
    rdfType = lambda space,rspace: dfType(space if rspace is None else rspace)
    return lambda space, rspace=None:[\
        "numpy",\
        ["dune/fem/function/adaptivefunction.hh","dune/fem/operator/linear/spoperator.hh"] +\
              space.cppIncludes,\
        dfType(space),\
        "Dune::Fem::SparseRowLinearOperator< " + dfType(space) + "," +\
        rdfType(space,rspace) + "," + "Dune::Fem::SparseRowMatrix<" + space.field + ",int>" ">",\
        solvers.femsolver,\
        "as_numpy"
    ]

def eigen():
    try:
        checkconfiguration.preprocessorAssert([ ("#if HAVE_EIGEN","Eigen package is not available") ])
    except checkconfiguration.ConfigurationError as err:
        print("configuration error while creating a discrete function with storage=eigen exiting...")
        print("You need to install the `eigen` package and reconfigure dune-py")
        print("adding -DEigen3_DIR='Path-to-eigen` to the CMAKE_FLAGS")
        raise

    dfType = lambda space: "Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< " +\
                           space.cppTypeName + ", Dune::Fem::EigenVector< " + space.field + " > > >"
    rdfType = lambda space,rspace: dfType(space if rspace is None else rspace)
    return lambda space,rspace=None:[\
        "eigen",\
        ["dune/fem/function/vectorfunction/managedvectorfunction.hh",\
                "dune/fem/storage/eigenvector.hh",\
                "dune/fem/operator/linear/eigenoperator.hh"] +\
              space.cppIncludes,\
        dfType(space),\
        "Dune::Fem::EigenLinearOperator< " + dfType(space) + "," + rdfType(space,rspace) + ">",\
        solvers.eigensolver,\
        "as_numpy"
    ]

def istl():
    dfType = lambda space: "Dune::Fem::ISTLBlockVectorDiscreteFunction< " + space.cppTypeName + " >"
    rdfType = lambda space,rspace: dfType(space if rspace is None else rspace)
    return lambda space,rspace=None:[\
        "istl",\
        ["dune/fem/function/blockvectorfunction.hh", "dune/fem/operator/linear/istloperator.hh"] +\
              space.cppIncludes,\
        dfType(space),\
        "Dune::Fem::ISTLLinearOperator< " + dfType(space) + "," + rdfType(space,rspace) + ">",
        solvers.istlsolver,\
        "as_istl"
    ]

def petsc():
    try:
        assertCMakeHave("HAVE_PETSC")
        dfType = lambda space: "Dune::Fem::PetscDiscreteFunction< " + space.cppTypeName + " >"
        rdfType = lambda space,rspace: dfType(space if rspace is None else rspace)
        def equalSpaces(space,rspace):
            # if not space==rspace and not rspace is None:
            #     raise NotImplementedError("Operator with petsc storage only with equal domain and range spaces implemented")
            return "Dune::Fem::PetscLinearOperator< " + dfType(space) + "," + rdfType(space,rspace) + ">"
        try:
            import petsc4py
            return lambda space, rspace=None:[\
                "petsc",\
                ["dune/fem/function/petscdiscretefunction.hh", "dune/fem/operator/linear/petscoperator.hh"] +\
                      [os.path.dirname(petsc4py.__file__)+"/include/petsc4py/petsc4py.h"] +\
                      space.cppIncludes,\
                dfType(space),\
                equalSpaces(space,rspace),\
                solvers.petscsolver,\
                "as_petsc"
            ]
        except:
            return lambda space,rspace=None:[\
                "petsc",\
                ["dune/fem/function/petscdiscretefunction.hh", "dune/fem/operator/linear/petscoperator.hh"] +\
                      space.cppIncludes,\
                dfType(space),\
                equalSpaces(space,rspace),\
                solvers.petscsolver,\
                "as_petsc"
            ]
    except ConfigurationError:
        raise ConfigurationError("petsc has not been found during configuration of dune - please add the path to petsc to the DUNE_CMAKE_FLAGS")
def petscadapt():
    try:
        assertCMakeHave("HAVE_PETSC")
        dfType = lambda space: "Dune::Fem::AdaptiveDiscreteFunction< " + space.cppTypeName + " >"
        rdfType = lambda space,rspace: dfType(space if rspace is None else rspace)
        def equalSpaces(space,rspace):
            if not space==rspace and not rspace is None:
                raise NotImplementedError("Operator with petsc storage only with equal domain and range spaces implemented")
            return "Dune::Fem::PetscLinearOperator< " + dfType(space) + "," + rdfType(space,rspace) + ">"
        return lambda space,rspace=None:[\
            "petscadapt",\
            ["dune/fem/function/adaptivefunction.hh", "dune/fem/operator/linear/petscoperator.hh"] +\
                  space.cppIncludes,\
            dfType(space),\
            equalSpaces(space,rspace),\
            solvers.petscsolver,\
            "as_numpy"
        ]
    except ConfigurationError:
        raise ConfigurationError("petsc has not been found during configuration of dune - please add the path to petsc to the DUNE_CMAKE_FLAGS")
