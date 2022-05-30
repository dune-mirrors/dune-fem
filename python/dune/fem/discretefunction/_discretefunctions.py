from __future__ import absolute_import, division, print_function, unicode_literals

import sys,os
import logging
logger = logging.getLogger(__name__)

from dune.common.checkconfiguration import assertCMakeHave, ConfigurationError
import dune.common.checkconfiguration as checkconfiguration
import dune.generator
from dune.deprecate import deprecated
from dune.common.utility import isString

from . import _solvers as solvers

def _storage( dfStorage="numpy", solverStorage=None ):

    assert isString(dfStorage)
    # if not selected, use same as discrete function
    if solverStorage is None:
        solverStorage = dfStorage
    else:
        # for petsc only the same storage for both works
        assert isString(solverStorage)
        if dfStorage == "petsc":
            assert solverStorage == "petsc"

    # make sure PETSc was found if selected
    if dfStorage == "petsc":
        try:
            assertCMakeHave("HAVE_PETSC")
        except ConfigurationError:
            raise ConfigurationError("petsc has not been found during configuration of dune - please add the path to petsc to the DUNE_CMAKE_FLAGS")

    # get storage options for discrete functions
    if dfStorage == "numpy":
        dfType  = lambda space: "Dune::Fem::AdaptiveDiscreteFunction< " + space.cppTypeName + " >"
        headers = ["dune/fem/function/adaptivefunction.hh"]
    elif dfStorage == "istl":
        dfType  = lambda space: "Dune::Fem::ISTLBlockVectorDiscreteFunction< " + space.cppTypeName + " >"
        headers = ["dune/fem/function/blockvectorfunction.hh"]
    elif dfStorage == "petsc":
        dfType  = lambda space: "Dune::Fem::PetscDiscreteFunction< " + space.cppTypeName + " >"
        headers = ["dune/fem/function/petscdiscretefunction.hh"]
    else:
        raise ValueError(f"_storage: wrong discrete function storage {dfStorage}. Valid are 'numpy','petsc','istl'")

    rdfType = lambda space,rspace: dfType(space if rspace is None else rspace)

    # as_numpy, as_istl, as_petsc
    asStorage = "as_" + dfStorage

    if solverStorage == "numpy":
        return lambda space, rspace=None:[\
            dfStorage,\
            headers + ["dune/fem/operator/linear/spoperator.hh"] + space.cppIncludes,\
            dfType(space),\
            "Dune::Fem::SparseRowLinearOperator< " + dfType(space) + "," +\
            rdfType(space,rspace) + "," + "Dune::Fem::SparseRowMatrix<" + space.field + ",int>" ">",\
            solvers.femsolver,\
            asStorage
            ]
    elif solverStorage == "istl":
        return lambda space,rspace=None:[\
            dfStorage,\
            headers + ["dune/fem/operator/linear/istloperator.hh"] + space.cppIncludes,\
            dfType(space),\
            "Dune::Fem::ISTLLinearOperator< " + dfType(space) + "," + rdfType(space,rspace) + ">",
            solvers.istlsolver,\
            asStorage
            ]
    elif solverStorage == "petsc":
        petscheader = []
        try:
            import petsc4py
            petscheader += [os.path.dirname(petsc4py.__file__)+"/include/petsc4py/petsc4py.h"]
        except:
            pass

        def equalSpaces(space,rspace):
            return "Dune::Fem::PetscLinearOperator< " + dfType(space) + "," + rdfType(space,rspace) + ">"
        return lambda space, rspace=None:[\
            dfStorage,\
            headers + ["dune/fem/operator/linear/petscoperator.hh"] + petscheader + space.cppIncludes,\
            dfType(space),\
            equalSpaces(space,rspace),\
            solvers.petscsolver,\
            asStorage
        ]
    else:
        raise ValueError(f"_storage: wrong discrete function storage {dfStorage}. Valid are 'numpy','petsc','istl'")

def numpy():
    return _storage(dfStorage="numpy")

def istl():
    return _storage(dfStorage="istl")

def petsc():
    return _storage(dfStorage="petsc")

@deprecated(name="storage=petscadapt",msg="'storage=\"petscadapt\" for the space is deprecated, use `storage=\"numpy\" instead")
def petscadapt():
    return _storage(dfStorage="numpy", solverStorage="petsc")

@deprecated(name="storage=fem",msg="'storage=\"fem\" for the space is deprecated, use `storage=\"numpy\" instead")
def fem():
    return numpy()

@deprecated(name="storage=fem",msg="'storage=\"adaptive\" for the space is deprecated, use `storage=\"numpy\" instead")
def adaptive():
    return numpy()

@deprecated(name="storage=eigen",msg="'storage=\"eigen\" for the space is deprecated, use `storage=\"numpy\" instead")
def eigen():
    return numpy()


#def eigen():
#    try:
#        checkconfiguration.preprocessorAssert([ ("#if HAVE_EIGEN","Eigen package is not available") ])
#    except checkconfiguration.ConfigurationError as err:
#        print("configuration error while creating a discrete function with storage=eigen exiting...")
#        print("You need to install the `eigen` package and reconfigure dune-py")
#        print("adding -DEigen3_DIR='Path-to-eigen` to the CMAKE_FLAGS")
#        raise
#
#    dfType = lambda space: "Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< " +\
#                           space.cppTypeName + ", Dune::Fem::EigenVector< " + space.field + " > > >"
#    rdfType = lambda space,rspace: dfType(space if rspace is None else rspace)
#    return lambda space,rspace=None:[\
#        "eigen",\
#        ["dune/fem/function/vectorfunction/managedvectorfunction.hh",\
#                "dune/fem/storage/eigenvector.hh",\
#                "dune/fem/operator/linear/eigenoperator.hh"] +\
#              space.cppIncludes,\
#        dfType(space),\
#        "Dune::Fem::EigenLinearOperator< " + dfType(space) + "," + rdfType(space,rspace) + ">",\
#        solvers.eigensolver,\
#        "as_numpy"
#    ]
