from __future__ import absolute_import, division, print_function, unicode_literals

import logging
logger = logging.getLogger(__name__)

def solvers(includes, storage, operator):
    includes += storage.includes + storage.linopIncludes +["dune/fempy/parameter.hh"]
    typeName = operator(storage.type, storage.linopType)
    return includes, typeName


def femsolver(storage,solverType=None):
    includes = ["dune/fem/solver/krylovinverseoperators.hh"]

    operator = lambda df,_: "Dune::Fem::KrylovInverseOperator< " + df + " >"

    includes, typeName = solvers(includes,storage,operator)
    parameter = {"linear.method":solverType} if solverType is not None else {}
    return "fem",includes,typeName, parameter

def istlsolver(storage,solverType=None):
    includes = ["dune/fem/solver/istlinverseoperators.hh"]

    operator = lambda df,_: "Dune::Fem::ISTLInverseOperator< " + df + " >"

    includes, typeName = solvers(includes,storage,operator)
    parameter = {"linear.method":solverType} if solverType is not None else {}
    return "istl",includes,typeName, parameter

def suitesparsesolver(storage,solverType="umfpack"):
    includes = ["dune/fem/solver/ldlsolver.hh", "dune/fem/solver/spqrsolver.hh", "dune/fem/solver/umfpacksolver.hh"]

    if solverType == "ldl":
        operator = lambda df,linop: "Dune::Fem::LDLInverseOperator< "    +df+", typename "+linop+"::MatrixType" + " >"
    elif solverType == "spqr_symmetric":
        operator = lambda df,linop: "Dune::Fem::SPQRInverseOperator< "   +df+", true, typename "+linop+"::MatrixType" + " >"
    elif solverType == "spqr_nonsymmetric":
        operator = lambda df,linop: "Dune::Fem::SPQRInverseOperator< "   +df+", false, typename "+linop+"::MatrixType" + " >"
    elif solverType == "umfpack":
        operator = lambda df,linop: "Dune::Fem::UMFPACKInverseOperator< "+df+", typename "+linop+"::MatrixType" + " >"
    else:
        raise ValueError("wrong krylov solver - only ldl,spqr_symmetric,spqr_nonsymmetric,umfpack available")

    includes, typeName = solvers(includes,storage,operator)
    return "suitesparse",includes,typeName,{}

def eigensolver(storage,solverType="bicgstab"):
    includes = ["dune/fem/solver/eigen.hh"]

    if solverType == "cg":
        operator = lambda df,_: "Dune::Fem::EigenCGInverseOperator< " + df + " >"
    elif solverType == "bicgstab":
        operator = lambda df,_: "Dune::Fem::EigenBiCGStabInverseOperator< " + df + " >"
    else:
        raise ValueError("wrong krylov solver - only cg,bicgstab available")

    includes, typeName = solvers(includes,storage,operator)
    return "eigen",includes,typeName, {}

def viennaclsolver(storage,solverType="gmres"):
    includes = ["dune/fem/solver/viennacl.hh"]

    if solverType == "cg":
        operator = lambda df,_: "Dune::Fem::ViennaCLCGInverseOperator< " + df + " >"
    elif solverType == "gmres":
        operator = lambda df,_: "Dune::Fem::ViennalCLGMResInverseOperator< " + df + " >"
    elif solverType == "bicgstab":
        operator = lambda df,_: "Dune::Fem::ViennalCLBiCGStabInverseOperator< " + df + " >"
    else:
        raise ValueError("wrong krylov solver - only cg,gmres,bicgstab available")

    includes, typeName = solvers(includes,storage,operator)
    return "viennacl",includes,typeName, {}

def petscsolver(storage,solverType=None):
    includes = ["dune/fem/solver/petscinverseoperators.hh"]

    operator = lambda df,_: "Dune::Fem::PetscInverseOperator< " + df + " >"

    includes, typeName = solvers(includes,storage,operator)
    parameter = {"nonlinear.linear.method":solverType} if solverType is not None else {}
    return "petsc",includes,typeName, parameter

def amgxsolver(storage,solverType="gmres"):
    print("AMGX sovler used ")
    includes = ["dune/fem/solver/amgxsolver.hh"]

    operator = lambda df,_: "Dune::Fem::AMGXInverseOperator< " + df + " >"

    includes, typeName = solvers(includes,storage,operator)
    parameter = {"nonlinear.linear.method":solverType} if solverType is not None else {}
    return "amgx",includes,typeName, parameter
