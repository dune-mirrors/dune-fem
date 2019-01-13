from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

def solvers(includes, storage, operator):
    _, dfIncludes, dfTypeName, linearOperatorType, _, _ = storage
    includes += dfIncludes + ["dune/fempy/parameter.hh"]
    typeName = operator(dfTypeName,linearOperatorType)
    return includes, typeName


def femsolver(storage,solverType=None):
    includes = ["dune/fem/solver/krylovinverseoperators.hh"]

    operator = lambda df,_: "Dune::Fem::KrylovInverseOperator< " + df + " >"

    includes, typeName = solvers(includes,storage,operator)
    parameter = {"fem.solver.newton.linear.krylovmethod":solverType} if solverType is not None else {}
    return "fem",includes,typeName, parameter

def istlsolver(storage,solverType=None):
    includes = ["dune/fem/solver/istlinverseoperators.hh"]

    operator = lambda df,_: "Dune::Fem::ISTLInverseOperator< " + df + " >"

    includes, typeName = solvers(includes,storage,operator)
    parameter = {"fem.solver.newton.linear.krylovmethod":solverType} if solverType is not None else {}
    return "istl",includes,typeName, parameter

def suitesparsesolver(storage,solverType="umfpack"):
    includes = ["dune/fem/solver/ldlsolver.hh", "dune/fem/solver/spqrsolver.hh", "dune/fem/solver/umfpacksolver.hh"]

    if solverType == "ldl":
        operator = lambda df,inop: "Dune::Fem::LDLOp<" + ",".join([df,linop]) + " >"
    elif solverType == "spqr_symmetric":
        operator = lambda df,linop: "Dune::Fem::SPQROp< " + ",".join([df,linop,"true"]) + " >"
    elif solverType == "spqr_nonsymmetric":
        operator = lambda df,linop: "Dune::Fem::SPQROp< " + ",".join([df,linop,"false"]) + " >"
    elif solverType == "umfpack":
        operator = lambda df,linop: "Dune::Fem::UMFPACKOp< " + ",".join([df,linop]) + " >"
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
    includes = ["dune/fem/solver/petscsolver.hh"]

    operator = lambda df,_: "Dune::Fem::PetscInverseOperator< " + df + " >"

    includes, typeName = solvers(includes,storage,operator)
    parameter = {"petsc.kspsolver.method":solverType} if solverType is not None else {}
    return "petsc",includes,typeName, parameter
