#ifndef DUNE_FEMPY_PY_TESTPETSC4PY_HH
#define DUNE_FEMPY_PY_TESTPETSC4PY_HH

#include <dune/fempy/pybind11/pybind11.hh>

#ifdef PETSC4PY_H // will be set it petsc4py.h was included (so import_petsc4py exists and the python module as well)
  #include <petsc.h>
  namespace Dune
  {
    namespace FemPy
    {
      namespace detail
      {
        void testPetsc4PyCompatibility()
        {
          if (import_petsc4py() != 0)
            throw pybind11::import_error("Failure t import petsc4py");
          pybind11::dict petscInfo =
                    pybind11::module::import("petsc4py").attr("get_config")();
          std::string petsc4pyDir = petscInfo["PETSC_DIR"].cast<std::string>();
          if (petsc4pyDir != PETSC_DIR)
          {
            throw pybind11::import_error(
  "petsc4py uses Petsc from "+petsc4pyDir+
  " dune was configured with Petsc from "+PETSC_DIR+" - the folders need to be the same.\n"+
  "To fix this you can reinstall petsc4py"+
  " setting the environment variable PETSC_DIR="+PETSC_DIR+".\n"+
  "Note that the version of petsc4py must match the install petsc, e.g., you can try\n"+
  "     pip uninstall petsc petsc4py\n"+
  "     PETSC_DIR="+PETSC_DIR+" pip install petsc4py=="+std::to_string(PETSC_VERSION_MAJOR)+"."+std::to_string(PETSC_VERSION_MINOR)+
  "\n");
          }
        }
      }
    }
  }
#else // #ifdef PETSC4PY_H
  // make headercheck happy - should not get here normally
  namespace Dune
  {
    namespace FemPy
    {
      namespace detail
      {
        void testPetsc4PyCompatibility()
        {
          throw std::runtime_error("this is just for headercheck - should not get here normally");
        }
      }
    }
  }

#endif // #ifdef PETSC4PY_H
#endif
