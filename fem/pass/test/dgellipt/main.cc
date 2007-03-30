#include <iostream>
#include <fstream>

#include <config.h>

#include <dune/common/stdstreams.cc>
#include <dune/common/misc.hh>
#include <dune/common/timer.hh>

#include <dune/grid/io/file/dgfparser/gridtype.hh>
#include <dune/common/mpihelper.hh>

using namespace Dune;

const int dim = dimworld;
const int ncomp = 1;

// include file with the description of the convection diffusion problem
#include "discretization.hh"

using namespace LDGExample;

template <typename Field, class GridImp, int dimR, int polOrd=0 >
struct DescriptionTraits
{
  enum { dim      = GridImp::dimension };
  enum { dimworld = GridImp::dimensionworld };

  enum { dimRange = dimR };

  typedef Field   FieldType;
  typedef GridImp GridType;

  // the model 
  typedef ModelParam  <FieldType,GridType,dimRange>  ModelParamType;
  typedef Model       <ModelParamType>  ModelType;

  // the discretisation 
  typedef DiscrParam  <ModelType,polOrd>  DiscrParamType;
};

int main (int argc, char **argv)
{
  // this method calls MPI_Init, if MPI is enabled
  MPIHelper::instance(argc,argv);

  // error message if called without parameter file
  if(argc < 2)
  {
    fprintf(stderr,"usage: %s <parameter file>\n",argv[0]);
  }
  
  // read data from the parameter file
  const char * paramname = "parameter";
  if(argc == 2) paramname = argv[1];

  std::string paramfile ( paramname );

  typedef DescriptionTraits <double,GridType,ncomp,3> DescrType;
  typedef DescrType :: ModelType ModelType;
  typedef DescrType :: DiscrParamType DiscrParamType;

  ModelType model(paramfile);

  Timer timer;
  simul<DiscrParamType>(model,paramfile);
  std::cout << "CPUtime: " << timer.elapsed() << " sec!\n";

  return 0;
}

