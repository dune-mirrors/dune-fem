//************************************************************
//
//  (C) written and directed by Robert Kloefkorn 
//
//************************************************************
#include <iostream>
#include <vector>
#include <assert.h>
#include <string.h>

#include <config.h>
#include <dune/common/stdstreams.cc>

using namespace Dune;


#define LARGE 1.0E308

static const int dim = DUNE_PROBLEM_DIM; 

//#include <dune/fem/space/lagrangebase.hh>
#include <dune/fem/discretefunction/dfadapt.hh>
#include <dune/common/stack.hh>

//#include <../../../space/dgspace/dgleafindexset.hh>
#include <dune/grid/common/gridpart.hh>

#include <dune/fem/io/file/grapedataio.hh>


typedef double REAL;

#include <dune/grid/io/visual/grapedatadisplay.hh>
#include <dune/grid/io/visual/combinedgrapedisplay.hh>

#include <dune/fem/visual/grape/datadisp/printhelp.cc>

#include "models.hh"
#include "stuff.cc"
#include "robstuff.cc"

typedef GridType                                        GR_GridType;
typedef DofManager<GR_GridType>                         GR_DofManagerType;
typedef DofManagerFactory <GR_DofManagerType>           GR_DofManagerFactoryType;
typedef DgType::GridPartType                            GR_GridPartType;
typedef GR_GridPartType::IndexSetType GR_IndexSetType;

typedef OutputType GR_DiscFuncType;
typedef OutputType GR_DiscFuncType;
typedef GrapeDataDisplay<GR_GridType > GrapeDispType;

void addError(GrapeDispType& disp,GR_GridType& grid,double time,
	      DgType::DestinationType& Uh) {
  typedef DgType::DestinationType::DiscreteFunctionSpaceType SpaceType;
  SpaceType* space = &SpaceType::instance(grid);
  DgType::DestinationType* lsg = new DgType::DestinationType("lsg",*space);
  lsg->set(0);
  InitialDataType problem(0.0,1,false);
  initialize(problem,*lsg,time);
  disp.addData(*lsg,"Lsg.",time);
  DgType::DestinationType* err = new DgType::DestinationType("err",*space);
  err->assign(*lsg);
  err->addScaled(Uh,-1.);
  disp.addData(*err,"error",time);

  typedef SpaceType::GridPartType GridPartType;
  enum { dimD = GridType :: dimension };
  typedef FunctionSpace < double , double, dimD , 1 > ScalarFSType;
  typedef DiscontinuousGalerkinSpace<ScalarFSType, GridPartType, 0>
    ConstDiscSType;
  ConstDiscSType* sspace = &ConstDiscSType::instance(grid);
  typedef DFAdapt<ConstDiscSType> ConstDiscFSType;
  ConstDiscFSType* l1err = new ConstDiscFSType("l1err",*sspace);
  l1err->set(1.);
  L1Error<DgType::DestinationType> L1err;
  L1err.norm(problem,Uh,time,*l1err);  
  disp.addData(*l1err,"l1-err",time);  
  
}

#include <dune/fem/visual/grape/datadisp/readtupledata.cc>
#include <dune/fem/visual/grape/datadisp/readtupparams.cc> 
#include <dune/fem/visual/grape/datadisp/readfile.cc>



int main(int argc, char **argv)
{

  if (argc < 2)
  {
    print_help("datadisp");
    return(0);
  }   

  if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "-help"))
  {
      print_help("dunedisp");
      return(0);
  }

  if(argc > 2) 
    return readParameterList(argc,argv);
  else 
    return readParameterFile(argc,argv);
 
}
