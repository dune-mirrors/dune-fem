/***********************************************************************************************
 Discretization is the interface class used to define the discretization of the problem
***********************************************************************************************/
#ifndef __LDGEXAMPLE_DISCRETIZATION_HH__
#define __LDGEXAMPLE_DISCRETIZATION_HH__


/* include definition of the physical problem (coefficient functions, etc.) */
#include "model.hh"

// choose discrete function type
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/discretefunction/dfadapt.hh>

#include <dune/grid/common/gridpart.hh>

// ascii parser for parameter files 
#include <dune/fem/io/file/asciiparser.hh>

// grape data io 
#include <dune/fem/io/file/grapedataio.hh>

// if grape was configured then include headers 
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#include <dune/grid/io/visual/combinedgrapedisplay.hh>
#endif

#include "discretemodels.hh"

#include <dune/fem/operator/inverseoperators.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>

#include <dune/fem/pass/dgelliptpass.hh>

// definition of L2Error 
#include <dune/fem/misc/l2error.hh>

#include "upwindflux.cc"

using namespace Dune;

namespace LDGExample { 

template <class ModelImpType, int polOrd=0 >
struct DiscrParam
{
  typedef ModelImpType  ModelType;

  enum { dimRange = ModelType::dimRange };
 
  typedef typename ModelType::GridType             GridType;

  enum { dim      = GridType::dimension };
  enum { dimworld = GridType::dimensionworld };
  enum { polyOrder = polOrd };

  typedef typename ModelType::FieldType            FieldType;
  typedef typename ModelType::FuncSpaceType        FuncSpaceType;

  typedef LeafGridPart < GridType > GridPartType ;
};


// The actual operator
template <class GradientModelType, class LaplaceModelType> 
class MySpaceOperator :
 public Operator<
    double,
    double,
    typename LaplaceModelType::Traits::DestinationType,
    typename LaplaceModelType::Traits::DestinationType>
{
  typedef typename LaplaceModelType :: Traits Traits;
public:
  typedef MySpaceOperator<GradientModelType,LaplaceModelType> ThisType;
  enum { polOrd = LaplaceModelType :: polynomialOrder };
  
  typedef typename Traits:: DestinationType DestinationType;
  typedef typename Traits:: DestinationType DiscreteFunctionType;

  typedef typename Traits::GridPartType GridPartType ;
  typedef typename GridPartType :: Traits :: GridType GridType;

  typedef DofManager<GridType> DofManagerType; 
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

  typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename GradientModelType::Traits::DiscreteFunctionSpaceType GradDiscreteFunctionSpaceType;

  typedef StartPass<DiscreteFunctionType> Pass0Type;
  // note, the destination type of the pass 0 is the argument type of pass 1
  typedef LocalDGElliptGradientPass<GradientModelType , Pass0Type> GradPassType;
  typedef LocalDGElliptPass<LaplaceModelType, GradPassType> LastPassType;

  typedef GradDiscreteFunctionSpaceType GradSpaceType;
  typedef DiscreteFunctionSpaceType LastSpaceType;

  typedef typename Traits :: FunctionSpaceType FuncSpaceType;
  typedef typename FuncSpaceType::RangeType RangeType;

  //! the exact solution to the problem for EOC calculation 
  class ExactSolution : public Function < FuncSpaceType , ExactSolution >
  {
    typedef typename FuncSpaceType::RangeType RangeType;
    typedef typename FuncSpaceType::RangeFieldType RangeFieldType;
    typedef typename FuncSpaceType::DomainType DomainType;
  public:
    ExactSolution (FuncSpaceType &f) : Function < FuncSpaceType , ExactSolution > ( f ) {}

    //! u(x,y,z) = (x-x^2)*(y-y^2)*(z-z^2)
    void evaluate (const DomainType & x , RangeType & ret) const
    {
      ret = exactSolution( &x[0] );
    }
    void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
    {
      evaluate ( x , ret );
    }
  };

  //! the exact solution to the problem for EOC calculation 
  template <class FuncSPCType>
  class ExactGradient : public Function < FuncSPCType , ExactGradient<
                        FuncSPCType > >
  {
    typedef typename FuncSPCType::RangeType RangeType;
    typedef typename FuncSPCType::RangeFieldType RangeFieldType;
    typedef typename FuncSPCType::DomainType DomainType;
  public:
    ExactGradient (const FuncSPCType &f) 
      : Function < FuncSPCType , ExactGradient< FuncSPCType > > ( f ) {}

    void evaluate (const DomainType & x , RangeType & ret) const
    {
      exactGradient( &x[0], &ret[0] );
    }

    void evaluate (const DomainType & x , double time , RangeType & ret) const
    {
      evaluate ( x , ret );
    }
  };

  MySpaceOperator ( GridType & grid , 
                    GradientModelType & gm, 
                    LaplaceModelType & lpm , double eps, bool verbose , int steps = 2 ) 
    : grid_(grid)
    , dm_(DofManagerFactoryType::getDofManager(grid_))
    , gridPart_(grid_)
    , gradSpace_(gridPart_)
    , lastSpace_(gridPart_)
    , pass0_()
    , pass1_( gm , pass0_, gradSpace_)
    , lastPass_( lpm, pass1_, lastSpace_ , eps , 3 , verbose )
    , steps_(steps)
  {
    std::cout << "Created SpaceOperator \n";
    std::cout.flush();
  }

  void operator()(const DestinationType& arg, DestinationType& dest) const 
  {
    const_cast<ThisType&> (*this).apply(arg,dest);
  }
  
  // apply space discretisation 
  void apply(const DestinationType& arg, DestinationType& dest)
  {
    DestinationType& Arg = const_cast<DestinationType&> (arg);

    typedef typename GradientModelType::Traits::DestinationType  GradFuncType; 
    typedef typename GradientModelType::Traits::DiscreteFunctionSpaceType GradientDiscreteFunctionSpaceType; 
    typedef typename GradientDiscreteFunctionSpaceType :: RangeType  GradRangeType;

    std::vector<RangeType> error(steps_);
    std::vector<GradRangeType> gradError(steps_);
    
    for(int i=0; i<steps_; ++i)
    {
      if(i > 0)
      {
        // refineGlobal is defined in description.hh
        grid_.globalRefine(refStepsForHalf);
        dm_.resize();
      }

      Arg.set(2.5);
      FuncSpaceType sp; 
      ExactSolution exact(sp); 

      lastPass_(arg,dest);

      {
        GradFuncType grad("gradient",gradSpace_);

        typedef typename GradFuncType :: FunctionSpaceType
          DFSpaceType; 
        typedef typename DFSpaceType :: RangeType RType; 

        lastPass_.evalGradient(dest,grad);

        L2Error < GradFuncType > l2errGrad;
        ExactGradient< GradientDiscreteFunctionSpaceType >
          exactGrad(gradSpace_);

        gradError[i] = l2errGrad.norm(exactGrad , grad);

      }

      //GrapeDataDisplay < GridType > grape( gridPart_.grid() ); 
      //grape.dataDisplay( dest , false );

      L2Error < DestinationType > l2err;
      // pol ord for calculation the error chould by higher than 
      // pol for evaluation the basefunctions 
      error[i] = l2err.norm(exact , dest);

      for(int k=0; k<RangeType::dimension; k++)
        std::cout << "\nError["<<i<<"] : " << error[i][k] << "\n";
      
      for(int k=0; k<GradRangeType::dimension; k++)
        std::cout << "GradError["<<i<<"] : " << gradError[i][k] << "\n";
    
      if( i > 0 )
      {
        for(int k=0; k<RangeType::dimension; k++)
        {
          double eoc = log( error[i-1][k]/error[i][k]) / M_LN2;
          std::cout << "\nEOC["<<i <<"] = " << eoc << " \n";
        }
        
        for(int k=0; k<GradRangeType::dimension; k++)
        {
          double eoc = log( gradError[i-1][k]/gradError[i][k]) / M_LN2;
          std::cout << "Grad EOC["<<i <<"] = " << eoc << " \n";
        }
      }
    }

    GrapeDataDisplay < GridType > grape( gridPart_.grid() ); 
    grape.dataDisplay( dest );
  }

  template <class TimeProviderType>
  void timeProvider(TimeProviderType & tp )
  {
    pass1_.timeProvider(&tp);
  }

  DestinationType * createDestinationFct (std::string name) 
  {
    return new DestinationType (name , lastSpace_ );
  }

private:
  GridType & grid_;
  DofManagerType & dm_;

  // we use the same index set and grid part for all spaces
  GridPartType gridPart_;

  mutable GradSpaceType gradSpace_;
  mutable LastSpaceType lastSpace_;

  mutable Pass0Type pass0_;
  mutable GradPassType pass1_;
  mutable LastPassType lastPass_;

  const int steps_; 
};

template <class DiscrType> 
void simul(typename DiscrType::ModelType & model, std::string paramFile) 
{
  typedef typename DiscrType::GridType                     GridType;
  typedef typename DiscrType::ModelType             ModelType;
  enum { polOrd = DiscrType::polyOrder };

  // choice of fluxes 
  typedef UpwindFlux<ModelType> NumericalFluxType;

  typedef LaplaceDiscreteModel < ModelType, NumericalFluxType, polOrd > LaplaceModelType;
  typedef GradientDiscreteModel < ModelType, NumericalFluxType, polOrd > GradientModelType;
  
  typedef MySpaceOperator <  GradientModelType, LaplaceModelType > 
                SpaceOperatorType; 
  typedef typename SpaceOperatorType :: DestinationType DestinationType;
   
  // initialize grid
  char dummyfile [4096];
  const char * paramfile = paramFile.c_str();
  readParameter(paramfile,"Grid",dummyfile);
  std::string macroGridName(dummyfile);

  int level;
  readParameter(paramfile,"StartLevel",level);

  GridPtr<GridType> gridptr(macroGridName); 
  GridType & grid = *gridptr;
  grid.globalRefine(refStepsForHalf*level);
  std::cout << "Grid size = " << grid.size(0) << "\n";

  double eps = 1e-20;
  readParameter(paramfile,"CGeps",eps);
  int verbose = 0;
  readParameter(paramfile,"verbose",verbose);

  int precon = 0;
  readParameter(paramfile,"Preconditioning",precon);
  bool preConditioning = (precon == 1) ? true : false;

  int eocsteps = 2;
  readParameter(paramfile,"EOCSteps",eocsteps);

  NumericalFluxType numericalFlux(model);
  
  LaplaceModelType lpm(model, numericalFlux, preConditioning);
  GradientModelType gm(model, numericalFlux, preConditioning);

  SpaceOperatorType spaceOp(grid , gm, lpm , eps, (verbose == 0) ? false : true, eocsteps );
  
  //! storage for the discrete solution and its update
  DestinationType *solution = spaceOp.createDestinationFct("solution");
  DestinationType *tmpRhs = spaceOp.createDestinationFct("tmpRhs");

  // initial data != 0
  solution->set(0.5);
  spaceOp(*tmpRhs,*solution);

  delete solution; 
  delete tmpRhs;
}

} // end namespace LDGExample
#endif
