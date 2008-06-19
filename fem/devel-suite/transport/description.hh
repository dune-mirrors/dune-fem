#ifndef __DESCRIPTION_HH__
#define __DESCRIPTION_HH__

// Dune headers
#include <dune/fem/dfadapt.hh>
#include <dune/fem/discfuncarray.hh>
#include <dune/fem/transfer/adaptoperator.hh>
#include <dune/fem/common/boundary.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/lagrangebase.hh>
#include <dune/fem/dofmanager.hh>

// Local includes
#include "../../misc/fvfemalt.cc"
#include "../../misc/fem.cc"
#include "../../misc/finitevolume.hh"
#include "../../misc/identity.hh"

// * Use this for NumericalFluxType
using namespace Adi;

using namespace Dune;

template <typename Field, 
          int dim, 
          int dimworld, 
          template <int,int> class GridImp, 
          int dimRange, 
          template <int,int,class> class ParameterImp , 
          int polOrd=0 >
struct Descr
{
  typedef GridImp< dim, dimworld> GridType;
  typedef FunctionSpace < Field, Field , dim , dimRange > FuncSpace;
  typedef DofManager<GridType> DofManagerType;
  typedef LeafGridPart<GridType> GridPartType;
  typedef LagrangeDiscreteFunctionSpace<FuncSpace,
                                        GridPartType,
                                        polOrd> FuncSpaceType ;
  typedef DFAdapt<FuncSpaceType> DiscreteFunctionType;
  typedef ParameterImp < dim , dimworld , FuncSpace > ProblemType;
  typedef ConstantLinearAdvectionUpwindFlux<FuncSpace> NumericalFluxType;
  typedef FieldVector<Field, dimworld> NormalType;
  typedef FieldVector<Field, dimRange> ValueType;
  typedef Identity<DiscreteFunctionType> IdentityType;
};
  
//! Problem Description, describes the problem physically and programatically
template <class DescrType> 
class TransportDescription 
{
  public:
  // useful typedefs 
  typedef typename DescrType::GridType GridType;
  typedef typename DescrType::FuncSpace FuncSpace;
  typedef typename DescrType::GridPartType GridPartType;
  typedef typename DescrType::FuncSpaceType FuncSpaceType;
  typedef typename DescrType::DofManagerType DofManagerType;
  typedef typename DescrType::ProblemType ProblemType;
  
  typedef typename DescrType::DiscreteFunctionType DiscreteFunctionType; 
  typedef typename DescrType::NumericalFluxType NumericalFluxType;

  typedef typename DescrType::NormalType NormalType;
  typedef typename DescrType::ValueType ValueType;
  typedef typename DescrType::IdentityType IdentityType;

  // Problem description typdef
  typedef TransportDescription < DescrType > MyType;


  typedef ScalarFV< NumericalFluxType,
                    DiscreteFunctionType> FVType;

  // * Those terms are not used here, but still exist in fvfem.cc
  //typedef ScalarDiffusion < MyType , DiscreteFunctionType > DiffType;
  //typedef Reaction < MyType, DiscreteFunctionType> ReacType;
  
  typedef Mapping<typename FuncSpace::DomainFieldType,
                  typename FuncSpace::RangeFieldType,
                  DiscreteFunctionType,
                  DiscreteFunctionType> MappingType;
  typedef DiscreteOperator<FVType, DiscreteFunctionType> SpaceOperatorType;

  typedef AdaptOperator<GridType,RestProlOperatorFV<DiscreteFunctionType>, DofManagerType > ADOperatorType;

private: 
  //! internal types  
  static const bool allLevels = true;

  FVType *localFV_;
  
  GridType *grid_;
  GridPartType *gridPart_;
  FuncSpaceType* funcSpace_;
  DofManagerType *dm_;
  
  DiscreteFunctionType *solution_;
  DiscreteFunctionType *rhs_;
  DiscreteFunctionType *errorFunc_;

  AdaptMapping adaptMap_;
  ADOperatorType *adop_;
  RestProlOperatorFV<DiscreteFunctionType>* newRP_;
 
  SpaceOperatorType* D0_;
  IdentityType D1_;
  
  int maxLevel_;
  int problem_;
  int order_; 

  int eocsteps_;
  int startlevel_;
  int endlevel_;

  char eocfile_ [1024];
  char datafile_[1024];
  char gridfile_[1024];
  char dmfile_  [1024];
  
  double *error_; 

  double cfl_;
  
  double startTime_;
  double endTime_;
    
  ProblemType * tp_;

  NumericalFluxType* numFlux_;

  double tolerance_;
  double eta_;

  int adaptive_;
  
public:
  ~TransportDescription () {
    if (grid_) delete grid_;
  }
  
  //! read Problem Description from file 
  TransportDescription ( const char * paramfile )
  {
    readParameter(paramfile,"MaxLevel",maxLevel_);
    readParameter(paramfile,"Problem",problem_);
    readParameter(paramfile,"Order",order_);
    readParameter(paramfile,"StartLevel",startlevel_);
    readParameter(paramfile,"EndLevel",endlevel_);
    readParameter(paramfile,"EOCSteps",eocsteps_);
     
    readParameter(paramfile,"CFL",cfl_);
    readParameter(paramfile,"StartTime",startTime_);
    readParameter(paramfile,"EndTime",endTime_);
        
    readParameter(paramfile,"EOCOut",eocfile_);
    readParameter(paramfile,"DataOut",datafile_);
    readParameter(paramfile,"GridOut",gridfile_);
    readParameter(paramfile,"DMOut",dmfile_);

    int adapt;
    readParameter(paramfile,"Adaptive",adapt);
    adaptive_ = adapt;
    readParameter(paramfile,"TOL",tolerance_);
    readParameter(paramfile,"eta",eta_);

    if(endlevel_ > maxLevel_) maxLevel_ = endlevel_;
    error_ = new double [maxLevel_];
#if SGRID
    // this leads to the same number of points for SGrid and AlbertGrid
    int n[DIM];
    double h[DIM];
    for(int i=0; i<DIM; i++) 
    {
      n[i] = 1; h[i] = 1.0;
    }

    grid_ = new GridType ((int *) &n, (double *) &h );
    for(int i=0;i<maxLevel_; i++)
      grid_->globalRefine(1);
#endif

#if AGRID 
    char gridName[1024];
    readParameter(paramfile,"Grid",gridName);

    // false means that we dont use the LevelIndex 
    grid_ = new GridType (gridName);
    
    for(int i=0;i<maxLevel_; i++)
      grid_->globalRefine ( 1 );
#endif
    
#if BGRID 
    char gridName[1024];
    readParameter(paramfile,"Grid",gridName);

    // false means that we dont use the LevelIndex 
    grid_ = new GridType ( gridName );
    
    //for(int i=0;i<maxLevel_; i++) grid_->globalRefine ( 1 );
    grid_->globalRefine ( maxLevel_ );

#endif
    
    std::cout << "Number of Elements: " 
              << grid_->size(grid_->maxLevel(),0)<< "\n";
    std::cout << "Number of Points  : " 
              << grid_->size(grid_->maxLevel(),DIM) << "\n";

    dm_ = &DofManagerFactory<DofManagerType>::getDofManager ( *grid_ );

    BoundaryManager<FuncSpaceType> bc(true);

    gridPart_ = new GridPartType(*grid_);
    funcSpace_ = new FuncSpaceType( *gridPart_);
    
    // the problem 
    tp_ = new ProblemType ( *funcSpace_ ); 
      
    // the unkowns  
    solution_ = new DiscreteFunctionType ("sol", *funcSpace_ );    
    rhs_      = new DiscreteFunctionType ("rhs", *funcSpace_ );
    errorFunc_ = new DiscreteFunctionType("error", *funcSpace_);

    // adaptation
    if (adaptive()) {
      typename FuncSpaceType::IteratorType it = funcSpace_->begin();
      switch(it->geometry().type()) {
      case simplex:
        newRP_ =
          new RestProlOperatorFV<DiscreteFunctionType>(*solution_, triangle);
          newRP_->setFatherChildWeight(0.5);
          break;
      case triangle:
        newRP_ = 
          new RestProlOperatorFV<DiscreteFunctionType>(*solution_, triangle);
          newRP_->setFatherChildWeight(0.5);
        break;
      case tetrahedron:
        newRP_ = new 
          RestProlOperatorFV<DiscreteFunctionType>(*solution_, tetrahedron);
          // for AlbertaGrid and ALUGrid this is different 
          if( grid_->type() == AlbertaGrid_Id) newRP_->setFatherChildWeight(0.5);
          else newRP_->setFatherChildWeight(0.125);
        break;
      case hexahedron:
        newRP_= new 
          RestProlOperatorFV<DiscreteFunctionType>(*solution_, hexahedron);
          newRP_->setFatherChildWeight(0.125);
        break;
      default:
        std::cerr << "Wrong element type" << std::endl;
        std::abort();
      }

      assert( newRP_ );

      typedef AdaptOperator<
        GridType,
        RestProlOperatorFV<DiscreteFunctionType>,
        DofManagerType
        > ADOperatorType;

      adop_ = new ADOperatorType ( *grid_ , *dm_, *newRP_ );
      adaptMap_ = (*adop_);

    }
  
    // the schemes
    NormalType velo;
    tp_->v(velo, velo); // * First paramater is only dummy here
    NumericalFluxType* numFlux_ = new NumericalFluxType(velo);
    localFV_ = new FVType(*numFlux_, 
                          bc, 
                          *funcSpace_ , 
                          errorFunc_, 
                          adaptive(), 
                          true );

    D0_ = new SpaceOperatorType(*localFV_, true);
    
    solution_->clear ();
    rhs_->clear ();

    std::cout << "Constructor finished \n";
  }

  bool adaptive() const {
    return adaptive_ > 0;
  }

  ProblemType & parameter () const 
  {
    return *tp_;
  }

  int eocSteps () const 
  {
    return eocsteps_;
  }
  
  int step () const 
  {
#if SGRID
    return 1;
#else
    return 2;
#endif
  }

  double refineTol () const 
  {
    return tolerance_;
  }
  
  double coarsenTol () const 
  {
    return eta_ * refineTol();
  }
  
  void addLevel (int step)  
  {
    startlevel_ += step;
    endlevel_   += step;
  }
  
  int startLevel () const 
  {
    return startlevel_;
  }
  
  int endLevel () const 
  {
    return endlevel_;
  }

  void setTime ( double time )
  {
    tp_->setTime(time);
  }
  
  double cfl () const 
  {
    return cfl_;
  } 
  
  // i.e. 0.0
  double startTime () const 
  {
    return startTime_;
  } 
 
  // i.e. T 
  double endTime () const 
  {
    return endTime_;
  } 
  
  // return maxlevel
  int level ()
  {
    return maxLevel_;
  }

  // clear solution and rhs 
  void clear ()
  {
    solution_->clear ();
    rhs_->clear ();
  };

  // return reference to solution 
  DiscreteFunctionType& solution ()
  {
    return *solution_;
  };

  //! allocate a vector of appropriate format; 
  //! memory for vector is managed extern
  DiscreteFunctionType* temporary ()
  {
    return ( new DiscreteFunctionType ("tmp", *funcSpace_ ) );
  };

  //! build right hand side, does not allocate b!
  DiscreteFunctionType & rhs (double time)//, DiscreteFunctionType & b)
  {
    // is done with setup of matrix 
    rhs_->clear ();
    return (*rhs_);
  };

  DiscreteFunctionType& initialData() 
  {
    L2Projection<DiscreteFunctionType> pro;
    //for(int level=startlevel_; level<=endlevel_; level += step())
    //{
    int level = grid_->maxLevel();
    solution_->clear();
    pro.template lumpi<1> (level, tp_->initialData() , *solution_);
    //}
    return *solution_;
  };

  bool explicitSpaceOperatorExists() const
  {
    return true;
  } 
  
  // return space operator D_0
  MappingType& D0() {
    return *D0_;
  }

  // return time operator D_1
  MappingType& D1() {
    return D1_;
  }

  void adapt()
  {
    //funcSpace_->adapt(*solution_,*rhs_);
    if(adaptive_)
    {
      adaptMap_.adapt();
      std::cout << "Number of global Elements: " << grid_->global_size(0) << "\n";
    }
    L1Norm l1n;
    std::cout << l1n.template norm<1>(grid_->maxLevel(), *solution_)  << " L1 Norm! \n";
  }

  GridType& grid() {
    return *grid_;
  }
  
  // print data 
  void printData(double time , int timestep)
  {
    grid_->write(xdr , gridfile_ ,time,timestep);
    solution_->write(xdr,datafile_,timestep);
    rhs_->write(xdr, "data/rhs" , timestep);
    dm_->write(xdr, dmfile_,timestep);
  }
  
  double eocCalc (double time, int norm)
  {
    std::cout << "Calc L"<< norm << " EOC \n";
    double *eoc = new double [maxLevel_];
    if(endlevel_ - startlevel_ >= step())
    {
      if(norm == 2)
      {
        for(int i=startlevel_; i<=endlevel_; i += step())
          error_[i] = L2error(i,time);
      }
      else 
      {
        for(int i=startlevel_; i<=endlevel_; i += step())
          error_[i] = error(i,time);
      }

      for(int i=startlevel_; i<=endlevel_-step(); i += step())
        eoc[i] = log ( error_[i] / error_[i+step()] ) / M_LN2;
      
      for(int i=startlevel_; i<=endlevel_-step(); i += step())
        std::cout << "Level " << i << " | EOC " << eoc[i] << " | Error "
          << error_[i] << "\n";
    }

    return eoc[startlevel_];
  }
  
  // calc error if axact solution exists 
  double error (int level, double time) 
  {
    L1Error l1err;
    if(tp_->hasExactSolution())
      return l1err.template error<3> ( level, tp_->exactSolution() ,*solution_, time);
    else 
      return -1.0;
  }
  
  // calc error if axact solution exists 
  double L2error (int level, double time) 
  {
    L2Error<DiscreteFunctionType> l2err;
    if(tp_->hasExactSolution())
      return l2err.template norm<3> ( level, tp_->exactSolution() ,*solution_, time);
    else 
      return -1.0;
  }
  
};

#endif
