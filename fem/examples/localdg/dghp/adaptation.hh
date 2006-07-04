/***********************************************************************************************
 
   Sourcefile:  adaptation.cc

   Titel:       grid and function adaptation due to error indicator

   Decription:  classes: Adaptation


***********************************************************************************************/


#ifndef __DUNE_ADAPTATION_HH__
#define __DUNE_ADAPTATION_HH__

// include grid parts  
#include <dune/grid/common/gridpart.hh> 

// include restricion, prolongation and adaptation operator classes for discrete functions
#include <dune/fem/space/dgspace/dgadaptoperator.hh>

// include discrete functions
#include <dune/fem/space/dgspace.hh>

#include <dune/fem/discretefunction/dfadapt.hh>
#include <dune/fem/discretefunction/adaptivefunction.hh>
#include <dune/fem/quadrature/cachequad.hh>

namespace Dune 
{


// class that holds the parameter for time discretization, like dt, time, ...
class TimeDiscrParam
{
public:
  //! constructor
  TimeDiscrParam (double dt = 1.0, double time = 0.0, int n = 0) : dt_(dt), time_(time), n_(n) {};

  //! destructor 
  ~TimeDiscrParam () {};

public:
  //! get and set time step size
  double getTimeStepSize () {return dt_;};

  void setTimeStepSize (double dt) { dt_ = dt; return;};

  //! get and set time step size for next step
  double getNewTimeStepSize () {return dtNew_;};

  void setNewTimeStepSize (double dt) { dtNew_ = dt; return;};


  //! get and set current time
  double getTime () {return time_;};

  void setTime (double time) {time_ = time; return;};

  //! get and set time step number
  int getTimeStepNumber () {return n_;};

  void setTimeStepNumber (int n) {n_ = n; return;};

private:
  double dt_;
  double dtNew_;
  double time_;
  int n_;
};




// class for the organization of the adaptation prozess
template <class DiscFuncType, class TimeDiscrParamType>
class Adaptation{

public:
  // useful enums and typedefs

  // discrete function type of adaptive functions
  typedef DiscFuncType                                      DiscreteFunctionType;
  typedef typename DiscreteFunctionType::FunctionSpaceType  DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::GridType      GridType;
  typedef typename DiscreteFunctionSpaceType::GridPartType  GridPartType;
  
  typedef DofManager<GridType>                              DofManagerType;
 

  enum { dim      = GridType::dimension };
  enum { dimworld = GridType::dimensionworld };


  // initialize functionspace, etc., for the indicator function
  typedef FunctionSpace < double, double, dim, 1 >          IndicatorFuncSpaceType;

  typedef DiscontinuousGalerkinSpace<IndicatorFuncSpaceType, 
			   GridPartType, 0, CachingStorage> IndicatorDiscFuncSpaceType;

  typedef DFAdapt < IndicatorDiscFuncSpaceType >            IndicatorDiscreteFunctionType; 

  typedef typename IndicatorDiscreteFunctionType::LocalFunctionType IndicatorLocalFuncType;


  // initialize restricion, prolongation and adaptation operator for discrete functions
  typedef RestProlOperator<DiscreteFunctionType>            RestProlType; 
  typedef RestProlOperator<IndicatorDiscreteFunctionType>   IndicatorRestProlType; 


  typedef AdaptOperator<GridType,RestProlType>              AdaptOperatorType;
  typedef AdaptOperator<GridType,IndicatorRestProlType>     IndicatorAdaptOperatorType;

  typedef AdaptMapping                                      AdaptMappingType;




public:
  //! constructor
  Adaptation (GridPartType &gridPart, TimeDiscrParamType &param,  const char * paramfile = 0) 
    : gridPart_(gridPart) , grid_(const_cast<GridType &>(gridPart_.grid()))
    , timeDiscrParam_(param)
  {
    // initialize discrete function space
    discFuncSpace_ = new IndicatorDiscFuncSpaceType( gridPart_ );

    // initialize indicator functions
    indicator_ = new IndicatorDiscreteFunctionType ("indicator", *discFuncSpace_);
    indicator_->set(0.0);


    // parameter die noch eingelesen werden muessen
    if(paramfile){
      readParameter(paramfile,"TOL",globalTolerance_);
      readParameter(paramfile,"Theta",coarsenTheta_);
      readParameter(paramfile,"EndTime",endTime_);
    }
    else{
      // set default values
      globalTolerance_ = 100000.;
      coarsenTheta_ = 0.1;
      endTime_ = 1.0;
    }

  };

  //! destructor
  ~Adaptation () {};

public: 

  void clearIndicator(){
    indicator_->set(0.0);
    return;
  }

  
  template<class EntityType>
  void addToLocalIndicator(EntityType &en, double val){
    //std::cout << "  addToLocalInd:  " << val << std::endl;
    IndicatorLocalFuncType localIndicator = indicator_->localFunction( en);
    localIndicator[0] += val;	
    return;
  };

  template<class EntityType>
  void setLocalIndicator(EntityType &en, double val){
    IndicatorLocalFuncType localIndicator = indicator_->localFunction( en);
    localIndicator[0] = val;
    return;
  };


  template<class EntityType>
  double getLocalIndicator(EntityType &en){
    IndicatorLocalFuncType localIndicator = indicator_->localFunction( en);
    return localIndicator[0];
  };


  double getGlobalEstimator(){    
    typedef typename IndicatorDiscFuncSpaceType::IteratorType IteratorType;
    
    double sum = 0.0;

    IteratorType endit = discFuncSpace_->end();
    for (IteratorType it = discFuncSpace_->begin(); it != endit; ++it)
    {
      sum += getLocalIndicator(*it);
    }

    return sum;
  };


  int numberOfElements (){
    typedef typename IndicatorDiscFuncSpaceType::IteratorType IteratorType;

    int num = 0;

    IteratorType endit = discFuncSpace_->end();
    for (IteratorType it = discFuncSpace_->begin(); it != endit; ++it)
    {
      num++;
    }
    
    std::cout << "Num El: "<< num << std::endl;

    return num;
  };

  double getLocalInTimeTolerance (){
    double dt = timeDiscrParam_.getTimeStepSize();
   
    localInTimeTolerance_ = globalTolerance_ * dt / endTime_;

    // inserted for testing!!
    localInTimeTolerance_ = 500.0;

    return localInTimeTolerance_;
  };


  double getLocalTolerance (){

    getLocalInTimeTolerance ();
    localTolerance_ = localInTimeTolerance_ / numberOfElements();

    // inserted for testing!!
    // localTolerance_ = 0.5;

    return localTolerance_;
  };


  void markEntities (){
    typedef typename IndicatorDiscFuncSpaceType::IteratorType IteratorType;
    
    getLocalTolerance();

    IteratorType endit = discFuncSpace_->end();
    for (IteratorType it = discFuncSpace_->begin(); it != endit; ++it)
    {
      std::cout << " ind  tol: " << getLocalIndicator(*it) << "  " <<  localTolerance_ << std::endl;
      if( (it->level() < 15) && (getLocalIndicator(*it) > localTolerance_) )
        grid_.mark(1, it);
      else if ( (it->level() > 4) && (getLocalIndicator(*it) < coarsenTheta_ * localTolerance_) )
      	grid_.mark(-1, it);
      else
	grid_.mark(0, it);
    }

    return;
    
  };


  void adapt(){
    //std::cout << " Adaptation called !! " << std::endl;
    
    adaptMapping_.adapt();
  };

  template <class DiscFuncImp>
  void calcIndicator(DiscFuncImp &func){
    typedef typename DiscFuncImp::DiscreteFunctionSpaceType FunctionSpaceType;
    typedef typename DiscFuncType::LocalFunctionType        LocalFunctionType;
    typedef typename FunctionSpaceType::GridType            GridType;
    typedef typename FunctionSpaceType::RangeType           RangeType;
    typedef typename GridType::template Codim<0>::Entity    EntityType;
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType; 
    typedef typename FunctionSpaceType::IteratorType        IteratorType;
    typedef typename EntityType::IntersectionIterator           IntersectionIteratorType;

    typedef CachingQuadrature<GridType,0>                   VolumeQuadratureType;

    const FunctionSpaceType & functionSpace_= func.getFunctionSpace();

    // get iterator from space 
    IteratorType it = functionSpace_.begin();
    IteratorType endit = functionSpace_.end();

    RangeType insRet;
    RangeType outRet;
    double val;

    // run through grid and apply the local operator
    for( ; it != endit; ++it ){

      EntityType & inside = const_cast<EntityType &> ( *it );
      
      VolumeQuadratureType volQuad(*it,0);

      IntersectionIteratorType endnit = it->iend();
      for(IntersectionIteratorType nit = it->ibegin(); nit != endnit; ++nit){

	if( nit.neighbor() ){

	  EntityPointerType ep = nit.outside();
	  EntityType & outside = const_cast<EntityType &> ( *ep);

	  LocalFunctionType insFunct  = func.localFunction ( inside  ); 
	  LocalFunctionType outFunct  = func.localFunction ( outside ); 

	  insFunct.evaluate (inside,  volQuad , 0, insRet);
	  outFunct.evaluate (outside, volQuad , 0, outRet);

	  outRet[0] = 0.0;


	  val = (insRet[0] - outRet[0])*(insRet[0] - outRet[0]);

	  addToLocalIndicator(*it, val);
	
	}
      }

    }


  }


  template <class DiscFuncImp>
  void addAdaptiveFunction(DiscFuncImp *func){
    typedef RestProlOperator<DiscFuncImp>                    RestProlImp; 
 
//    typedef AdaptOperator<GridType,RestProlImp,DofManagerType> AdaptOpImp;
//  last template parameter is superfluous
    typedef AdaptOperator<GridType,RestProlImp>              AdaptOpImp;
    
    RestProlImp *restProl;
    AdaptOpImp  *adaptOp;

    GeometryType eltype = (*(grid_.template lbegin<0>(0))).geometry().type();

    restProl = new RestProlImp ( *func);

    typedef DofManagerFactory<DofManagerType> DofManagerFactoryType; 
 
    adaptOp  = new AdaptOpImp( grid_  , *restProl );

    adaptMapping_ = (*adaptOp);
  }

 template <class DiscFuncImp>
 void addAdaptiveFunction(DiscFuncImp *func, DiscFuncImp *func2){
    typedef RestProlOperator<DiscFuncImp>                    RestProlImp; 
 
//    typedef AdaptOperator<GridType,RestProlImp,DofManagerType> AdaptOpImp;
    typedef AdaptOperator<GridType,RestProlImp> AdaptOpImp;
 
    
    RestProlImp *restProl;
    AdaptOpImp  *adaptOp;
    
    RestProlImp *restProl2;
    AdaptOpImp  *adaptOp2;

    GeometryType eltype = (*(grid_.template lbegin<0>(0))).geometry().type();


    restProl = new RestProlImp ( *func , eltype );

     restProl2 = new RestProlImp ( *func2 , eltype );

    typedef DofManagerFactory<DofManagerType> DofManagerFactoryType; 
   
    adaptOp  = new AdaptOpImp( grid_ , *restProl );
    adaptOp2  = new AdaptOpImp( grid_ , *restProl2 );


    adaptMapping_ = (*adaptOp) + (*adaptOp);

  };

  //! export indicator function
  IndicatorDiscreteFunctionType& indicator () {  
    return *indicator_;
  };


private:
  //! grid part, has grid and ind set 
  GridPartType  &gridPart_;

  //! the grid 
  GridType & grid_;

  //! function space for discrete solution
  IndicatorDiscFuncSpaceType *discFuncSpace_;


  //! indicator function
  IndicatorDiscreteFunctionType *indicator_;


  //! parameters for adaptation
  double globalTolerance_;
  double localInTimeTolerance_;
  double localTolerance_;
  double coarsenTheta_;

  //! timestep size in time discretization parameters und endTime 
  TimeDiscrParamType & timeDiscrParam_;
  double endTime_;


  //! colloect all adaptive mappings for adaption
  AdaptMappingType adaptMapping_;

 };

} // end namespace 

#endif
