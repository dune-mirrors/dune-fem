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

#include "localfunction.hh"

namespace Dune 
{
template <class DiscFuncType> 
struct CoarseningError {
  typedef typename DiscFuncType::LocalFunctionType LDiscFuncType;
  typedef typename LDiscFuncType::BaseFunctionSetType BaseSetType;
  typedef typename DiscFuncType::RangeType RangeType;
  typedef typename DiscFuncType::DomainType DomainType;
  typedef typename DiscFuncType::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::GridType GridType;
  typedef typename GridType::template Codim<0>::Entity EntityType;
  typedef typename GridType::template Codim<0>::Entity::Geometry GeometryType;
  static bool compute(const DiscFuncType& df,
		      const EntityType& father,
		      RangeType& error) {
    error = 0.;
    const FunctionSpaceType& space = df.getFunctionSpace();  
    const BaseSetType& bset = space.getBaseFunctionSet(father);
    int polOrd = space.polynomOrder();
    int quadOrd = 2 * polOrd + 2;
    CachingQuadrature <GridType , 0 > quad(father,quadOrd);
    LocalFuncHelper<LDiscFuncType> ldffather(space);
    ldffather.init(father);
    ldffather.assign(0);
    RangeType u(0.0);
    RangeType v(0.0);
    RangeType phi(0.0);
    if(!father.isLeaf()) {
      typedef typename EntityType::HierarchicIterator HierarchicIterator;
      { // set function in father
	HierarchicIterator it = father.hbegin( father.level() + 1 );
	HierarchicIterator endit = father.hend  ( father.level() + 1 );
	double invdet = 1./father.geometry().integrationElement(DomainType(0.));
	for( ; it != endit; ++it) {
	  if (!(it->isLeaf())) return false;
	  LDiscFuncType ldfson = df.localFunction(*it);
	  const GeometryType& geometryInFather = it->geometryInFather();
	  for(int qP = 0; qP < quad.nop(); ++qP) {
	    double det = quad.weight(qP)
	      * it->geometry().integrationElement(quad.point(qP)) * invdet;
	    ldfson.evaluate(*it,quad,qP,u);
	    DomainType x = geometryInFather.global(quad.point(qP));
	    for (int i=0;i<ldffather.numDofs();++i) {
	      bset.eval(i,quad,qP,phi);
	      ldffather[i] += det * (u * phi);
	    }
	  }
	}
      }
      { // compute error
 	HierarchicIterator it = father.hbegin( father.level() + 1 );
	HierarchicIterator endit = father.hend  ( father.level() + 1 );
	for( ; it != endit; ++it) {
	  LDiscFuncType ldfson = df.localFunction(*it);
	  const GeometryType& geometryInFather = it->geometryInFather();
	  for(int qP = 0; qP < quad.nop(); ++qP) {
	    double det = quad.weight(qP)
	      * it->geometry().integrationElement(quad.point(qP));
	    ldfson.evaluate(*it,quad,qP,u);
	    DomainType x = geometryInFather.global(quad.point(qP));
	    ldffather.evaluateLocal(father,x,polOrd,v);
	    v -= u;
	    error += det*v.one_norm();
	  }
	}
      }
    }
    return true;
  }
};
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
  typedef typename DiscFuncType::RangeType                  RangeType;
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
  Adaptation (GridPartType &gridPart, TimeDiscrParamType &param,  const double tol,double T) 
    : gridPart_(gridPart) , grid_(const_cast<GridType &>(gridPart_.grid()))
    , timeDiscrParam_(param), globalTolerance_(tol)
  {
    // initialize discrete function space
    discFuncSpace_ = new IndicatorDiscFuncSpaceType( gridPart_ );

    // initialize indicator functions
    indicator_ = new IndicatorDiscreteFunctionType ("indicator", *discFuncSpace_);
    indicator_->set(0.0);

    // set default values
    initialTheta_ = 0.01;
    coarsenTheta_ = 0.4;
    endTime_ = T;
    alphaSigSet_ = 0.01;
    maxLevFlag_ = 1;

  }
  TimeDiscrParamType &param() {
    return timeDiscrParam_;
  }
  
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


  inline double numberOfElements (){
    /*std::cout << "    Num El: "<< indicator_->size() << std::endl;*/
    return double(grid_.size(0));
  };


  int calcSignificantElements (){
    typedef typename IndicatorDiscFuncSpaceType::IteratorType IteratorType;
    
    int maxlevel = grid_.maxLevel();
    
    numSigSet_ = 1;
    numMaxLev_ = 0;
    tolSigSetInv_ = 0.0;
    tolSigSet_ = 0.0;
    tolMaxLevSet_ = 0.0;

    IteratorType endit = discFuncSpace_->end();
    for (IteratorType it = discFuncSpace_->begin(); it != endit; ++it)
    {     
      double help = getLocalIndicator(*it);
      if( help < localTolerance_ * alphaSigSet_){
	      tolSigSetInv_ += help;
      }
      if ( it->level() == maxlevel){
	      tolMaxLevSet_ += help;
	      numMaxLev_++;
      }
      else{
	      tolSigSet_ += help;
	      numSigSet_++;
      }
    }
    if( (numSigSet_ > 1) && (tolMaxLevSet_ + tolSigSetInv_) < (localInTimeTolerance_ * 0.5) ){
      std::cout << "    CalcSigSet: AVOID MAXLEVEL REFINEMENT ... " << "\n";
      maxLevFlag_ = 0;
      //numSigSet_ += numMaxLev_;
      tolSigSetInv_ += tolMaxLevSet_;
    }
    else{
      maxLevFlag_ = 1;
      numSigSet_ += numMaxLev_;
      tolSigSet_ += tolMaxLevSet_;
    }
    // maxLevFlag_ = 1;
    return numSigSet_;
    double num = double(numberOfElements());
    maxLevAlpha_=10.;
    maxLevBeta_=(num-maxLevAlpha_*double(numMaxLev_))/(num-double(numMaxLev_));
    std::cout << "    ALPHA: " << " " << num << " " << numMaxLev_ << " " << maxLevAlpha_ << " " << maxLevBeta_ << std::endl;
    if (maxLevBeta_<0.25) {
      std::cout << "     No MaxLevel FIX!!!!!!!!" << std::endl;
      maxLevBeta_ = maxLevAlpha_ = 1.0;
    }
    assert(maxLevBeta_>0);
    return numSigSet_;
    
    /*
    typedef typename IndicatorDiscreteFunctionType::DofIteratorType DofIteratorType;

    numSigSet_ = 1;
    tolSigSetInv_ = 0.0;
    tolSigSet_ = 0.0;
    
    DofIteratorType endit = indicator_->dend();
    for (DofIteratorType it = indicator_->dbegin(); it != endit; ++it)
    {
       double help = *it;
       if( help < localTolerance_ * alphaSigSet_)
	 tolSigSetInv_ += help;
       else{
	 tolSigSet_ += help;
	 numSigSet_++;	 
       }
    }
    return numSigSet_;   
    */
  };

  double getLocalInTimeTolerance (){
    double dt = timeDiscrParam_.getTimeStepSize();
   
    localInTimeTolerance_ =  (1. - initialTheta_) * globalTolerance_ * dt / endTime_;

    return localInTimeTolerance_;
  };

  inline double getInitTolerance (){
    localInTimeTolerance_ = initialTheta_ * globalTolerance_;
    return (localInTimeTolerance_);
  };


  double getLocalTolerance (){

    getLocalInTimeTolerance ();
    localToleranceOrig_ = localInTimeTolerance_ / numberOfElements();
    localTolerance_ = localToleranceOrig_;

    int num = calcSignificantElements();

    std::cout << "   Estimator = " <<  tolSigSetInv_ + tolSigSet_ 
	      << "   Tol_Sig = " << (localInTimeTolerance_ - tolSigSetInv_) 
	      << "   Tol = " << localInTimeTolerance_ 
	      << "   Num Sig El: " << num 
	      << "   Num El: " <<  numberOfElements() << "\n";
   
    // localTolerance_  = (localInTimeTolerance_ - tolSigSetInv_)/double(num);
    localTolerance_ = localToleranceOrig_;

    return localTolerance_;
  };

  double getLocalCrsTolerance (){

    getLocalInTimeTolerance ();
    localToleranceOrig_ = localInTimeTolerance_ / numberOfElements();
    localTolerance_ = localToleranceOrig_;

    int num = calcSignificantElements();

    std::cout << "   Estimator = " <<  tolSigSetInv_ + tolSigSet_ 
	      << "   Tol_Sig = " << (localInTimeTolerance_ - tolSigSetInv_) 
	      << "   Tol = " << localInTimeTolerance_ 
	      << "   Num Sig El: " << num 
	      << "   Num El: " <<  numberOfElements() << "\n";
   
    localTolerance_  = (localInTimeTolerance_ - tolSigSetInv_)/double(num);
    // localTolerance_ = localToleranceOrig_;

    return localTolerance_;
  };

  double getInitLocalTolerance (){

    getInitTolerance ();
    localToleranceOrig_ = localInTimeTolerance_ / numberOfElements();

    int num = calcSignificantElements();

    std::cout << "   Estimator = " <<  tolSigSetInv_ + tolSigSet_ 
	      << "   Tol_Sig = " << (localInTimeTolerance_ - tolSigSetInv_) 
	      << "   Tol = " << localInTimeTolerance_ 
	      << "   Num Sig El: " << num 
	      << "   Num El: " <<  numberOfElements() << "\n";
    
    //localTolerance_  = (localInTimeTolerance_ - tolSigSetInv_)/double(num);
    localTolerance_ = localToleranceOrig_;

    return localTolerance_;
  };

  void markNeighbours (){
    typedef typename IndicatorDiscFuncSpaceType::IteratorType IteratorType;
    typedef typename IndicatorDiscFuncSpaceType::GridType      GridType;
    typedef typename GridType::template Codim<0>::Entity      EntityType;
    typedef typename EntityType::IntersectionIterator         IntersectionIteratorType;
    
    int maxlevel = grid_.maxLevel();

    IteratorType endit = discFuncSpace_->end();
    for (IteratorType it = discFuncSpace_->begin(); it != endit; ++it)
    {     
      if( getLocalIndicator(*it) > localTolerance_ ){
	IntersectionIteratorType endnit = it->iend();
	for(IntersectionIteratorType nit = it->ibegin(); nit != endnit; ++nit){
	  if( nit.neighbor() && (nit.outside().level() >= it->level()  ) ){
	    grid_.mark(1,nit.outside());
	  }
	}
      }
    }

    return;
    
  };  

  void markEntities (){
    typedef typename IndicatorDiscFuncSpaceType::IteratorType IteratorType;
    
    getLocalTolerance();

    IteratorType endit = discFuncSpace_->end();
    for (IteratorType it = discFuncSpace_->begin(); it != endit; ++it)
    {
      // std::cout << " ind  tol: " << getLocalIndicator(*it) << "  " <<  localTolerance_ << std::endl;
      if( getLocalIndicator(*it) > localTolerance_ )
        grid_.mark(1, it);
      else if ( (getLocalIndicator(*it) < coarsenTheta_ * localTolerance_) )
      	grid_.mark(-1, it);
      else
	grid_.mark(0, it);
    }

    return;
    
  };

  void markInitRefineEntities (){
    typedef typename IndicatorDiscFuncSpaceType::IteratorType IteratorType;
 
    int maxlevel = grid_.maxLevel();
   
    getInitLocalTolerance();

    IteratorType endit = discFuncSpace_->end();
    for (IteratorType it = discFuncSpace_->begin(); it != endit; ++it)
    {
      //std::cout << " ind  tol: " << getLocalIndicator(*it) << "  " <<  localTolerance_ << std::endl;
      double help = getLocalIndicator(*it);
      if( (help > localTolerance_) ){
	if ( maxLevFlag_)
	  grid_.mark(1, it);
	else if (it->level() < maxlevel)
	  grid_.mark(1, it);
      }
      else
	grid_.mark(0, it);
    }

    return;
    
  };

  void markRefineEntities (){
    typedef typename IndicatorDiscFuncSpaceType::IteratorType IteratorType;
 
    int maxlevel = grid_.maxLevel();
   
    getLocalTolerance();
    std::cout << "    LocalRef: MARK FOR REFINEMENT ...   MaxLevel = " << maxlevel << "\n";

    IteratorType endit = discFuncSpace_->end();
    for (IteratorType it = discFuncSpace_->begin(); it != endit; ++it)
    {
      //std::cout << " ind  tol: " << getLocalIndicator(*it) << "  " <<  localTolerance_ << std::endl;
      double help = getLocalIndicator(*it);
      if( (help > localTolerance_) ){
	      if ( maxLevFlag_)
	        grid_.mark(1, it);
	      else if (it->level() < maxlevel)
	        grid_.mark(1, it);
      }
      /*
      if (it->level() == maxlevel && help > localTolerance_*maxLevAlpha_)
        grid_.mark(1,it);
      else if (it->level()<maxlevel && help > localTolerance_*maxLevBeta_)
        grid_.mark(1,it);
      else
	      grid_.mark(0, it);
        */
    }
   
    return;
    
  };


  void markRefineEntities_Old (){
    typedef typename IndicatorDiscFuncSpaceType::IteratorType IteratorType;
     
    double tolMaxlevel = 0.0;
    double tolRest = 0.0;
    int maxlevel = grid_.maxLevel();
    int num = 0;

    getLocalTolerance();
    std::cout << "    LocalRef: MARK FOR REFINEMENT ...   MaxLevel = " << maxlevel << "\n";

    IteratorType endit = discFuncSpace_->end();
    for (IteratorType it = discFuncSpace_->begin(); it != endit; ++it)
    {
      //std::cout << " ind  tol: " << getLocalIndicator(*it) << "  " <<  localTolerance_ << std::endl;
      double help = getLocalIndicator(*it);
      if( (help > localTolerance_) ){
	      if ( it->level() == maxlevel){
	        tolMaxlevel += help;
	        num++;
	      }
	      else
	        tolRest += help;
        grid_.mark(1, it);
      }
      else
	      grid_.mark(0, it);
    }

    if( ( (tolMaxlevel + tolSigSetInv_ + (tolSigSet_ - tolMaxlevel)*0.25) < localInTimeTolerance_) ){
      std::cout << "    LocalRef: AVOID MAXLEVEL REFINEMENT ... " << "\n";
      IteratorType endit = discFuncSpace_->end();
      for (IteratorType it = discFuncSpace_->begin(); it != endit; ++it)
	    {
	      if ( it->level() == maxlevel)
	        grid_.mark(0, it);
	      } 
      }
    // else
    //  markNeighbours ();
    return;
    
  };


  template <class DFType,class RhoType>
  void markCoarsenEntities (DFType& df,RhoType& rho) {
    typedef typename IndicatorDiscFuncSpaceType::IteratorType IteratorType;
    int marked=0,wouldhave=0;
    getLocalCrsTolerance();

    IteratorType endit = discFuncSpace_->end();
    for (IteratorType it = discFuncSpace_->begin(); it != endit; ++it)
    {
      double localInd = getLocalIndicator(*it);
      if ( it->level()>0 && (localInd < coarsenTheta_ * localTolerance_) ) {
	RangeType cerror(0);
	bool leafs = CoarseningError<DFType>::compute(df,*(it->father()),cerror);
	if (leafs) {
	  double cInd = cerror.one_norm() * rho.localFunction(*it)[0] ;
	  if (localInd+cInd < coarsenTheta_ * localTolerance_) {
	    grid_.mark(-1, it);
	    marked++;
	  }
	  else {
	    grid_.mark(0, it);
	    wouldhave++;
	  }
	}
      }
      else
	grid_.mark(0, it);
    }
    std::cout << "     HCOARSEN: marked = " << marked << " from = " << wouldhave+marked << std::endl;
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
    typedef AdaptOperator<GridType,RestProlImp>              AdaptOpImp;
    
    RestProlImp *restProl;
    AdaptOpImp  *adaptOp;

    restProl = new RestProlImp ( *func);
 
    adaptOp  = new AdaptOpImp( grid_  , *restProl );

    adaptMapping_ = (*adaptOp);
  }


  template <class DiscFuncImp,class DiscFuncImp2>
 void addAdaptiveFunction(DiscFuncImp *func, DiscFuncImp2 *func2){
    typedef RestProlOperator<DiscFuncImp>                RestProlImp; 
    typedef RestProlOperator<DiscFuncImp2>               RestProlImp2; 
    typedef AdaptOperator<GridType,RestProlImp>          AdaptOpImp;
    typedef AdaptOperator<GridType,RestProlImp2>         AdaptOpImp2;
 
    
    RestProlImp *restProl;
    AdaptOpImp  *adaptOp;
    
    RestProlImp2 *restProl2;
    AdaptOpImp2  *adaptOp2;


    restProl  = new RestProlImp  ( *func  );
    restProl2 = new RestProlImp2 ( *func2 );

   
    adaptOp   = new AdaptOpImp  ( grid_ , *restProl );
    adaptOp2  = new AdaptOpImp2 ( grid_ , *restProl2 );


    adaptMapping_ = (*adaptOp) + (*adaptOp2);

  };


 template <class DiscFuncImp,  class DiscFuncImp2, class DiscFuncImp3,
	   class DiscFuncImp4, class DiscFuncImp5, class DiscFuncImp6, class DiscFuncImp7, class DiscFuncImp8>
 void addAdaptiveFunction(DiscFuncImp *func, DiscFuncImp2 *func2, DiscFuncImp3 *func3, 
			  DiscFuncImp4 *func4, DiscFuncImp5 *func5, DiscFuncImp6 *func6, DiscFuncImp7 *func7, DiscFuncImp8 *func8){
  
    typedef RestProlOperator<IndicatorDiscreteFunctionType>   RestProlImp0;   
    typedef RestProlOperator<DiscFuncImp>                     RestProlImp;   
    typedef RestProlOperator<DiscFuncImp2>                    RestProlImp2; 
    typedef RestProlOperator<DiscFuncImp3>                    RestProlImp3; 
    typedef RestProlOperator<DiscFuncImp4>                    RestProlImp4; 
    typedef RestProlOperator<DiscFuncImp5>                    RestProlImp5; 
    typedef RestProlOperator<DiscFuncImp6>                    RestProlImp6; 
    typedef RestProlOperator<DiscFuncImp7>                    RestProlImp7; 
    typedef RestProlOperator<DiscFuncImp8>                    RestProlImp8; 
 
    typedef AdaptOperator<GridType,RestProlImp0> AdaptOpImp0;
    typedef AdaptOperator<GridType,RestProlImp> AdaptOpImp;
    typedef AdaptOperator<GridType,RestProlImp2> AdaptOpImp2;
    typedef AdaptOperator<GridType,RestProlImp3> AdaptOpImp3;
    typedef AdaptOperator<GridType,RestProlImp4> AdaptOpImp4;
    typedef AdaptOperator<GridType,RestProlImp5> AdaptOpImp5;
    typedef AdaptOperator<GridType,RestProlImp6> AdaptOpImp6;
    typedef AdaptOperator<GridType,RestProlImp7> AdaptOpImp7;
    typedef AdaptOperator<GridType,RestProlImp8> AdaptOpImp8;
 
    
    RestProlImp0 *restProl0;
    RestProlImp  *restProl;
    RestProlImp2 *restProl2;
    RestProlImp3 *restProl3;
    RestProlImp4 *restProl4;
    RestProlImp5 *restProl5;
    RestProlImp6 *restProl6;
    RestProlImp7 *restProl7;
    RestProlImp8 *restProl8;


    AdaptOpImp0  *adaptOp0;
    AdaptOpImp   *adaptOp;
    AdaptOpImp2  *adaptOp2;
    AdaptOpImp3  *adaptOp3;
    AdaptOpImp4  *adaptOp4;
    AdaptOpImp5  *adaptOp5;
    AdaptOpImp6  *adaptOp6;
    AdaptOpImp7  *adaptOp7;
    AdaptOpImp8  *adaptOp8;


    restProl0 = new RestProlImp0 ( *indicator_  );
    restProl  = new RestProlImp  ( *func  );
    restProl2 = new RestProlImp2 ( *func2 );
    restProl3 = new RestProlImp3 ( *func3 );
    restProl4 = new RestProlImp4 ( *func4 );
    restProl5 = new RestProlImp5 ( *func5 );
    restProl6 = new RestProlImp6 ( *func6 );
    restProl7 = new RestProlImp7 ( *func7 );
    restProl8 = new RestProlImp8 ( *func8 );
   
    adaptOp0  = new AdaptOpImp0( grid_ , *restProl0  );
    adaptOp   = new AdaptOpImp ( grid_ , *restProl  );
    adaptOp2  = new AdaptOpImp2( grid_ , *restProl2 );
    adaptOp3  = new AdaptOpImp3( grid_ , *restProl3 );
    adaptOp4  = new AdaptOpImp4( grid_ , *restProl4 );
    adaptOp5  = new AdaptOpImp5( grid_ , *restProl5 );
    adaptOp6  = new AdaptOpImp6( grid_ , *restProl6 );
    adaptOp7  = new AdaptOpImp7( grid_ , *restProl7 );
    adaptOp8  = new AdaptOpImp8( grid_ , *restProl8 );



    adaptMapping_ = (*adaptOp0 ) 
      + (*adaptOp ) + (*adaptOp2) + (*adaptOp3) + (*adaptOp4) 
      + (*adaptOp5) + (*adaptOp6) + (*adaptOp7) + (*adaptOp8);

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
  double localToleranceOrig_;
  double localTolerance_;
  double coarsenTheta_;
  double initialTheta_;

  double alphaSigSet_;
  double tolSigSet_;
  double tolSigSetInv_;
  double tolMaxLevSet_;
  int    numSigSet_;
  int    numMaxLev_;
  int    maxLevFlag_;
  double maxLevAlpha_;
  double maxLevBeta_;
    
  //! timestep size in time discretization parameters und endTime 
  TimeDiscrParamType & timeDiscrParam_;
  double endTime_;


  //! colloect all adaptive mappings for adaption
  AdaptMappingType adaptMapping_;

 };

} // end namespace 

#endif
