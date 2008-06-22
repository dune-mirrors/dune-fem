#ifndef DUNE_ADAPTMANAGER_HH
#define DUNE_ADAPTMANAGER_HH

//- local includes 
#include <dune/common/timer.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/operator/common/objpointer.hh>

#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/loadbalancer.hh>

namespace Dune{

/** @addtogroup Adaptation Adaptation 
    Here the interfaces and algorithms for adapatation of a grid are
    described and implemented. 

    The strategy for restrict/prolong of data is chosen through the
    parameter
     \parametername \c fem.adaptation.method \n
                    method used for 
                    prolongation/restriction of data during grid
                    refinement;
                    defaults to generic \n
    Values are:
     - 0 == none: no adaptation is performed
     - 1 == generic: works on all grids 
     - 2 == callback: only AlbertaGrid and ALUGrid
     .
     
    \remarks 
    The Interface of Adaptation Manager is described by the class
    AdaptationManagerInterface.
    @{
 **/

/** \brief AdaptationManagerInterface class. 
 
 This Class is the result of a combination of different
 AdaptationOperators. It is the same principle as with Mapping. 
*/ 
class AdaptationManagerInterface : virtual public LoadBalancerInterface 
{
public:
  //! \brief default constructor 
  AdaptationManagerInterface () : am_ (0) {}

  //! \brief destructor 
  virtual ~AdaptationManagerInterface () {}
  
  /** \brief on call of this method the internal adaptation operator is
      called. 
  */
  virtual void adapt ()  
  {
    //std::cout << "called AdaptationManagerInterface::adapt()" << std::endl;
    if(am_) am_->adapt();  
    else 
    {
      std::cerr << "WARNING: AdaptationManagerInterface::adapt: no adaptation manager assigned! \n";
    }
  }

  /** \brief returns true if adaptation manager as adaptation method different to NONE
     \return \b true if adaptation method is not NONE, \b false otherwise
  */
  virtual bool adaptive () const  
  { 
    return (am_) ? (am_->adaptive()) : false; 
  } 

  /** \brief returns name of adaptation method 
     \return name of adaptation method 
  */
  virtual const char * methodName() const 
  {
    return (am_) ? (am_->methodName()) : "unknown method";
  }
    
  /** \brief Assignment operator, pointer to adaptation manager is stored 
      \return reference to this (i.e. *this) 
  */
  AdaptationManagerInterface & operator = (const AdaptationManagerInterface & am)
  {
      /** \todo This const-casting seems strange to me! */
    am_ = const_cast<AdaptationManagerInterface *> (&am);
    return (*this);
  }

  /** @copydoc LoadBalancerInterface::loadBalance */
  virtual bool loadBalance () 
  { 
    return (am_) ? (am_->loadBalance()) : false; 
  }

  /** @copydoc LoadBalancerInterface::balanceCounter */
  virtual int balanceCounter () const 
  { 
    return (am_) ? (am_->balanceCounter()) : 0; 
  }

  /** \brief time that last adaptation cycle took */
  virtual double adaptationTime () const 
  {
    return 0.0;
  }
 
private: 
  //! pointer to adaptation manager 
  AdaptationManagerInterface* am_;
};

//! deprecated typedef 
typedef AdaptationManagerInterface AdaptMapping;

/** \brief AdaptationMethod is a simple adaptation method reader class. */
template <class GridType>
class AdaptationMethod : virtual public AdaptationManagerInterface 
{  
public:  
  //! type of adaptation method 
  enum AdaptationMethodType { none = 0, //!< no adaptation is performed 
                              generic = 1, //!< a generic restriction and prolongation algorithm is used  
                              callback = 2 //!< the callback mechanism from AlbertaGrid and ALUGrid is used 
                            };
  
public:
  /** \brief constructor of AdaptationMethod 
     \param grid Grid that adaptation method is read for  
     \param paramFile optional parameter file which contains 
        the following two lines:
        # 0 == none, 1 == generic, 2 == call back (only AlbertaGrid and ALUGrid)  
        AdaptationMethod: 1 # default value 

  */   
  AdaptationMethod(const GridType & grid,
                   std::string paramFile = "") 
    : adaptationMethod_(generic)
  {
    const bool output = (grid.comm().rank() == 0);
    int am = 1;
    if( paramFile != "")
    {
      readParameter(paramFile,"AdaptationMethod",am, output);
    }
    else {
      am = Parameter::getValidValue<int>("fem.adaptation.method",am,
                      ValidateInterval<int,true,true>(0,2));
    }
    if(am == 2) adaptationMethod_ = callback;
    else if(am == 1) adaptationMethod_ = generic;
    else adaptationMethod_ = none;

    // for structred grid adaptation is disabled 
    if(! Capabilities::IsUnstructured<GridType>::v ) 
    {
      std::cerr << "WARNING: AdaptationMethod: adaptation disabled for structured grid! \n";
      adaptationMethod_ = none;
    }
      
    if( output )
    {
      std::cout << "Created AdaptationMethod: adaptation method = " << methodName() << std::endl;
    }
  }

  //! virtual destructor 
  virtual ~AdaptationMethod () {}

  /** @copydoc AdaptationManagerInterface::methodName */
  virtual const char * methodName() const 
  {
    switch (adaptationMethod_) {
      case generic: return "generic";
      case callback: return "callback";              
      case none: return "no adaptation";
      default:  return "unknown method";
    }
  }

  /** @copydoc AdaptationManagerInterface::adaptive */
  virtual bool adaptive () const { return adaptationMethod_ != none; }
  
protected:  
  //! method identifier 
  AdaptationMethodType adaptationMethod_;
};

/*! This could be seen as a hack, but is not
  With this method we define the class CombineRestProl which is only 
  for combining of local grid operations without using virtual methods.
  The interface of the two defined methods of the class (restrictLocal and
  prolongLocal) is given by the implementation (see below ) and 
  has to be the same for all local operators you want to combine.
*/
#define PARAM_CLASSNAME CombinedRestProl 
#define PARAM_FUNC_1 restrictLocal 
#define PARAM_FUNC_2 prolongLocal 
#define PARAM_FUNC_3 calcFatherChildWeight 
#define PARAM_FUNC_4 addToCommunicator 
#include <dune/fem/operator/common/combine.inc>

/*! \brief This class manages the adaptation process. 
 If the method adapt is called, then the grid is adapted and also 
 all the data belonging to the given dof manager will be rearranged 
 for data set where it is necessary to keep the data.
 */
template <class GridType, class RestProlOperatorImp >
class AdaptationManager 
: public Dune :: AdaptationMethod< GridType >,
  public ObjPointerStorage 
{  
  typedef Dune :: AdaptationMethod< GridType > BaseType;
  typedef typename BaseType :: AdaptationMethodType AdaptationMethodType; 
  
  template <class AdaptManager, class GridImp, bool isGoodGrid> 
  struct AdaptationMethod
  {
    template <class DofManagerImp, class RPOpImp>
    static void adapt(const AdaptManager& am, GridImp & grid, 
                      DofManagerImp& dm , RPOpImp& rpop,
                      AdaptationMethodType adaptMethod) 
    {
      // use generic adapt method 
      if( adaptMethod == BaseType :: generic ) 
      {
        am.template genericAdapt<All_Partition> ();
        return ;
      }
      
      // use grid call back adapt method 
      if( adaptMethod == BaseType :: callback ) 
      {
        grid.adapt(dm,rpop); 
        return ;
      }
    }
  };
  
  template <class AdaptManager, class GridImp> 
  struct AdaptationMethod<AdaptManager,GridImp,false>
  {
    template <class DofManagerImp, class RPOpImp>
    static void adapt(const AdaptManager& am, GridImp & grid, 
                      DofManagerImp& dm , RPOpImp& rpop,
                      AdaptationMethodType adaptMethod) 
    {
      // use generic adapt method 
      if(adaptMethod != BaseType :: none ) 
      {
        am.template genericAdapt<All_Partition> ();
        return ;
      }
    }
  };
  
  typedef AdaptationManager<GridType,RestProlOperatorImp> ThisType;
  typedef DofManager< GridType > DofManagerType; 
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;


public:
  typedef typename GridType :: Traits :: LocalIdSet LocalIdSet;

  /** \brief constructor of AdaptationManager 
     \param grid Grid that adaptation is done for 
     \param rpOp restriction and prlongation operator that describes how the 
      user data is projected to other grid levels
     \param paramFile optional parameter file which contains 
        the following two lines:
        # 0 == none, 1 == generic, 2 == call back (only AlbertaGrid and ALUGrid)  
        AdaptationMethod: 1 # default value 
  */   
  AdaptationManager (GridType & grid, RestProlOperatorImp & rpOp,
      std::string paramFile = "") 
    : BaseType(grid,paramFile)
    , grid_(grid) 
    , dm_ ( DofManagerFactoryType::getDofManager(grid_) )
    , rpOp_ (rpOp) 
    , adaptTime_(0.0)
  {
  }

  //! destructor 
  virtual ~AdaptationManager () {}

  /*! 
   Add to AdaptationManagers means that the RestProlOperators will be combined.
   See DiscreteOperatorImp.
   */
  template <class RestProlOperatorType> 
  AdaptationManager<GridType,
  CombinedRestProl <RestProlOperatorImp,RestProlOperatorType> > & 
  operator + (const AdaptationManager<GridType,RestProlOperatorType> &op)
  {
    //std::cout << "Operator + of AdaptationManager\n";
    typedef AdaptationManager<GridType,RestProlOperatorType> CopyType;
    typedef CombinedRestProl <RestProlOperatorImp,RestProlOperatorType> COType;
     
    COType *newRPOp = new COType ( rpOp_  , const_cast<CopyType &> (op).getRestProlOp() );

    typedef AdaptationManager < GridType, COType > OPType;
   
    OPType *adaptOp = new OPType ( grid_ , *newRPOp );    

    // memorize this new generated object because is represents this
    // operator and is deleted if this operator is deleted
    saveObjPointer( adaptOp , newRPOp );
   
    //std::cout << "Created " << adaptOp << " \n";
    return *adaptOp;
  }
 
  //! no public method, but has to be public, because all AdaptationManagers 
  //! must be able to call this method and the template parameters are
  //! allways different 
  RestProlOperatorImp & getRestProlOp ()  
  {
    return rpOp_;
  }
  
  /** \brief 
     according to adaption method parameter 
     the adaption procedure is done, 
     0 == no adaptation
     1 == generic adaption 
     2 == grid call back adaptation (only in AlbertaGrid and ALUGrid)
  */
  virtual void adapt ()
  {
    // get stopwatch 
    Timer timer; 

    AdaptationMethod<ThisType,GridType,
      Conversion<GridType,HasHierarchicIndexSet>::exists>::
        adapt(*this,grid_,dm_,rpOp_,this->adaptationMethod_);
    
    // take time 
    adaptTime_ = timer.elapsed();
  }

  //! \brief default load balancing method which does nothing
  virtual bool loadBalance () 
  { 
    return false; 
  }

  //! default load balancing counter is zero 
  virtual int balanceCounter () const 
  { 
    return 0; 
  }

  /** @copydoc AdaptationManagerInterface::adaptationTime */
  virtual double adaptationTime() const 
  {
    return adaptTime_;
  }

private:  
  /** \brief generic adaptation procedure
     adapt defines the grid walkthrough before and after grid adaptation.
     Note that the LocalOperator can be an combined Operator 
     Domain and Range are defined through class Operator
  */
  template <PartitionIteratorType pitype>
  void genericAdapt () const 
  {
    // call pre-adapt, returns true if at least 
    // one element is marked for coarsening 
    bool restr = grid_.preAdapt();  

    if(restr)
    {
      // get macro grid iterator 
      typedef typename GridType::template Codim<0>::
        template Partition<pitype> :: LevelIterator LevelIterator;

      // make a hierarchical to insert all elements 
      // that are father of elements that might be coarsened 
      {
        // get macro iterator 
        LevelIterator endit  = grid_.template lend<0,pitype>   ( 0 );
        for(LevelIterator it = grid_.template lbegin<0,pitype> ( 0 );
            it != endit; ++it )
        {
          hierarchicRestrict( *it , dm_.indexSetRPop() );
        }
      }

      // now resize memory 
      dm_.resizeForRestrict();

      // now project all data to fathers 
      {
        // get macro iterator 
        LevelIterator endit  = grid_.template lend<0,pitype>   ( 0 );
        for(LevelIterator it = grid_.template lbegin<0,pitype> ( 0 );
            it != endit; ++it )
        {
          hierarchicRestrict( *it , rpOp_ );
        }
      }
    }
    
    // adapt grid due to preset markers
    // returns true if at least one element was refined 
    bool ref = grid_.adapt();

    if(ref)
    {
      // resizes the index sets and resizes the memory
      dm_.resize();
      typedef typename DofManagerType :: IndexSetRestrictProlongType IndexSetRPType;
      typedef CombinedRestProl <IndexSetRPType,RestProlOperatorImp> COType;
      COType tmpop ( dm_.indexSetRPop() , rpOp_ );
      
      typedef typename GridType::template Codim<0>::
        template Partition<pitype> :: LevelIterator LevelIterator;

      // make run through grid to project data 
      LevelIterator endit = grid_.template lend<0,pitype> ( 0 );
      for(LevelIterator it = grid_.template lbegin<0,pitype> ( 0 );
          it != endit; ++it )
      {
        hierarchicProlong( *it , tmpop );
      }
    }

    // if grid was coarsend or refined, do dof compress 
    if(restr || ref)
    {
      // compress index sets and data 
      dm_.compress();
    }

    // do cleanup 
    grid_.postAdapt();
  }
  
private:
  //! make hierarchic walk trough for restriction 
  template <class EntityType, class RestrictOperatorType  >
  bool hierarchicRestrict ( EntityType& en, RestrictOperatorType & restop ) const 
  {
    if(!en.isLeaf())
    {
      // true means we are going to restrict data 
      bool doRestrict = true;
      
      // if the children have children then we have to go deeper 
      const int childLevel = en.level() + 1;
      typedef typename EntityType::HierarchicIterator HierarchicIterator; 
      
      // check all children first 
      {
        HierarchicIterator endit  = en.hend  ( childLevel );
        for(HierarchicIterator it = en.hbegin( childLevel ); it != endit; ++it)
        {
          doRestrict &= hierarchicRestrict( *it , restop );
        }
      }

      // if doRestrict is still true, restrict data 
      if(doRestrict)
      {
        // true for first child, otherwise false 
        bool initialize = true;
        HierarchicIterator endit  = en.hend  ( childLevel );
        for(HierarchicIterator it = en.hbegin( childLevel ); it != endit; ++it)
        {
          restop.restrictLocal( en , *it , initialize);     
          initialize = false;
        }
      }
    }
    // if all children return mightBeCoarsened,
    // then doRestrict on father remains true 
    return en.mightBeCoarsened();
  }

  template <class EntityType, class ProlongOperatorType >
  void hierarchicProlong ( EntityType &en, ProlongOperatorType & prolop ) const 
  {
    typedef typename EntityType::HierarchicIterator HierarchicIterator; 
    typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType; 
    
    // NOTE: initialize not working here
    // because we call hierarchically  
    
    // first call on this element 
    bool initialize = true;
    
    HierarchicIterator endit  = en.hend  ( grid_.maxLevel() );
    for(HierarchicIterator it = en.hbegin( grid_.maxLevel() ); 
        it != endit; ++it)
    {
      // should only get here on non-leaf entities 
      assert( !en.isLeaf() );

      // don't call this method on ghosts 
      assert( en.partitionType() != GhostEntity );
      
      EntityType & son = *it; 
      if( son.wasRefined() )
      {
        EntityPointerType vati = son.father();
        prolop.prolongLocal( *vati , son , initialize ); 
        initialize = false;
      }
    }
  }
 
protected:  
  //! corresponding grid 
  mutable GridType & grid_;

  //! DofManager corresponding to grid
  mutable DofManagerType & dm_;
  
  //! Restriction and Prolongation Operator 
  mutable RestProlOperatorImp & rpOp_;

  //! time that adaptation took 
  double adaptTime_;
};

/*! \brief This class manages the adaptation process including a load
  balancing after the adaptation step. 
*/
template <class GridType, class RestProlOperatorImp >
class AdaptationLoadBalanceManager :
  public AdaptationManager<GridType,RestProlOperatorImp> ,
  public LoadBalancer<GridType> 
{  
  typedef AdaptationManager<GridType,RestProlOperatorImp> BaseType; 
  typedef LoadBalancer<GridType> Base2Type;

  mutable CommunicationManagerList commList_;
  double balanceTime_ ;

  // do not copy 
  AdaptationLoadBalanceManager(const AdaptationLoadBalanceManager&);
public:  
  AdaptationLoadBalanceManager(
      GridType & grid, RestProlOperatorImp & rpOp, std::string paramFile = "", 
      int balanceCounter = 0 ) 
    : BaseType(grid,rpOp,paramFile) 
    , Base2Type( grid , rpOp , paramFile , balanceCounter )
    , commList_(rpOp)
  {
  }

  /** @copydoc LoadBalancerInterface::loadBalance */
  virtual bool loadBalance () 
  {
    // call load balance 
    return Base2Type :: loadBalance();
  }

  /** @copydoc LoadBalancerInterface::balanceCounter */ 
  virtual int balanceCounter () const 
  { 
    return Base2Type :: balanceCounter ();
  }

  /** @copydoc LoadBalancerInterface::loadBalanceTime */
  virtual double loadBalanceTime() const 
  {
    return balanceTime_;
  }
  
  /** @copydoc AdaptationManagerInterface::adapt */ 
  virtual void adapt () 
  {
    // adapt grid 
    BaseType :: adapt ();

    // if adaptation is enabled 
    if( this->adaptive() )
    {
      // get stopwatch 
      Timer timer; 
    
      // do load balancing 
      loadBalance ();

      // exchange all modified data 
      commList_.exchange();

      // get time  
      this->balanceTime_ = timer.elapsed();
    }
  }
};
/** @} end documentation group */

} // end namespace Dune 
#endif
