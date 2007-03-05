#ifndef DUNE_ADAPTMANAGER_HH
#define DUNE_ADAPTMANAGER_HH

//- Dune includes 
#include <dune/common/array.hh>

//- local includes 
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/operator/common/objpointer.hh>
namespace Dune{

/** @defgroup Adaptation Adaptation
 * Concept for local grid refinement
 * @{
 **/


  
/** @defgroup AdaptationManagerImp Adaptation Manager
    @{
**/

/*! This could be seen as a hack, but is not
  With this method we define the class CombineRestProl which is only 
  for combining of local grid operations without using virtual methods.
  The interface of the two defined methods of the class (restrictLocal and
  prolongLocal) is given by the implementation (see below ) and 
  has to be the same for all local operators you want to combine 
*/
#define PARAM_CLASSNAME CombinedRestProl 
#define PARAM_FUNC_1 restrictLocal 
#define PARAM_FUNC_2 prolongLocal 
#define PARAM_FUNC_3 calcFatherChildWeight 
#include <dune/fem/operator/common/combine.hh>

/*! Combination of different AdaptationManagers

 This Class is the result of a combination of different
 AdaptationOperators. It is the same principle as with Mapping and
 DiscreteOperatorImp. 
*/ 
class AdaptationManagerInterface 
{
public:
  //! default constructor 
  AdaptationManagerInterface () : am_ (0) {}

  //! destructor 
  virtual ~AdaptationManagerInterface () {}
  
  //! all adaptation operators have this method which adapts the
  //! corresponding grid and organizes the restriction prolongation process
  //! of the underlying function spaces
  virtual void adapt () const 
  {
    //std::cout << "called AdaptationManagerInterface::adapt()" << std::endl;
    if(am_) am_->adapt();  
    else 
    {
      std::cerr << "WARNING: adapt! \n";
    }
  };

  //! returns true if adaptation manager does something during adapt call
  virtual bool adaptive () const  
  { 
    return (am_) ? (am_->adaptive()) : false; 
  } 

  virtual const char * methodName() const 
  {
    return (am_) ? (am_->methodName()) : "unknown method";
  }
    
  //! Assignement operator
  AdaptationManagerInterface & operator = (const AdaptationManagerInterface & am)
  {
      /** \todo This const-casting seems strange to me! */
    am_ = const_cast<AdaptationManagerInterface *> (&am);
    return (*this);
  }
private: 
  AdaptationManagerInterface *am_;
};

typedef AdaptationManagerInterface AdaptMapping;

/*! \brief This class manages the adaptation process. 
 If the method adapt is called, then the grid is adapted and also 
 all the data belonging to the given dof manager will be rearranged 
 for data set where it is necessary to keep the data.
 */
template <class GridType, class RestProlOperatorImp >
class AdaptationManager :
  public AdaptationManagerInterface , public ObjPointerStorage 
{  
  enum AdaptationMethodType { none = 0, generic = 1, callback = 2 };
  
  template <class AdaptManager, class GridImp, bool isGoodGrid> 
  struct AdaptationMethod
  {
    template <class DofManagerImp, class RPOpImp>
    static void adapt(const AdaptManager& am, GridImp & grid, 
                      DofManagerImp& dm , RPOpImp& rpop,
                      AdaptationMethodType adaptMethod) 
    {
      // use generic adapt method 
      if( adaptMethod == generic ) 
      {
        am.template genericAdapt<All_Partition> ();
        return ;
      }
      
      // use grid call back adapt method 
      if( adaptMethod == callback ) 
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
      if(adaptMethod != none ) 
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
  const char * methodName() const 
  {
    switch (adaptationMethod_) {
      case generic: return "generic";
      case callback: return "callback";              
      case none: return "no adaptation";
      default:  return "unknown method";
    }
  }

  typedef typename GridType :: Traits :: LocalIdSet LocalIdSet;
public:
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
    : grid_(grid) 
    , dm_ ( DofManagerFactoryType::getDofManager(grid_) )
    , rpOp_ (rpOp) 
    , adaptationMethod_(generic)
  {
    if( paramFile != "")
    {
      int am = 1;
      readParameter(paramFile,"AdaptationMethod",am);
      if(am == 2) adaptationMethod_ = callback;
      else if(am == 1) adaptationMethod_ = generic;
      else adaptationMethod_ = none;
    }

    // for structred grid adaptation is disabled 
    if(! Capabilities::IsUnstructured<GridType>::v ) 
    {
      std::cerr << "WARNING: AdaptationManager: adaptation disabled for structured grid! \n";
      adaptationMethod_ = none;
    }
      
    std::cout << "Created AdaptationManager: adaptation method = " << methodName() << std::endl;
  }

  virtual ~AdaptationManager () {}

  //! returns true if adaptation is enabled 
  bool adaptive () const { return adaptationMethod_ != none; }

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
  
  //! according to adaption method parameter 
  //! the adaption procedure is done, 
  //! 0 == no adaptation
  //! 1 == generic adaption 
  //! 2 == grid call back adaptation (only in AlbertaGrid and ALUGrid)
  void adapt () const 
  {
    AdaptationMethod<ThisType,GridType,
      Conversion<GridType,HasHierarchicIndexSet>::exists>::
        adapt(*this,grid_,dm_,rpOp_,adaptationMethod_);
  }
 
private:  
  //! generic adaptation procedure
  //! adapt defines the grid walkthrough before and after grid adaptation.
  //! Note that the LocalOperator can be an combined Operator 
  //! Domain and Range are defined through class Operator
  template <PartitionIteratorType pitype>
  void genericAdapt () const 
  {
    // call pre-adapt, returns true if at least 
    // one element is marked for coarsening 
    bool restr = grid_.preAdapt();  

    if(restr)
    {
      dm_.resizeForRestrict();
      
      typedef typename DofManagerType :: IndexSetRestrictProlongType IndexSetRPType;
      typedef CombinedRestProl <IndexSetRPType,RestProlOperatorImp> COType;
      COType tmpop ( dm_.indexSetRPop() , rpOp_ );

      typedef typename GridType::template Codim<0>::
        template Partition<pitype> :: LevelIterator LevelIterator;

      // make a hierarchical run through grid 
      {
        // get macro iterator 
        LevelIterator endit  = grid_.template lend<0,pitype>   ( 0 );
        for(LevelIterator it = grid_.template lbegin<0,pitype> ( 0 );
            it != endit; ++it )
        {
          hierarchicRestrict( *it , tmpop );
        }
      }
    }
    
    // adapt grid due to preset markers
    // returns true if at least one element was refined 
    bool ref = grid_.adapt();

    if(ref)
    {
      dm_.resize();
      typedef typename DofManagerType :: IndexSetRestrictProlongType IndexSetRPType;
      typedef CombinedRestProl <IndexSetRPType,RestProlOperatorImp> COType;
      COType tmpop ( dm_.indexSetRPop() , rpOp_ );
      
      typedef typename GridType::template Codim<0>::
        template Partition<pitype> :: LevelIterator LevelIterator;

      // make run through grid 
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
      dm_.dofCompress();
    }

    // do cleanup 
    grid_.postAdapt();
  }
  
private:
  // make hierarchic walk trough 
  template <class EntityType, class RestrictOperatorType  >
  bool hierarchicRestrict ( EntityType &en, RestrictOperatorType & restop ) const 
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
 
  //! corresponding grid 
  mutable GridType & grid_;

  //! DofManager corresponding to grid
  mutable DofManagerType & dm_;
  
  //! Restriction and Prolongation Operator 
  mutable RestProlOperatorImp & rpOp_;

  AdaptationMethodType adaptationMethod_;
};
/** @} end documentation group */
}

#endif
