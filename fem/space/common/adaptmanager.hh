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
  typedef AdaptationManager<GridType,RestProlOperatorImp> MyType;
  typedef DofManager< GridType > DofManagerType; 
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;
public:
  typedef typename GridType :: Traits :: LocalIdSet LocalIdSet;
  //! create DiscreteOperator with a LocalOperator 
  AdaptationManager (GridType & grid, RestProlOperatorImp & rpOp) 
    : grid_(grid) 
    , dm_ ( DofManagerFactoryType::getDofManager(grid_) )
    , rpOp_ (rpOp) 
  {}

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
  
  //! adapt defines the grid walkthrough before and after grid adaptation.
  //! Note that the LocalOperator can be an combined Operator 
  //! Domain and Range are defined through class Operator
  void adapt () const 
  {
    bool restr = grid_.preAdapt();  

    if(restr)
    {
      dm_.resizeForRestrict();
      
      typedef typename DofManagerType :: IndexSetRestrictProlongType IndexSetRPType;
      typedef CombinedRestProl <IndexSetRPType,RestProlOperatorImp> COType;
      COType tmpop ( dm_.indexSetRPop() , rpOp_ );
      //typedef CombineInterface<RestrictProlongPair,IndexSetRPType&,RestProlOperatorImp&> COType;
      //COType tmpop ( dm_.indexSetRPop() , rpOp_ );

      typedef typename GridType::template Codim<0>::LevelIterator LevelIterator;

      // make run through grid 
      for(int l=0; l<grid_.maxLevel(); l++)
      {
        LevelIterator endit  = grid_.template lend<0>   ( l );
        for(LevelIterator it = grid_.template lbegin<0> ( l );
              it != endit; ++it )
        {
          hierarchicRestrict( *it , tmpop );
        }
      }
    }
    
    bool ref = grid_.adapt();

    if(ref)
    {
      dm_.resize();
      typedef typename DofManagerType :: IndexSetRestrictProlongType IndexSetRPType;
      typedef CombinedRestProl <IndexSetRPType,RestProlOperatorImp> COType;
      COType tmpop ( dm_.indexSetRPop() , rpOp_ );
      
      typedef typename GridType::template Codim<0>::LevelIterator LevelIterator;

      // make run through grid 
      LevelIterator endit = grid_.template lend<0> ( 0 );
      for(LevelIterator it = grid_.template lbegin<0> ( 0 );
          it != endit; ++it )
      {
        hierarchicProlong( *it , tmpop );
      }
    }

    // if grid was coarsend or refined, do dof compress 
    if(restr || ref)
      dm_.dofCompress();

    // do cleanup 
    grid_.postAdapt();
  }
  
private:
  // make hierarchic walk trough 
  template <class EntityType, class RestrictOperatorType  >
  void hierarchicRestrict ( EntityType &en, RestrictOperatorType & restop ) const 
  {
    if(!en.isLeaf())
    {
      typedef typename EntityType::HierarchicIterator HierarchicIterator; 
      HierarchicIterator it    = en.hbegin( en.level() + 1 );

      // if the children have children then we have to go deeper 
      HierarchicIterator endit = en.hend  ( en.level() + 1 );
     
      // ok because we checked en.isLeaf 
      if(!it->isLeaf()) return;
      
      // true for first child, otherwise false 
      bool initialize = true;
      
      for( ; it != endit; ++it)
      {
        EntityType & son = *it;
        if( son.mightBeCoarsened() )
        {
          restop.restrictLocal( en , son, initialize);     
          initialize = false;
        }
      }
    }
  }

  template <class EntityType, class ProlongOperatorType >
  void hierarchicProlong ( EntityType &en, ProlongOperatorType & prolop ) const 
  {
    typedef typename EntityType::HierarchicIterator HierarchicIterator; 
    
    bool initialize = true;
    
    HierarchicIterator endit  = en.hend  ( grid_.maxLevel() );
    for(HierarchicIterator it = en.hbegin( grid_.maxLevel() ); 
        it != endit; ++it)
    {
      assert( !en.isLeaf() );
      EntityType & son = *it; 
      if( son.wasRefined() )
      {
        prolop.prolongLocal( *(son.father()), son , initialize );     
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
};
/** @} end documentation group */
}

#endif
