#ifndef DUNE_FVSPACEMAPPER_HH
#define DUNE_FVSPACEMAPPER_HH

//- Dune includes 
#include <dune/fem/space/common/dofmapperinterface.hh>
#include <dune/fem/space/dgspace/dgmapper.hh>

namespace Dune
{
   
  /** \brief 
   The FVSpaceMapper maps local to global dofs for the FVSpace. 
    */
  template< class GridPartImp, int polOrd, int dimrange >
  class FiniteVolumeMapper;



  template< class GridPartImp, int polOrd, int dimrange >
  struct FiniteVolumeMapperTraits
  {
    typedef GridPartImp GridPartType;

    typedef typename GridPartType :: IndexSetType IndexSetType;

    typedef typename GridPartType :: template Codim< 0 > :: IteratorType :: Entity
      EntityType;

    typedef DGDofMapIterator DofMapIteratorType;

    typedef FiniteVolumeMapper< GridPartType, polOrd, dimrange >
      DofMapperType;
  };


  //*****************************************************************
  //
  // specialisation for polynom order 0 and arbitray dimrange 
  //
  //*****************************************************************
  template< class GridPartImp, int dimrange >
  class FiniteVolumeMapper< GridPartImp, 0, dimrange >
  : public DofMapperDefault< FiniteVolumeMapperTraits< GridPartImp, 0, dimrange > > 
  {
  public:
    typedef FiniteVolumeMapperTraits< GridPartImp, 0, dimrange > Traits;

    typedef typename Traits :: IndexSetType IndexSetType;

    typedef typename Traits :: EntityType EntityType;

    typedef typename Traits :: DofMapIteratorType DofMapIteratorType;

  protected:
    enum { numCodims = IndexSetType::ncodim };

    const IndexSetType &indexSet_;
      
  public:
    inline FiniteVolumeMapper ( const IndexSetType &indexSet,
                                int numDofs ) 
    : indexSet_( indexSet )
    {
      assert( numDofs == dimrange );
    }

    //! return size of function space, here number of elements 
    int size () const
    {
      return indexSet_.size( 0 ) * dimrange;
    }

     /** \copydoc Dune::DofMapperInterface::begin(const EntityType &entity) const */
    inline DofMapIteratorType begin ( const EntityType &entity ) const
    {
      const int baseIndex = indexSet_.index( entity ) * dimrange;
      typename DofMapIteratorType :: BeginIterator type;
      return DofMapIteratorType( type, baseIndex );
    }
    
    /** \copydoc Dune::DofMapperInterface::end(const EntityType &entity) const */
    inline DofMapIteratorType end ( const EntityType &entity ) const
    {
      typename DofMapIteratorType :: EndIterator type; 
      return DofMapIteratorType( type, dimrange );
    }

    /** \copydoc Dune::DofMapperInterface::mapToGlobal(const EntityType &entity,int localNum) const */
    inline int mapToGlobal ( const EntityType &entity, int localNum ) const
    {
      return indexSet_.index( entity ) * dimrange + localNum;
    }

    /** \copydoc Dune::DofMapperInterface::mapToGlobal(const EntityType &entity,int localNum) const 
        \note This method returns zero since the current FV
        implementation only supports p = 0! 
    */
    template <class EntityImp>
    inline int mapToGlobal ( const EntityImp &entity, int localNum ) const
    {
      return 0; 
    }

    //! return old index for hole 
    int oldIndex ( const int hole, int ) const
    {   
      // corresponding number of set is newn 
      const int newn  = static_cast<int> (hole/dimrange);
      // local number of dof is local 
      const int local = (hole % dimrange); 
      return (dimrange * indexSet_.oldIndex(newn,0)) + local;
    }

    //! return new index for hole 
    int newIndex (const int hole, int ) const
    {
      // corresponding number of set is newn 
      const int newn = static_cast<int> (hole / dimrange);
      // local number of dof is local 
      const int local = (hole % dimrange); 
      return (dimrange * indexSet_.newIndex(newn,0)) + local;
    }

    //! return numbber of exsiting hole 
    int numberOfHoles ( int ) const
    {   
      // this index set works only for codim = 0 at the moment
      return dimrange * indexSet_.numberOfHoles(0);
    }
    
    // is called once and calcs the insertion points too
    int newSize() const 
    {
      return this->size();
    }

    //! return number of dof per element 
    int numDofs () const 
    {
      return dimrange;
    }

    //! return the sets needsCompress 
    bool needsCompress () const { return indexSet_.needsCompress(); }
  };

#if 0
  template <class IndexSetImp>
  class FiniteVolumeMapper<IndexSetImp,0,1>
  : public DofMapperDefault < FiniteVolumeMapper <IndexSetImp,0,1> > 
  {
    // corresp. index set 
    const IndexSetImp & indexSet_;
  public:
    typedef IndexSetImp IndexSetType;
    
    FiniteVolumeMapper ( const IndexSetType  & is , int numDofs ) : indexSet_ (is) {}

    // we have virtual function ==> virtual destructor 
    virtual ~FiniteVolumeMapper () {}

    //! return size of function space, here number of elements 
    int size () const
    {
      return indexSet_.size(0);
    }
    
    //! map Entity an local Dof number to global Dof number 
    //! for polOrd = 0
    template <class EntityType>
    int mapToGlobal (EntityType &en, int localNum ) const
    {
      return indexSet_.template index<0> (en,localNum);
    }

    //! return old index of hole 
    int oldIndex (const int hole, int ) const
    {   
      return indexSet_.oldIndex(hole,0);
    }

    //! return new index of hole 
    int newIndex (const int hole, int ) const
    {
      return indexSet_.newIndex(hole,0);
    }

    //! return number of holes 
    int numberOfHoles ( int ) const
    {   
      // this index set works only for codim = 0 at the moment
      return indexSet_.numberOfHoles(0);
    }
    
    // is called once and calcs the insertion points too
    int newSize() const 
    {
      return this->size();
    }

    //! return number of dof per entity, here this method returns 1
    int numDofs () const 
    {
      return 1;
    }

    //! return the sets needsCompress 
    bool needsCompress () const { return indexSet_.needsCompress(); }
  };
#endif

} // end namespace Dune 

#endif
