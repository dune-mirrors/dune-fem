#ifndef DUNE_FVSPACEMAPPER_HH
#define DUNE_FVSPACEMAPPER_HH

//- Dune includes 
#include <dune/fem/space/common/dofmapper.hh>
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

  private:
    typedef DofMapperDefault< Traits > BaseType;

  public:
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

     /** \copydoc Dune::DofMapper::begin(const EntityType &entity) const */
    inline DofMapIteratorType begin ( const EntityType &entity ) const
    {
      const int baseIndex = indexSet_.index( entity ) * dimrange;
      typename DofMapIteratorType :: BeginIterator type;
      return DofMapIteratorType( type, baseIndex );
    }
    
    /** \copydoc Dune::DofMapper::end(const EntityType &entity) const */
    inline DofMapIteratorType end ( const EntityType &entity ) const
    {
      typename DofMapIteratorType :: EndIterator type; 
      return DofMapIteratorType( type, dimrange );
    }

    /** \copydoc Dune::DofMapper::mapToGlobal(const EntityType &entity,int localNum) const */
    inline int mapToGlobal ( const EntityType &entity, int localNum ) const
    {
      return indexSet_.index( entity ) * dimrange + localNum;
    }
    
    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    int maxNumDofs () const 
    {
      return dimrange;
    }

    /** \copydoc Dune::DofMapper::oldIndex(const int hole,const int block) const */
    int oldIndex ( const int hole, const int block ) const
    {
      assert( block == 0 );

      const int ishole = hole / dimrange;
      const int local = hole % dimrange;
      return dimrange * indexSet_.oldIndex( ishole, 0) + local;
    }

    /** \copydoc Dune::DofMapper::newIndex(const int hole,const int block) const */
    int newIndex ( const int hole, const int block ) const
    {
      assert( block == 0 );

      const int ishole = hole / dimrange;
      const int local = hole % dimrange;
      return dimrange * indexSet_.newIndex( ishole, 0) + local;
    }

    /** \copydoc Dune::DofMapper::numberOfHoles(const int block) const */
    int numberOfHoles ( const int block ) const
    {
      assert( block == 0 );
      return dimrange * indexSet_.numberOfHoles( 0 );
    }
    
    // is called once and calcs the insertion points too
    int newSize() const 
    {
      return this->size();
    }

    //! return the sets consecutive 
    bool consecutive () const
    {
      return BaseType::checkConsecutive(indexSet_);
    }
  };

} // end namespace Dune 

#endif
