#ifndef DUNE_FEM_REDUCEDBASISSPACE_MAPPER_HH
#define DUNE_FEM_REDUCEDBASISSPACE_MAPPER_HH

#include <dune/fem/space/common/dofmapperinterface.hh>

namespace Dune
{

  template< class GridPartImp, class BaseFunctionListImp >
  class ReducedBasisMapper;



  template< class GridPartImp, class BaseFunctionListImp >
  struct ReducedBasisMapperTraits
  {
    typedef GridPartImp GridPartType;
    
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType :: Entity
      EntityType;
    
    typedef BaseFunctionListImp BaseFunctionListType;

    typedef ReducedBasisMapper< GridPartType, BaseFunctionListType >
      DofMapperType;
    
    typedef DefaultDofMapIterator< GridPartImp, DofMapperType >
      DofMapIteratorType;
  };



  /** \class ReducedBasisMapper
   *  \brief provides the mapper for the reduced basis space 
   *
   *  This mapper just performs the identity mapping.
   */
  template< class GridPartImp, class BaseFunctionListImp >
  class ReducedBasisMapper
  : public DofMapperDefault
    < ReducedBasisMapperTraits< GridPartImp, BaseFunctionListImp > >
  {
  public:
    typedef ReducedBasisMapperTraits< GridPartImp, BaseFunctionListImp >
      Traits;
   
    typedef typename Traits :: BaseFunctionListType BaseFunctionListType;
    
    typedef typename Traits :: GridPartType GridPartType;
    typedef typename Traits :: EntityType EntityType;
    typedef typename Traits :: DofMapIteratorType DofMapIteratorType;

  private:
    typedef ReducedBasisMapper< GridPartType, BaseFunctionListType > ThisType;
    typedef DofMapperDefault< Traits > BaseType;

  protected:
    const BaseFunctionListType &baseFunctionList_;

  public:
    inline explicit ReducedBasisMapper ( const BaseFunctionListType &baseFunctionList )
    : baseFunctionList_( baseFunctionList )
    {
    }

    /** \copydoc Dune::DofMapperInterface::begin(const EntityType &entity) const */
    inline DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: beginIterator, entity, *this );
    }
    
    /** \copydoc Dune::DofMapperInterface::end(const EntityType &entity) const */
    inline DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: endIterator, entity, *this );
    }

    /** \copydoc Dune::DofMapperInterface::mapToGlobal(const EntityType &entity,int localDof) const */
    int mapToGlobal ( const EntityType &entity, int localDof ) const
    {
      return localDof;
    }

    bool needsCompress () const
    {
      return false;
    }

    int newIndex ( const int hole, const int block ) const
    {
      return -1;
    }

    int newSize () const
    {
      return size();
    }

    int numberOfHoles ( const int block ) const
    {
      return 0;
    }

    int numDofs () const
    {
      return size();
    }

    int oldIndex ( const int hole, const int block ) const
    {
      return -1;
    }

    int size () const
    {
      return baseFunctionList_.size();
    }
  };
  
};

#endif
