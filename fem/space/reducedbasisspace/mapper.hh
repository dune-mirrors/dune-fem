#ifndef DUNE_FEM_REDUCEDBASISSPACE_MAPPER_HH
#define DUNE_FEM_REDUCEDBASISSPACE_MAPPER_HH

#include <dune/fem/space/common/dofmapperinterface.hh>

namespace Dune
{

/*======================================================================*/
/*!
 *  \class ReducedBasisMapper
 *  \brief provides the mapper for the reduced basis space 
 *
 *  This class just performs the identity mapping.
 *
 */
/*======================================================================*/
  template< class BaseFunctionListImp >
  class ReducedBasisMapper
  : public DofMapperDefault< ReducedBasisMapper< BaseFunctionListImp > >
  {
  public:
    typedef BaseFunctionListImp BaseFunctionListType;

  private:
    typedef ReducedBasisMapper< BaseFunctionListType > ThisType;

  protected:
    const BaseFunctionListType &baseFunctionList_;

  public:
    inline explicit ReducedBasisMapper ( const BaseFunctionListType &baseFunctionList )
    : baseFunctionList_( baseFunctionList )
    {
    }

    template< class EntityType >
    int mapToGlobal ( const EntityType &entity, int localDof )
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
