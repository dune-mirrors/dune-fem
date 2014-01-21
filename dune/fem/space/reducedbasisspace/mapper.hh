#ifndef DUNE_FEM_REDUCEDBASISSPACE_MAPPER_HH
#define DUNE_FEM_REDUCEDBASISSPACE_MAPPER_HH

#include <dune/fem/space/mapper/dofmapper.hh>

namespace Dune
{

  namespace Fem 
  {

    template< class GridPart, class BaseFunctionList >
    class ReducedBasisMapper;


    template< class GridPart, class BaseFunctionList >
    struct ReducedBasisMapperTraits
    {
      typedef ReducedBasisMapper< GridPart, BaseFunctionList > DofMapperType;

      typedef GridPart GridPartType;
      typedef BaseFunctionList BaseFunctionListType;
      
      typedef typename GridPartType::template Codim< 0 >::EntityType ElementType;

      typedef int SizeType;
    };


    /** \class ReducedBasisMapper
     *  \brief provides the mapper for the reduced basis space 
     *
     *  This mapper just performs the identity mapping.
     */
    template< class GridPart, class BaseFunctionList >
    class ReducedBasisMapper
    : public DofMapper< ReducedBasisMapperTraits< GridPart, BaseFunctionList > >
    {
      typedef ReducedBasisMapper< GridPart, BaseFunctionList > ThisType;
      typedef DofMapper< ReducedBasisMapperTraits< GridPart, BaseFunctionList > > BaseType;

    public:
      typedef typename BaseType::Traits Traits;
      typedef typename Traits::SizeType SizeType;

      typedef typename BaseType::ElementType ElementType;
     
      typedef typename Traits::GridPartType GridPartType;
      typedef typename Traits::BaseFunctionListType BaseFunctionListType;

    public:
      explicit ReducedBasisMapper ( const BaseFunctionListType &baseFunctionList )
      : baseFunctionList_( baseFunctionList )
      {}

      /** \copydoc Dune::Fem::DofMapper::mapEach */
      template< class Functor >
      void mapEach ( const ElementType &element, Functor f ) const
      {
        const int n = size();
        for( int i = 0; i < n; ++i )
          f( i, i );
      }

      /** \copydoc Dune::Fem::DofMapper::mapToGlobal(const ElementType &entity,const int localDof) const */
      int mapToGlobal ( const ElementType &entity, const int localDof ) const
      {
        return localDof;
      }

      /** \copydoc Dune::Fem::DofMapper::mapEachEntityDof(const Entity &entity,Functor f) const */
      template< class Entity, class Functor >
      void mapEachEntityDof ( const Entity &entity, Functor f ) const
      {
        DUNE_THROW( NotImplemented, "ReducedBasisSpace cannot map entity DoFs." );
      }

      /** \copydoc Dune::Fem::DofMapper::numDofs(const ElementType &element) const */
      int numDofs ( const ElementType &element ) const
      {
        return size();
      }

      /** \copydoc Dune::Fem::DofMapper::maxNumDofs() const */
      int maxNumDofs () const
      {
        return size();
      }

      /** \copydoc Dune::Fem::DofMapper::numEntityDofs(const Entity &entity) const */
      template< class Entity >
      int numEntityDofs ( const Entity &entity ) const
      {
        DUNE_THROW( NotImplemented, "ReducedBasisSpace cannot map entity DoFs." );
        return 0;
      }
     
      /** \copydoc Dune::Fem::DofMapper::size() const */
      SizeType size () const
      {
        return baseFunctionList_.size();
      }

    protected:
      const BaseFunctionListType &baseFunctionList_;
    };

  } // namespace Fem
  
} // namespace Dune

#endif // #ifndef DUNE_FEM_REDUCEDBASISSPACE_MAPPER_HH
