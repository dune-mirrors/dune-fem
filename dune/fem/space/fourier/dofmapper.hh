#ifndef DUNE_FEM_SPACE_FOURIER_DOFMAPPER_HH
#define DUNE_FEM_SPACE_FOURIER_DOFMAPPER_HH

#include <cstddef>

#include <dune/fem/space/mapper/dofmapper.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------

    template< class GridPart, int order > class FourierDofMapper;



    // FourierDofMapperTraits
    // ----------------------

    template< class GridPart, int order >
    struct FourierDofMapperTraits
    {
      typedef FourierDofMapper< GridPart, order > DofMapperType;

      typedef GridPart GridPartType;
      typedef typename GridPartType::template Codim< 0 >::EntityType ElementType;

      typedef std::size_t SizeType;
    };



    // FourierDofMapper
    // ----------------

    template< class GridPart, int order >
    class FourierDofMapper
    : public AdaptiveDofMapper< FourierDofMapperTraits< GridPart, order > >
    {
      typedef FourierDofMapper< GridPart, order > ThisType;
      typedef AdaptiveDofMapper< FourierDofMapperTraits< GridPart, order > > BaseType;

    public:
      typedef typename BaseType::Traits Traits;
      typedef typename BaseType::ElementType ElementType;
      typedef typename BaseType::GlobalKeyType GlobalKeyType;
      typedef typename BaseType::SizeType SizeType;


      /////////////////////////////////
      // DofMapper interface methods //
      /////////////////////////////////

      /** @copydoc Dune::Fem::DofMapper::size */
      static SizeType size () { return 1; }

      /** @copydoc Dune::Fem::DofMapper::contains */
      static bool contains ( int codim ) { return false; }

      /** @copydoc Dune::Fem::DofMapper::fixedDataSize */
      static bool fixedDataSize ( int codim ) { return true; }

      /** @copydoc Dune::Fem::DofMapper::mapEach */
      template< class Functor >
      static void mapEach ( const ElementType &element, Functor f )
      {
        f( 0, 0 );
      }

      /** @copydoc Dune::Fem::DofMapper::mapEachEntityDof */
      template< class Entity, class Functor >
      static void mapEachEntityDof ( const Entity &entity, Functor f )
      {}

      /** @copydoc Dune::Fem::DofMapper::maxNumDofs */
      static SizeType maxNumDofs () { return size(); }

      /** @copydoc Dune::Fem::DofMapper::numDofs */
      static SizeType numDofs ( const ElementType &element ) { return size(); }

      /** @copydoc Dune::Fem::DofMapper::numEntityDofs */
      template< class Entity >
      static SizeType numEntityDofs ( const Entity &entity )
      {
        return 0;
      }


      //////////////////////////////////////////
      // AdaptviveDofMapper interface methods //
      //////////////////////////////////////////

      /** @copydoc Dune::Fem::AdaptiveDofMapper::numberOfHoles */
      SizeType numberOfHoles ( const int block ) const { return SizeType( 0 ); }

      /** @copydoc Dune::Fem::AdaptiveDofMapper::oldIndex */
      GlobalKeyType oldIndex ( const int hole, const int block ) const
      {
        DUNE_THROW( Dune::NotImplemented, "Method oldIndex() not implemented yet" );
      }

      /** @copydoc Dune::Fem::AdaptiveDofMapper::newIndex */
      GlobalKeyType newIndex ( const int hole, const int block ) const
      {
        DUNE_THROW( Dune::NotImplemented, "Method newIndex() not implemented yet" );
      }

      /** @copydoc Dune::Fem::AdaptiveDofMapper::consecutive */
      bool consecutive () const { return true; }

      /** @copydoc Dune::Fem::AdaptiveDofMapper::oldOffSet */
      SizeType oldOffSet ( const int block ) const { return SizeType( 0 ); }

      /** @copydoc Dune::Fem::AdaptiveDofMapper::offSet */
      SizeType offSet ( const int block ) const { return SizeType( 0 ); }

      /** @copydoc Dune::Fem::AdaptiveDofMapper::numBlocks */
      SizeType numBlocks () const { return SizeType( 1 ); }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FOURIER_DOFMAPPER_HH
