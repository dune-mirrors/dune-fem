#ifndef DUNE_FEM_SPACE_RANNACHERTUREK_DOFMAPPERCODE_HH
#define DUNE_FEM_SPACE_RANNACHERTUREK_DOFMAPPERCODE_HH

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>

// dune-fem includes
#include <dune/fem/space/mapper/code.hh>
#include <dune/fem/space/mapper/compile.hh>
#include <dune/fem/space/mapper/indexsetdofmapper.hh>


namespace Dune
{

  namespace Fem
  {

    // RannacherTurekBlockMapperSingletonKey
    // -------------------------------------

    template< class GridPart >
    class RannacherTurekBlockMapperSingletonKey
    {
      typedef RannacherTurekBlockMapperSingletonKey< GridPart > ThisType;

    public:
      typedef GridPart GridPartType;

      RannacherTurekBlockMapperSingletonKey ( const GridPartType &gridPart )
      : gridPart_( gridPart )
      {}

      const GridPartType &gridPart () const { return gridPart_; }

      bool operator== ( const ThisType &other ) const
      {
        return ( &gridPart_ == &(other.gridPart_) );
      }

      bool operator!= ( const ThisType &other ) const
      {
        return !( *this == other );
      }

    private:
      const GridPartType &gridPart_;
    };



    // RannacherTurekDofMapperCodeFactory
    // ----------------------------------

    template< class LocalCoefficients >
    struct RannacherTurekDofMapperCodeFactory
    {
      typedef LocalCoefficients LocalCoefficientsType;

      RannacherTurekDofMapperCodeFactory ( const LocalCoefficients &localCoefficients = LocalCoefficients() )
      : localCoefficients_( localCoefficients )
      {}

      template< class RefElement,
                std::enable_if_t< std::is_same< std::decay_t< decltype( std::declval< const RefElement & >().size( 0 ) ) >, int >::value, int > = 0,
                std::enable_if_t< std::is_same< std::decay_t< decltype( std::declval< const RefElement & >().type( 0, 0 ) ) >, GeometryType >::value, int > = 0 >
      Dune::Fem::DofMapperCode operator() ( const RefElement &refElement ) const
      {
        if ( refElement.type( 0, 0 ).isCube() )
          return Dune::Fem::compile( refElement, localCoefficients() );
        else
          return Dune::Fem::DofMapperCode();
      }

      const LocalCoefficientsType &localCoefficients () const
      {
        return localCoefficients_;
      }

    private:
      LocalCoefficientsType localCoefficients_;
    };



    // RannacherTurekBlockMapperFactory
    // --------------------------------

    template< class GridPart, class LocalCoefficients >
    struct RannacherTurekBlockMapperFactory
    {
      typedef GridPart GridPartType;
      typedef LocalCoefficients LocalCoefficientsType;

      typedef Dune::Fem::IndexSetDofMapper< GridPartType > BlockMapperType;

      static BlockMapperType *createObject ( const RannacherTurekBlockMapperSingletonKey< GridPartType > &key )
      {
        return createObject( key.gridPart() );
      }

      static BlockMapperType *createObject ( const GridPart &gridPart )
      {
        RannacherTurekDofMapperCodeFactory< LocalCoefficientsType > codeFactory;
        return new BlockMapperType( gridPart, codeFactory );
      }

      static void deleteObject ( BlockMapperType *blockMapper )
      {
        delete blockMapper;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_RANNACHERTUREK_DOFMAPPERCODE_HH
