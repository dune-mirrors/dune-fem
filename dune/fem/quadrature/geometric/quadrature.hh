#ifndef DUNE_FEM_QUADRATURE_GEOMETRIC_QUADRATURE_HH
#define DUNE_FEM_QUADRATURE_GEOMETRIC_QUADRATURE_HH

#include <cstddef>

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

#include <dune/fem/quadrature/quadrature.hh>

namespace Dune
{

  namespace Fem
  {

    // GeometricQuadrature
    // -------------------

    template< class Field, int mydim, int dim, class Implementation >
    class GeometricQuadrature
    {
      using ThisType = GeometricQuadrature< Field, mydim, dim, Implementation >;

    public:
      /** \brief field type */
      using FieldType = Field;

      /** \brief mydimension */
      static const int mydimension = mydim;
      /** \brief dimension */
      static const int dimension = dim;

      /** \brief coordinate type */
      using CoordinateType = Dune::FieldVector< FieldType, dimension >;
      /** \brief local coordinate type */
      using LocalCoordinateType = Dune::FieldVector< FieldType, mydimension >;

      /** \brief quadrature point wrapper type */
      using QuadraturePointWrapperType = Dune::Fem::QuadraturePointWrapper< Implementation >;

    protected:
#ifndef DOXYGEN
      GeometricQuadrature () = default;
#endif // #ifndef DOXYGEN

    public:
      /** \name Copying and assignment
       *  \{
       */

      /** \brief copy constructor */
      GeometricQuadrature ( const ThisType & ) = default;

      /** \brief move constructor */
      GeometricQuadrature ( ThisType && ) = default;

      /** \brief assignment operator */
      GeometricQuadrature &operator= ( const ThisType & ) = default;

      /** \brief move assignment operator */
      GeometricQuadrature &operator= ( ThisType && ) = default;

      /** \} */

      /** \name Query methods
       *  \{
       */

      /** \brief return geometry type */
      Dune::GeometryType type () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().type() );
        return impl().type();
      }

      /** \brief return order */
      int order () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().order() );
        return impl().order();
      }

      /** \} */

      /** \name Quadrature points and weights
       *  \{
       */

      /** \brief return number of quadrature points */
      std::size_t nop () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().nop() );
        return impl().nop();
      }

      /** \brief return coordinates of \f$i\f$-th quadrature point */
      const CoordinateType &point ( std::size_t i ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().point( i ) );
        return impl().point( i );
      }

      /** \brief return local coordinates of \f$i\f$-th quadrature point */
      const LocalCoordinateType &localPoint ( std::size_t i ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().localPoint( i ) );
        return impl().localPoint( i );
      }

      /** \brief return quadrature weight */
      FieldType weight ( std::size_t i ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().weight( i ) );
        return impl().weight( i );
      }

      /** \} */

      /** \name Quadrature point wrapper
       *  \{
       */

      /** \brief return quadrature point wrapper */
      const QuadraturePointWrapperType operator[] ( std::size_t i ) const
      {
        return QuadraturePointWrapperType( impl(), i );
      }

      /** \} */

    protected:
      const Implementation &impl () const
      {
        return static_cast< const Implementation & >( *this );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_QUADRATURE_GEOMETRIC_QUADRATURE_HH
