#ifndef DUNE_FEM_QUADRATURE_GEOMETRIC_GEOMETRY_HH
#define DUNE_FEM_QUADRATURE_GEOMETRIC_GEOMETRY_HH

#include <cstddef>

#include <functional>

#include <dune/geometry/type.hh>

#include "quadrature.hh"

namespace Dune
{

  namespace Fem
  {

    // GeometryQuadrature
    // ------------------

    template< class QuadratureRule >
    class GeometryQuadrature
    : public GeometricQuadrature< typename QuadratureRule::CoordType, QuadratureRule::d, QuadratureRule::d, GeometryQuadrature< QuadratureRule > >
    {
      using ThisType = GeometryQuadrature< QuadratureRule >;
      using BaseType = GeometricQuadrature< typename QuadratureRule::CoordType, QuadratureRule::d, QuadratureRule::d, GeometryQuadrature< QuadratureRule > >;

    public:
      /** \brief type of <tt>Dune::Geometry::QuadratureRule</tt> */
      using QuadratureRuleType = QuadratureRule;

      /** \copydoc Dune::Fem::GeometricQuadrature::FieldType */
      using FieldType = typename BaseType::FieldType;

      /** \copydoc Dune::Fem::GeometricQuadrature::CoordinateType */
      using CoordinateType = typename BaseType::CoordinateType;
      /** \copydoc Dune::Fem::GeometricQuadrature::LocalCoordinateType */
      using LocalCoordinateType = typename BaseType::LocalCoordinateType;

      /** \name Construction
       *  \{
       */

      explicit GeometryQuadrature ( const QuadratureRuleType &quadratureRule )
        : quadratureRule_( quadratureRule )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      /** \brief copy constructor */
      GeometryQuadrature ( const ThisType & ) = default;

      /** \brief move constructor */
      GeometryQuadrature ( ThisType && ) = default;

      /** \brief assignment operator */
      GeometryQuadrature &operator= ( const ThisType & ) = default;

      /** \brief move assignment operator */
      GeometryQuadrature &operator= ( ThisType && ) = default;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::GeometricQuadrature::type */
      Dune::GeometryType type () const { return quadratureRule().type(); }

      /** \copydoc Dune::Fem::GeometricQuadrature::order */
      int order () const { return quadratureRule().order(); }

      /** \copydoc Dune::Fem::GeometricQuadrature::nop */
      std::size_t nop () const { return quadratureRule().size(); }

      /** \copydoc Dune::Fem::GeometricQuadrature::point */
      const CoordinateType &point ( std::size_t i ) const
      {
        return quadratureRule()[ i ].position();
      }

      /** \copydoc Dune::Fem::GeometricQuadrature::localPoint */
      const LocalCoordinateType &localPoint ( std::size_t i ) const
      {
        return point( i );
      }

      /** \copydoc Dune::Fem::GeometricQuadrature::weight */
      FieldType weight ( std::size_t i ) const
      {
        return quadratureRule()[ i ].weight();
      }

      /** \} */

    private:
      const QuadratureRuleType &quadratureRule () const { return quadratureRule_.get(); }

      std::reference_wrapper< const QuadratureRuleType > quadratureRule_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_QUADRATURE_GEOMETRIC_GEOMETRY_HH
