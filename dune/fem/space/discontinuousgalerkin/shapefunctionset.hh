#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SHAPEFUNCTIONSET_HH

// C++ includes
#include <cstdlib>

// dune-common includes
#include <dune/common/static_assert.hh>

// dune-geometry includes
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/shapefunctionset/capabilities.hh>

// local includes 
#include "orthonormalbase_mod.hh"

/**
  @file
  @author Christoph Gersbacher
  @brief Provides orthonormal shape function set
*/


namespace Dune
{

  namespace Fem
  {

    // Forward declaration
    // -------------------

    template< class FunctionSpace, int polOrder >
    class OrthonormalShapeFunctionSet;



    namespace ShapeFunctionSetCapabilities
    {

      // Template specialization of hasStaticSize
      // ----------------------------------------

      template< class FunctionSpace, int polOrder >
      class hasStaticSize< OrthonormalShapeFunctionSet< FunctionSpace, polOrder > > 
      {
        template< int order, int dimDomain >
        struct NumShapeFunctions;

        template< int order >
        struct NumShapeFunctions< order, 1 >
        {
          static const int v = order + 1;
        };

        template< int order >
        struct NumShapeFunctions< order, 2 >
        {
          static const int v = (order + 2) * (order + 1) / 2;
        };

        template< int order >
        struct NumShapeFunctions< order, 3 >
        {
          static const int v = ((order+1)*(order+2)*(2*order+3)/6
                                 + (order+1)*(order+2)/2)/2;
        };

      public:
        static const bool v = true;
        static const int size = NumShapeFunctions< polOrder, FunctionSpace::dimDomain >::v;
      };

    } // namespace ShapeFunctionSetCapabilities



    // OrthonormalShapeFunctionHelper
    // ------------------------------

    template< class FunctionSpace >
    struct OrthonormalShapeFunctionHelper 
    {
      dune_static_assert( (FunctionSpace::dimRange == 1),
                          "FunctionSpace must be scalar (i.e., dimRange = 1)." );

      typedef FunctionSpace FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

    protected:
      // line
      typedef Dune::GenericGeometry::Prism< Dune::GenericGeometry::Point > Line;

      // quadrilateral
      typedef Dune::GenericGeometry::Prism< Line > Quadrilateral;
      // triangle
      typedef Dune::GenericGeometry::Pyramid< Line > Triangle;

      // hexahedron
      typedef Dune::GenericGeometry::Pyramid< Quadrilateral > Hexahedron;
      // prism
      typedef Dune::GenericGeometry::Prism< Triangle > Prism;
      // pyramid
      typedef Dune::GenericGeometry::Prism< Quadrilateral > Pyramid;
      // tetrahedron
      typedef Dune::GenericGeometry::Pyramid< Triangle > Tetrahedron;

      // make unique topology id
      template< unsigned int topologyId >
      struct UniqueId
      {
        static const unsigned int v = topologyId | (unsigned int)Dune::GenericGeometry::prismConstruction;
      };

    public:
      template< class Topology >
      struct Evaluate;

      template< class Topology >
      struct Jacobian;

      template< class Topology >
      struct Hessian;
    };



    // Implementation of OrthonormalShapeFunctionHelper::Evaluate
    // ----------------------------------------------------------

    template< class FunctionSpace >
    template< class Topology >
    struct OrthonormalShapeFunctionHelper< FunctionSpace >::Evaluate
    {
      template< class Functor >
      static void apply ( const DomainType &x, Functor functor )
      {
        const unsigned int id = UniqueId< Topology::id >::v;
        for( std::size_t i = 0; i < size; ++i )
          functor( i, evaluate( integral_constant< unsigned int, id >(), x ) );
      }

    protected:
      static RangeType evaluate ( const integral_constant< unsigned int, Line::id >(), const DomainType &x )
      {
        return OrthonormalBase_1D::eval_line( baseNum, &x[ 0 ] );
      }
      static RangeType evaluate ( const integral_constant< unsigned int, Quadrilateral::id >(), const DomainType &x )
      {
        return OrthonormalBase_2D::eval_quadrilateral_2d( baseNum, &x[ 0 ] );
      }
      static RangeType evaluate ( const integral_constant< unsigned int, Triangle::id >(), const DomainType &x )
      {
        return OrthonormalBase_2D::eval_triangle_2d( baseNum, &x[ 0 ] );
      }
      static RangeType evaluate ( const integral_constant< unsigned int, Hexahedron::id >(), const DomainType &x )
      {
        return OrthonormalBase_3D::eval_hexahedron_3d( baseNum, &x[ 0 ] );
      }
      static RangeType evaluate ( const integral_constant< unsigned int, Prism::id >(), const DomainType &x )
      {
        return OrthonormalBase_3D::eval_prism_3d( baseNum, &x[ 0 ] );
      }
      static RangeType evaluate ( const integral_constant< unsigned int, Pyramid::id >(), const DomainType &x )
      {
        return OrthonormalBase_3D::eval_pyramid_3d( baseNum, &x[ 0 ] );
      }
      static RangeType evaluate ( const integral_constant< unsigned int, Tetrahedron::id >(), const DomainType &x )
      {
        return OrthonormalBase_3D::eval_tetrahedron_3d( baseNum, &x[ 0 ] );
      }
    };



    // OrthonormalShapeFunctionSet
    // ---------------------------
    
    template< class FunctionSpace, int polOrder >
    class OrthonormalShapeFunctionSet
    {
      dune_static_assert( (FunctionSpace::dimRange == 1),
                          "FunctionSpace must be scalar (i.e., dimRange = 1)." );

    public:
      //! \brief function space type
      typedef typename FunctionSpace FunctionSpaceType;

      //! \brief domain type
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      OrthonormalShapeFunctionSet ( const GeometryType &type )
      : topologyId_( type.id() )
      {}

      /** @copydoc Dune::Fem::ShapeFunctionSet::size */
      std::size_t size () const
      {
        return ShapeFunctionSetCapabilities::hasStaticSize< ThisType >::size;
      }

      /** @copydoc Dune::Fem::ShapeFunctionSet::evaluateEach */
      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        GenericGeometry::IfTopology< ShapeFunctionSetHelperType::Evaluate, dimension >
          ::apply( topologyId_, coordinate( x ), functor );
      }

      /** @copydoc Dune::Fem::ShapeFunctionSet::jacobianEach */
      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        GenericGeometry::IfTopology< ShapeFunctionSetHelperType::Jacobian, dimension >
          ::apply( topologyId_, coordinate( x ), functor );
      }

      /** @copydoc Dune::Fem::ShapeFunctionSet::hessianEach */
      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        GenericGeometry::IfTopology< ShapeFunctionSetHelperType::Hessian, dimension >
          ::apply( topologyId_, coordinate( x ), functor );
      }

    private:
      unsigned int topologyId_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SHAPEFUNCTIONSET_HH
