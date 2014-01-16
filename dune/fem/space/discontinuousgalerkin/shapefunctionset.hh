#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SHAPEFUNCTIONSET_HH

// C++ includes
#include <cassert>
#include <cstdlib>

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/static_assert.hh>

// dune-geometry includes
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/type.hh>

// local includes 
#include "orthonormalbase_1d.hh"
#include "orthonormalbase_2d.hh"
#include "orthonormalbase_3d.hh"

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



    // OrthonormalShapeFunctionSetSize
    // -------------------------------

    template< class FunctionSpace, int polOrder >
    class OrthonormalShapeFunctionSetSize
    {
      dune_static_assert( (FunctionSpace::dimDomain <= 3),
                          "Shape function set only implemented up to dimension 3." );

      template< int order, int dimension >
      struct ShapeFunctionSetSize;

      template< int order >
      struct ShapeFunctionSetSize< order, 1 >
      {
        static const std::size_t v = order + 1;
      };

      template< int order >
      struct ShapeFunctionSetSize< order, 2 >
      {
        static const std::size_t v = (order + 2) * (order + 1) / 2;
      };

      template< int order >
      struct ShapeFunctionSetSize< order, 3 >
      {
        static const size_t v = ((order+1)*(order+2)*(2*order+3)/6
                               + (order+1)*(order+2)/2)/2;
      };

    public:
      static const std::size_t v = ShapeFunctionSetSize< polOrder, FunctionSpace::dimDomain >::v;
    };



    // OrthonormalShapeFunctionHelper
    // ------------------------------

    template< class FunctionSpace, int polOrder >
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
      // pyramid
      typedef Dune::GenericGeometry::Prism< Quadrilateral > Hexahedron;
      // hexahedron
      typedef Dune::GenericGeometry::Pyramid< Quadrilateral > Pyramid;
      // prism
      typedef Dune::GenericGeometry::Prism< Triangle > Prism;
      // tetrahedron
      typedef Dune::GenericGeometry::Pyramid< Triangle > Tetrahedron;

      // mapping from topology to basic geometry type (see list above)
      template< class Topology >
      class BasicGeometryType 
      {
        template< unsigned int topologyId >
        struct UniqueId
        {
          static const unsigned int v = topologyId | (unsigned int)Dune::GenericGeometry::prismConstruction;
        };

      public:
        typedef typename Dune::GenericGeometry::Topology< UniqueId< Topology::id >::v, Topology::dimension >::type Type;
      };

      // instantiate a basic geometry type
      template< class Topology >
      static typename BasicGeometryType< Topology >::Type basicGeometryType ()
      {
        return typename BasicGeometryType< Topology >::Type();
      }

    public:
      template< class Topology >
      struct EvaluateEach;

      template< class Topology >
      struct JacobianEach;

      template< class Topology >
      struct HessianEach;
    };



    // Implementation of OrthonormalShapeFunctionHelper::EvaluateEach
    // --------------------------------------------------------------

    template< class FunctionSpace, int polOrder >
    template< class Topology >
    struct OrthonormalShapeFunctionHelper< FunctionSpace, polOrder >::EvaluateEach
    {
      template< class Functor >
      static void apply ( const DomainType &x, Functor functor )
      {
        const std::size_t size = OrthonormalShapeFunctionSetSize< FunctionSpace, polOrder >::v;
        for( std::size_t i = 0; i < size; ++i )
        {
          assert( x.size() == Topology::dimension );
          functor( i, evaluate( basicGeometryType< Topology>(), i, x ) );
        }
      }


    protected:
      typedef typename DomainType :: value_type  DomainField;
      typedef typename RangeType  :: value_type  RangeField;
      typedef OrthonormalBase_1D< DomainField, RangeField > OrthonormalBase1d;
      typedef OrthonormalBase_2D< DomainField, RangeField > OrthonormalBase2d;
      typedef OrthonormalBase_3D< DomainField, RangeField > OrthonormalBase3d;

      static RangeType evaluate ( const Line &, std::size_t i, const DomainType &x )
      {
        return OrthonormalBase1d::eval_line( i, &x[ 0 ] );
      }
      static RangeType evaluate ( const Quadrilateral &, std::size_t i, const DomainType &x )
      {
        return OrthonormalBase2d::eval_quadrilateral_2d( i, &x[ 0 ] );
      }
      static RangeType evaluate ( const Triangle &, std::size_t i, const DomainType &x )
      {
        return OrthonormalBase2d::eval_triangle_2d( i, &x[ 0 ] );
      }
      static RangeType evaluate ( const Pyramid &, std::size_t i, const DomainType &x )
      {
        return OrthonormalBase3d::eval_pyramid_3d( i, &x[ 0 ] );
      }
      static RangeType evaluate ( const Hexahedron &, std::size_t i, const DomainType &x )
      {
        return OrthonormalBase3d::eval_hexahedron_3d( i, &x[ 0 ] );
      }
      static RangeType evaluate ( const Prism &, std::size_t i, const DomainType &x )
      {
        return OrthonormalBase3d::eval_prism_3d( i, &x[ 0 ] );
      }
      static RangeType evaluate ( const Tetrahedron &, std::size_t i, const DomainType &x )
      {
        return OrthonormalBase3d::eval_tetrahedron_3d( i, &x[ 0 ] );
      }
    };



    // Implementation of OrthonormalShapeFunctionHelper::JacobianEach
    // --------------------------------------------------------------

    template< class FunctionSpace, int polOrder >
    template< class Topology >
    struct OrthonormalShapeFunctionHelper< FunctionSpace, polOrder >::JacobianEach
    {
      template< class Functor >
      static void apply ( const DomainType &x, Functor functor )
      {
        JacobianRangeType jacobian;
        const std::size_t size = OrthonormalShapeFunctionSetSize< FunctionSpace, polOrder >::v;
        for( std::size_t i = 0; i < size; ++i )
        {
          assert( x.size() == Topology::dimension );
          evaluate( basicGeometryType< Topology >(), i, x, jacobian );
          functor( i, jacobian );
        }
      }

    protected:
      typedef typename DomainType :: value_type       DomainField;
      typedef typename JacobianRangeType:: value_type RangeField;
      typedef OrthonormalBase_1D< DomainField, RangeField > OrthonormalBase1d;
      typedef OrthonormalBase_2D< DomainField, RangeField > OrthonormalBase2d;
      typedef OrthonormalBase_3D< DomainField, RangeField > OrthonormalBase3d;

      static void evaluate ( const Line &, std::size_t i, const DomainType &x, JacobianRangeType &jacobian )
      {
        OrthonormalBase1d::grad_line( i , &x[ 0 ], &jacobian[ 0 ][ 0 ] );
      }
      static void evaluate ( const Quadrilateral &, std::size_t i, const DomainType &x, JacobianRangeType &jacobian )
      {
        OrthonormalBase2d::grad_quadrilateral_2d( i , &x[ 0 ], &jacobian[ 0 ][ 0 ] );
      }
      static void evaluate ( const Triangle &, std::size_t i, const DomainType &x, JacobianRangeType &jacobian )
      {
        OrthonormalBase2d::grad_triangle_2d( i , &x[ 0 ], &jacobian[ 0 ][ 0 ] );
      }
      static void evaluate ( const Pyramid &, std::size_t i, const DomainType &x, JacobianRangeType &jacobian )
      {
        OrthonormalBase3d::grad_pyramid_3d( i , &x[ 0 ], &jacobian[ 0 ][ 0 ] );
      }
      static void evaluate ( const Hexahedron &, std::size_t i, const DomainType &x, JacobianRangeType &jacobian )
      {
        OrthonormalBase3d::grad_hexahedron_3d( i , &x[ 0 ], &jacobian[ 0 ][ 0 ] );
      }
      static void evaluate ( const Prism &, std::size_t i, const DomainType &x, JacobianRangeType &jacobian )
      {
        OrthonormalBase3d::grad_prism_3d( i , &x[ 0 ], &jacobian[ 0 ][ 0 ] );
      }
      static void evaluate ( const Tetrahedron &, std::size_t i, const DomainType &x, JacobianRangeType &jacobian )
      {
        OrthonormalBase3d::grad_tetrahedron_3d( i , &x[ 0 ], &jacobian[ 0 ][ 0 ] );
      }
    };

    // Implementation of OrthonormalShapeFunctionHelper::HessianEach
    // -------------------------------------------------------------

    template< class FunctionSpace, int polOrder >
    template< class Topology >
    struct OrthonormalShapeFunctionHelper< FunctionSpace, polOrder >::HessianEach
    {
      template< class Functor >
      static void apply ( const DomainType &x, Functor functor )
      {
        typedef typename DomainType :: value_type Field;
        Field hess[3] = { 0, 0, 0 };
        HessianRangeType hessian;
        const std::size_t size = OrthonormalShapeFunctionSetSize< FunctionSpace, polOrder >::v;
        for( std::size_t i = 0; i < size; ++i )
        {
          assert( x.size() == Topology::dimension );
          evaluate( basicGeometryType< Topology >(), i, x, hess );
          for( unsigned int j = 0; j < FunctionSpace::dimDomain;  ++j )
            for( unsigned int k = 0; k < FunctionSpace::dimDomain; ++k )
              hessian[ 0 ][ j ][ k ] = hess[ j + k ];

          functor( i, hessian );
        }
      }

    protected:
      typedef typename DomainType       :: value_type  DomainField;
      typedef typename HessianRangeType :: value_type  RangeField;
      typedef OrthonormalBase_2D< DomainField, RangeField > OrthonormalBase2d;

      static void evaluate ( const Line &, std::size_t i, const DomainType &x, double (&hessian)[3] )
      {
        DUNE_THROW( NotImplemented, "On orthonormal shape function set HessianAll() is only implemented for triangles" );
      }
      static void evaluate ( const Quadrilateral &, std::size_t i, const DomainType &x, double (&hessian)[3] )
      {
        DUNE_THROW( NotImplemented, "On orthonormal shape function set HessianAll() is only implemented for triangles" );
      }
      static void evaluate ( const Triangle &, std::size_t i, const DomainType &x, double (&hessian)[3] )
      {
        OrthonormalBase2d::hess_triangle_2d( i , &x[ 0 ], hessian );
      }
      static void evaluate ( const Pyramid &, std::size_t i, const DomainType &x, double (&hessian)[3] )
      {
        DUNE_THROW( NotImplemented, "On orthonormal shape function set HessianAll() is only implemented for triangles" );
      }
      static void evaluate ( const Hexahedron &, std::size_t i, const DomainType &x, double (&hessian)[3] )
      {
        DUNE_THROW( NotImplemented, "On orthonormal shape function set HessianAll() is only implemented for triangles" );
      }
      static void evaluate ( const Prism &, std::size_t i, const DomainType &x, double (&hessian)[3] )
      {
        DUNE_THROW( NotImplemented, "On orthonormal shape function set HessianAll() is only implemented for triangles" );
      }
      static void evaluate ( const Tetrahedron &, std::size_t i, const DomainType &x, double (&hessian)[3] )
      {
        DUNE_THROW( NotImplemented, "On orthonormal shape function set HessianAll() is only implemented for triangles" );
      }
    };


    // OrthonormalShapeFunctionSet
    // ---------------------------
    
    template< class FunctionSpace, int polOrder >
    class OrthonormalShapeFunctionSet
    {
      dune_static_assert( (FunctionSpace::dimRange == 1),
                          "FunctionSpace must be scalar (i.e., dimRange = 1)." );

      // this type
      typedef OrthonormalShapeFunctionSet< FunctionSpace, polOrder > ThisType;

      // helper class
      typedef OrthonormalShapeFunctionHelper< FunctionSpace, polOrder > ShapeFunctionSetHelperType;

    public:
      //! \brief function space type
      typedef FunctionSpace FunctionSpaceType;

      //! \brief dimension
      static const int dimension = FunctionSpaceType::dimDomain;

      //! \brief domain type
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      explicit OrthonormalShapeFunctionSet ( const GeometryType &type )
      : topologyId_( type.id() )
      {}

      /** @copydoc Dune::Fem::ShapeFunctionSet::order*/
      int order () const { return polOrder; }


      /** @copydoc Dune::Fem::ShapeFunctionSet::size */
      std::size_t size () const
      {
        return OrthonormalShapeFunctionSetSize< FunctionSpace, polOrder >::v;
      }

      /** @copydoc Dune::Fem::ShapeFunctionSet::evaluateEach */
      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        Dune::GenericGeometry::IfTopology< ShapeFunctionSetHelperType::template EvaluateEach, dimension >
          ::apply( topologyId_, coordinate( x ), functor );
      }

      /** @copydoc Dune::Fem::ShapeFunctionSet::jacobianEach */
      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        Dune::GenericGeometry::IfTopology< ShapeFunctionSetHelperType::template JacobianEach, dimension >
          ::apply( topologyId_, coordinate( x ), functor );
      }

      /** @copydoc Dune::Fem::ShapeFunctionSet::hessianEach */
      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        Dune::GenericGeometry::IfTopology< ShapeFunctionSetHelperType::template HessianEach, dimension >
          ::apply( topologyId_, coordinate( x ), functor );
      }

    private:
      unsigned int topologyId_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SHAPEFUNCTIONSET_HH
