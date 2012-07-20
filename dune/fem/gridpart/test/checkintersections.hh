#ifndef DUNE_FEM_GRIDPART_TEST_CHECKINTERSECTIONS_HH
#define DUNE_FEM_GRIDPART_TEST_CHECKINTERSECTIONS_HH

//- dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/stdstreams.hh>

//- dune-geometry includes
#include <dune/geometry/quadraturerules/gaussquadrature.hh>
#include <dune/geometry/referenceelements.hh>

//- dune-grid includes
#include <dune/grid/test/checkgeometry.cc>

//- dune-fem includes
#include <dune/fem/gridpart/test/failure.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>

namespace Dune
{
  namespace Fem
  {
    template< class GridPartType, class FailureHandler >
    struct CheckIntersections
    {
      /** \brief type of intersection iterator */
      typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
      /** \brief type of intersection */
      typedef typename GridPartType::IntersectionType IntersectionType;
      /** \brief loal geometry type */
      typedef typename IntersectionType::Geometry GeometryType;
      /** \brief loal geometry type */
      typedef typename IntersectionType::LocalGeometry LocalGeometryType;
      /** \brief entity pointer type */
      typedef typename IntersectionType::EntityPointer EntityPointerType;
      /** \brief entity type */
      typedef typename EntityPointerType::Entity EntityType;
    
      /** \brief failure type */
      struct AssignmentFailure;
      /** \brief failure type */
      struct SumNormalFailure;

    private:
      // static tests
      dune_static_assert( ( is_same< typename IntersectionType::ctype,
                                     typename GridPartType::ctype 
                                   >::value 
                          ), "Conflicting types for Intersection" );

      dune_static_assert( ( static_cast< int>( IntersectionType::dimension )
                              == static_cast< int >( GridPartType::dimension )
                          ), "IntersectionIterator has wrong dimension" );
      dune_static_assert( ( static_cast< int>( IntersectionType::dimensionworld )
                              == static_cast< int >( GridPartType::dimensionworld )
                          ), "IntersectionIterator has wrong dimensionworld" );

      // check intersection iterator assignment
      template< class IntersectionIterator >
      static void checkIntersectionIteratorAssignment( const IntersectionIterator &it, 
                                                       FailureHandler &failureHandler );

      // check intersection iterator assignment
      static void checkOuterNormals ( const EntityType &entity,
                                      const GridPartType &gridPart,
                                      FailureHandler &failureHandler );

      // check local geometries 
      static void checkLocalGeometries ( const IntersectionType &intersection,
                                         const EntityType &entity,
                                         FailureHandler &failureHandler );

    public:
      static void check ( const GridPartType &gridPart, FailureHandler &failureHandler )
      {
        typedef typename GridPartType::template Codim< 0 >::IteratorType IteratorType;
        typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

        const IteratorType end = gridPart.template end< 0>();
        for( IteratorType it = gridPart.template begin< 0 >(); it != end; ++it )
        {
          const EntityType &entity = *it;

          typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
          const IntersectionIteratorType iend = gridPart.iend( entity );
          for( IntersectionIteratorType iit = gridPart.ibegin( entity ); iit != iend; ++iit )
          {
            const IntersectionType &intersection = *iit;

            // create intersection caching quadrature to check twists 
            IntersectionQuadrature< CachingQuadrature< GridPartType, 1 >, true > 
              inter( gridPart, intersection, 2 );

            checkIntersectionIteratorAssignment( iit, failureHandler );
            checkLocalGeometries( intersection, entity, failureHandler );
          }

          checkOuterNormals( entity, gridPart, failureHandler );
          if( SumNormalFailure::failed() )
            failureHandler( SumNormalFailure::instance() );
        }
      }
    };



    // Implementation of CheckIntersections
    // ------------------------------------

    template< class GridPartType, class FailureHandler >
    template< class IntersectionIterator >
    inline void CheckIntersections< GridPartType, FailureHandler >
      ::checkIntersectionIteratorAssignment( const IntersectionIterator &it, FailureHandler &failureHandler )
    {
      // type of intersection iterator
      typedef IntersectionIterator IntersectionIteratorType;

      AssignmentFailure failure;
      // assert assignment
      IntersectionIteratorType it2 = it;
      if( it != it2 )
        failureHandler( failure );
      if( it->inside() != it2->inside() )
        failureHandler( failure );
      if( it->neighbor() )
        if( !it->neighbor() || it->outside() != it2->outside() )
          failureHandler( failure );

      // now increment second iterator
      ++it2;
      if( it == it2 )
        failureHandler( failure );
    }


    template< class GridPartType, class FailureHandler >
    inline void CheckIntersections< GridPartType, FailureHandler >
      ::checkOuterNormals ( const EntityType &entity,
                            const GridPartType &gridPart,
                            FailureHandler &failureHandler )
    {
      // global coordinate type
      typedef typename GeometryType::GlobalCoordinate GlobalCoordinateType;
      // local coordinate type
      typedef typename GeometryType::LocalCoordinate LocalCoordinateType;
      // single coordinate type
      typedef typename GridPartType::ctype ctype;

      // get Gauss quadrature points
      typedef Dune::GenericGeometry::GaussPoints< ctype > Points;
      typedef Dune::GenericGeometry::GenericQuadratureFactory< 
          LocalGeometryType::mydimension, ctype, Points 
        > QuadratureFactory;
      const typename QuadratureFactory::Object &quad = *QuadratureFactory::create( entity.type(), 3 );

      // compute sum of all integration outer normals
      GlobalCoordinateType sum( 0 );
      const IntersectionIteratorType end = gridPart.iend( entity );
      for( IntersectionIteratorType it = gridPart.ibegin( entity ); it != end; ++it )
      {
        for ( size_t i = 0; i < quad.size(); ++i )
        {
          const LocalCoordinateType &point = quad[ i ].position();
          GlobalCoordinateType integrationOuterNormal = it->integrationOuterNormal( point );
          sum.axpy( quad[ i ].weight(), integrationOuterNormal );
        }
      }

      // check wheter error is within fixed tolerance
      double error = sum.two_norm();
      const double tolerance = SumNormalFailure::tolerance;
      if( (sum.two_norm() > tolerance ) && ( entity.partitionType() != Dune::GhostEntity) )
        SumNormalFailure::add( entity, error );
    }


    template< class GridPartType, class FailureHandler >
    inline void CheckIntersections< GridPartType, FailureHandler >
      ::checkLocalGeometries ( const IntersectionType &intersection,
                               const EntityType &entity,
                               FailureHandler &failureHandler )
    {
      LocalGeometryType geometryInInside = intersection.geometryInInside();
      geometryInInside.type();
      Dune::checkLocalGeometry( geometryInInside, entity.type(), "geometryInInside" );

      if( intersection.neighbor() )
      {
        LocalGeometryType geometryInOutside = intersection.geometryInOutside();
        Dune::checkLocalGeometry( geometryInOutside, entity.type(), "geometryInOutside" );
      }
    }



    // Implementation of Failures used in CheckIntersections 
    // -----------------------------------------------------

    template< class GridPartType, class FailureHandler >
    struct CheckIntersections< GridPartType, FailureHandler >::AssignmentFailure
    : public Failure
    {
      virtual void writeTo ( std::ostream &out ) const
      {
        out <<  __FILE__
          << ":" << __LINE__ << ": Failure :"
          << "Assignment failure";
      }
    };



    template< class GridPartType, class FailureHandler >
    struct CheckIntersections< GridPartType, FailureHandler >::SumNormalFailure
    : public Failure
    {
      /** \brief entity type */
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      /** \brief tolerance */
      static const double tolerance;

    protected:
      /** \brief constructor */
      SumNormalFailure ()
      : fails_( 0 )
      {}

    public:
      /** \brief return instance of singleton */
      static SumNormalFailure &instance ()
      {
        static SumNormalFailure instance_;
        return instance_;
      }

      /** \brief return true, if actual failure */
      static bool failed ()
      {
        return ( instance().fails_ > 0 );
      }

      /** \brief add failed entity */
      static void add ( const EntityType &entity, double error )
      {
        instance().fails_ += 1;
      }

      /** \brief write message to stream */
      virtual void writeTo ( std::ostream &out ) const
      {
        assert( failed() );
        out <<  __FILE__
          << ":" << __LINE__ << ": Failure :"
          << "Summation of normals fails for " 
          << instance().fails_ << " entities in grid part";
      }

    private:
      int fails_;
    };



    template< class GridPartType, class FailureHandler >
    const double CheckIntersections< GridPartType, FailureHandler >::SumNormalFailure::tolerance = 1e-8;

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_TEST_CHECKINTERSECTIONS_HH
