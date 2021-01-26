#ifndef DUNE_FEM_LINESEGMENTSAMPLER_HH
#define DUNE_FEM_LINESEGMENTSAMPLER_HH

#include <limits>
#include <vector>
#include <cmath>
#include <type_traits>

#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/function/localfunction/const.hh>

namespace Dune
{

  namespace Fem
  {

    // LineSegmentSampler
    // ------------------

    /** \class LineSegmentSampler
     *  \brief samples values of a discrete function along a given line segment
     *
     *  The class LineSegmentSampler provides a method for sampling the values
     *  of given discrete function along an arbitrary line contained in some
     *  GridPart. The sampling points are always equidistant and include the
     *  line segment's end points.
     *
     *  \note The grid is required to be flat, i.e., grid dimension and world
     *        dimension must coincide.
     *
     *  \tparam  GridPart  type of grid part to sample on
     */
    template< class GridPart >
    struct LineSegmentSampler
    {
      typedef GridPart GridPartType;

      typedef typename GridPartType::GridType::ctype DomainFieldType;

      static const int dimDomain = GridPartType::dimensionworld;
      static const int dimGrid = GridPartType::dimension;

      typedef FieldVector< DomainFieldType, dimDomain > DomainType;
      typedef FieldVector< DomainFieldType, dimGrid > LocalDomainType;

    private:
      static_assert( dimDomain == dimGrid, "LineSegmentSampler supports only flat grids." );

      template< class Vector >
      struct Reduce;

      typedef Dune::ReferenceElements< DomainFieldType, dimGrid > ReferenceElements;

    public:
      /** \brief constructor
       *
       *  \param[in]  gridPart  the grid part to sample over
       *  \param[in]  left      left end point of the line segment
       *  \param[in]  right     right end point of the line segment
       *
       *  \note Actually, left and right can be exchanged. The order only
       *        influences the parametrization of the line segment.
       *  \note The entire line segment must be a subset of the grid part.
       */
      LineSegmentSampler ( const GridPart &gridPart, const DomainType &left, const DomainType &right )
      : gridPart_( gridPart ), left_( left ), right_( right )
      {}

      /** \brief sample a given function
       *
       *  The operator() actually samples the values of a given grid function.
       *
       *  \param[in]   f        grid function to sample
       *  \param[out]  samples  std::vector receiving the samples
       *
       *  \note The number of sampling points is determined from the size of @a
       *        samples, which may not be less than 2.
       */
      template< class GridFunction >
      void operator() ( const GridFunction &f, std::vector< typename GridFunction::RangeType > &samples ) const;

      /** \brief returns sampling points
       *
       *  The operator() actually samples the values of a given grid function.
       *
       *  \param[out]  points  std::vector receiving the points
       *
       *  \note The number of sampling points is determined from the size of @a
       *        points, which may not be less than 2.
       */
      void samplePoints( std::vector< DomainType > &points ) const;

      /** @brief obtain grid part on which the LineSegmentSampler works */
      const GridPart &gridPart () const { return gridPart_; }

    private:
      const GridPart &gridPart_;
      DomainType left_, right_;
    };



    // LineSegmentSampler::Reduce
    // --------------------------

    template< class GridPart >
    template< class Vector >
    struct LineSegmentSampler< GridPart >::Reduce
    {
      Vector operator() ( const Vector &a, const Vector &b ) const
      {
        Vector c;
        for( int k = 0; k < Vector::dimension; ++k )
          c[ k ] = std::min( a[ k ], b[ k ] );
        return c;
      }
    };



    // Implementation of LineSegmentSampler
    // ------------------------------------

    template< class GridPart >
    template< class GridFunction >
    inline void LineSegmentSampler< GridPart >
      ::operator() ( const GridFunction &f, std::vector< typename GridFunction::RangeType > &samples ) const
    {
      typedef typename GridPartType::template Codim< 0 >::IteratorType IteratorType;
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      const int numSamples = samples.size();
      if( numSamples < 2 )
        DUNE_THROW( InvalidStateException, "LineSegmentSampler cannot sample less than 2 points." );
      DomainType ds = right_ - left_;
      ds /= DomainFieldType( numSamples - 1 );

      typedef typename GridFunction::RangeFieldType RangeFieldType;
      const RangeFieldType invalid
        = std::numeric_limits< RangeFieldType >::quiet_NaN();
      for( int i = 0; i < numSamples; ++i )
        samples[ i ] = typename GridFunction::RangeType( invalid );

      ConstLocalFunction< GridFunction > lf( f );
      const IteratorType end = gridPart().template end< 0 >();
      for( IteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
      {
        const EntityType &entity = *it;
        const typename EntityType::Geometry &geometry = entity.geometry();

        DomainFieldType lower = std::numeric_limits< DomainFieldType >::min();
        DomainFieldType upper = std::numeric_limits< DomainFieldType >::max();

        auto &refElement = ReferenceElements::general( geometry.type() );
        const int numFaces = refElement.size( 1 );
        for( int face = 0; face < numFaces; ++face )
        {
          const LocalDomainType &center = refElement.position( face, 1 );

          DomainType normal;
          const LocalDomainType &refNormal = refElement.integrationOuterNormal( face );
          geometry.jacobianInverseTransposed( center ).mv( refNormal, normal );

          const DomainFieldType nds = normal * ds;
          const DomainFieldType ncl = normal * (geometry.global( center ) - left_);
          if( std::abs( nds ) > 1e-8 )
          {
            // ds is not parallel to the face
            const DomainFieldType alpha = ncl / nds;
            if( nds > 0 )
              upper = std::min( upper, alpha );
            else
              lower = std::max( lower, alpha );
          }
          else if( ncl < -1e-8 )
          {
            // ds is parallel to the face and the line lies on the outside
            lower = std::numeric_limits< DomainFieldType >::max();
            upper = std::numeric_limits< DomainFieldType >::min();
          }
        }

        if( lower <= upper )
        {
          lf.bind( entity );
          const int i_upper = std::min( int( std::floor( upper + 1e-8 ) ), numSamples - 1 );
          const int i_lower = std::max( int( std::ceil( lower - 1e-8 ) ), 0 );
          for( int i = i_lower; i <= i_upper; ++i )
          {
            DomainType xi = left_;
            xi.axpy( DomainFieldType( i ), ds );
            lf.evaluate( geometry.local( xi ), samples[ i ] );
          }
          lf.unbind();
        }
      }

      typedef Reduce< typename GridFunction::RangeType > Op;
      gridPart().comm().template allreduce< Op >( &(samples[ 0 ]), numSamples );

      bool valid = true;

      // only use isnan when field type is a floating point
      if constexpr ( std::is_floating_point< RangeFieldType >::value )
      {
        for( int i = 0; i < numSamples; ++i )
        {
          for( int d=0; d<GridFunction::RangeType::dimension; ++d )
          {
            valid &= ! std::isnan( samples[ i ][ d ] );
          }
        }
      }

      if( !valid )
        DUNE_THROW( InvalidStateException, "LineSegmentSampler could not find all samples." );
    }


    template< class GridPart >
    inline void LineSegmentSampler< GridPart >
      :: samplePoints ( std::vector< DomainType > &points ) const
    {
      const int numSamples = points.size();
      if( numSamples < 2 )
        DUNE_THROW( InvalidStateException, "LineSegmentSampler cannot sample less than 2 points." );
      DomainType ds = right_ - left_;
      ds /= DomainFieldType( numSamples - 1 );

      for( int i = 0; i < numSamples; ++i )
      {
        points[ i ] = left_;
        points[ i ].axpy( DomainFieldType( i ), ds );
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_LINESEGMENTSAMPLER_HH
