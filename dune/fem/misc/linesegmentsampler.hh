#ifndef DUNE_FEM_LINESEGMENTSAMPLER_HH
#define DUNE_FEM_LINESEGMENTSAMPLER_HH

#include <limits>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/grid/common/genericreferenceelements.hh>

namespace Dune
{

  namespace Fem
  {

    /** @brief class LineSegmentSampler samples values of discrete function along a given line
     *         
     *  Class LineSegmentSampler provides method for sampling the values
     *  of given discrete function along arbitrary line intersecting the @sa GridPart.
     *  The constructor takes grid part reference and two points (@a left and @a right) in 
     *  n-dimensional world space in which the grid part resides. 
     *  The method operator(@a f,@a samples) where @a f is discrete function
     *  and @a samples is a standard vector of length @a m>1, calculates values of @a f in points
     *  \f$ left + i * (right-left)/(m-1) \f$ for \f$ i=1,...,m-1 \f$.
    */
    template< class GridPart >
    struct LineSegmentSampler
    {
      typedef GridPart GridPartType;

      typedef typename GridPartType::GridType::ctype DomainFieldType;

      static const int dimDomain = GridPartType::GridType::dimensionworld;
      static const int dimGrid = GridPartType::GridType::dimension;

      typedef FieldVector< DomainFieldType, dimDomain  > DomainType;
      typedef FieldVector< DomainFieldType, dimGrid > LocalDomainType;

    private:
      dune_static_assert( dimDomain == dimGrid, "LineSegmentSampler supports only flat grids." );

      typedef GenericReferenceElement< DomainFieldType, dimGrid > ReferenceElement;
      typedef GenericReferenceElements< DomainFieldType, dimGrid > ReferenceElements;

    public:
      LineSegmentSampler ( const GridPart &gridPart, const DomainType &left, const DomainType &right )
      : gridPart_( gridPart ), left_( left ), right_( right )
      {}

      /** @brief operator() samples the values of provided discrete function along the line
       *  @param f a @sa DiscreteFunction
       *  @return samples a vector where the sampled values of function f are returned
       *
       *  Note that @a samples implicitly defines the number of sampling points. This number
       *  is simply given as a length of vector @a samples
      */
      template< class GridFunction >
      void operator() ( const GridFunction &f, std::vector< typename GridFunction::RangeType > &samples ) const;

      /** @brief returns reference to grid part
      */
      const GridPart &gridPart () const { return gridPart_; }

    private:
      const GridPart &gridPart_;
      DomainType left_, right_;
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

      typedef typename GridFunction::LocalFunctionType LocalFunctionType;

      const int numSamples = samples.size();
      if( numSamples < 2 )
        DUNE_THROW( InvalidStateException, "LineSegmentSampler cannot sample less than 2 points." );
      DomainType ds = right_ - left_;
      ds /= DomainFieldType( numSamples - 1 );

      const typename GridFunction::RangeFieldType nan
        = std::numeric_limits< typename GridFunction::RangeFieldType >::quiet_NaN();
      for( int i = 0; i < numSamples; ++i )
        samples[ i ] = typename GridFunction::RangeType( nan );

      const IteratorType end = gridPart().template end< 0 >();
      for( IteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
      {
        const EntityType &entity = *it;
        const typename EntityType::Geometry &geometry = entity.geometry();

        DomainFieldType lower = std::numeric_limits< DomainFieldType >::min();
        DomainFieldType upper = std::numeric_limits< DomainFieldType >::max();

        const ReferenceElement &refElement = ReferenceElements::general( geometry.type() );
        const int numFaces = refElement.size( 1 );
        for( int face = 0; face < numFaces; ++face )
        {
          const LocalDomainType &center = refElement.position( face, 1 );

          DomainType normal;
          const LocalDomainType &refNormal = refElement.volumeOuterNormal( face );
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
          const LocalFunctionType &lf = f.localFunction( entity );
          const int i_upper = std::min( int( std::floor( upper + 1e-8 ) ), numSamples - 1 );
          const int i_lower = std::max( int( std::ceil( lower - 1e-8 ) ), 0 );
          for( int i = i_lower; i <= i_upper; ++i )
          {
            DomainType xi = left_;
            xi.axpy( DomainFieldType( i ), ds );
            lf.evaluate( geometry.local( xi ), samples[ i ] );
          }
        }
      }

      bool valid = true;
      for( int i = 0; i < numSamples; ++i )
        valid &= (samples[ i ] == samples[ i ]);
      if( !valid )
        DUNE_THROW( InvalidStateException, "LineSegmentSampler could not find all samples." );
    }

  }

}

#endif // #ifndef DUNE_FEM_LINESEGMENTSAMPLER_HH
