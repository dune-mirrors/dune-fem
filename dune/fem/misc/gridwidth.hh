#ifndef DUNE_FEM_GRIDWIDTH_HH
#define DUNE_FEM_GRIDWIDTH_HH

//- system includes
#include <limits>

//- Dune includes
#include <dune/common/fvector.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/common/binaryfunctions.hh>

namespace Dune
{

  namespace Fem
  {

    template <class GridPartImp, class MinMax = Min<double> >
    class GridWidthProvider;

    //! utility functions for calculating the maximum grid width
    class GridWidth
    {
    protected:
      template< class GridPartType, class EntityType,
                class ElemGeoInfoType, class FaceGeoInfoType >
      static inline double
      maxEdgeWidth ( const GridPartType & gridPart,
                     const EntityType& entity,
                     const ElemGeoInfoType &geoInfo,
                     const FaceGeoInfoType &faceGeoInfo )
      {
        typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
        double faceVol = std::numeric_limits<double>::max() ;
        int indexInInside = -1;
        double currVol = std::numeric_limits<double>::min() ;
        double refFaceVol = std::numeric_limits<double>::min() ;

        const IntersectionIteratorType endit = gridPart.iend( entity );
        for( IntersectionIteratorType it = gridPart.ibegin( entity );
             it != endit; ++it)
        {
          typedef typename IntersectionIteratorType::Intersection Intersection;
          const Intersection &intersection = *it;

          typedef typename Intersection::Geometry LocalGeometryType;
          const LocalGeometryType &interGeo = intersection.geometry();

          const int localIndex = intersection.indexInInside();
          // calculate face Volume also for non-conforming intersections
          if( indexInInside != localIndex )
          {
            if( indexInInside >= 0 )
              faceVol = std::min( faceVol, currVol*refFaceVol );

            refFaceVol = faceGeoInfo.referenceVolume( intersection.type() );
            indexInInside = localIndex;
            currVol = 0.0;
          }

          currVol += interGeo.volume();
        }

        // do check for last set of intersections
        if( indexInInside >= 0 )
          faceVol = std::min( faceVol, currVol*refFaceVol );

        const double elemVol = entity.geometry().volume()
                             * geoInfo.referenceVolume( entity.type() );
        return elemVol / faceVol;
      }

    public:
      template <class MinMax>
      struct MinMaxInit ;

      template < class T >
      struct MinMaxInit< Min< T > >
      {
        static T init ()
        {
          return std::numeric_limits< T >::max() ;
        }
      };

      template < class T >
      struct MinMaxInit< Max< T > >
      {
        static T init ()
        {
          return std::numeric_limits< T >::min() ;
        }
      };

      /** \brief calculate minimal grid width as
          \f$h_{E} :=\frac{|E|}{h^{E}_{m}}\f$ with \f$h^{E}_{m} := \min_{e \in \mathcal{I}_{E}} |e| \f$.

          \param gridPart set of element the minmal with is calculated for
          \param communicate do global communication to get minimal width for all processes
          (default = true)
          \param minmax Function which switches between the overall minimun and maximum
      */
      template <class GridPartType, class MinMax>
      static inline
      double calcGridWidth (const GridPartType & gridPart,
                            const bool communicate,
                            const MinMax& minmax
                           )
      {
        typedef typename GridPartType :: GridType GridType;
        typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;

        // get geo infor for elements
        typedef AllGeomTypes< typename GridPartType :: IndexSetType, GridType > GeomInfoType;
        GeomInfoType geoInfo( gridPart.indexSet() );

        // get geo infor for faces
        typedef GeometryInformation< GridType , 1 > FaceGeometryInformationType;
        FaceGeometryInformationType faceGeoInfo( geoInfo.geomTypes(1) );

        // initialize with big value
        double width = MinMaxInit< MinMax > :: init() ;

        // cartesian case
        if(  Dune::Capabilities::isCartesian<GridType>::v &&
           ! Capabilities::isLocallyAdaptive<GridType>::v  )
        {
          // here we only need to check one element
          IteratorType it = gridPart.template begin<0> ();
          if( it != gridPart.template end<0> () )
          {
            width = minmax( maxEdgeWidth(gridPart,
                                         *it,
                                         geoInfo,
                                         faceGeoInfo ),
                            width );
          }
        }
        else
        {
          // unstructured case
          const IteratorType endit = gridPart.template end<0> ();
          for(IteratorType it = gridPart.template begin<0> ();
              it != endit; ++it )
          {
            width = minmax( maxEdgeWidth(gridPart,
                                         *it,
                                         geoInfo,
                                         faceGeoInfo ),
                            width );
          }
        }

        // calculate global minimum
        if( communicate )
        {
          double w = width ;
          gridPart.comm().template allreduce<MinMax> ( &w, &width, 1 );
        }

        return width;
      }
      template <class GridPartType>
      static inline
      double calcGridWidth (const GridPartType & gridPart,
                            const bool communicate = true)
      {
        return calcGridWidth ( gridPart, communicate, Min<double>() );
      }
    };

    //! utility functions for calculating the maximum grid width
    template <class GridType, class MinMax >
    class GridWidthProvider
    {
      typedef GridWidthProvider < GridType, MinMax > ThisType;
    public:
      //! type of singleton provider
      typedef SingletonList< const GridType* , ThisType > ProviderType;

    protected:
      typedef DofManager<GridType> DofManagerType;

      typedef LeafGridPart< GridType > GridPartType;

      const GridType& grid_;
      const DofManagerType& dm_;

      GridPartType gridPart_;

      mutable double width_;
      mutable int sequence_;

      GridWidthProvider( const ThisType& );
    public:
      //! constructor taking grid part
      GridWidthProvider(const GridType* grid)
        : grid_( *grid )
        , dm_( DofManagerType :: instance( grid_ ))
        , gridPart_( const_cast<GridType& > (grid_) )
        , width_(-1.0)
        , sequence_(-1)
      {
      }

      //! return characteristic grid width
      double gridWidth () const
      {
        calcWidths();
        return width_;
      }

    protected:
      void calcWidths() const
      {
#ifndef NDEBUG
        // make sure grid width calculation is done on every process
        int check = (dm_.sequence() != sequence_) ? 1 : 0;
        int willCalc = gridPart_.comm().min( check );
        assert( check == willCalc );
#endif

        if( dm_.sequence() != sequence_ )
        {
          // calculate grid width
          width_ = GridWidth :: calcGridWidth ( gridPart_ , true );

          assert( width_ > 0 );
          sequence_ = dm_.sequence();
        }
      }
    };

  } // namespace Fem

} // namespace Dune

#endif
