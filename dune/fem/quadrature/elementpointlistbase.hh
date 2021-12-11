#ifndef DUNE_FEM_ELEMENTPOINTLISTBASE_HH
#define DUNE_FEM_ELEMENTPOINTLISTBASE_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/quadrature/quadrature.hh>

#include <dune/fem/gridpart/common/capabilities.hh>

namespace Dune
{

  namespace Fem
  {

    /** \brief ElementPointListBase */
    template< class GridPartImp, int codim, class IntegrationTraits >
    class ElementPointListBase;


    template< class GridPartImp, class IntegrationTraits >
    class ElementPointListBase< GridPartImp, 0, IntegrationTraits >
    {
      typedef ElementPointListBase< GridPartImp, 0, IntegrationTraits > This;

    public:
      //! type of the grid partition
      typedef GridPartImp GridPartType;

      //! inside and outside flags
      enum Side { INSIDE, OUTSIDE };

      //! codimension of the integration point list
      static const int codimension = 0;

      //! coordinate type
      typedef typename GridPartType::ctype RealType;

      //! dimension of the grid
      static const int dimension = GridPartType::dimension;

      //! type of the integration point list
      typedef typename IntegrationTraits::IntegrationPointListType IntegrationPointListType;

      typedef typename IntegrationTraits::CoordinateType CoordinateType;
      typedef typename IntegrationPointListType::CoordinateType LocalCoordinateType;

      typedef typename IntegrationPointListType :: QuadratureKeyType  QuadratureKeyType;

      /** \brief constructor
       *
       *  \param[in]  geometry  geometry type, the quadrature lives on
       *  \param[in]  order     desired minimal order of the quadrature
       */
      ElementPointListBase ( const GeometryType &geometry, const QuadratureKeyType& quadKey )
      : quad_( geometry, quadKey )
      {}

      /** \brief constructor
       *
       *  \param[in]  geometry  geometry type, the quadrature lives on
       *  \param[in]  order     desired minimal order of the quadrature
       */
      ElementPointListBase ( const IntegrationPointListType& ipList )
      : quad_( ipList )
      {}

      /** \copydoc Dune::Fem::IntegrationPointList::nop */
      size_t nop () const
      {
        return quadImp().nop();
      }

      /** \brief obtain local coordinates of i-th integration point
       *
       *  This method returns a reference to the local coordinates of the i-th
       *  integration point for 0 <= i < nop(). Here, local coordinates means
       *  coordinates with respect to the reference element of the subentity.
       *
       *  \param[in]  i  number of the integration point, 0 <= i < nop()
       *
       *  \returns reference to i-th integration point
       */
      const LocalCoordinateType &localPoint( size_t i ) const
      {
        return quadImp().point( i );
      }

      /** \copydoc Dune::Fem::IntegrationPointList::id
       */
      size_t id () const
      {
        return quadImp().id();
      }

      /** \copydoc Dune::Fem::IntegrationPointList::order
       */
      int order () const
      {
        return quadImp().order();
      }

      /** \copydoc Dune::Fem::IntegrationPointList::geometry
       */
      GeometryType geometry () const
      {
        return quadImp().geometryType();
      }

      /** \copydoc Dune::Fem::IntegrationPointList::geometry
       */
      GeometryType type () const
      {
        return quadImp().geometryType();
      }

      /** \copydoc Dune::Fem::IntegrationPointList::geometry
       */
      GeometryType geometryType () const
      {
        return quadImp().geometryType();
      }

      /** \brief obtain GeometryType of the corresponding codim-0 the integration
       *         point list belongs to
       *
       *  An element integration point list can return the coordinates of integration
       *  points with resepct to the codim-0 reference element and the reference
       *  element corresponding to the subentity the quadrature actually lives on.
       *  This method returns the geometry of the codim-0 entity.
       *
       *  \note Calling this method yields a virtual function call, so do not
       *        call this method unnecessarily.
       *
       *  \returns GeometryType for this integration point list
       */
      GeometryType elementGeometry () const
      {
        return quadImp().geometry();
      }

      /** \brief convenience implementation for Dune::Fem::CachingInterface */
      size_t cachingPoint( const size_t quadraturePoint ) const
      {
        return quadraturePoint;
      }

      /** \brief convenience implementation for Dune::Fem::CachingInterface */
      size_t localCachingPoint( const size_t quadraturePoint ) const
      {
        return quadraturePoint;
      }

      /** \brief convenience implementation for Dune::Fem::CachingInterface */
      static constexpr bool twisted () { return false; }

      /** \brief convenience implementation for Dune::Fem::CachingInterface */
      inline int twistId () const { return 0; }

      int localFaceIndex () const
      {
        return 0;
      }

      inline int nCachingPoints () const { return nop(); }
      inline int cachingPointStart () const { return 0; }

    protected:
      /** \brief obtain the actual implementation of the quadrature
       *
       *  \note This method may only be used in derived classes.
       *
       *  \returns a reference to the actual implementation of the quadrature
       */
      const IntegrationPointListType &quadImp () const
      {
        return quad_;
      }

    protected:
      IntegrationPointListType quad_;
    };



    /** \copydoc ElementIntegrationPointList */
    template< class GridPartImp, int codim, class IntegrationTraits >
    class ElementPointListBase
    {
      typedef ElementPointListBase< GridPartImp, codim, IntegrationTraits > This;

    public:
      //! type of the grid partition
      typedef GridPartImp GridPartType;

      //! inside and outside flags
      enum Side { INSIDE, OUTSIDE };

      //! codimension of the element integration point list
      static const int codimension = codim;

      //! coordinate type
      typedef typename GridPartType::ctype RealType;

      //! dimension of the grid
      static const int dimension = GridPartType::dimension;

      //! type of the integration point list
      typedef typename IntegrationTraits::IntegrationPointListType IntegrationPointListType;

      typedef typename IntegrationTraits::CoordinateType CoordinateType;
      typedef typename IntegrationPointListType::CoordinateType LocalCoordinateType;

      typedef typename IntegrationPointListType :: QuadratureKeyType  QuadratureKeyType;

      /** \brief constructor
       *
       *  \param[in]  elementGeo      geometry type of the element
       *  \param[in]  faceGeo         geometry type of the subentity
       *  \param[in]  localFaceIndex  index of the subentity
       *  \param[in]  order           desired order of the quadrature
       */
      ElementPointListBase ( const GeometryType &elementGeo,
                             const GeometryType &faceGeo,
                             const int localFaceIndex,
                             const QuadratureKeyType& quadKey )
      : quad_( faceGeo, quadKey ),
        elementGeometry_( elementGeo ),
        localFaceIndex_( localFaceIndex )
      {}

      /** \brief constructor
       *
       *  \param[in]  elementGeo      geometry type of the element
       *  \param[in]  localFaceIndex  index of the subentity
       *  \param[in]  order           desired order of the quadrature
       */
      ElementPointListBase ( const GeometryType &elementGeo,
                             const int localFaceIndex,
                             const QuadratureKeyType& quadKey )
      : quad_( getFaceGeometry( elementGeo, localFaceIndex ), quadKey ),
        elementGeometry_( elementGeo ),
        localFaceIndex_( localFaceIndex )
      {}

      /** \copydoc Dune::Fem::IntegrationPointList::nop
       */
      size_t nop () const
      {
        return quadImp().nop();
      }

      /** \brief obtain local coordinates of i-th integration point
       *
       *  This method returns a reference to the local coordinates of the i-th
       *  integration point for 0 <= i < nop(). Here, local coordinates means
       *  coordinates with respect to the reference element of the subentity.
       *
       *  \param[in]  i  number of the integration point, 0 <= i < nop()
       *
       *  \returns reference to i-th integration point
       */
      const LocalCoordinateType &localPoint ( size_t i ) const
      {
        return quad_.point( i );
      }

      /** \copydoc Dune::Fem::IntegrationPointList::id
       */
      size_t id () const
      {
        return quadImp().id();
      }

      /** \copydoc Dune::Fem::IntegrationPointList::order
       */
      int order () const
      {
        return quadImp().order();
      }

      /** \brief obtain GeometryType for this integration point list
       */
      GeometryType geometry () const
      {
        return quadImp().geo();
      }

      /** \brief obtain GeometryType of the corresponding codim-0 the integration
       *         point list belongs to
       *
       *  An element integration point list can return the coordinates of integration
       *  points with resepct to the codim-0 reference element and the reference
       *  element corresponding to the subentity the quadrature actually lives on.
       *  This method returns the geometry of the codim-0 entity.
       *
       *  \note Calling this method yields a virtual function call, so do not
       *        call this method unnecessarily.
       *
       *  \returns GeometryType for this integration point list
       */
      GeometryType elementGeometry () const
      {
        return elementGeometry_;
      }

      size_t cachingPoint( const size_t quadraturePoint ) const
      {
        return quadraturePoint;
      }

      size_t localCachingPoint( const size_t quadraturePoint ) const
      {
        return quadraturePoint;
      }

      /** \brief convenience implementation for Dune::Fem::CachingInterface */
      static constexpr bool twisted () { return false; }

      /** \brief convenience implementation for Dune::Fem::CachingInterface */
      inline int twistId () const { return 0; }

      inline int nCachingPoints () const { return nop(); }
      inline int cachingPointStart () const { return 0; }

      int localFaceIndex () const
      {
        return localFaceIndex_;
      }

    protected:
      /** \brief obtain the actual implementation of the quadrature
       *
       *  \note This method may only be used in derived classes.
       *
       *  \returns a reference to the actual implementation of the quadrature
       */
      const IntegrationPointListType &quadImp() const
      {
        return quad_;
      }

      static GeometryType
      getFaceGeometry ( const GeometryType &elementGeo, const int face )
      {
        // for cube and simplex geom types the dim-1 geom type
        // is also cube or simplex

        static const bool isCube =
              GridPartCapabilities::hasSingleGeometryType< GridPartType >::v &&
              GridPartCapabilities::hasSingleGeometryType< GridPartType >::topologyId == Dune::GeometryTypes::cube(dimension).id();

        static const bool isSimplex =
              GridPartCapabilities::hasSingleGeometryType< GridPartType >::v &&
              GridPartCapabilities::hasSingleGeometryType< GridPartType >::topologyId == Dune::GeometryTypes::simplex(dimension).id();

        if( isCube || isSimplex )
        {
          assert( elementGeo.dim() == dimension );
          if( isCube )
          {
            return Dune::GeometryTypes::cube( dimension-1 );
          }
          else
          {
            assert( isSimplex );
            return Dune::GeometryTypes::simplex( dimension-1 );
          }
        }
        else if( elementGeo.isNone() )
        {
          // if cell geometry is none and dim is 2 then the
          // face is a normal edge which is of type cube
          if( elementGeo.dim() == 2 )
          {
            return Dune::GeometryTypes::cube( 1 );
          }
          else
          {
            return Dune::GeometryTypes::none( elementGeo.dim()-1 );
          }
        }
        else // use reference element to determine type
        {
          assert( ! elementGeo.isNone() );
          typedef Dune::ReferenceElements< RealType, dimension > RefElements;
          return RefElements::general( elementGeo ).type( face, codimension );
        }
      }

    private:
      IntegrationPointListType quad_;
      GeometryType elementGeometry_;
      int localFaceIndex_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ELEMENTPOINTLISTBASE_HH
