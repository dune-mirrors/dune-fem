#ifndef DUNE_FEM_ELEMENTPOINTLISTBASE_HH
#define DUNE_FEM_ELEMENTPOINTLISTBASE_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/quadrature/quadrature.hh>

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

      /** \brief constructor
       *  
       *  \param[in]  geometry  geometry type, the quadrature lives on
       *  \param[in]  order     desired minimal order of the quadrature
       */
      ElementPointListBase ( const GeometryType &geometry, int order )
      : quad_( geometry, order )
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

      size_t cachingPoint( const size_t quadraturePoint ) const 
      {
        return quadraturePoint; 
      }

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

    private:
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
      typedef typename GridPartType::GridType::ctype RealType;

      //! dimension of the grid
      static const int dimension = GridPartType::GridType::dimension;
      
      //! type of the integration point list 
      typedef typename IntegrationTraits::IntegrationPointListType IntegrationPointListType;

      typedef typename IntegrationTraits::CoordinateType CoordinateType;
      typedef typename IntegrationPointListType::CoordinateType LocalCoordinateType;

      //! the type of the quadrature point 
      typedef QuadraturePointWrapper< This > QuadraturePointWrapperType;

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
                             const int order )
      : quad_( faceGeo, order ),
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
                             const int order )
      : quad_( getFaceGeometry( elementGeo, localFaceIndex ), order ),
        elementGeometry_( elementGeo ),
        localFaceIndex_( localFaceIndex )
      {}

      const QuadraturePointWrapperType operator[] ( size_t i ) const
      {
        return QuadraturePointWrapperType( *this, i );
      }
     
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

      int localFaceIndex () const
      {
        return localFaceIndex_;
      }

      static GeometryType
      getFaceGeometry ( const GeometryType &elementGeo, const int face )
      {
        typedef Dune::ReferenceElements< RealType, dimension > RefElements;
        return RefElements::general( elementGeo ).type( face, codimension );
      }

    private:
      IntegrationPointListType quad_;
      GeometryType elementGeometry_;
      int localFaceIndex_;
    };

	} // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ELEMENTPOINTLISTBASE_HH
