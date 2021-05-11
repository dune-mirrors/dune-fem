#ifndef DUNE_FEM_FEMQUADRATURES_HH
#define DUNE_FEM_FEMQUADRATURES_HH

#include <dune/geometry/type.hh>

#include <dune/fem/quadrature/quadratureimp.hh>

// quadrature storage classes
#include "gausspoints.hh"
#include "pyramidpoints.hh"
#include "pardgsimplexquadrature.hh"

namespace Dune
{

  namespace Fem
  {

    struct SimplexMaxOrder
    {
      // uses implementation from parDG
      enum { maxOrder1 = 39, maxOrder2 = 13, maxOrder3 = 12 };

      static int maxOrder( const int dim )
      {
        if( dim == 1 )
          return maxOrder1 ;
        else if( dim == 2 )
          return maxOrder2 ;
        else if( dim == 3 )
          return maxOrder3 ;
        else
        {
          DUNE_THROW(NotImplemented,"SimplexMaxOrder::maxOrder: wrong dimension");
          return -1;
        }
      }

    };

    /*  \class SimplexQuadrature
     *  \ingroup Quadrature
     *  \brief generic quadrature class for simplices
     *
     *  SimplexQuadrature implements the geometry-specific part of the quadrature
     *  and initialises the vector quadrature points and weights.
     *
     *  \note The UG quadrature rules are used here.
     */
    template< class FieldImp, int dim >
    class SimplexQuadrature
    : public QuadratureImp< FieldImp, dim >
    {
    public:
      typedef FieldImp FieldType;

    private:
      typedef SimplexQuadrature< FieldType, dim > ThisType;
      typedef QuadratureImp< FieldType, dim > BaseType;

    public:
      /** \copydoc Dune::Fem::QuadratureImp::CoordinateType */
      typedef typename BaseType :: CoordinateType CoordinateType;

    protected:
      int order_;

    public:
      /** \brief constructor filling the list of points and weights
       *
       *  \param[in]  gemoetry  geometry type for which a quadrature is desired
       *  \param[in]  order     desired order (provided by the user)
       *  \param[in]  id        unique identifier (provided by QuadratureProvider)
       */
      SimplexQuadrature( const GeometryType& geometry, int order, size_t id );

      /** \copydoc Dune::Fem::QuadratureImp::geometry
       */
      virtual GeometryType geometryType () const
      {
        return Dune::GeometryTypes::simplex(dim);
      }

      /** \copydoc Dune::Fem::QuadratureImp::order
       */
      virtual int order () const
      {
        return order_;
      }

      //! maximal order of available quadratures
      static size_t maxOrder ()
      {
        if( dim == 1 )
          return SimplexMaxOrder::maxOrder1;
        if( dim == 2 )
          return SimplexMaxOrder::maxOrder2;
        if( dim == 3 )
          return SimplexMaxOrder::maxOrder3;
        DUNE_THROW( NotImplemented, "SimplexQuadratures from dim > 3 not implemented." );
      }
    };



    /*  \class CubeQuadrature
     *  \ingroup Quadrature
     *  \brief generic quadrature class for cubes
     *
     *  CubeQuadrature implements the geometry-specific part of the quadrature
     *  and initialises the vector quadrature points and weights.
     *
     *  \note The quadrature uses the 1d gauss points (and their tensorial
     *        product) as quadrature points
     */
    template< class FieldImp, int dim >
    class CubeQuadrature
    : public QuadratureImp< FieldImp, dim >
    {
    public:
      typedef FieldImp FieldType;

    private:
      typedef CubeQuadrature< FieldType, dim > ThisType;
      typedef QuadratureImp< FieldType, dim > BaseType;

    public:
      /** \copydoc Dune::Fem::QuadratureImp::CoordinateType */
      typedef typename BaseType :: CoordinateType CoordinateType;

    protected:
      int order_;

    public:
      /** \brief constructor filling the list of points and weights
       *
       *  \param[in]  gemoetry  geometry type for which a quadrature is desired
       *  \param[in]  order     desired order (provided by the user)
       *  \param[in]  id        unique identifier (provided by QuadratureProvider)
       */
      CubeQuadrature( const GeometryType &geometry, int order, size_t id );

      /** \copydoc Dune::Fem::QuadratureImp::geometry */
      virtual GeometryType geometryType () const
      {
        return Dune::GeometryTypes::cube(dim);
      }

      /** \copydoc Dune::Fem::QuadratureImp::order */
      virtual int order () const
      {
        return order_;
      }

      /** \brief maximal order of available quadratures */
      static size_t maxOrder ()
      {
        return GaussPts :: highestOrder;
      }
    };



    /*  \class PrismQuadrature
     *  \ingroup Quadrature
     *  \brief quadrature class for prisms
     *
     *  PrismQuadrature implements the geometry-specific part of the quadrature
     *  and initialises the vector quadrature points and weights.
     *
     *  \note The HD stuff is used here, but needs some rework since only one
     *        rule is provided.
     */
    template< class FieldImp >
    class PrismQuadrature
    : public QuadratureImp< FieldImp, 3 >
    {
    public:
      typedef FieldImp FieldType;

    private:
      typedef PrismQuadrature< FieldType > ThisType;
      typedef QuadratureImp< FieldType, 3 > BaseType;

    public:
      /** \copydoc Dune::Fem::QuadratureImp::CoordinateType */
      typedef typename BaseType :: CoordinateType CoordinateType;

    private:
      int order_;

    public:
      /** \brief constructor filling the list of points and weights
       *
       *  \param[in]  gemoetry  geometry type for which a quadrature is desired
       *  \param[in]  order     desired order (provided by the user)
       *  \param[in]  id        unique identifier (provided by QuadratureProvider)
       */
      PrismQuadrature( const GeometryType &geometry, int order, size_t id );

      /** \copydoc Dune::Fem::QuadratureImp::geometry */
      virtual GeometryType geometryType () const
      {
        return Dune::GeometryTypes::prism;
      }

      /** \copydoc Dune::Fem::QuadratureImp::order */
      virtual int order () const
      {
        return order_;
      }

      /** \brief maximal order of available quadratures */
      static size_t maxOrder ()
      {
        return SimplexMaxOrder::maxOrder2;
      }
    };



    /*  \class PyramidQuadrature
     *  \ingroup Quadrature
     *  \brief quadrature class for pyramids
     *
     *  PyramidQuadrature implements the geometry-specific part of the quadrature
     *  and initialises the vector quadrature points and weights.
     *
     *  \note The HD stuff is used here, but needs some rework since only one
     *        rule is provided.
     */
    template< class FieldImp >
    class PyramidQuadrature
    : public QuadratureImp< FieldImp, 3 >
    {
    public:
      typedef FieldImp FieldType;

    private:
      typedef PyramidQuadrature< FieldType > ThisType;
      typedef QuadratureImp< FieldType, 3 > BaseType;

    public:
      /** \copydoc Dune::Fem::QuadratureImp::CoordinateType */
      typedef typename BaseType :: CoordinateType CoordinateType;

    private:
      int order_;

    public:
      /** \brief constructor filling the list of points and weights
       *
       *  \param[in]  gemoetry  geometry type for which a quadrature is desired
       *  \param[in]  order     desired order (provided by the user)
       *  \param[in]  id        unique identifier (provided by QuadratureProvider)
       */
      PyramidQuadrature( const GeometryType &geometry, int order, size_t id );

      /** \copydoc Dune::Fem::QuadratureImp::geometry */
      virtual GeometryType geometryType () const
      {
        return Dune::GeometryTypes::pyramid;
      }

      /** \copydoc Dune::Fem::QuadratureImp::order */
      virtual int order () const
      {
        return order_;
      }

      /** \brief maximal order of available quadratures */
      static size_t maxOrder ()
      {
        return PyramidPoints :: highest_order;
      }
    };


    /*  \class PolygonQuadrature
     *  \ingroup Quadrature
     *  \brief quadrature class for polygons
     *
     */
    template< class FieldImp, int dim >
    class PolyhedronQuadrature
    : public QuadratureImp< FieldImp, dim >
    {
    public:
      typedef FieldImp FieldType;

    private:
      typedef PolyhedronQuadrature< FieldType, dim > ThisType;
      typedef QuadratureImp< FieldType, dim > BaseType;

    public:
      /** \copydoc Dune::Fem::QuadratureImp::CoordinateType */
      typedef typename BaseType :: CoordinateType CoordinateType;

      using BaseType::addQuadraturePoint;

    private:
      GeometryType geometryType_;
      int order_;

      static const unsigned int topologyId = -1;

    protected:
      using BaseType :: points_;
      using BaseType :: weights_;

    public:
      /** \brief constructor filling the list of points and weights
       *
       *  \param[in]  gemoetry  geometry type for which a quadrature is desired
       *  \param[in]  order     desired order (provided by the user)
       *  \param[in]  id        unique identifier (provided by QuadratureProvider)
       */
      PolyhedronQuadrature( const GeometryType &geometry, int order, size_t id );

      /** \copydoc Dune::Fem::QuadratureImp::geometry */
      virtual GeometryType geometryType () const
      {
        return geometryType_;
      }

      /** \copydoc Dune::Fem::QuadratureImp::order */
      virtual int order () const
      {
        return order_;
      }

      /** \brief maximal order of available quadratures */
      static size_t maxOrder ()
      {
        return SimplexMaxOrder::maxOrder( dim );
      }

      void reset( const int order, const int nop )
      {
        order_ = order;
        points_.clear();
        points_.reserve( nop );
        weights_.clear();
        weights_.reserve( nop );
      }
    };



  } // end namespace Fem

} // end namespace Dune

#include "femquadratures_inline.hh"

#endif // #ifndef DUNE_FEM_FEMQUADRATURES_HH
