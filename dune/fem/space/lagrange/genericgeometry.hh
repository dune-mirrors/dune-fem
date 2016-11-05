#ifndef DUNE_FEM_SPACE_LAGRANGE_GENERICGEOMETRY_HH
#define DUNE_FEM_SPACE_LAGRANGE_GENERICGEOMETRY_HH

#include <type_traits>

// dune-common includes
#include <dune/common/fvector.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/misc/metaprogramming.hh>


namespace Dune
{

  namespace Fem
  {

    /** \addtogroup Impl
     *
     *  Generic geometries are a way of constructing new geometries out of
     *  existing ones. This allows for code depending on the geometry to be
     *  written generically. For example, the LagrangeBaseFunction is written
     *  in this generic way. Therefore, you could create a LagrangeBaseFunction
     *  for a 7 dimensional simplex (and no new code is needed).
     *
     *  Generic geometries are created out of simpler ones by the following
     *  rules:
     *  - A point (denoted by \f$p\f$) is a 0-dimensional generic geometry
     *    (PointGeometry). Its reference element is the origin of the
     *    coordinate system.
     *  - Let \f$g\f$ be a \f$d\f$-dimensional generic geometry. Then the pyramid
     *    over \f$g\f$ (denoted by \f$g^\cdot\f$) with the reference geometry
     *    \f[
     *    g^\cdot = \bigl\lbrace (s x,s) \in {R}^{d+1}
     *              \bigm\vert s \in [0,1], x \in g \bigr\rbrace
     *    \f]
     *    is also a generic geometry (PyramidGeometry).
     *  - Let \f$g_1\f$, \f$g_2\f$ be generic geometries. The their product
     *    \f$g_1 \times g_2\f$ (with the obvious reference geometry) is also a
     *    generic geometry (ProductGeometry).
     *  .
     *
     *  Consider the following examples of generic geometries:
     *  - the line \f$p^\cdot\f$,
     *  - the tetrahedron \f$p^{\cdot\cdot\cdot}\f$,
     *  - the cube (3-dimensioal) \f$p^\cdot \times p^\cdot \times p^\cdot\f$,
     *  - the prism \f$p^{\cdot\cdot} \times p^\cdot\f$,
     *  - the 4-sided pyramid \f$(p^\cdot \times p^\cdot)^\cdot\f$.
     *  .
     *  \note All reference geometries defined in dune-grid can be expressed as
     *        generic geometries.
     */

    /** \class PointGeometry
     *  \ingroup Impl
     *  \brief generic geometry modelling a single point
     */
    class PointGeometry
    {
    public:
      /** \brief dimension of the geometry object */
      static const unsigned int dimension = 0;

      template< unsigned int codim >
      class Codim
      {
        static_assert( (codim <= dimension), "Codimension must be less or equal to dimension." );

      public:
        static const unsigned int numSubEntities = ((codim == 0) ? 1 : 0);
      };

      /** \brief number of subentites of a given codimension */
      static unsigned int numSubEntities ( unsigned int codim )
      {
        return ((codim == 0) ? 1 : 0);
      }
    };



    /** \class PyramidGeometry
     *  \ingroup Impl
     *  \brief generic geometry modelling a pyramid over a base geometry
     */
    template< class BaseGeometry >
    class PyramidGeometry
    {
    public:
      /** \brief type of base geometry */
      typedef BaseGeometry BaseGeometryType;

      /** \brief dimension of the geometry object */
      static const unsigned int dimension = BaseGeometryType::dimension + 1;

      template< unsigned int codim >
      class Codim
      {
        static_assert( (codim <= dimension), "Codimension must be less or equal to dimension." );

      public:
        static const unsigned int numSubEntities
          = MetaIf< (codim > 0),
                    MetaPlus< MetaInt< Protect< BaseGeometryType::template Codim, codim-1, PointGeometry::template Codim< 0 >, 0 >::numSubEntities >,
                              MetaInt< Protect< BaseGeometryType::template Codim, codim, PointGeometry::template Codim< 0 >, dimension >::numSubEntities > >,
                    MetaInt< 1 > >::value;
      };

      /** \brief number of subentites of a given codimension */
      static unsigned int numSubEntities ( unsigned int codim )
      {
        if( codim > 0 )
        {
          const unsigned int sameCodimCount = BaseGeometryType::numSubEntities( codim-1 );
          if( codim < dimension )
            return sameCodimCount + BaseGeometryType::numSubEntities( codim );
          else
            return (codim == dimension ? sameCodimCount+1 : 0);
        }
        else
          return 1;
      }
    };



    /** \class ProductGeometry
     *  \ingroup Impl
     *  \brief generic geometry modelling the product of two base geometries
     */
    template< class FirstGeometry, class SecondGeometry >
    class ProductGeometry
    {
    public:
      /** \brief type of the first base geometry */
      typedef FirstGeometry FirstGeometryType;
      /** \brief type of the second base geometry */
      typedef SecondGeometry SecondGeometryType;

      /** \brief dimension of the geometry object */
      static const unsigned int dimension = FirstGeometryType::dimension + SecondGeometryType::dimension;

      template< unsigned int codim >
      class Codim
      {
        static_assert( (codim <= dimension), "Codimension must be less or equal to dimension." );

        template< unsigned int i >
        struct NumSubEntities
        : public MetaInt< FirstGeometryType::template Codim< codim-i >::numSubEntities * SecondGeometryType::template Codim< i >::numSubEntities >
        {};

      public:
        static const unsigned int numSubEntities = Loop< MetaPlus, NumSubEntities, codim >::value;
      };

      /** \brief number of subentites of a given codimension */
      static unsigned int numSubEntities ( unsigned int codim )
      {
        unsigned int cnt = 0;
        for( unsigned int i = 0; i <= codim; ++i )
          cnt += FirstGeometryType::numSubEntities( codim - i ) * SecondGeometryType::numSubEntities( i );
        return cnt;
      }
    };



    template< unsigned int id, unsigned int dim >
    class GeometryWrapper
    {
      static_assert( (id < (1 << dim)), "id too large." );

      static const bool isPrism = ((id >> (dim-1)) != 0);

      typedef GeometryWrapper< (id & ~(1 << (dim-1))), dim-1 > DimensionReductionType;

      template< bool >
      struct Prism
      {
        typedef GeometryWrapper< (id & 1), 1 > LineGeometryType;
        typedef ProductGeometry< typename DimensionReductionType::ImplType, typename LineGeometryType::ImplType >
          ImplType;
      };

      template< bool >
      struct Pyramid
      {
        typedef PyramidGeometry< typename DimensionReductionType::ImplType >
          ImplType;
      };

    public:
      static const unsigned int dimension = dim;

      typedef typename std::conditional< isPrism, Prism< true >, Pyramid< false > >::type::ImplType
        ImplType;
    };

    template< unsigned int id >
    class GeometryWrapper< id, 1 >
    {
      static_assert( (id < 2), "id too large." );

    public:
      static const unsigned int dimension = 1;

      typedef PyramidGeometry< PointGeometry > ImplType;
    };

    template< unsigned int id >
    class GeometryWrapper< id, 0 >
    {
      static_assert( (id < 1), "id too large." );

    public:
      static const unsigned int dimension = 0;

      typedef PointGeometry ImplType;
    };



    // Local Coordinates
    // -----------------

    template< class Geometry, class Field, unsigned int offset = 0 >
    class LocalCoordinate;



    template< class Field, unsigned int offset >
    class LocalCoordinate< PointGeometry, Field, offset >
    {
      typedef LocalCoordinate< PointGeometry, Field, offset > ThisType;

    public:
      typedef PointGeometry GeometryType;

      static const unsigned int dimension = GeometryType::dimension;

      typedef Field FieldType;

      inline LocalCoordinate ()
      {}

      template< int sz >
      explicit LocalCoordinate ( const FieldVector< FieldType, sz > &x )
      {
        static_assert( (sz >= offset + dimension), "Invalid vector size" );
      }

      ThisType &operator= ( const FieldType s )
      {
        return *this;
      }

      template< int sz >
      ThisType &operator= ( const FieldVector< FieldType, sz > &x )
      {
        static_assert( (sz >= offset + dimension), "Invalid vector size" );
        return *this;
      }

      ThisType &operator= ( const ThisType &v )
      {
        return *this;
      }

      ThisType &operator*= ( const FieldType s )
      {
        return *this;
      }

      ThisType &operator+= ( const ThisType &v )
      {
        return *this;
      }

      ThisType &operator-= ( const ThisType &v )
      {
        return *this;
      }

      const FieldType &operator[] ( const unsigned int i ) const
      {
        DUNE_THROW( RangeError, "LocalCoordinate: No such index." );
      }

      FieldType &operator[] ( const unsigned int i )
      {
        DUNE_THROW( RangeError, "LocalCoordinate: No such index." );
      }
    };



    template< class BaseGeometry, class Field, unsigned int offset >
    class LocalCoordinate< PyramidGeometry< BaseGeometry >, Field, offset >
    {
      typedef LocalCoordinate< PyramidGeometry< BaseGeometry >, Field, offset > ThisType;

    public:
      typedef BaseGeometry BaseGeometryType;

      typedef PyramidGeometry< BaseGeometryType > GeometryType;

      static const unsigned int dimension = GeometryType::dimension;

      typedef Field FieldType;

      typedef LocalCoordinate< BaseGeometry, FieldType, offset > BaseCoordinateType;

      static const unsigned int index = offset + BaseGeometryType::dimension;

      LocalCoordinate ()
      {}

      template< int sz >
      explicit LocalCoordinate ( const FieldVector< FieldType, sz > &x )
      : myCoordinate_( x[ index ] ),
        baseCoordinate_( x )
      {
        static_assert( (sz >= offset + dimension), "Invalid vector size" );
      }

      ThisType &operator= ( const FieldType s )
      {
        myCoordinate_ = s;
        baseCoordinate_ = s;
        return *this;
      }

      template< int sz >
      ThisType &operator= ( const FieldVector< FieldType, sz > &x )
      {
        static_assert( (sz >= offset + dimension), "Invalid vector size" );

        myCoordinate_ = x[ index ];
        baseCoordinate_ = x;
        return *this;
      }

      ThisType &operator= ( const ThisType &v )
      {
        myCoordinate_ = v.myCoordinate_;
        baseCoordinate_ = v.baseCoordinate_;
        return *this;
      }

      ThisType &operator*= ( const FieldType s )
      {
        myCoordinate_ *= s;
        baseCoordinate_ *= s;
        return *this;
      }

      ThisType &operator+= ( const ThisType &v )
      {
        myCoordinate_ += v.myCoordinate_;
        baseCoordinate_ += v.baseCoordinate_;
        return *this;
      }

      ThisType &operator-= ( const ThisType &v )
      {
        myCoordinate_ -= v.myCoordinate_;
        baseCoordinate_ -= v.baseCoordinate_;
        return *this;
      }

      const FieldType &operator[] ( const unsigned int i ) const
      {
        if( i == index )
          return myCoordinate_;
        else
          return baseCoordinate_[ i ];
      }

      FieldType &operator[] ( const unsigned int i )
      {
        if( i == index )
          return myCoordinate_;
        else
          return baseCoordinate_[ i ];
      }

      const FieldType &operator* () const
      {
        return myCoordinate_;
      }

      FieldType &operator* ()
      {
        return myCoordinate_;
      }

      const BaseCoordinateType &base () const
      {
        return baseCoordinate_;
      }

      BaseCoordinateType &base ()
      {
        return baseCoordinate_;
      }

    private:
      FieldType myCoordinate_;
      BaseCoordinateType baseCoordinate_;
    };



    template< class FirstGeometry, class SecondGeometry, class Field, unsigned int offset >
    class LocalCoordinate< ProductGeometry< FirstGeometry, SecondGeometry >, Field, offset >
    {
      typedef LocalCoordinate< ProductGeometry< FirstGeometry, SecondGeometry >, Field, offset > ThisType;

    public:
      typedef FirstGeometry FirstGeometryType;
      typedef SecondGeometry SecondGeometryType;
      typedef ProductGeometry< FirstGeometryType, SecondGeometryType > GeometryType;

      static const unsigned int dimension = GeometryType::dimension;

      typedef Field FieldType;

    protected:
      static const unsigned int firstOffset = offset;
      static const unsigned int secondOffset = offset + FirstGeometryType::dimension;

    public:
      typedef LocalCoordinate< FirstGeometryType, FieldType, firstOffset > FirstCoordinateType;
      typedef LocalCoordinate< SecondGeometryType, FieldType, secondOffset > SecondCoordinateType;

      LocalCoordinate ()
      {}

      template< int sz >
      explicit LocalCoordinate ( const FieldVector< FieldType, sz > &x )
      : firstCoordinate_( x ),
        secondCoordinate_( x )
      {
        static_assert( (sz >= offset + dimension), "Invalid vector size" );
      }

      ThisType &operator= ( const FieldType s )
      {
        firstCoordinate_ = s;
        secondCoordinate_ = s;
        return *this;
      }

      template< int sz >
      ThisType &operator= ( const FieldVector< FieldType, sz > &x )
      {
        static_assert( (sz >= offset + dimension), "Invalid vector size" );

        firstCoordinate_ = x;
        secondCoordinate_ = x;
        return *this;
      }

      ThisType &operator= ( const ThisType &v )
      {
        firstCoordinate_ = v;
        secondCoordinate_ = v;
        return *this;
      }

      ThisType &operator*= ( const FieldType s )
      {
        firstCoordinate_ *= s;
        secondCoordinate_ *= s;
        return *this;
      }

      ThisType &operator+= ( const ThisType &v )
      {
        firstCoordinate_ += v;
        secondCoordinate_ += v;
        return *this;
      }

      ThisType &operator-= ( const ThisType &v )
      {
        firstCoordinate_ -= v;
        secondCoordinate_ -= v;
        return *this;
      }

      const FieldType &operator[] ( const unsigned int i ) const
      {
        if( i < secondOffset )
          return firstCoordinate_[ i ];
        else
          return secondCoordinate_[ i ];
      }

      FieldType &operator[] ( const unsigned int i )
      {
        if( i < secondOffset )
          return firstCoordinate_[ i ];
        else
          return secondCoordinate_[ i ];
      }

      const FirstCoordinateType &first () const
      {
        return firstCoordinate_;
      }

      FirstCoordinateType &first ()
      {
        return firstCoordinate_;
      }

      const SecondCoordinateType &second () const
      {
        return secondCoordinate_;
      }

      SecondCoordinateType &second ()
      {
        return secondCoordinate_;
      }

    private:
      FirstCoordinateType firstCoordinate_;
      SecondCoordinateType secondCoordinate_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_GENERICGEOMETRY_HH
