#ifndef DUNE_GENERICGEOMETRY_HH
#define DUNE_GENERICGEOMETRY_HH

#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>
#include <dune/common/static_assert.hh>

#include <dune/fem/misc/metaprogramming.hh>

namespace Dune
{

  /** \addtogroup GenericGeometry
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
   *  \ingroup GenericGeometry
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
      dune_static_assert( (codim <= dimension), "Codimension must be less or equal to dimension." );

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
   *  \ingroup GenericGeometry
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
      dune_static_assert( (codim <= dimension), "Codimension must be less or equal to dimension." );

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
   *  \ingroup GenericGeometry
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
      dune_static_assert( (codim <= dimension), "Codimension must be less or equal to dimension." );

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



  template< GeometryType::BasicType type, unsigned int dim >
  class GeometryWrapper;



  template< unsigned int dim >
  class GeometryWrapper< GeometryType::simplex, dim >
  {
    typedef GeometryWrapper< GeometryType::simplex, dim > ThisType;
    typedef GeometryWrapper< GeometryType::simplex, dim-1 > DimensionReductionType;
    
    // disallow use for dim > 3
    dune_static_assert( (dim <= 3), "Dimension must be less or equal to 3." );

  public:
    static const unsigned int dimension = dim;

    typedef PyramidGeometry< typename DimensionReductionType::GenericGeometryType > GenericGeometryType;
    
    static void duneSubEntity ( const unsigned int codim, unsigned int &subEntity )
    {}
  };


  
  template<>
  class GeometryWrapper< GeometryType::simplex, 0 >
  {
    typedef GeometryWrapper< GeometryType::simplex, 0 > ThisType;

  public:
    static const unsigned int dimension = 0;

    typedef PointGeometry GenericGeometryType;
    
    static void duneSubEntity ( const unsigned int codim, unsigned int &subEntity )
    {}
  };


  
  template<>
  class GeometryWrapper< GeometryType::simplex, 2 >
  {
    typedef GeometryWrapper< GeometryType::simplex, 2 > ThisType;
    typedef GeometryWrapper< GeometryType::simplex, 1 > DimensionReductionType;

  public:
    static const unsigned int dimension = 2;

    typedef PyramidGeometry< DimensionReductionType::GenericGeometryType > GenericGeometryType;

    static void duneSubEntity ( const unsigned int codim, unsigned int &subEntity )
    {
      subEntity = (codim == 1 ? 2 - subEntity : subEntity);
    }
  };
  

  
  template<>
  class GeometryWrapper< GeometryType::simplex, 3 >
  {
    typedef GeometryWrapper< GeometryType::simplex, 3 > ThisType;
    typedef GeometryWrapper< GeometryType::simplex, 2 > DimensionReductionType;

  public:
    static const unsigned int dimension = 3;

    typedef PyramidGeometry< DimensionReductionType::GenericGeometryType > GenericGeometryType;

    static void duneSubEntity ( const unsigned int codim, unsigned int &subEntity )
    {
      if( codim == 1 )
        subEntity = 3 - subEntity;
      else if( codim == 2 )
      {
        if( subEntity == 1 )
          subEntity = 2;
        else if( subEntity == 2 )
          subEntity = 1;
      }
    }
  };



  template< unsigned int dim >
  class GeometryWrapper< GeometryType::cube, dim >
  {
    typedef GeometryWrapper< GeometryType::cube, dim > ThisType;
    typedef GeometryWrapper< GeometryType::simplex, 1 > LineGeometryType;
    typedef GeometryWrapper< GeometryType::cube, dim-1 > DimensionReductionType;

  public:
    static const unsigned int dimension = dim;

    // disallow use for dim > 3
    dune_static_assert( (dim <= 3), "Dimension must be less or equal to 3." );

    typedef ProductGeometry< typename DimensionReductionType::GenericGeometryType, typename LineGeometryType::GenericGeometryType >
      GenericGeometryType;

    static void duneSubEntity ( const unsigned int codim, unsigned int &subEntity )
    {}
  };


  
  template<>
  class GeometryWrapper< GeometryType::cube, 0 >
  {
    typedef GeometryWrapper< GeometryType::cube, 0 > ThisType;
    typedef GeometryWrapper< GeometryType::simplex, 1 > LineGeometryType;

  public:
    static const unsigned int dimension = 0;
    
    typedef PointGeometry GenericGeometryType;
    
    static void duneSubEntity ( const unsigned int codim, unsigned int &subEntity )
    {}
  };


  
  template<>
  class GeometryWrapper< GeometryType::cube, 3 >
  {
    typedef GeometryWrapper< GeometryType::cube, 3 > ThisType;
    typedef GeometryWrapper< GeometryType::simplex, 1 > LineGeometryType;
    typedef GeometryWrapper< GeometryType::cube, 2 > DimensionReductionType;

  public:
    static const unsigned int dimension = 3;

    typedef ProductGeometry< DimensionReductionType::GenericGeometryType, LineGeometryType::GenericGeometryType >
      GenericGeometryType;

    static void duneSubEntity ( const unsigned int codim, unsigned int &subEntity )
    {
      if( codim != 2 )
        return;

      if( subEntity == 6 )
        subEntity = 8;
      else if( subEntity == 7 )
        subEntity = 9;
      else if( subEntity == 8 )
        subEntity = 6;
      else if( subEntity == 9 )
        subEntity = 7;
    }
  };



  template<>
  class GeometryWrapper< GeometryType::pyramid, 3 >
  {
    typedef GeometryWrapper< GeometryType::pyramid, 3 > ThisType;
    typedef GeometryWrapper< GeometryType::cube, 2 > BaseGeometryType;

  public:
    static const unsigned int dimension = 3;

    typedef PyramidGeometry< BaseGeometryType::GenericGeometryType > GenericGeometryType;

    static void duneSubEntity ( const unsigned int codim, unsigned int &subEntity )
    {}
  };



  template<>
  class GeometryWrapper< GeometryType::prism, 3 >
  {
    typedef GeometryWrapper< GeometryType::prism, 3 > ThisType;
    typedef GeometryWrapper< GeometryType::simplex, 2 > FirstGeometryType;
    typedef GeometryWrapper< GeometryType::simplex, 1 > SecondGeometryType;

  public:
    static const unsigned int dimension = 3;
    
    typedef ProductGeometry< FirstGeometryType::GenericGeometryType, SecondGeometryType::GenericGeometryType >
      GenericGeometryType;

    static void duneSubEntity ( const unsigned int codim, unsigned int &subEntity )
    {}
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
      dune_static_assert( (sz >= offset + dimension), "Invalid vector size" );
    }

    ThisType &operator= ( const FieldType s )
    {
      return *this;
    }

    template< int sz >
    ThisType &operator= ( const FieldVector< FieldType, sz > &x )
    {
      dune_static_assert( (sz >= offset + dimension), "Invalid vector size" );
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
      dune_static_assert( (sz >= offset + dimension), "Invalid vector size" );
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
      dune_static_assert( (sz >= offset + dimension), "Invalid vector size" );
      
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
      dune_static_assert( (sz >= offset + dimension), "Invalid vector size" );
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
      dune_static_assert( (sz >= offset + dimension), "Invalid vector size" );
      
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

}

#endif
