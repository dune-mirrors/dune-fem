#ifndef DUNE_GENERICGEOMETRY_HH
#define DUNE_GENERICGEOMETRY_HH

#include <dune/common/misc.hh>
#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>

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
    enum
    {
      /** \brief dimension of the geometry object */
      dimension = 0
    };

    template< unsigned int codim >
    class Codim
    {
    private:
      CompileTimeChecker< (codim <= dimension) >
        __CODIM_MUST_BE_LESS_EQUAL_TO_DIMENSION__;

    public:
      enum { numSubEntities = ((codim == 0) ? 1 : 0) };
    };

    /** \brief number of subentites of a given codimension */
    inline static unsigned int numSubEntities ( unsigned int codim )
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

    enum
    {
      /** \brief dimension of the geometry object */
      dimension = BaseGeometryType :: dimension + 1
    };

    template< unsigned int codim >
    class Codim
    {
    private:
      CompileTimeChecker< (codim <= dimension) >
        __CODIM_MUST_BE_LESS_EQUAL_TO_DIMENSION__;

    public:
      enum
      {
        numSubEntities
          = MetaIf
            <
              (codim > 0),
              MetaPlus
              <
                MetaInt
                <
                  Protect< BaseGeometryType :: template Codim, codim - 1,
                           PointGeometry :: template Codim< 0 >, 0
                         > :: numSubEntities
                >,
                MetaInt
                <
                  Protect< BaseGeometryType :: template Codim, codim,
                           PointGeometry :: template Codim< 0 >, dimension
                         > :: numSubEntities
                >
              >,
              MetaInt< 1 >
            > :: value
      };
    };

    /** \brief number of subentites of a given codimension */
    inline static unsigned int numSubEntities ( unsigned int codim )
    {
      if( codim > 0 )
      {
        const unsigned int sameCodimCount
          = BaseGeometryType :: numSubEntities( codim - 1 );
        if( codim < dimension )
          return sameCodimCount + BaseGeometryType :: numSubEntities( codim );
        else
          return (codim == dimension ? sameCodimCount + 1 : 0);
      } else
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

    enum
    {
      /** \brief dimension of the geometry object */
      dimension = FirstGeometryType :: dimension
                  + SecondGeometryType :: dimension
    };

    template< unsigned int codim >
    class Codim
    {
    private:
      CompileTimeChecker< (codim <= dimension) >
        __CODIM_MUST_BE_LESS_EQUAL_TO_DIMENSION__;

    private:
      template< unsigned int i >
      struct NumSubEntities
      : public MetaInt
        < FirstGeometryType :: template Codim< codim - i > :: numSubEntities
          * SecondGeometryType :: template Codim< i > :: numSubEntities
        >
      {
      };

    public:
      enum { numSubEntities = Loop< MetaPlus, NumSubEntities, codim > :: value };
    };

    /** \brief number of subentites of a given codimension */
    inline static unsigned int numSubEntities ( unsigned int codim )
    {
      enum { firstDimension = FirstGeometryType :: dimension };
        
      unsigned int cnt = 0;
      for( unsigned int i = 0; i <= codim; ++i ) {
        cnt += FirstGeometryType :: numSubEntities( codim - i )
             * SecondGeometryType :: numSubEntities( i );
      }
      return cnt;
    }
  };



  template< GeometryType :: BasicType type, unsigned int dim >
  class GeometryWrapper;



  template< unsigned int dim >
  class GeometryWrapper< GeometryType :: simplex, dim >
  {
  public:
    enum { dimension = dim };

  private:
    typedef GeometryWrapper< GeometryType :: simplex, dimension > ThisType;
    typedef GeometryWrapper< GeometryType :: simplex,  dimension - 1 >
      DimensionReductionType;
    
    // disallow use for dim > 3
    CompileTimeChecker< (dim <= 3) > __assert_dimension__;


  public:
    typedef PyramidGeometry< typename DimensionReductionType :: GenericGeometryType >
      GenericGeometryType;
    
    static inline void duneSubEntity ( const unsigned int codim,
                                       unsigned int &subEntity )
    {
    }
  };


  
  template<>
  class GeometryWrapper< GeometryType :: simplex, 0 >
  {
  public:
    enum { dimension = 0 };

  private:
    typedef GeometryWrapper< GeometryType :: simplex,  dimension > ThisType;

  public:
    typedef PointGeometry GenericGeometryType;
    
    static inline void duneSubEntity ( const unsigned int codim,
                                       unsigned int &subEntity )
    {
    }
  };


  
  template<>
  class GeometryWrapper< GeometryType :: simplex, 2 >
  {
  public:
    enum { dimension = 2 };

  private:
    typedef GeometryWrapper< GeometryType :: simplex, dimension > ThisType;
    typedef GeometryWrapper< GeometryType :: simplex,  dimension - 1 >
      DimensionReductionType;

  public:
    typedef PyramidGeometry< DimensionReductionType :: GenericGeometryType >
      GenericGeometryType;

    static inline void duneSubEntity ( const unsigned int codim,
                                       unsigned int &subEntity )
    {
      if( codim == 1 )
        subEntity = 2 - subEntity;
    }
  };
  

  
  template<>
  class GeometryWrapper< GeometryType :: simplex, 3 >
  {
  public:
    enum { dimension = 3 };

  private:
    typedef GeometryWrapper< GeometryType :: simplex, dimension > ThisType;
    typedef GeometryWrapper< GeometryType :: simplex,  dimension - 1 >
      DimensionReductionType;

  public:
    typedef PyramidGeometry< DimensionReductionType :: GenericGeometryType >
      GenericGeometryType;

    static inline void duneSubEntity ( const unsigned int codim,
                                       unsigned int &subEntity )
    {
      if( codim == 1 )
        subEntity = 3 - subEntity;
      else if( codim == 2 ) {
        if( subEntity == 1 )
          subEntity = 2;
        else if( subEntity == 2 )
          subEntity = 1;
      }
    }
  };



  template< unsigned int dim >
  class GeometryWrapper< GeometryType :: cube, dim >
  {
  public:
    enum { dimension = dim };

  private:
    typedef GeometryWrapper< GeometryType :: cube, dimension > ThisType;
    typedef GeometryWrapper< GeometryType :: simplex, 1 > LineGeometryType;
    typedef GeometryWrapper< GeometryType :: cube, dimension - 1 >
      DimensionReductionType;

    // disallow use for dim > 3
    CompileTimeChecker< (dim <= 3) > __assert_dimension__;

  public:
    typedef ProductGeometry< typename DimensionReductionType :: GenericGeometryType,
                             typename LineGeometryType :: GenericGeometryType >
      GenericGeometryType;

    static inline void duneSubEntity ( const unsigned int codim,
                                       unsigned int &subEntity )
    {
    }
  };


  
  template<>
  class GeometryWrapper< GeometryType :: cube, 0 >
  {
  public:
    enum { dimension = 0 };
    
  private:
    typedef GeometryWrapper< GeometryType :: cube, dimension > ThisType;
    typedef GeometryWrapper< GeometryType :: simplex, 1 > LineGeometryType;

  public:
    typedef PointGeometry GenericGeometryType;
    
    static inline void duneSubEntity ( const unsigned int codim,
                                       unsigned int &subEntity )
    {
    }
  };


  
  template<>
  class GeometryWrapper< GeometryType :: cube, 3 >
  {
  public:
    enum { dimension = 3 };

  private:
    typedef GeometryWrapper< GeometryType :: cube, dimension > ThisType;
    typedef GeometryWrapper< GeometryType :: simplex, 1 > LineGeometryType;
    typedef GeometryWrapper< GeometryType :: cube, dimension - 1 >
      DimensionReductionType;

  public:
    typedef ProductGeometry< DimensionReductionType :: GenericGeometryType,
                             LineGeometryType :: GenericGeometryType >
      GenericGeometryType;

    static inline void duneSubEntity ( const unsigned int codim,
                                       unsigned int &subEntity )
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
  class GeometryWrapper< GeometryType :: pyramid, 3 >
  {
  public:
    enum { dimension = 3 };

  private:
    typedef GeometryWrapper< GeometryType :: pyramid, 3 > ThisType;
    typedef GeometryWrapper< GeometryType :: cube, 2 > BaseGeometryType;

  public:
    typedef PyramidGeometry< BaseGeometryType :: GenericGeometryType >
      GenericGeometryType;

    static inline void duneSubEntity ( const unsigned int codim,
                                       unsigned int &subEntity )
    {
    }
  };



  template<>
  class GeometryWrapper< GeometryType :: prism, 3 >
  {
  public:
    enum { dimension = 3 };
    
  private:
    typedef GeometryWrapper< GeometryType :: prism, 3 > ThisType;
    typedef GeometryWrapper< GeometryType :: simplex, 2 >  FirstGeometryType;
    typedef GeometryWrapper< GeometryType :: simplex, 1 > SecondGeometryType;

  public:
    typedef ProductGeometry< FirstGeometryType :: GenericGeometryType,
                             SecondGeometryType :: GenericGeometryType >
      GenericGeometryType;

    static inline void duneSubEntity ( const unsigned int codim,
                                       unsigned int &subEntity )
    {
    }
  };



  // Local Coordinates
  // -----------------
  
  template< class Geometry, class Field, unsigned int offset = 0 >
  class LocalCoordinate;


  
  template< class Field, unsigned int offset >
  class LocalCoordinate< PointGeometry, Field, offset >
  {
  public:
    typedef PointGeometry GeometryType;

    enum { dimension = GeometryType :: dimension };

    typedef Field FieldType;

  private:
    typedef LocalCoordinate< GeometryType, FieldType, offset > ThisType;
    
  public:
    inline LocalCoordinate ()
    {}
    
    template< int sz >
    inline explicit LocalCoordinate ( const FieldVector< FieldType, sz > &x )
    {
      typedef CompileTimeChecker< (sz >= offset + dimension) >
        __CHECK_VECTOR_SIZE__;
    }

    inline ThisType &operator= ( const FieldType s )
    {
      return *this;
    }

    template< int sz >
    inline ThisType &operator= ( const FieldVector< FieldType, sz > &x )
    {
      typedef CompileTimeChecker< (sz >= offset + dimension) >
        __CHECK_VECTOR_SIZE__;
      return *this;
    }
    
    inline ThisType &operator= ( const ThisType &v )
    {
      return *this;
    }

    inline ThisType &operator*= ( const FieldType s )
    {
      return *this;
    }

    inline ThisType &operator+= ( const ThisType &v )
    {
      return *this;
    }

    inline ThisType &operator-= ( const ThisType &v )
    {
      return *this;
    }

    inline const FieldType &operator[] ( const unsigned int i ) const
    {
      DUNE_THROW( RangeError, "LocalCoordinate: No such index." );
    }

    inline FieldType &operator[] ( const unsigned int i )
    {
      DUNE_THROW( RangeError, "LocalCoordinate: No such index." );
    }
  };



  template< class BaseGeometry, class Field, unsigned int offset >
  class LocalCoordinate< PyramidGeometry< BaseGeometry >, Field, offset >
  {
  public:
    typedef BaseGeometry BaseGeometryType;
    
    typedef PyramidGeometry< BaseGeometryType > GeometryType;

    enum { dimension = GeometryType :: dimension };

    typedef Field FieldType;
    
    typedef LocalCoordinate< BaseGeometry, FieldType, offset >
      BaseCoordinateType;

    enum { index = offset + BaseGeometryType :: dimension };
    
  private:
    typedef LocalCoordinate< GeometryType, FieldType, offset > ThisType;

  private:
    FieldType myCoordinate_;
    BaseCoordinateType baseCoordinate_;

  public:
    inline LocalCoordinate ()
    {}
    
    template< int sz >
    inline explicit LocalCoordinate ( const FieldVector< FieldType, sz > &x )
    : myCoordinate_( x[ index ] ),
      baseCoordinate_( x ) 
    {
      typedef CompileTimeChecker< (sz >= offset + dimension) >
        __CHECK_VECTOR_SIZE__;
    }
    
    inline ThisType &operator= ( const FieldType s )
    {
      myCoordinate_ = s;
      baseCoordinate_ = s;
      return *this;
    }
    
    template< int sz >
    inline ThisType &operator= ( const FieldVector< FieldType, sz > &x )
    {
      typedef CompileTimeChecker< (sz >= offset + dimension) >
        __CHECK_VECTOR_SIZE__;
      
      myCoordinate_ = x[ index ];
      baseCoordinate_ = x;
      return *this;
    }

    inline ThisType &operator= ( const ThisType &v )
    {
      myCoordinate_ = v.myCoordinate_;
      baseCoordinate_ = v.baseCoordinate_;
      return *this;
    }

    inline ThisType &operator*= ( const FieldType s )
    {
      myCoordinate_ *= s;
      baseCoordinate_ *= s;
      return *this;
    }
    
    inline ThisType &operator+= ( const ThisType &v )
    {
      myCoordinate_ += v.myCoordinate_;
      baseCoordinate_ += v.baseCoordinate_;
      return *this;
    }
    
    inline ThisType &operator-= ( const ThisType &v )
    {
      myCoordinate_ -= v.myCoordinate_;
      baseCoordinate_ -= v.baseCoordinate_;
      return *this;
    }
    
    inline const FieldType &operator[] ( const unsigned int i ) const
    {
      if( i == index )
        return myCoordinate_;
      else
        return baseCoordinate_[ i ];
    }

    inline FieldType &operator[] ( const unsigned int i )
    {
      if( i == index )
        return myCoordinate_;
      else
        return baseCoordinate_[ i ];
    }

    inline const FieldType &operator* () const
    {
      return myCoordinate_;
    }

    inline FieldType &operator* ()
    {
      return myCoordinate_;
    }

    inline const BaseCoordinateType &base () const
    {
      return baseCoordinate_;
    }

    inline BaseCoordinateType &base ()
    {
      return baseCoordinate_;
    }
  };



  template< class FirstGeometry, class SecondGeometry, class Field,
            unsigned int offset >
  class LocalCoordinate
    < ProductGeometry< FirstGeometry, SecondGeometry >, Field, offset >
  {
  public:
    typedef FirstGeometry FirstGeometryType;
    typedef SecondGeometry SecondGeometryType;
    typedef ProductGeometry< FirstGeometryType, SecondGeometryType >
      GeometryType;

    enum { dimension = GeometryType :: dimension };

    typedef Field FieldType;

  protected:
    enum {
      firstOffset = offset,
      secondOffset = offset + FirstGeometryType :: dimension
    };

  public:
    typedef LocalCoordinate< FirstGeometryType, FieldType, firstOffset >
      FirstCoordinateType;
    typedef LocalCoordinate< SecondGeometryType, FieldType, secondOffset >
      SecondCoordinateType;

  private:
    typedef LocalCoordinate< GeometryType, FieldType, offset > ThisType;

  private:
    FirstCoordinateType firstCoordinate_;
    SecondCoordinateType secondCoordinate_;

  public:
    inline LocalCoordinate ()
    {}
    
    template< int sz >
    inline explicit LocalCoordinate ( const FieldVector< FieldType, sz > &x )
    : firstCoordinate_( x ),
      secondCoordinate_( x )
    {
      typedef CompileTimeChecker< (sz >= offset + dimension) >
        __CHECK_VECTOR_SIZE__;
    }

    inline ThisType &operator= ( const FieldType s )
    {
      firstCoordinate_ = s;
      secondCoordinate_ = s;
      return *this;
    }
    
    template< int sz >
    inline ThisType &operator= ( const FieldVector< FieldType, sz > &x )
    {
      typedef CompileTimeChecker< (sz >= offset + dimension) >
        __CHECK_VECTOR_SIZE__;
      
      firstCoordinate_ = x;
      secondCoordinate_ = x;
      return *this;
    }

    inline ThisType &operator= ( const ThisType &v )
    {
      firstCoordinate_ = v;
      secondCoordinate_ = v;
      return *this;
    }

    inline ThisType &operator*= ( const FieldType s )
    {
      firstCoordinate_ *= s;
      secondCoordinate_ *= s;
      return *this;
    }
    
    inline ThisType &operator+= ( const ThisType &v )
    {
      firstCoordinate_ += v;
      secondCoordinate_ += v;
      return *this;
    }
    
    inline ThisType &operator-= ( const ThisType &v )
    {
      firstCoordinate_ -= v;
      secondCoordinate_ -= v;
      return *this;
    }

    inline const FieldType &operator[] ( const unsigned int i ) const
    {
      if( i < secondOffset )
        return firstCoordinate_[ i ];
      else
        return secondCoordinate_[ i ];
    }

    inline FieldType &operator[] ( const unsigned int i )
    {
      if( i < secondOffset )
        return firstCoordinate_[ i ];
      else
        return secondCoordinate_[ i ];
    }

    inline const FirstCoordinateType &first () const
    {
      return firstCoordinate_;
    }

    inline FirstCoordinateType &first ()
    {
      return firstCoordinate_;
    }

    inline const SecondCoordinateType &second () const
    {
      return secondCoordinate_;
    }

    inline SecondCoordinateType &second ()
    {
      return secondCoordinate_;
    }
  };

}

#endif
