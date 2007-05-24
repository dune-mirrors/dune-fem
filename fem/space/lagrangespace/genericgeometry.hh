#ifndef DUNE_GENERICGEOMETRY_HH
#define DUNE_GENERICGEOMETRY_HH

#include <dune/common/misc.hh>
#include <dune/common/geometrytype.hh>

namespace Dune
{

  class PointGeometry
  {
  public:
    enum { dimension = 0 };

    inline static unsigned int numSubEntities ( unsigned int codim )
    {
      return ((codim == 0) ? 1 : 0);
    }
  };



  template< class BaseGeometry >
  class PyramidGeometry
  {
  public:
    typedef BaseGeometry BaseGeometryType;

    enum { dimension = BaseGeometryType :: dimension + 1 };

    inline static unsigned int numSubEntities ( unsigned int codim )
    {
      if( codim > 0 ) {
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



  template< class FirstGeometry, class SecondGeometry >
  class ProductGeometry
  {
  public:
    typedef FirstGeometry FirstGeometryType;
    typedef SecondGeometry SecondGeometryType;

    enum { dimension = FirstGeometryType :: dimension
                     + SecondGeometryType :: dimension };

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



  template< class GeometryType, class FieldVectorImp, unsigned int offset = 0 >
  class LocalCoordinate;


  
  template< class FieldVectorImp, unsigned int offset >
  class LocalCoordinate< PointGeometry, FieldVectorImp, offset >
  {
  public:
    typedef PointGeometry GeometryType;

    enum { dimension = GeometryType :: dimension };

    typedef FieldVectorImp FieldVectorType;
    typedef typename FieldVectorType :: field_type FieldType;

  private:
    typedef LocalCoordinate< GeometryType, FieldVectorType, offset > ThisType;
    
  public:
    inline LocalCoordinate ( FieldVectorType &x )
    {
    }

    inline ThisType& operator= ( FieldType s )
    {
      return *this;
    }
  };

  template< class BaseGeometryType, class FieldVectorImp, unsigned int offset >
  class LocalCoordinate< PyramidGeometry< BaseGeometryType >, FieldVectorImp, offset >
  {
  public:
    typedef PyramidGeometry< BaseGeometryType > GeometryType;

    enum { dimension = GeometryType :: dimension };

    typedef FieldVectorImp FieldVectorType;
    typedef typename FieldVectorType :: field_type FieldType;

    typedef LocalCoordinate< BaseGeometryType, FieldVectorType, offset >
      BaseCoordinateType;

    enum { index = offset + BaseGeometryType :: dimension };
    
  private:
    typedef LocalCoordinate< GeometryType, FieldVectorType, offset > ThisType;

  private:
    FieldType &myCoordinate_;
    BaseCoordinateType baseCoordinate_;

  public:
    inline LocalCoordinate ( FieldVectorType &x )
    : myCoordinate_( x[ index ] ),
      baseCoordinate_( x ) 
    {
    }
    
    inline ThisType& operator= ( FieldType s )
    {
      myCoordinate_ = s;
      baseCoordinate_ = s;
      return *this;
    }

    inline const FieldType& operator*() const
    {
      return myCoordinate_;
    }

    inline FieldType& operator*()
    {
      return myCoordinate_;
    }

    inline const BaseCoordinateType& base () const
    {
      return baseCoordinate_;
    }

    inline BaseCoordinateType& base ()
    {
      return baseCoordinate_;
    }
  };



  template< class FirstGeometryType,
            class SecondGeometryType,
            class FieldVectorImp,
            unsigned int offset >
  class LocalCoordinate< ProductGeometry< FirstGeometryType,
                                          SecondGeometryType >,
                         FieldVectorImp,
                         offset >
  {
  public:
    typedef ProductGeometry< FirstGeometryType, SecondGeometryType >
      GeometryType;

    enum { dimension = GeometryType :: dimension };

    typedef FieldVectorImp FieldVectorType;
    typedef typename FieldVectorType :: field_type FieldType;

    typedef LocalCoordinate< FirstGeometryType, FieldVectorType, offset >
      FirstCoordinateType;
    typedef LocalCoordinate< SecondGeometryType,
                             FieldVectorType,
                             offset + FirstGeometryType :: dimension >
      SecondCoordinateType;

  private:
    typedef LocalCoordinate< GeometryType, FieldVectorType, offset > ThisType;

  private:
    FirstCoordinateType firstCoordinate_;
    SecondCoordinateType secondCoordinate_;

  public:
    inline LocalCoordinate ( FieldVectorType &x )
    : firstCoordinate_( x ),
      secondCoordinate_( x )
    {
    }

    inline ThisType& operator= ( FieldType s )
    {
      firstCoordinate_ = s;
      secondCoordinate_ = s;
      return *this;
    }

    inline const FirstCoordinateType& first () const
    {
      return firstCoordinate_;
    }

    inline FirstCoordinateType& first ()
    {
      return firstCoordinate_;
    }

    inline const SecondCoordinateType& second () const
    {
      return secondCoordinate_;
    }

    inline SecondCoordinateType& second ()
    {
      return secondCoordinate_;
    }
  };
}

#endif
