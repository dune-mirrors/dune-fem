#ifndef DUNE_LAGRANGESPACE_GENERICLAGRANGEPOINTS_HH
#define DUNE_LAGRANGESPACE_GENERICLAGRANGEPOINTS_HH

#include <iostream>

#include <dune/config.h>
#include <dune/common/fvector.hh>

#include "genericgeometry.hh"

namespace Dune
{

  template< class GenericGeometryType, unsigned int order >
  class GenericLagrangePoint;

  template< class GeometryType, unsigned int order, unsigned int codim >
  class GenericLagrangePointCodim;
 

  
  template< unsigned int order >
  class GenericLagrangePoint< PointGeometry, order >
  {
  public:
    typedef PointGeometry GeometryType;
    enum { dimension = GeometryType :: dimension };
    typedef FieldVector< unsigned int, dimension > DofCoordinateType;
 
    enum { polynomialOrder = order };
    
    template< class GenericGeometryType, unsigned int porder >
    friend class GenericLagrangePoint;
   
    template< class FunctionSpaceType, class GeometryType, unsigned int porder >
    friend class GenericLagrangeBaseFunction;

  private:
    typedef GenericLagrangePoint< GeometryType, polynomialOrder > ThisType;
 
  public:
    enum { numLagrangePoints = 1 };

  protected:
    DofCoordinateType dofCoordinate_;
    LocalCoordinate< GeometryType, DofCoordinateType > localDofCoordinate_;

  public:
    template< unsigned int codim >
    class Codim
    {
    public:
      static inline unsigned int maxDofs ()
      {
        return ((codim == 0) ? 1 : 0);
      }
    };

  public:
    inline GenericLagrangePoint ( unsigned int index )
    : localDofCoordinate_( dofCoordinate_ )
    {
      dofCoordinate( index, localDofCoordinate_ );
    }
    
    inline GenericLagrangePoint ( const ThisType &point )
    : dofCoordinate_( point.dofCoordinate_ ),
      localDofCoordinate_( dofCoordinate_ )
    {
    }

    template< class LocalCoordinateType >
    static inline void dofSubEntity( LocalCoordinateType &coordinate,
                                     unsigned int &codim,
                                     unsigned int &subEntity )
    {
      codim = 0;
      subEntity = 0;
    }

    template< class LocalCoordinateType >
    static inline void dofSubEntity( LocalCoordinateType &coordinate,
                                     unsigned int &codim,
                                     unsigned int &subEntity,
                                     unsigned int &dofNumber )
    {
      codim = 0;
      subEntity = 0;
      dofNumber = 0;
    }

    inline void dofSubEntity ( unsigned int &codim,
                               unsigned int &subEntity )
    {
      dofSubEntity( localDofCoordinate_, codim, subEntity );
    }
    
    inline void dofSubEntity ( unsigned int &codim,
                               unsigned int &subEntity,
                               unsigned int &dofNumber )
    {
      dofSubEntity( localDofCoordinate_, codim, subEntity, dofNumber );
    }
   
    static inline unsigned int entityDofNumber ( unsigned int codim,
                                                 unsigned int subEntity,
                                                 unsigned int dof )
    {
      //assert( (codim == 0) && (subEntity == 0) && (dof == 0) );
      return 0;
    }
    
    template< class LocalCoordinateType >
    static inline unsigned int height ( LocalCoordinateType &coordinate )
    {
      return polynomialOrder;
    }

    inline unsigned int height ()
    {
      return height( localDofCoordinate_ );
    }

    template< class FieldType >
    inline void local ( FieldVector< FieldType, dimension > &coordinate ) const
    {
      const FieldType factor = 1 / (FieldType)polynomialOrder;
      for( int i = 0; i < dimension; ++i )
        coordinate[ i ] = factor * dofCoordinate_[ i ]; 
    }
    
    static inline unsigned int maxDofs ( unsigned int codim )
    {
      return ((codim == 0) ? 1 : 0);
    }


    static inline unsigned int numDofs ( unsigned int codim,
                                         unsigned int subEntity )
    {
      return ((codim == 0) ? 1 : 0);
    }

  protected:
    template< class LocalCoordinateType >
    static inline void dofCoordinate ( unsigned int index,
                                       LocalCoordinateType &coordinate )
    {
      assert( index <= numLagrangePoints );
      coordinate = 0;
    }
  };


  
  template< class BaseGeometryType >
  class GenericLagrangePoint< PyramidGeometry< BaseGeometryType >, 0 >
  {
  public:
    typedef PyramidGeometry< BaseGeometryType > GeometryType;
    enum { dimension = GeometryType :: dimension };
    typedef FieldVector< unsigned int, dimension > DofCoordinateType;
 
    enum { polynomialOrder = 0 };

    template< class GenericGeometryType, unsigned int porder >
    friend class GenericLagrangePoint;
 
    template< class FunctionSpaceType, class GeometryType, unsigned int porder >
    friend class GenericLagrangeBaseFunction;
  
  private:
    typedef GenericLagrangePoint< GeometryType, polynomialOrder > ThisType;
    
  public:
    enum { numLagrangePoints = 1 };

  protected:
    DofCoordinateType dofCoordinate_;
    LocalCoordinate< GeometryType, DofCoordinateType > localDofCoordinate_;

  public:
    template< unsigned int codim >
    class Codim
    {
    public:
      static inline unsigned int maxDofs ()
      {
        return ((codim == 0) ? 1 : 0);
      }
    
      static inline unsigned int maxDofsReduction ()
      {
        return ((codim == dimension) ? 1 : 0);
      }
    };
    
  public:
    inline GenericLagrangePoint ( unsigned int index )
    : localDofCoordinate_( dofCoordinate_ )
    {
      dofCoordinate( index, localDofCoordinate_ );
    }

    inline GenericLagrangePoint ( const ThisType &point )
    : dofCoordinate_( point.dofCoordinate_ ),
      localDofCoordinate_( dofCoordinate_ )
    {
    }
    
    template< class LocalCoordinateType >
    static inline void dofSubEntity ( LocalCoordinateType &coordinate,
                                      unsigned int &codim,
                                      unsigned int &subEntity )
    {
      codim = 0;
      subEntity = 0;
    }
    
    template< class LocalCoordinateType >
    static inline void dofSubEntity ( LocalCoordinateType &coordinate,
                                      unsigned int &codim,
                                      unsigned int &subEntity,
                                      unsigned int &dofNumber )
    {
      codim = 0;
      subEntity = 0;
      dofNumber = 0;
    }
    
  
    inline void dofSubEntity ( unsigned int &codim,
                               unsigned int &subEntity )
    {
      dofSubEntity( localDofCoordinate_, codim, subEntity );
    } 
    
    inline void dofSubEntity ( unsigned int &codim,
                               unsigned int &subEntity,
                               unsigned int &dofNumber )
    {
      dofSubEntity( localDofCoordinate_, codim, subEntity, dofNumber );
    } 
    
    template< class LocalCoordinateType >
    static inline void dofSubEntityReduction ( LocalCoordinateType &coordinate,
                                               unsigned int &codim,
                                               unsigned int &subEntity )
    {
      codim = dimension;
      subEntity = 0;
    }

    template< class LocalCoordinateType >
    static inline void dofSubEntityReduction ( LocalCoordinateType &coordinate,
                                               unsigned int &codim,
                                               unsigned int &subEntity,
                                               unsigned int &dofNumber )
    {
      codim = dimension;
      subEntity = 0;
      dofNumber = 0;
    }

    static inline unsigned int entityDofNumber ( unsigned int codim,
                                                 unsigned int subEntity,
                                                 unsigned int dof )
    {
      //assert( (codim == 0) && (subEntity == 0) && (dof == 0) );
      return 0;
    }

    static inline unsigned int entityDofNumberReduction ( unsigned int codim,
                                                          unsigned int subEntity,
                                                          unsigned int dof )
    {
      //assert( (codim == dimension) && (subEntity == 0) && (dof == 0) );
      return 0;
    }

    template< class LocalCoordinateType >
    static inline unsigned int height ( LocalCoordinateType &coordinate )
    {
      return polynomialOrder;
    }

    inline unsigned int height ()
    {
      return height( localDofCoordinate_ );
    }

    template< class FieldType >
    inline void local ( FieldVector< FieldType, dimension > &coordinate ) const
    {
      const FieldType factor = 1 / (FieldType)polynomialOrder;
      for( int i = 0; i < dimension; ++i )
        coordinate[ i ] = factor * dofCoordinate_[ i ]; 
    }
    
    static inline unsigned int maxDofs ( unsigned int codim )
    {
      return ((codim == 0) ? 1 : 0);
    }
    
    static inline unsigned int maxDofsReduction ( unsigned int codim )
    {
      return ((codim == dimension) ? 1 : 0);
    }

    static inline unsigned int numDofs ( unsigned int codim,
                                         unsigned int subEntity )
    {
      return ((codim == 0) ? 1 : 0);
    }

    static inline unsigned int numDofsReduction ( unsigned int codim,
                                                  unsigned int subEntity )
    {
      return ((codim == dimension) ? 1 : 0);
    }

  protected:
    template< class LocalCoordinateType >
    static inline void dofCoordinate ( unsigned int index,
                                       LocalCoordinateType &coordinate )
    {
      assert( index <= numLagrangePoints );
      coordinate = 0;
    }
  };



  template< class BaseGeometryType, unsigned int order >
  class GenericLagrangePoint< PyramidGeometry< BaseGeometryType >, order >
  {
  public:
    typedef PyramidGeometry< BaseGeometryType > GeometryType;
    enum { dimension = GeometryType :: dimension };
    typedef FieldVector< unsigned int, dimension > DofCoordinateType;
 
    enum { polynomialOrder = order };
    
    template< class GenericGeometryType, unsigned int porder >
    friend class GenericLagrangePoint;
    
    template< class FunctionSpaceType, class GeometryType, unsigned int porder >
    friend class GenericLagrangeBaseFunction;
    
    template< class GeometryType, unsigned int porder, unsigned int codim >
    friend class GenericLagrangePointCodim;

  private:
    typedef GenericLagrangePoint< GeometryType, polynomialOrder > ThisType;
 
    typedef GenericLagrangePoint< GeometryType, polynomialOrder - 1 >
      OrderReductionType;
    typedef GenericLagrangePoint< BaseGeometryType, polynomialOrder >
      DimensionReductionType;

  public:
    enum { numLagrangePoints = DimensionReductionType :: numLagrangePoints
                             + OrderReductionType :: numLagrangePoints };
   
  protected:
    DofCoordinateType dofCoordinate_;
    LocalCoordinate< GeometryType, DofCoordinateType > localDofCoordinate_;

  public:
    template< unsigned int codim >
    class Codim
    : public GenericLagrangePointCodim< GeometryType, polynomialOrder, codim >
    {
    };

  public:
    inline GenericLagrangePoint ( unsigned int index )
    : localDofCoordinate_( dofCoordinate_ )
    {
      dofCoordinate( index, localDofCoordinate_ );
    }
    
    inline GenericLagrangePoint ( const ThisType &point )
    : dofCoordinate_( point.dofCoordinate_ ),
      localDofCoordinate_( dofCoordinate_ )
    {
    }

    template< class LocalCoordinateType >
    static inline void dofSubEntity ( LocalCoordinateType &coordinate,
                                      unsigned int &codim,
                                      unsigned int &subEntity )
    {
      if( !useDimReduction( coordinate ) ) {
        --(*coordinate);
        OrderReductionType :: dofSubEntityReduction
          ( coordinate, codim, subEntity );
        ++(*coordinate);
        if( codim > 0 )
          subEntity += BaseGeometryType :: numSubEntities( codim - 1 );
      } else {
        DimensionReductionType :: dofSubEntity
          ( coordinate.base(), codim, subEntity );
        ++codim;
      }
    }
    
    template< class LocalCoordinateType >
    static inline void dofSubEntity ( LocalCoordinateType &coordinate,
                                      unsigned int &codim,
                                      unsigned int &subEntity,
                                      unsigned int &dofNumber )
    {
      if( !useDimReduction( coordinate ) ) {
        --(*coordinate);
        OrderReductionType :: template dofSubEntityReduction
          ( coordinate, codim, subEntity, dofNumber );
        ++(*coordinate);
        if( codim > 0 )
          subEntity += BaseGeometryType :: numSubEntities( codim - 1 );
      } else {
        DimensionReductionType :: dofSubEntity
          ( coordinate.base(), codim, subEntity, dofNumber );
        ++codim;
      }
    }

 
    inline void dofSubEntity ( unsigned int &codim,
                               unsigned int &subEntity )
    {
      dofSubEntity( localDofCoordinate_, codim, subEntity );
    }

    inline void dofSubEntity ( unsigned int &codim,
                               unsigned int &subEntity,
                               unsigned int &dofNumber )
    {
      dofSubEntity( localDofCoordinate_, codim, subEntity, dofNumber );
    }

    template< class LocalCoordinateType >
    static inline void dofSubEntityReduction ( LocalCoordinateType &coordinate,
                                               unsigned int &codim,
                                               unsigned int &subEntity )
    {
      if( !useDimReduction( coordinate ) ) {
        --(*coordinate);
        OrderReductionType :: dofSubEntityReduction
          ( coordinate, codim, subEntity );
        ++(*coordinate);
      } else
        DimensionReductionType :: dofSubEntity
          ( coordinate.base(), codim, subEntity );
    }
 
    template< class LocalCoordinateType >
    static inline void dofSubEntityReduction ( LocalCoordinateType &coordinate,
                                               unsigned int &codim,
                                               unsigned int &subEntity,
                                               unsigned int &dofNumber )
    {
      if( !useDimReduction( coordinate ) ) {
        --(*coordinate);
        OrderReductionType :: dofSubEntityReduction
          ( coordinate, codim, subEntity, dofNumber );
        ++(*coordinate);
        dofNumber += DimensionReductionType :: numDofs( codim, subEntity );
      } else
        DimensionReductionType :: dofSubEntity
          ( coordinate.base(), codim, subEntity, dofNumber );
    }
 
    static inline unsigned int entityDofNumber ( unsigned int codim,
                                                 unsigned int subEntity,
                                                 unsigned int dof )
    {
      if( codim == 0 )
        return DimensionReductionType :: numLagrangePoints
               + OrderReductionType :: entityDofNumberReduction
                   ( codim, subEntity, dof );

      const unsigned int baseSubEntities = BaseGeometryType :: numSubEntities( codim - 1 );
      if( subEntity >= baseSubEntities )
        return DimensionReductionType :: numLagrangePoints
               + OrderReductionType :: entityDofNumberReduction
                   ( codim, subEntity - baseSubEntities, dof );
      else
        return DimensionReductionType :: entityDofNumber( codim - 1, subEntity, dof );
    }

    static inline unsigned int entityDofNumberReduction ( unsigned int codim,
                                                          unsigned int subEntity,
                                                          unsigned int dof )
    {
      const unsigned int baseEntityDofs
        = DimensionReductionType :: numDofs( codim, subEntity );
      if( dof >=  baseEntityDofs )
        return DimensionReductionType :: numLagrangePoints
               + OrderReductionType :: entityDofNumberReduction
                   ( codim, subEntity, dof - baseEntityDofs );
      else
        return DimensionReductionType :: entityDofNumber( codim, subEntity, dof );
    }
    
    template< class LocalCoordinateType >
    static inline unsigned int height ( LocalCoordinateType &coordinate )
    {
      if( !useDimReduction( coordinate ) ) {
        --(*coordinate);
        unsigned int h = OrderReductionType :: height( coordinate );
        ++(*coordinate);
        return h;
      } else
        return DimensionReductionType :: height( coordinate.base() );
    }

    inline unsigned int height ()
    {
      return height( localDofCoordinate_ );
    }

    template< class FieldType >
    inline void local ( FieldVector< FieldType, dimension > &coordinate ) const
    {
      const FieldType factor = 1 / (FieldType)polynomialOrder;
      for( int i = 0; i < dimension; ++i )
        coordinate[ i ] = factor * dofCoordinate_[ i ]; 
    }
    
    static inline unsigned int maxDofs ( unsigned int codim )
    {
      const unsigned int maxOrderDofs
        = OrderReductionType :: maxDofsReduction( codim );

      if( codim == 0 )
        return maxOrderDofs;

      const unsigned int maxDimDofs
        = DimensionReductionType :: maxDofs( codim - 1 );

      return (maxDimDofs > maxOrderDofs) ? maxDimDofs : maxOrderDofs;
    }

    static inline unsigned int maxDofsReduction ( unsigned int codim )
    {
      return DimensionReductionType :: maxDofs( codim )
             + OrderReductionType :: maxDofsReduction( codim );
    }
   
    static inline unsigned int numDofs ( unsigned int codim,
                                         unsigned int subEntity )
    {
      if( codim == 0 )
        return OrderReductionType :: numDofsReduction( codim, subEntity );
      
      const unsigned int baseSubEntities
        = BaseGeometryType :: numSubEntities( codim - 1 );
      if( subEntity < baseSubEntities )
        return DimensionReductionType :: numDofs( codim - 1, subEntity );
      else
        return OrderReductionType
               :: numDofsReduction( codim, subEntity - baseSubEntities );
    }

    static inline unsigned int numDofsReduction ( unsigned int codim,
                                                  unsigned int subEntity )
    {
      return DimensionReductionType :: numDofs( codim, subEntity )
             + OrderReductionType :: numDofsReduction( codim, subEntity );
    }

    template< class LocalCoordinateType >
    static inline bool useDimReduction ( const LocalCoordinateType &coordinate )
    {
      return (*coordinate == 0);
    }

  protected:
    template< class LocalCoordinateType >
    static inline void dofCoordinate ( unsigned int index,
                                       LocalCoordinateType &coordinate )
    {
      assert( index <= numLagrangePoints );

      if( index < DimensionReductionType :: numLagrangePoints ) {
        (*coordinate) = 0;
        DimensionReductionType :: dofCoordinate( index, coordinate.base() );
      }
      else {
        const int orderIndex
          = index - DimensionReductionType :: numLagrangePoints;
        OrderReductionType :: dofCoordinate( orderIndex, coordinate );
        ++(*coordinate);
      }
    }
  };


 
  template< class BaseGeometryType, unsigned int order, unsigned int codim >
  class GenericLagrangePointCodim
    < PyramidGeometry< BaseGeometryType >, order, codim >
  {
  private:
    typedef PyramidGeometry< BaseGeometryType > GeometryType;

    typedef GenericLagrangePoint< GeometryType, order >
      LagrangePointType;
    
    typedef typename LagrangePointType :: DimensionReductionType
      DimensionReductionType;
    typedef typename LagrangePointType :: OrderReductionType
      OrderReductionType;
 
  public:
    static inline unsigned int maxDofs ()
    {
      const unsigned int maxOrderDofs
        = OrderReductionType :: template Codim< codim > :: maxDofsReduction();

      const unsigned int maxDimDofs
        = DimensionReductionType :: template Codim< codim - 1 > :: maxDofs();
      return (maxDimDofs > maxOrderDofs) ? maxDimDofs : maxOrderDofs;
    }

    static inline unsigned int maxDofsReduction ()
    {
      return DimensionReductionType :: template Codim< codim > :: maxDofs()
             + OrderReductionType :: template Codim< codim > :: maxDofsReduction();
    }
  };


  
  template< class BaseGeometryType, unsigned int order >
  class GenericLagrangePointCodim
    < PyramidGeometry< BaseGeometryType >, order, 0 >
  {
  private:
    typedef PyramidGeometry< BaseGeometryType > GeometryType;

    typedef GenericLagrangePoint< GeometryType, order >
      LagrangePointType;
    
    typedef typename LagrangePointType :: DimensionReductionType
      DimensionReductionType;
    typedef typename LagrangePointType :: OrderReductionType
      OrderReductionType;

    enum { codim = 0 };

  public:
    static inline unsigned int maxDofs ()
    {
      const unsigned int maxOrderDofs
        = OrderReductionType :: template Codim< codim > :: maxDofsReduction();

      return maxOrderDofs;
    }

    static inline unsigned int maxDofsReduction ()
    {
      return DimensionReductionType :: template Codim< codim > :: maxDofs()
             + OrderReductionType :: template Codim< codim > :: maxDofsReduction();
    }
  };


  
  template< class FirstGeometryType,
            class SecondGeometryType,
            unsigned int order >
  class GenericLagrangePoint< ProductGeometry< FirstGeometryType,
                                               SecondGeometryType >,
                              order >
  {
  public:
    typedef ProductGeometry< FirstGeometryType, SecondGeometryType >
      GeometryType;
    enum { dimension = GeometryType :: dimension };
    typedef FieldVector< unsigned int, dimension > DofCoordinateType;

    enum { polynomialOrder = order };

    template< class GenericGeometryType, unsigned int porder >
    friend class GenericLagrangePoint;
 
    template< class FunctionSpaceType, class GeometryType, unsigned int porder >
    friend class GenericLagrangeBaseFunction;

  private:
    typedef GenericLagrangePoint< GeometryType, polynomialOrder > ThisType;

    typedef GenericLagrangePoint< FirstGeometryType, polynomialOrder >
      FirstReductionType;
    typedef GenericLagrangePoint< SecondGeometryType, polynomialOrder >
      SecondReductionType;
 
  public:
   enum { numLagrangePoints = FirstReductionType :: numLagrangePoints
                             * SecondReductionType :: numLagrangePoints };
   
  protected:
    DofCoordinateType dofCoordinate_;
    LocalCoordinate< GeometryType, DofCoordinateType > localDofCoordinate_;

  private:
    template< unsigned int codim, unsigned int i >
    class CodimIterator
    {
    public:
      static inline unsigned int maxDofs ()
      {
        const unsigned int n
          = FirstReductionType :: template Codim< codim - i > :: maxDofs()
          * SecondReductionType :: template Codim< i > :: maxDofs();
      
        const unsigned int m = CodimIterator< codim, i-1 > :: maxDofs();
        return ((m > n) ? m : n);
      }
    };
    
    template< unsigned int codim >
    class CodimIterator< codim, 0 >
    {
    private:
      enum { i = 0 };

    public:
      static inline unsigned int maxDofs ()
      {
        const unsigned int n
          = FirstReductionType :: template Codim< codim - i > :: maxDofs()
          * SecondReductionType :: template Codim< i > :: maxDofs();
        
        return n;
      }
    };

  public:
    template< unsigned int codim >
    class Codim
    {
    public:
      static inline unsigned int maxDofs ()
      {
        return CodimIterator< codim, codim > :: maxDofs();
      }
    };

  public:
    inline GenericLagrangePoint ( unsigned int index )
    : localDofCoordinate_( dofCoordinate_ )
    {
      dofCoordinate( index, localDofCoordinate_ );
    }

    inline GenericLagrangePoint ( const ThisType &point )
    : dofCoordinate_( point.dofCoordinate_ ),
      localDofCoordinate_( dofCoordinate_ )
    {
    }

    template< class LocalCoordinateType >
    static inline void dofSubEntity( LocalCoordinateType &coordinate,
                                     unsigned int &codim,
                                     unsigned int &subEntity )
    {
      unsigned int firstCodim, secondCodim;
      unsigned int firstSubEntity, secondSubEntity;

      FirstReductionType :: dofSubEntity( coordinate.first(),
                                          firstCodim,
                                          firstSubEntity );
      SecondReductionType :: dofSubEntity( coordinate.second(),
                                           secondCodim,
                                           secondSubEntity );

      codim = firstCodim + secondCodim;

      subEntity = 0;
      for( unsigned int i = 0; i < secondCodim; ++i )
        subEntity += FirstGeometryType :: numSubEntities( codim - i  )
                   * SecondGeometryType :: numSubEntities( i );
      subEntity += firstSubEntity + secondSubEntity
                 * FirstGeometryType :: numSubEntities( firstCodim );
    }

    template< class LocalCoordinateType >
    static inline void dofSubEntity( LocalCoordinateType &coordinate,
                                     unsigned int &codim,
                                     unsigned int &subEntity,
                                     unsigned int &dofNumber )
    {
      unsigned int firstCodim, secondCodim;
      unsigned int firstSubEntity, secondSubEntity;
      unsigned int firstDofNumber, secondDofNumber;

      FirstReductionType :: dofSubEntity( coordinate.first(),
                                          firstCodim,
                                          firstSubEntity,
                                          firstDofNumber );
      SecondReductionType :: dofSubEntity( coordinate.second(),
                                           secondCodim,
                                           secondSubEntity,
                                           secondDofNumber );

      codim = firstCodim + secondCodim;

      subEntity = 0;
      for( unsigned int i = 0; i < secondCodim; ++i )
        subEntity += FirstGeometryType :: numSubEntities( codim - i  )
                   * SecondGeometryType :: numSubEntities( i );
      subEntity += firstSubEntity + secondSubEntity
                 * FirstGeometryType :: numSubEntities( firstCodim );

      dofNumber = firstDofNumber + secondDofNumber
                * FirstReductionType :: numDofs( firstCodim, firstSubEntity );
    }
   
    inline void dofSubEntity ( unsigned int &codim, unsigned int &subEntity )
    {
      dofSubEntity( localDofCoordinate_, codim, subEntity );
    } 
    
    inline void dofSubEntity ( unsigned int &codim,
                               unsigned int &subEntity,
                               unsigned int &dofNumber )
    {
      dofSubEntity( localDofCoordinate_, codim, subEntity, dofNumber );
    } 

    static inline unsigned int entityDofNumber ( unsigned int codim,
                                                 unsigned int subEntity,
                                                 unsigned int dofNumber )
    {
      unsigned int firstCodim = codim;
      unsigned int secondCodim = 0;
      for( ; secondCodim < codim; --firstCodim, ++secondCodim ) {
        const unsigned int num
          = FirstGeometryType :: numSubEntities( firstCodim )
          * SecondGeometryType :: numSubEntities( secondCodim );

        if( subEntity < num )
          break;
        subEntity -= num;
      }
      
      const unsigned int n = FirstGeometryType :: numSubEntities( firstCodim );
      const unsigned int firstSubEntity = subEntity % n;
      const unsigned int secondSubEntity = subEntity / n;

      const unsigned int m
        = FirstReductionType :: numDofs( firstCodim, firstSubEntity );
      const unsigned int firstDofNumber = dofNumber % m;
      const unsigned int secondDofNumber = dofNumber / m;
   
      const unsigned int firstEntityDofNumber
        = FirstReductionType :: entityDofNumber
            ( firstCodim, firstSubEntity, firstDofNumber );
      const unsigned int secondEntityDofNumber
        = SecondReductionType :: entityDofNumber
            ( secondCodim, secondSubEntity, secondDofNumber );
      
      return firstEntityDofNumber
             + secondEntityDofNumber * FirstReductionType :: numLagrangePoints;
    }

    template< class LocalCoordinateType >
    static inline unsigned int height ( LocalCoordinateType &coordinate )
    {
      const unsigned int firstHeight
        = FirstReductionType :: height( coordinate.first() );
      const unsigned int secondHeight
        = SecondReductionType :: height( coordinate.second() );
        
      return ((firstHeight < secondHeight) ? firstHeight : secondHeight);
    }

    inline unsigned int height ()
    {
      return height( localDofCoordinate_ );
    }

    template< class FieldType >
    inline void local ( FieldVector< FieldType, dimension > &coordinate ) const
    {
      const FieldType factor = 1 / (FieldType)polynomialOrder;
      for( int i = 0; i < dimension; ++i )
        coordinate[ i ] = factor * dofCoordinate_[ i ]; 
    }
    
    static inline unsigned int maxDofs ( unsigned int codim )
    {
      unsigned int max = 0;
      for( unsigned int i = 0; i <= codim; ++i ) {
        const unsigned int n
          = FirstReductionType :: maxDofs( codim - i )
          * SecondReductionType :: maxDofs( i );
        max = (max >= n) ? max : n;
      }
      return max;
    }

    static inline unsigned int numDofs ( unsigned int codim,
                                         unsigned int subEntity )
    {
      unsigned int firstCodim = codim;
      unsigned int secondCodim = 0;
      for( ; secondCodim < codim; --firstCodim, ++secondCodim ) {
        const unsigned int num
          = FirstGeometryType :: numSubEntities( firstCodim )
          * SecondGeometryType :: numSubEntities( secondCodim );

        if( subEntity < num )
          break;
        subEntity -= num;
      }
      
      const unsigned int n = FirstGeometryType :: numSubEntities( firstCodim );
      const unsigned int firstSubEntity = subEntity % n;
      const unsigned int secondSubEntity = subEntity / n;
     
      return FirstReductionType :: numDofs( firstCodim, firstSubEntity )
             * SecondReductionType :: numDofs( secondCodim, secondSubEntity );
    }

  protected:
    template< class LocalCoordinateType >
    static inline void dofCoordinate ( unsigned int index,
                                       LocalCoordinateType &coordinate )
    {
      assert( index <= numLagrangePoints );

      const unsigned int firstIndex
        = index % FirstReductionType :: numLagrangePoints;
      const unsigned int secondIndex
        = index / FirstReductionType :: numLagrangePoints;

      FirstReductionType :: dofCoordinate( firstIndex, coordinate.first() );
      SecondReductionType :: dofCoordinate( secondIndex, coordinate.second() );
    }
  };
}

#endif
