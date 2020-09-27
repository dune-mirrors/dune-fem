#ifndef DUNE_FEM_SPACE_LAGRANGE_GENERICBASEFUNCTIONS_HH
#define DUNE_FEM_SPACE_LAGRANGE_GENERICBASEFUNCTIONS_HH

#include <dune/common/fvector.hh>

#include "genericlagrangepoints.hh"

namespace Dune
{

  namespace Fem
  {

    template< class FunctionSpace, class GeometryType, unsigned int order >
    class GenericLagrangeBaseFunction;


    template< class FunctionSpace, unsigned int order >
    class GenericLagrangeBaseFunction< FunctionSpace, PointGeometry, order >
    {
      typedef GenericLagrangeBaseFunction< FunctionSpace, PointGeometry, order > ThisType;

      static_assert( (FunctionSpace::dimRange == 1), "FunctionSpace must be scalar." );

    public:
      typedef FunctionSpace FunctionSpaceType;

      typedef PointGeometry GeometryType;

      static constexpr unsigned int polynomialOrder = order;

      typedef GenericLagrangePoint< GeometryType, polynomialOrder >
        LagrangePointType;
      static const unsigned int numBaseFunctions = LagrangePointType :: numLagrangePoints;

      typedef typename FunctionSpaceType :: DomainType DomainType;
      typedef typename FunctionSpaceType :: RangeType RangeType;

      typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;

    private:
      const LagrangePointType lagrangePoint_;

    public:
      explicit GenericLagrangeBaseFunction ( unsigned int baseNum )
      : lagrangePoint_( baseNum )
      {}

      template< class LocalDofCoordinateType, class LocalCoordinateType >
      inline static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                                    const FieldVector< int, 0 > &diffVariable,
                                    DomainFieldType factor,
                                    const LocalCoordinateType &x,
                                    RangeType &phi )
      {
        phi[ 0 ] = 1;
      }

      template< class LocalDofCoordinateType, class LocalCoordinateType >
      inline static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                                    const FieldVector< int, 1 > &diffVariable,
                                    DomainFieldType factor,
                                    const LocalCoordinateType &x,
                                    RangeType &phi )
      {
        phi[ 0 ] = 0;
      }

      template< class LocalDofCoordinateType, class LocalCoordinateType >
      inline static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                                    const FieldVector< int, 2 > &diffVariable,
                                    DomainFieldType factor,
                                    const LocalCoordinateType &x,
                                    RangeType &phi )
      {
        phi[ 0 ] = 0;
      }

      template< int diffOrder >
      inline void
      evaluate ( const FieldVector< int, diffOrder > &diffVariable,
                 const DomainType &x,
                 RangeType &phi ) const
      {
        const LocalCoordinate< GeometryType, DomainFieldType > xlocal( x );
        LagrangePointType point( lagrangePoint_ );
        evaluate( point.dofCoordinate_, diffVariable, 1, xlocal, phi );
      }
    };


    template< class FunctionSpace, class BaseGeometryType >
    class GenericLagrangeBaseFunction< FunctionSpace, PyramidGeometry< BaseGeometryType >, 0 >
    {
      typedef GenericLagrangeBaseFunction< FunctionSpace, PyramidGeometry< BaseGeometryType >, 0 > ThisType;

      static_assert( (FunctionSpace::dimRange == 1), "FunctionSpace must be scalar." );

    public:
      typedef FunctionSpace FunctionSpaceType;

      typedef PyramidGeometry< BaseGeometryType > GeometryType;

      enum { polynomialOrder = 0 };

      typedef GenericLagrangePoint< GeometryType, polynomialOrder > LagrangePointType;
      static const unsigned int numBaseFunctions = LagrangePointType :: numLagrangePoints;


      typedef typename FunctionSpaceType :: DomainType DomainType;
      typedef typename FunctionSpaceType :: RangeType RangeType;

      typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;

    private:
      const LagrangePointType lagrangePoint_;

    public:
      explicit GenericLagrangeBaseFunction ( unsigned int baseNum )
      : lagrangePoint_( baseNum )
      {}

      template< unsigned int porder, class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 0 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        phi[ 0 ] = 1;
      }

      template< class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 0 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        return evaluate< polynomialOrder >( dofCoordinate, diffVariable, factor, x, phi );
      }

      template< unsigned int porder, class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 1 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        phi[ 0 ] = 0;
      }

      template< class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 1 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        return evaluate< polynomialOrder >( dofCoordinate, diffVariable, factor, x, phi );
      }

      template< unsigned int porder, class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 2 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        phi[ 0 ] = 0;
      }

      template< class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 2 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        return evaluate< polynomialOrder >( dofCoordinate, diffVariable, factor, x, phi );
      }

      template< int diffOrder >
      void evaluate ( const FieldVector< int, diffOrder > &diffVariable,
                      const DomainType &x,
                      RangeType &phi ) const
      {
        const LocalCoordinate< GeometryType, DomainFieldType > xlocal( x );
        LagrangePointType point( lagrangePoint_ );
        evaluate( point.dofCoordinate_, diffVariable, 1, xlocal, phi );
      }
    };


    template< class FunctionSpace, class BaseGeometryType, unsigned int order >
    class GenericLagrangeBaseFunction< FunctionSpace, PyramidGeometry< BaseGeometryType >, order >
    {
      typedef GenericLagrangeBaseFunction< FunctionSpace, PyramidGeometry< BaseGeometryType >, order > ThisType;

      static_assert( (FunctionSpace::dimRange == 1), "FunctionSpace must be scalar." );

    public:
      typedef FunctionSpace FunctionSpaceType;

      typedef PyramidGeometry< BaseGeometryType > GeometryType;

      static constexpr unsigned int polynomialOrder = order;

      typedef GenericLagrangePoint< GeometryType, polynomialOrder >
        LagrangePointType;
      static const unsigned int numBaseFunctions = LagrangePointType :: numLagrangePoints;


      typedef typename FunctionSpaceType :: DomainType DomainType;
      typedef typename FunctionSpaceType :: RangeType RangeType;

      typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;

    private:
      typedef GenericLagrangeBaseFunction
        < FunctionSpaceType, BaseGeometryType, polynomialOrder >
        DimensionReductionType;
      typedef GenericLagrangeBaseFunction
        < FunctionSpaceType, GeometryType, polynomialOrder - 1 >
        OrderReductionType;

    private:
      const LagrangePointType lagrangePoint_;

    public:
      explicit GenericLagrangeBaseFunction ( unsigned int baseNum )
      : lagrangePoint_( baseNum )
      {}

      template< unsigned int porder, class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 0 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        const DomainFieldType divisor = DomainFieldType( 1 ) / ((DomainFieldType)polynomialOrder);
        const DomainFieldType myfactor = porder * divisor;
        const RealType myshift = (porder - polynomialOrder) * divisor;

        if( LagrangePointType :: useDimReduction( dofCoordinate ) )
        {
          DimensionReductionType::evaluate( dofCoordinate.base(), diffVariable, myfactor * factor, x.base(), phi );

          const unsigned int height
            = LagrangePointType :: height( dofCoordinate );
          for( unsigned int i = 0; i < height; ++i )
          {
            ++(*dofCoordinate);
            RangeType psi;
            evaluate< porder >( dofCoordinate, diffVariable, factor, x, psi );
            phi -= psi;
          }
          (*dofCoordinate) -= height;
        }
        else
        {
          --(*dofCoordinate);
          OrderReductionType::template evaluate< porder >( dofCoordinate, diffVariable, factor, x, phi );
          ++(*dofCoordinate);
          phi[ 0 ] *= (polynomialOrder / ((RealType)(*dofCoordinate)))
                    * (myfactor * factor * (*x) - myshift);
        }
      }

      template< class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 0 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        return evaluate< polynomialOrder >( dofCoordinate, diffVariable, factor, x, phi );
      }

      template< unsigned int porder, class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 1 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        const DomainFieldType divisor = polynomialOrder;
        const DomainFieldType myfactor = porder / divisor;
        const RealType myshift = (porder - polynomialOrder) / divisor;

        FieldVector< int, 0 > dv;

        if( LagrangePointType :: useDimReduction( dofCoordinate ) )
        {
          if( (unsigned int)diffVariable[ 0 ] != LocalDofCoordinateType::index )
            DimensionReductionType::evaluate( dofCoordinate.base(), diffVariable, myfactor * factor, x.base(), phi );
          else
            phi = 0;

          const unsigned int height
            = LagrangePointType :: height( dofCoordinate );
          for( unsigned int i = 0; i < height; ++i )
          {
            ++(*dofCoordinate);
            RangeType psi;
            evaluate< porder >( dofCoordinate, diffVariable, factor, x, psi );
            phi -= psi;
          }
          (*dofCoordinate) -= height;
        }
        else
        {
          --(*dofCoordinate);
          OrderReductionType::template evaluate< porder >( dofCoordinate, diffVariable, factor, x, phi );
          phi *= (myfactor * factor * (*x) - myshift);

          if( (unsigned int)diffVariable[ 0 ] == LocalDofCoordinateType::index )
          {
            RangeType psi;
            OrderReductionType::template evaluate< porder >( dofCoordinate, dv, factor, x, psi );
            phi.axpy( myfactor * factor, psi );
          }
          ++(*dofCoordinate);
          phi *= (polynomialOrder) / ((RealType)(*dofCoordinate));
        }
      }

      template< class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 1 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        return evaluate< polynomialOrder >( dofCoordinate, diffVariable, factor, x, phi );
      }

      template< unsigned int porder, class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 2 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        const DomainFieldType divisor = polynomialOrder;
        const DomainFieldType myfactor = porder / divisor;
        const RealType myshift = (porder - polynomialOrder) / divisor;

        FieldVector< int, 1 > dv0( diffVariable[ 0 ] );
        FieldVector< int, 1 > dv1( diffVariable[ 1 ] );

        if( LagrangePointType :: useDimReduction( dofCoordinate ) )
        {
          if( ((unsigned int)diffVariable[ 0 ] != LocalDofCoordinateType::index)
              && ((unsigned int)diffVariable[ 1 ] != LocalDofCoordinateType::index) )
          DimensionReductionType::evaluate( dofCoordinate.base(), diffVariable, myfactor * factor, x.base(), phi );
          else
            phi = 0;

          const unsigned int height = LagrangePointType::height( dofCoordinate );
          for( unsigned int i = 0; i < height; ++i )
          {
            ++(*dofCoordinate);
            RangeType psi;
            evaluate< porder >( dofCoordinate, diffVariable, factor, x, psi );
            phi -= psi;
          }
          (*dofCoordinate) -= height;
        }
        else
        {
          RangeType psi;
          --(*dofCoordinate);
          OrderReductionType::template evaluate< porder >( dofCoordinate, diffVariable, factor, x, phi );
          phi *= (myfactor * factor * (*x) - myshift);

          if( (unsigned int)diffVariable[ 0 ] == LocalDofCoordinateType::index )
          {
            OrderReductionType::template evaluate< porder >( dofCoordinate, dv1, factor, x, psi );
            phi.axpy( myfactor * factor, psi );
          }

          if( (unsigned int)diffVariable[ 1 ] == LocalDofCoordinateType::index )
          {
            OrderReductionType::template evaluate< porder >( dofCoordinate, dv0, factor, x, psi );
            phi.axpy( myfactor * factor, psi );
          }

          ++(*dofCoordinate);
          phi *= (polynomialOrder) / ((RealType)(*dofCoordinate));
        }
      }

      template< class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 2 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        return evaluate< polynomialOrder >( dofCoordinate, diffVariable, factor, x, phi );
      }

      template< int diffOrder >
      void evaluate ( const FieldVector< int, diffOrder > &diffVariable,
                      const DomainType &x,
                      RangeType &phi ) const
      {
        const LocalCoordinate< GeometryType, DomainFieldType > xlocal( x );
        LagrangePointType point( lagrangePoint_ );
        evaluate( point.dofCoordinate_, diffVariable, 1, xlocal, phi );
      }
    };


    template< class FunctionSpace, class FirstGeometryType, class SecondGeometryType, unsigned int order >
    class GenericLagrangeBaseFunction< FunctionSpace, ProductGeometry< FirstGeometryType, SecondGeometryType >, order >
    {
      typedef GenericLagrangeBaseFunction< FunctionSpace, ProductGeometry< FirstGeometryType, SecondGeometryType >, order > ThisType;

      static_assert( (FunctionSpace::dimRange == 1), "FunctionSpace must be scalar." );

    public:
      typedef FunctionSpace FunctionSpaceType;

      typedef ProductGeometry< FirstGeometryType, SecondGeometryType > GeometryType;

      static constexpr unsigned int polynomialOrder = order;

      typedef GenericLagrangePoint< GeometryType, polynomialOrder > LagrangePointType;
      static const unsigned int numBaseFunctions = LagrangePointType :: numLagrangePoints;

      typedef typename FunctionSpaceType :: DomainType DomainType;
      typedef typename FunctionSpaceType :: RangeType RangeType;

      typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;

    private:
      typedef GenericLagrangeBaseFunction
        < FunctionSpaceType, FirstGeometryType, polynomialOrder >
        FirstReductionType;
      typedef GenericLagrangeBaseFunction
        < FunctionSpaceType, SecondGeometryType, polynomialOrder >
        SecondReductionType;

    private:
      const LagrangePointType lagrangePoint_;

    public:
      explicit GenericLagrangeBaseFunction ( unsigned int baseNum )
      : lagrangePoint_( baseNum )
      {}

      template< class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 0 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        RangeType psi;
        FirstReductionType::evaluate( dofCoordinate.first(), diffVariable, factor, x.first(), phi );
        SecondReductionType::evaluate( dofCoordinate.second(), diffVariable, factor, x.second(), psi );
        phi[ 0 ] *= psi[ 0 ];
      }

      template< class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 1 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        FieldVector< int, 0 > dv;
        RangeType psi1, psi2;

        FirstReductionType::evaluate( dofCoordinate.first(), diffVariable, factor, x.first(), psi1 );
        SecondReductionType::evaluate( dofCoordinate.second(), dv, factor, x.second(), psi2 );
        phi[ 0 ] = psi1[ 0 ] * psi2[ 0 ];

        FirstReductionType::evaluate( dofCoordinate.first(), dv, factor, x.first(), psi1 );
        SecondReductionType::evaluate( dofCoordinate.second(), diffVariable, factor, x.second(), psi2 );
        phi[ 0 ] += psi1[ 0 ] * psi2[ 0 ];
      }

      template< class LocalDofCoordinateType, class LocalCoordinateType >
      static void evaluate ( LocalDofCoordinateType &dofCoordinate,
                             const FieldVector< int, 2 > &diffVariable,
                             DomainFieldType factor,
                             const LocalCoordinateType &x,
                             RangeType &phi )
      {
        FieldVector< int, 0 > dv;
        FieldVector< int, 1 > dv0( diffVariable[ 0 ] );
        FieldVector< int, 1 > dv1( diffVariable[ 1 ] );
        RangeType psi1, psi2;

        FirstReductionType::evaluate( dofCoordinate.first(), diffVariable, factor, x.first(), psi1 );
        SecondReductionType::evaluate( dofCoordinate.second(), dv, factor, x.second(), psi2 );
        phi[ 0 ] = psi1[ 0 ] * psi2[ 0 ];

        FirstReductionType::evaluate( dofCoordinate.first(), dv0, factor, x.first(), psi1 );
        SecondReductionType::evaluate( dofCoordinate.second(), dv1, factor, x.second(), psi2 );
        phi[ 0 ] += psi1[ 0 ] * psi2[ 0 ];

        FirstReductionType::evaluate( dofCoordinate.first(), dv1, factor, x.first(), psi1 );
        SecondReductionType::evaluate( dofCoordinate.second(), dv0, factor, x.second(), psi2 );
        phi[ 0 ] += psi1[ 0 ] * psi2[ 0 ];

        FirstReductionType::evaluate( dofCoordinate.first(), dv, factor, x.first(), psi1 );
        SecondReductionType::evaluate( dofCoordinate.second(), diffVariable, factor, x.second(), psi2 );
        phi[ 0 ] += psi1[ 0 ] * psi2[ 0 ];
      }

      template< int diffOrder >
      void evaluate ( const FieldVector< int, diffOrder > &diffVariable,
                      const DomainType &x,
                      RangeType &phi ) const
      {
        const LocalCoordinate< GeometryType, DomainFieldType > xlocal( x );
        LagrangePointType point( lagrangePoint_ );
        evaluate( point.dofCoordinate_, diffVariable, 1, xlocal, phi );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_GENERICBASEFUNCTIONS_HH
