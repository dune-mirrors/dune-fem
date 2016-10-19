#ifndef DUNE_FEM_SCHEMES_INTEGRANDS_HH
#define DUNE_FEM_SCHEMES_INTEGRANDS_HH

#include <algorithm>
#include <tuple>
#include <utility>

#include <dune/common/ftraits.hh>

namespace Dune
{

  namespace Fem
  {

    // DiffusionModelIntegrands
    // ------------------------

    template< class Model >
    struct DiffusionModelIntegrands
    {
      typedef typename Model::GridPartType GridPartType;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
      typedef typename GridPartType::IntersectionType IntersectionType;

      typedef typename Model::RangeType RangeType;
      typedef typename Model::JacobianRangeType JacobianRangeType;

      typedef std::tuple< RangeType, JacobianRangeType > ValueType;

      explicit DiffusionModelIntegrands ( const Model &model ) : model_( &model ) {}

      bool init ( const EntityType &entity ) { return model().init( entity ); }

      bool init ( const IntersectionType &intersection )
      {
        return (intersection.boundary() && model().hasNeumanBoundary() && model().init( intersection.inside() ));
      }

      template< class Point >
      ValueType interior ( const Point &x, const ValueType &u ) const
      {
        RangeType source( 0 );
        model().source( x, std::get< 0 >( u ), std::get< 1 >( u ), source );
        JacobianRangeType dFlux( 0 );
        model().diffusiveFlux( x, std::get< 0 >( u ), std::get< 1 >( u ), dFlux );
        return std::make_tuple( source, dFlux );
      }

      template< class Point >
      auto linearizedInterior ( const Point &x, const ValueType &u ) const
      {
        return [ this, x, u ] ( const ValueType &phi ) {
            RangeType source( 0 );
            model().linSource( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), source );
            JacobianRangeType dFlux( 0 );
            model().linDiffusiveFlux( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), dFlux );
            return std::make_tuple( source, dFlux );
          };
      }

      template< class Point >
      ValueType boundary ( const Point &x, const ValueType &u ) const
      {
        RangeType alpha( 0 );
        model().alpha( x, std::get< 0 >( u ), alpha );
        return std::make_tuple( alpha, 0 );
      }

      template< class Point >
      auto linearizedBoundary ( const Point &x, const ValueType &u ) const
      {
        return [ this, x, u ] ( const ValueType &phi ) {
            RangeType alpha( 0 );
            model().linAlpha( std::get< 0 >( u ), x, std::get< 0 >( phi ), alpha );
            return std::make_tuple( alpha, 0 );
          };
      }

      const Model &model () const { return *model_; }

    private:
      const Model *model_;
    };



    // DGDiffusionModelIntegrands
    // --------------------------

    template< class Model >
    struct DGDiffusionModelIntegrands
    {
      typedef typename Model::GridPartType GridPartType;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
      typedef typename GridPartType::IntersectionType IntersectionType;

      typedef typename Model::RangeType RangeType;
      typedef typename Model::JacobianRangeType JacobianRangeType;

      typedef typename FieldTraits< RangeType >::field_type RangeFieldType;

      typedef std::tuple< RangeType, JacobianRangeType > ValueType;

      DGDiffusionModelIntegrands ( const Model &model, RangeFieldType penalty )
        : model_( &model ), penalty_( penalty )
      {}

      bool init ( const EntityType &entity )
      {
        intersection_ = nullptr;
        return model().init( entity );
      }

      bool init ( const IntersectionType &intersection )
      {
        intersection_ = &intersection;
        if( intersection.boundary() )
        {
          const EntityType inside = intersection.inside();
          beta_ = penalty_ * intersection.geometry().volume() / inside.geometry().volume();
          return (model().hasNeumanBoundary() && model().init( inside ));
        }
        else if( intersection.neighbor() )
        {
          const auto volIn = intersection.inside().geometry().volume();
          const auto volOut = intersection.outside().geometry().volume();
          beta_ = penalty_ * intersection.geometry().volume() / std::min( volIn, volOut );
          return true;
        }
      }

      template< class Point >
      ValueType interior ( const Point &x, const ValueType &u ) const
      {
        RangeType source( 0 );
        model().source( x, std::get< 0 >( u ), std::get< 1 >( u ), source );
        JacobianRangeType dFlux( 0 );
        model().diffusiveFlux( x, std::get< 0 >( u ), std::get< 1 >( u ), dFlux );
        return ValueType( source, dFlux );
      }

      template< class Point >
      auto linearizedInterior ( const Point &x, const ValueType &u ) const
      {
        return [ this, x, u ] ( const ValueType &phi ) {
            RangeType source( 0 );
            model().linSource( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), source );
            JacobianRangeType dFlux( 0 );
            model().linDiffusiveFlux( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), dFlux );
            return ValueType( source, dFlux );
          };
      }

      template< class Point >
      ValueType boundary ( const Point &x, const ValueType &u ) const
      {
        RangeType alpha( 0 );
        model().alpha( x, std::get< 0 >( u ), alpha );
        return ValueType( alpha, 0 );
      }

      template< class Point >
      auto linearizedBoundary ( const Point &x, const ValueType &u ) const
      {
        return [ this, x, u ] ( const ValueType &phi ) {
            RangeType alpha( 0 );
            model().linAlpha( std::get< 0 >( u ), x, std::get< 0 >( phi ), alpha );
            return ValueType( alpha, 0 );
          };
      }

      template< class Point >
      std::pair< ValueType, ValueType > skeleton ( const Point &xIn, const ValueType &uIn, const Point &xOut, const ValueType &uOut ) const
      {
        const EntityType inside = intersection().inside();
        const EntityType outside = intersection().outside();

        const RangeFieldType half = RangeFieldType( 1 ) / RangeFieldType( 2 );
        const auto normal = intersection().unitOuterNormal( xIn.localPosition() );

        ValueType uJump( 0, 0 );
        std::get< 0 >( uJump ) = std::get< 0 >( uOut ) - std::get< 0 >( uIn );
        for( int i = 0; i < RangeType::dimension; ++i )
          std::get< 1 >( uJump )[ i ].axpy( std::get< 0 >( uJump )[ i ], normal );

        model().init( outside );
        JacobianRangeType dFluxOut( 0 ), dFluxPrimeOut( 0 );
        model().diffusiveFlux( xOut, std::get< 0 >( uOut ), std::get< 1 >( uJump ), dFluxPrimeOut );
        model().diffusiveFlux( xOut, std::get< 0 >( uOut ), 0, dFluxOut );
        dFluxPrimeOut -= dFluxOut;
        dFluxOut = 0;
        model().diffusiveFlux( xOut, std::get< 0 >( uOut ), std::get< 1 >( uOut ), dFluxOut );

        model().init( inside );
        JacobianRangeType dFluxIn( 0 ), dFluxPrimeIn( 0 );
        model().diffusiveFlux( xIn, std::get< 0 >( uIn ), std::get< 1 >( uJump ), dFluxPrimeIn );
        model().diffusiveFlux( xIn, std::get< 0 >( uIn ), 0, dFluxIn );
        dFluxPrimeIn -= dFluxIn;
        dFluxIn = 0;
        model().diffusiveFlux( xIn, std::get< 0 >( uIn ), std::get< 1 >( uIn ), dFluxIn );

        RangeType int0 = std::get< 0 >( uJump );
        int0 *= beta_;
        dFluxIn += dFluxOut;
        dFluxIn.usmv( -half, normal, int0 );

        dFluxPrimeIn *= -half;
        dFluxPrimeOut *= -half;

        return std::make_pair( ValueType( -int0, dFluxPrimeIn ), ValueType( int0, dFluxPrimeOut ) );
      }

      template< class Point >
      auto linearizedSkeleton ( const Point &xIn, const ValueType &uIn, const Point &xOut, const ValueType &uOut ) const
      {
        const auto normal = intersection().unitOuterNormal( xIn.localPosition() );

        ValueType uJump( 0, 0 );
        std::get< 0 >( uJump ) = std::get< 0 >( uOut ) - std::get< 0 >( uIn );
        for( int i = 0; i < RangeType::dimension; ++i )
          std::get< 1 >( uJump )[ i ].axpy( std::get< 0 >( uJump )[ i ], normal );

        auto intIn = [ this, xIn, uIn, xOut, uOut, normal, uJump ] ( const ValueType &phiIn ) {
          const EntityType inside = intersection().inside();
          const EntityType outside = intersection().outside();

          const RangeFieldType half = RangeFieldType( 1 ) / RangeFieldType( 2 );

          ValueType phiJump( 0, 0 );
          std::get< 0 >( phiJump ) -= std::get< 0 >( phiIn );
          for( int i = 0; i < RangeType::dimension; ++i )
            std::get< 1 >( phiJump )[ i ].axpy( std::get< 0 >( phiJump )[ i ], normal );

          model().init( outside );
          JacobianRangeType dFluxPrimeOut( 0 );
          model().linDiffusiveFlux( std::get< 0 >( uOut ), std::get< 1 >( uJump ), xOut, 0, std::get< 1 >( phiJump ), dFluxPrimeOut );

          model().init( inside );
          JacobianRangeType dFluxIn( 0 ), dFluxPrimeIn( 0 );
          model().linDiffusiveFlux( std::get< 0 >( uIn ), std::get< 1 >( uJump ), xIn, std::get< 0 >( phiIn ), std::get< 1 >( phiJump ), dFluxPrimeIn );
          model().linDiffusiveFlux( std::get< 0 >( uIn ), 0, xIn, std::get< 0 >( phiIn ), 0, dFluxIn );
          dFluxPrimeIn -= dFluxIn;
          dFluxIn = 0;
          model().linDiffusiveFlux( std::get< 0 >( uIn ), std::get< 1 >( uIn ), xIn, std::get< 0 >( phiIn ), std::get< 1 >( phiIn ), dFluxIn );

          RangeType int0 = std::get< 0 >( phiJump );
          int0 *= beta_;
          dFluxIn.usmv( -half, normal, int0 );

          dFluxPrimeIn *= -half;
          dFluxPrimeOut *= -half;

          return std::make_pair( ValueType( -int0, dFluxPrimeIn ), ValueType( int0, dFluxPrimeOut ) );
        };

        auto intOut = [ this, xIn, uIn, xOut, uOut, normal, uJump ] ( const ValueType &phiOut ) {
          const EntityType inside = intersection().inside();
          const EntityType outside = intersection().outside();

          const RangeFieldType half = RangeFieldType( 1 ) / RangeFieldType( 2 );

          ValueType phiJump( 0, 0 );
          std::get< 0 >( phiJump ) = std::get< 0 >( phiOut );
          for( int i = 0; i < RangeType::dimension; ++i )
            std::get< 1 >( phiJump )[ i ].axpy( std::get< 0 >( phiJump )[ i ], normal );

          model().init( outside );
          JacobianRangeType dFluxOut( 0 ), dFluxPrimeOut( 0 );
          model().linDiffusiveFlux( std::get< 0 >( uOut ), std::get< 1 >( uJump ), xOut, std::get< 0 >( phiOut ), std::get< 1 >( phiJump ), dFluxPrimeOut );
          model().linDiffusiveFlux( std::get< 0 >( uOut ), 0, xOut, std::get< 0 >( phiOut ), 0, dFluxOut );
          dFluxPrimeOut -= dFluxOut;
          dFluxOut = 0;
          model().linDiffusiveFlux( std::get< 0 >( uOut ), std::get< 1 >( uOut ), xOut, std::get< 0 >( phiOut ), std::get< 1 >( phiOut ), dFluxOut );

          model().init( inside );
          JacobianRangeType dFluxPrimeIn( 0 );
          model().linDiffusiveFlux( std::get< 0 >( uIn ), std::get< 1 >( uJump ), xIn, 0, std::get< 1 >( phiJump ), dFluxPrimeIn );

          RangeType int0 = std::get< 0 >( phiJump );
          int0 *= beta_;
          dFluxOut.usmv( -half, normal, int0 );

          dFluxPrimeIn *= -half;
          dFluxPrimeOut *= -half;

          return std::make_pair( ValueType( -int0, dFluxPrimeIn ), ValueType( int0, dFluxPrimeOut ) );
        };

        return std::make_pair( intIn, intOut );
      }

      const Model &model () const { return *model_; }

    private:
      const IntersectionType &intersection () const { assert( intersection_ ); return *intersection_; }

      const Model *model_;
      RangeFieldType penalty_;
      const IntersectionType *intersection_ = nullptr;
      RangeFieldType beta_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_INTEGRANDS_HH
