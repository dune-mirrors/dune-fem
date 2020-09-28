#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LOCALRESTRICTPROLONG_HH_
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LOCALRESTRICTPROLONG_HH_

// dune-fem includes
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>

// local includes
#include "declaration.hh"
#include "localdgmassmatrix.hh"


namespace Dune
{

  namespace Fem
  {

    /** @ingroup RestrictProlongImpl
        @{
    **/


    // DiscontinuousGalerkinLocalRestrictProlong
    // -----------------------------------------

    template< class DiscreteFunctionSpace, bool applyInverse >
    class DiscontinuousGalerkinLocalRestrictProlong
    {
      typedef DiscontinuousGalerkinLocalRestrictProlong< DiscreteFunctionSpace, applyInverse > ThisType;

    public:
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

      typedef CachingQuadrature< GridPartType, 0 > QuadratureType;

      typedef LocalMassMatrix< DiscreteFunctionSpaceType, QuadratureType > LocalMassMatrixType;

      DiscontinuousGalerkinLocalRestrictProlong ( const DiscreteFunctionSpaceType& space )
      : localMassMatrix_( space, space.order() * 2 ),
        weight_( -1 ),
        temp_( space )
      {}

      void setFatherChildWeight ( const DomainFieldType &weight )
      {
        weight_ = weight;
      }

      //! restrict data to father
      template< class LFFather, class LFSon, class LocalGeometry >
      void restrictLocal ( LFFather &lfFather, const LFSon &lfSon,
                           const LocalGeometry &geometryInFather, bool initialize ) const
      {
        typedef ConstantLocalRestrictProlong< DiscreteFunctionSpaceType > ConstantLocalRestrictProlongType;
        const DomainFieldType weight = (weight_ < DomainFieldType( 0 ) ? ConstantLocalRestrictProlongType::calcWeight( lfFather.entity(), lfSon.entity() ) : weight_);

        assert( weight > 0.0 );

        if( initialize )
          lfFather.clear();

        if( applyInverse )
        {
          temp_.init( lfFather.entity() );
          temp_.clear();
        }

        typedef typename LFSon :: EntityType  EntityType ;
        typedef typename EntityType :: Geometry   Geometry;
        const EntityType& sonEntity = lfSon.entity();
        const Geometry sonGeo = sonEntity.geometry();

        QuadratureType quad( sonEntity, 2*lfFather.order()+1 );
        const int nop = quad.nop();
        for( int qp = 0; qp < nop; ++qp )
        {
          RangeFieldType quadWeight = quad.weight( qp );

          // in case of non-orthonormal basis we have to
          // apply the integration element and the
          // inverse mass matrix later
          if( applyInverse )
          {
            quadWeight *= sonGeo.integrationElement( quad.point(qp) );
          }
          else
            quadWeight *= weight ;

          RangeType value;
          lfSon.evaluate( quad[ qp ], value );
          value *= quadWeight;

          if( applyInverse )
            temp_.axpy( geometryInFather.global( quad.point( qp ) ), value );
          else
            lfFather.axpy( geometryInFather.global( quad.point( qp ) ), value );
        }

        if( applyInverse )
        {
          localMassMatrix_.applyInverse( temp_ );
          lfFather += temp_;
        }
      }
      template< class LFFather >
      void restrictFinalize ( LFFather &lfFather ) const
      {}

      template< class LFFather, class LFSon, class LocalGeometry >
      void prolongLocal ( const LFFather &lfFather, LFSon &lfSon,
                          const LocalGeometry &geometryInFather, bool initialize ) const
      {
        lfSon.clear();

        typedef typename LFSon :: EntityType  EntityType ;
        typedef typename EntityType :: Geometry   Geometry;
        const EntityType& sonEntity = lfSon.entity();
        const Geometry sonGeo = sonEntity.geometry();

        QuadratureType quad( sonEntity, 2*lfSon.order()+1 );
        const int nop = quad.nop();
        for( int qp = 0; qp < nop; ++qp )
        {
          RangeFieldType quadWeight = quad.weight( qp );

          // in case of non-orthonormal basis we have to
          // apply the integration element and the
          // inverse mass matrix later
          if( applyInverse )
          {
            quadWeight *= sonGeo.integrationElement( quad.point(qp) );
          }

          RangeType value;
          lfFather.evaluate( geometryInFather.global( quad.point( qp ) ), value );
          value *= quadWeight;
          lfSon.axpy( quad[ qp ], value );
        }

        if( applyInverse )
        {
          localMassMatrix_.applyInverse( sonEntity, lfSon );
        }
      }

      bool needCommunication () const { return true; }

    protected:
      LocalMassMatrixType localMassMatrix_;
      DomainFieldType weight_;
      mutable TemporaryLocalFunction< DiscreteFunctionSpace > temp_;
    };



    // DefaultLocalRestrictProlong for DiscontinuousGalerkinSpace
    // ----------------------------------------------------------

    template< class FunctionSpaceImp, class GridPartImp, int polOrd, class StorageImp >
    class DefaultLocalRestrictProlong< DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > >
    : public DiscontinuousGalerkinLocalRestrictProlong< DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp >, false >
    {
    public:
      typedef DiscontinuousGalerkinLocalRestrictProlong< DiscontinuousGalerkinSpace<
        FunctionSpaceImp, GridPartImp, polOrd, StorageImp >, false  >  BaseType;
      DefaultLocalRestrictProlong ( const DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > & space )
        : BaseType( space )
      {}
    };

    template< class FunctionSpaceImp, class GridPartImp, class StorageImp >
    class DefaultLocalRestrictProlong< DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
    : public ConstantLocalRestrictProlong< DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
    {
    public:
      DefaultLocalRestrictProlong ( const DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > & )
      {}
    };



    // DefaultLocalRestrictProlong for LegendreDiscontinuousGalerkinSpace
    // ------------------------------------------------------------------

    template< class FunctionSpaceImp, class GridPartImp, int polOrd, class StorageImp >
    class DefaultLocalRestrictProlong< LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > >
    : public DiscontinuousGalerkinLocalRestrictProlong< LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp >, false >
    {
    public:
      typedef DiscontinuousGalerkinLocalRestrictProlong<
        LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd,  StorageImp >, false > BaseType;
      DefaultLocalRestrictProlong ( const LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > & space )
        : BaseType( space )
      {}
    };

    template< class FunctionSpaceImp, class GridPartImp, class StorageImp >
    class DefaultLocalRestrictProlong< LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
    : public ConstantLocalRestrictProlong< LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
    {
    public:
      DefaultLocalRestrictProlong ( const LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > & )
      {}
    };



    // DefaultLocalRestrictProlong for HierarchicLegendreDiscontinuousGalerkinSpace
    // ----------------------------------------------------------------------------

    template< class FunctionSpaceImp, class GridPartImp, int polOrd, class StorageImp >
    class DefaultLocalRestrictProlong< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > >
    : public DiscontinuousGalerkinLocalRestrictProlong< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp >, false >
    {
    public:
      typedef DiscontinuousGalerkinLocalRestrictProlong<
        HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd,  StorageImp >, false > BaseType;
      DefaultLocalRestrictProlong ( const HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > & space )
        : BaseType( space )
      {}
    };

    template< class FunctionSpaceImp, class GridPartImp, class StorageImp >
    class DefaultLocalRestrictProlong< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
    : public ConstantLocalRestrictProlong< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
    {
    public:
      DefaultLocalRestrictProlong ( const HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > & )
      {}
    };



    // DefaultLocalRestrictProlong for LagrangeDiscontinuousGalerkinSpace
    // ------------------------------------------------------------------

    template< class FunctionSpaceImp, class GridPartImp, int polOrd, class StorageImp >
    class DefaultLocalRestrictProlong< LagrangeDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > >
    : public DiscontinuousGalerkinLocalRestrictProlong< LagrangeDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp >, true >
    {
    public:
      typedef DiscontinuousGalerkinLocalRestrictProlong< LagrangeDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp >, true >  BaseType ;
      DefaultLocalRestrictProlong ( const LagrangeDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > & space )
        : BaseType( space )
      {}
    };

    template< class FunctionSpaceImp, class GridPartImp, class StorageImp >
    class DefaultLocalRestrictProlong< LagrangeDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
    : public ConstantLocalRestrictProlong< LagrangeDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
    {
    public:
      DefaultLocalRestrictProlong ( const LagrangeDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > & )
      {}
    };

    ///@}

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LOCALRESTRICTPROLONG_HH_
