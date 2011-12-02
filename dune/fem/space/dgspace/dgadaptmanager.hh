#ifndef DUNE_DGADAPTMANAGERIMP_HH
#define DUNE_DGADAPTMANAGERIMP_HH

#include <dune/fem/function/localfunction/temporarylocalfunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/dgspace/dgspace.hh>

namespace Dune
{

  /** @ingroup RestrictProlongImpl
      @{
  **/


  // DiscontinuousGalerkinLocalRestrictProlong
  // -----------------------------------------

  template< class DiscreteFunctionSpace >
  class DiscontinuousGalerkinLocalRestrictProlong
  {
    typedef DiscontinuousGalerkinLocalRestrictProlong< DiscreteFunctionSpace > ThisType;

  public:
    typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;

    DiscontinuousGalerkinLocalRestrictProlong ()
    : weight_( -1 )
    {}

    void setFatherChildWeight ( const DomainFieldType &weight )
    {
      weight_ = weight;
    }

    //! restrict data to father 
    template< class FT, class ST, class LocalGeometry >
    void restrictLocal ( LocalFunction< FT > &lfFather, const LocalFunction< ST > &lfSon, 
                         const LocalGeometry &geometryInFather, bool initialize ) const
    {
      typedef ConstantLocalRestrictProlong< DiscreteFunctionSpaceType > ConstantLocalRestrictProlongType;
      const DomainFieldType weight = (weight_ < DomainFieldType( 0 ) ? ConstantLocalRestrictProlongType::calcWeight( lfFather.entity(), lfSon.entity() ) : weight_); 

      assert( weight > 0.0 );
      //assert( std::abs( geometryInFather.volume() - weight ) < 1e-8 );

      if( initialize )
        lfFather.clear();
      
      QuadratureType quad( lfSon.entity(), 2*lfFather.order()+1 );
      const int nop = quad.nop();
      for( int qp = 0; qp < nop; ++qp )
      {
        RangeType value;
        lfSon.evaluate( quad[ qp ], value );
        value *= weight * quad.weight( qp );
        lfFather.axpy( geometryInFather.global( quad.point( qp ) ), value );
      }
    }

    template< class FT, class ST, class LocalGeometry >
    void prolongLocal ( const LocalFunction< FT > &lfFather, LocalFunction< ST > &lfSon,
                        const LocalGeometry &geometryInFather, bool initialize ) const
    {
      lfSon.clear();

      QuadratureType quad( lfSon.entity(), 2*lfSon.order()+1 );
      const int nop = quad.nop();
      for( int qp = 0; qp < nop; ++qp )
      {
        RangeType value;
        lfFather.evaluate( geometryInFather.global( quad.point( qp ) ), value );
        value *= quad.weight( qp );
        lfSon.axpy( quad[ qp ], value );
      }
    }

    bool needCommunication () const { return true; }

  protected:
    DomainFieldType weight_;
  };



  // DefaultLocalRestrictProlong for DiscontinuousGalerkinSpace
  // ----------------------------------------------------------

  template< class FunctionSpaceImp, class GridPartImp, int polOrd, template< class > class StorageImp >
  struct DefaultLocalRestrictProlong< DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > >
  : public DiscontinuousGalerkinLocalRestrictProlong< DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > >
  {
    DefaultLocalRestrictProlong ( const DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > & )
    {}
  };

  template< class FunctionSpaceImp, class GridPartImp, template< class > class StorageImp >
  struct DefaultLocalRestrictProlong< DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
  : public ConstantLocalRestrictProlong< DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
  {
    DefaultLocalRestrictProlong ( const DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > & )
    {}
  };



  // DefaultLocalRestrictProlong for LegendreDiscontinuousGalerkinSpace
  // ------------------------------------------------------------------

  template< class FunctionSpaceImp, class GridPartImp, int polOrd, template< class > class StorageImp >
  struct DefaultLocalRestrictProlong< LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > >
  : public DiscontinuousGalerkinLocalRestrictProlong< LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > >
  {
    DefaultLocalRestrictProlong ( const LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > & )
    {}
  };

  template< class FunctionSpaceImp, class GridPartImp, template< class > class StorageImp >
  struct DefaultLocalRestrictProlong< LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
  : public ConstantLocalRestrictProlong< LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
  {
    DefaultLocalRestrictProlong ( const LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > & )
    {}
  };

  ///@}

} // namespace Dune

#endif // #ifndef DUNE_DGADAPTMANAGERIMP_HH
