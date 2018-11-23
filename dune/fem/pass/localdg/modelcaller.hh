#ifndef DUNE_FEM_PASS_LOCALDG_MODELCALLER_HH
#define DUNE_FEM_PASS_LOCALDG_MODELCALLER_HH

#if HAVE_DUNE_FEM_DG
#error "Outdated header, #include <dune/fem-dg/pass/modelcaller.hh> instead!"
#endif

#include <cstddef>
#include <vector>

#include <dune/fem/common/tupletypetraits.hh>
#include <dune/fem/common/tupleutility.hh>
#include <dune/fem/common/typeindexedtuple.hh>
#include <dune/fem/pass/common/localfunctiontuple.hh>
#include <dune/fem/pass/common/pass.hh>
#include <dune/fem/pass/common/selector.hh>

namespace Dune
{

  namespace Fem
  {

    // DGDiscreteModelCaller
    // ---------------------

    /** \brief   model caller for local DG pass
     *  \class   DGDiscreteModelCaller
     *  \ingroup PassHyp
     *
     *  \todo please doc me
     */
    template< class DiscreteModel, class Argument, class PassIds >
    class DGDiscreteModelCaller
    {
      typedef DGDiscreteModelCaller< DiscreteModel, Argument, PassIds > ThisType;

    public:
      // discrete model type
      typedef DiscreteModel DiscreteModelType;
      // argument type
      typedef Argument ArgumentType;

      /** \brief selector (tuple of integral constants) */
      typedef typename DiscreteModelType::Selector Selector;

      // entity type
      typedef typename DiscreteModelType::EntityType EntityType;
      // intersection type
      typedef typename DiscreteModelType::IntersectionType IntersectionType;

      // type of volume quadrature
      typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
      // type of face quadrature
      typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;

      // type of mass factor
      typedef typename DiscreteModelType::MassFactorType MassFactorType;

      // function space type
      typedef typename DiscreteModelType::FunctionSpaceType FunctionSpaceType;
      // range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      // jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    private:
      typedef Dune::MakeSubTuple< ArgumentType, typename Dune::FirstTypeIndexTuple< PassIds, Selector >::type > FilterType;
      typedef typename FilterType::type DiscreteFunctionPointerTupleType;

      typedef typename TupleTypeTraits< DiscreteFunctionPointerTupleType >::PointeeTupleType DiscreteFunctionTupleType;

    protected:
      typedef LocalFunctionTuple< DiscreteFunctionTupleType, EntityType > LocalFunctionTupleType;

      typedef typename LocalFunctionTupleType::RangeTupleType RangeTupleType;
      typedef typename LocalFunctionTupleType::JacobianRangeTupleType JacobianRangeTupleType;

    public:
      DGDiscreteModelCaller ( ArgumentType &argument, DiscreteModelType &discreteModel )
      : discreteModel_( discreteModel ),
        time_( 0. ),
        discreteFunctions_( FilterType::apply( argument ) ),
        localFunctionsInside_( Dune::DereferenceTuple< DiscreteFunctionPointerTupleType >::apply( discreteFunctions_ ) ),
        localFunctionsOutside_( Dune::DereferenceTuple< DiscreteFunctionPointerTupleType >::apply( discreteFunctions_ ) )
      {}

      // return true, if discrete model has flux
      bool hasFlux () const { return discreteModel().hasFlux(); }
      // return true, if discrete model has mass
      bool hasMass () const { return discreteModel().hasMass(); }
      // return true, if discrete model has source
      bool hasSource () const { return discreteModel().hasSource(); }

      // set time
      void setTime ( double time ) { time_ = time; }
      // return time
      double time () const { return time_; }

      // set inside local funtions to entity
      void setEntity ( const EntityType &entity )
      {
        localFunctionsInside_.init( entity );
        discreteModel().setEntity( entity );
      }

      // set outside local functions to entity
      void setNeighbor ( const EntityType &entity )
      {
        localFunctionsOutside_.init( entity );
        discreteModel().setNeighbor( entity );
      }

      // evaluate inside local functions in all quadrature points
      void setEntity ( const EntityType &entity, const VolumeQuadratureType &quadrature )
      {
        setEntity( entity );
        values_.resize( quadrature.nop() );
        jacobians_.resize( quadrature.nop() );
        localFunctionsInside_.evaluateQuadrature( quadrature, values_ );
        localFunctionsInside_.evaluateQuadrature( quadrature, jacobians_ );
      }

      // evaluate outside local functions in all quadrature points
      template< class QuadratureType >
      void setNeighbor ( const EntityType &neighbor,
                         const QuadratureType &inside,
                         const QuadratureType &outside )
      {
        // we assume that setEntity() was called in advance!
        valuesInside_.resize( inside.nop() );
        localFunctionsInside_.evaluateQuadrature( inside, valuesInside_ );

        setNeighbor( neighbor );

        valuesOutside_.resize( outside.nop() );
        localFunctionsOutside_.evaluateQuadrature( outside, valuesOutside_ );
      }

      // please doc me
      template< class QuadratureType >
      void setBoundary ( const EntityType &entity, const QuadratureType &quadrature )
      {
        valuesInside_.resize( quadrature.nop() );
        localFunctionsInside_.evaluateQuadrature( quadrature, valuesInside_ );
      }

      // evaluate analytical flux
      void analyticalFlux ( const EntityType &entity,
                            const VolumeQuadratureType &quadrature,
                            const int qp,
                            JacobianRangeType &flux )
      {
        assert( hasFlux() );
        discreteModel().analyticalFlux( entity, time(), quadrature.point( qp ), values_[ qp ], flux );
      }

      // evaluate analytical flux
      double source ( const EntityType &entity,
                      const VolumeQuadratureType &quadrature,
                      const int qp,
                      RangeType &source )
      {
        assert( hasSource() );
        return discreteModel().source( entity, time(), quadrature.point( qp ), values_[ qp ], jacobians_[ qp ], source );
      }


      // evaluate both analytical flux and source
      double analyticalFluxAndSource ( const EntityType &entity,
                                       const VolumeQuadratureType &quadrature,
                                       const int qp,
                                       JacobianRangeType &flux,
                                       RangeType &source )
      {
        // we may only assume that hasSource() == true, cf. pass.hh
        if( hasFlux() )
          analyticalFlux( entity, quadrature, qp, flux );
        return ThisType::source( entity, quadrature, qp, source );
      }

      // evaluate numerical flux
      template< class QuadratureType >
      double numericalFlux ( const IntersectionType &intersection,
                             const QuadratureType &inside,
                             const QuadratureType &outside,
                             const int qp,
                             RangeType &gLeft, RangeType &gRight )
      {
        return discreteModel().numericalFlux( intersection, time(), inside.localPoint( qp ), valuesInside_[ qp ], valuesOutside_[ qp ], gLeft, gRight );
      }

      // evaluate boundary flux
      double boundaryFlux ( const IntersectionType &intersection,
                            const FaceQuadratureType &quadrature,
                            const int qp,
                            RangeType &gLeft )
      {
        return discreteModel().boundaryFlux( intersection, time(), quadrature.localPoint( qp ), valuesInside_[ qp ], gLeft );
      }

      // evaluate mass
      void mass ( const EntityType &entity,
                  const VolumeQuadratureType &quadrature,
                  const int qp,
                  MassFactorType &massFactor )
      {
        discreteModel().mass( entity, time(), quadrature.point( qp ), values_[ qp ], massFactor );
      }

      DGDiscreteModelCaller ( const ThisType & ) = delete;
      ThisType operator= ( const ThisType & ) = delete;

    protected:
      DiscreteModelType &discreteModel () { return discreteModel_; }
      const DiscreteModelType &discreteModel () const { return discreteModel_; }

    private:
      DiscreteModelType &discreteModel_;
      double time_;
      DiscreteFunctionPointerTupleType discreteFunctions_;

    protected:
      LocalFunctionTupleType localFunctionsInside_, localFunctionsOutside_;
      std::vector< Dune::TypeIndexedTuple< RangeTupleType, Selector > > values_, valuesInside_, valuesOutside_;
      std::vector< Dune::TypeIndexedTuple< JacobianRangeTupleType, Selector > > jacobians_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_LOCALDG_MODELCALLER_HH
