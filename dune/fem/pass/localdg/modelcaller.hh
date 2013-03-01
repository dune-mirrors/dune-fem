#ifndef DUNE_FEM_PASS_LOCALDG_MODELCALLER_HH
#define DUNE_FEM_PASS_LOCALDG_MODELCALLER_HH

#include <cstddef>
#include <vector>

#include <dune/fem/pass/common/filter.hh>
#include <dune/fem/pass/common/localfunctiontuple.hh>
#include <dune/fem/pass/common/pass.hh>
#include <dune/fem/pass/common/pointertuple.hh>
#include <dune/fem/pass/common/selector.hh>
#include <dune/fem/pass/common/tupletypetraits.hh>

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
      typedef typename DiscreteModelType::VolumeQudratureType VolumeQuadratureType;
      // type of face quadrature
      typedef typename DiscreteModelType::FaceQudratureType FaceQuadratureType;

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
      typedef PointerTuple< typename FilterType::type > DiscreteFunctionPointerTupleType;

      typedef typename DiscreteFunctionPointerTupleType::ElementType DiscreteFunctionTupleType;
      typedef LocalFunctionTuple< DiscreteFunctionTupleType, EntityType > LocalFunctionTupleType;

      typedef typename LocalFunctionTupleType::RangeTupleType RangeTupleType;
      typedef typename LocalFunctionTupleType::JacobianRangeTupleType JacobianRangeTupleType;

    public:
      DGDiscreteModelCaller ( ArgumentType &argument, DiscreteModelType &discreteModel )
      : discreteModel_( discreteModel ),
        time_( 0. ),
        discreteFunctions_( FilterType::apply( argument ) ),
        localFunctionsInside_( *discreteFunctions_ ),
        localFunctionsOutside_( *discreteFunctions_ )
#ifndef NDEBUG
      //
#endif // #ifndef NDEBUG
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
        localFunctionsInside_.setEntity( entity );
        discreteModel().setEntity( entity );
      }

      // set outside local functions to entity 
      void setNeighbor ( const EntityType &entity )
      {
        localFunctionsOutside_.setEntity( entity );
        discreteModel().setNeighbor( entity );
      }

      // evaluate inside local functions in all quadrature points
      void setEntity ( const EntityType &entity, const VolumeQuadratureType &quadrature )
      {
        setEntity( entity );
        resize( values_, quadrature.nop() );
        localFunctionsInside_.evaluateQuadrature( quadrature, values_ );

#ifndef NDEBUG
        // ...
#endif // #ifndef NDEBUG
      }

      // evaluate outside local functions in all quadrature points
      template< class QuadratureType >
      void setNeighbor ( const EntityType &neighbor,
                         const QuadratureType &inner,
                         const QuadratureType &outer )
      {
        // we assume that setEntity() was called in advance!
        resize( valuesInside_, inner.nop() );
        localFunctionsInside_.evaluateQuadrature( inner, valuesInside_ );

        setNeighbor( neighbor );

        resize( valuesOutside_, outer.nop() );
        localFunctionsOutside_.evaluateQuadrature( outer, valuesOutside_ );

#ifndef NDEBUG
        // ...
#endif // #ifndef NDEBUG
      }

      // please doc me
      template< class QuadratureType >
      void setBoundary ( const EntityType &entity, const QuadratureType &quadrature )
      {
        resize( valuesInside_, quadrature.nop() );
        localFunctionsInside_.evaluateQuadrature( quadrature, valuesInside_ );

#ifndef NDEBUG
        // ...
#endif // #ifndef NDEBUG
      }

      // evaluate analytical flux
      void analyticalFlux ( const EntityType &entity,
                            const VolumeQuadratureType &quadrature,
                            const int qp,
                            JacobianRangeType &flux ) const
      {
#ifndef NDEBUG
        //
#endif // #ifndef NDEBUG

        assert( hasFlux() );
        discreteModel().analyticalFlux( entity, time(), quadrature.point( qp ), values_, flux );
      }

      // evaluate analytical flux
      void source ( const EntityType &entity,
                    const VolumeQuadratureType &quadrature,
                    const int qp,
                    RangeType &source ) const
      {
#ifndef NDEBUG
        //
#endif // #ifndef NDEBUG

        assert( hasSource() );
        discreteModel().source( entity, time(), quadrature.point( qp ), values_[ qp ], jacobians_[ qp ], source );
      }


      // evaluate both analytical flux and source
      double analyticalFluxAndSource ( const EntityType &entity,
                                       const VolumeQuadratureType &quadrature,
                                       const int qp,
                                       JacobianRangeType &flux,
                                       RangeType &source ) const
      {
        analyticalFlux( entity, quadrature, qp, flux );
        return source( entity, quadrature, qp, source );
      }

      // evaluate numerical flux
      template< class QuadratureType >
      double numericalFlux ( const IntersectionType &intersection,
                             const QuadratureType &inner,
                             const QuadratureType &outer,
                             const int qp,
                             RangeType &gLeft, RangeType &gRight ) const
      {
#ifndef NDEBUG
        // ...
#endif // #ifndef NDEBUG

        return discreteModel().numericalFlux( intersection, time(), inner.localPoint( qp ), valuesInside_[ qp ], valuesOutside_[ qp ], gLeft, gRight );
      }

      // evaluate boundary flux
      double boundaryFlux ( const IntersectionType &intersection,
                            const FaceQuadratureType &quadrature,
                            const int qp,
                            RangeType &gLeft ) const
      {
#ifndef NDEBUG
        // ...
#endif // #ifndef NDEBUG

        return discreteModel().boundaryFlux( intersection, time(), quadrature.localPoint( qp ), valuesInside_[ qp ], gLeft );
      }

      // evaluate mass
      void mass ( const EntityType &entity,
                  const VolumeQuadratureType &quadrature,
                  const int qp,
                  MassFactorType &massFactor ) const
      {
#ifndef NDEBUG
        // ...
#endif // #ifndef NDEBUG

        discreteModel().mass( entity, time(), quadrature.point( qp ), values_[ qp ], massFactor );
      }

    protected:
      DiscreteModelType &discreteModel () { return discreteModel_; }
      const DiscreteModelType &discreteModel () const { return discreteModel_; }

    private:
      // resize vector if necessary
      template< class Quadrature, class Vector >
      bool resize ( Vector &vector, const int size )
      {
        if( vector.size() >= size )
          return false;

        vector.resize( size );
        return true;
      }

      // forbid copying and assignment
      DGDiscreteModelCaller ( const ThisType & );
      ThisType operator= ( const ThisType & );

      DiscreteModelType &discreteModel_;
      double time_;
      DiscreteFunctionPointerTupleType discreteFunctions_;
      LocalFunctionTupleType localFunctionsInside_, localFunctionsOutside_;
      std::vector< RangeTupleType > values_, valuesInside_, valuesOutside_;
      std::vector< JacobianRangeTupleType > jacobians_;
#ifndef NDEBUG
      //
#endif // #ifndef NDEBUG
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_LOCALDG_MODELCALLER_HH
