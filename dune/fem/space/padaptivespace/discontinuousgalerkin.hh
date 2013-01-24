#ifndef DUNE_FEM_SPACE_PADAPTIVE_DISCONTINUOUSGALERKIN_HH
#define DUNE_FEM_SPACE_PADAPTIVE_DISCONTINUOUSGALERKIN_HH

#include <dune/fem/operator/projection/dgl2projection.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

#include "adaptmanager.hh"
#include "declaration.hh"
#include "generic.hh"
#include "lagrange.hh"
#include "mapper.hh"
#include "restrictprolong.hh"


namespace Dune
{

  namespace Fem 
  {

    /** \addtogroup PAdaptiveDGSpace
     *
     *  Provides access to base function sets for different element types in
     *  one grid and size of function space and maps from local to global dof
     *  number.
     *
     *  \note This space can only be used with special index sets. If you want
     *  to use the PAdaptiveDGSpace with an index set only
     *  supporting the index set interface you will have to use the
     *  IndexSetWrapper class to provide the required functionality.
     *
     *  \note For adaptive calculations one has to use index sets that are
     *  capable of adaption (i.e. the method adaptive returns true). See also
     *  AdaptiveLeafIndexSet.
     */

    // PAdaptiveDGSpaceTraits
    // ----------------------

    template< class FunctionSpace, class GridPart, unsigned int polOrder,
              template< class > class BaseFunctionStorage >
    struct PAdaptiveDGSpaceTraits 
      : public PAdaptiveLagrangeSpaceTraits
          < FunctionSpace, GridPart, polOrder, BaseFunctionStorage >
    {
      typedef PAdaptiveDGSpace
        < FunctionSpace, GridPart, polOrder, BaseFunctionStorage >
        DiscreteFunctionSpaceType;

      enum { localBlockSize = FunctionSpace :: dimRange };

      //! this is a continuous space 
      static const bool continuousSpace = false ;

      // mapper for block
      typedef PAdaptiveDGMapper< GridPart, polOrder > BlockMapperType;
      typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;
      
      /** \brief defines type of communication data handle for this type of space
       */
      template< class DiscreteFunction,
                class Operation = DFCommunicationOperation :: Copy >
      struct CommDataHandle
      {
        //! type of data handle 
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        //! type of operatation to perform on scatter 
        typedef Operation OperationType;
      };
    };



    // PAdaptiveDGSpace
    // ----------------

    /** \class   PAdaptiveDGSpace
     *  \ingroup PAdaptiveDGSpace
     *  \brief   adaptive DG discrete function space
     */
    template< class FunctionSpaceImp,
              class GridPartImp,
              int polOrder,
              template< class > class BaseFunctionStorageImp = CachingStorage >
    class PAdaptiveDGSpace
    : public GenericDiscreteFunctionSpace
             < PAdaptiveDGSpaceTraits< FunctionSpaceImp,
                                       GridPartImp,
                                       polOrder,
                                       BaseFunctionStorageImp > >
    {
    public:
      //! traits for the discrete function space
      typedef PAdaptiveDGSpaceTraits< FunctionSpaceImp,
                                      GridPartImp,
                                      polOrder,
                                      BaseFunctionStorageImp >
        Traits;
      //! type of the discrete function space
      typedef PAdaptiveDGSpace< FunctionSpaceImp,
                                GridPartImp,
                                polOrder,
                                BaseFunctionStorageImp >
              PAdaptiveDGSpaceType;

    private:
      typedef GenericDiscreteFunctionSpace< Traits > BaseType;
      typedef PAdaptiveDGSpaceType ThisType;
    public:

      typedef typename Traits :: GridPartType GridPartType;
      typedef typename Traits :: GridType GridType;
      typedef typename Traits :: IndexSetType IndexSetType;
      typedef typename Traits :: IteratorType IteratorType;
      typedef typename BaseType ::IntersectionType IntersectionType;

      //! maximum polynomial order of functions in this space
      enum { polynomialOrder = Traits :: polynomialOrder };
      
      //! type of compiled local key 
      typedef typename Traits :: CompiledLocalKeyType  CompiledLocalKeyType;

      // deprecated name 
      typedef CompiledLocalKeyType LagrangePointSetType;

      //! mapper used to implement mapToGlobal
      typedef typename Traits :: MapperType MapperType;

      //! mapper used to for block vector function 
      typedef typename Traits :: BlockMapperType BlockMapperType;

      //! size of local blocks
      enum { localBlockSize = Traits :: localBlockSize };

      //! dimension of a value
      enum { dimVal = 1 };

      //! type of DoF manager
      typedef DofManager< GridType > DofManagerType;

      using BaseType :: gridPart;
      using BaseType :: blockMapper;
      using BaseType :: compiledLocalKey;

      //! default communication interface 
      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;

      //! default communication direction 
      static const CommunicationDirection defaultDirection = ForwardCommunication;

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part for the Lagrange space
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       */
      explicit PAdaptiveDGSpace
        ( GridPartType &gridPart,
          const InterfaceType commInterface = defaultInterface,
          const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection )
      {
      }

      //! copy constructor needed for p-adaption 
      PAdaptiveDGSpace( const PAdaptiveDGSpace& other ) 
      : BaseType( other )
      {
      }

    protected:
      using BaseType :: dfList_;
      using BaseType :: searchFunction;

    public:
      /*! \brief add function to discrete function space for p-adaptation 
          (currently only supported by AdaptiveDiscreteFunction )
       */
      template <class DiscreteFunction> 
      void addFunction( DiscreteFunction& df ) const
      {
        assert( searchFunction( df ) == dfList_.end() );
        // select L2Porjection to be the LocalInterpolation 
        typedef typename BaseType :: template PAdaptiveDiscreteFunctionEntry< 
            DiscreteFunction, DGL2ProjectionImpl > RealEntryType ;
        typedef typename BaseType :: PAdaptiveDiscreteFunctionEntryInterface
          EntryInterface;

        EntryInterface* entry = new RealEntryType( df );
        assert( entry );
        dfList_.push_front( entry );
      }

      //! deprecated method 
      template< class EntityType >
      inline const CompiledLocalKeyType &lagrangePointSet( const EntityType &entity ) const
      {
        return compiledLocalKey( entity.type(),
                                 blockMapper().polynomOrder( entity ) );
      }

      //! deprecated method 
      inline const CompiledLocalKeyType &lagrangePointSet( const GeometryType type ) const
      {
        return compiledLocalKey( type, polynomialOrder );
      }

      //! deprecated method 
      inline const CompiledLocalKeyType &lagrangePointSet( const GeometryType type, const int order ) const
      {
        return compiledLocalKey( type, order );
      }

      using BaseType::continuous;
      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      inline bool continuous (const IntersectionType &intersection) const
      { 
        return false;
      }

    };

  } // namespace Fem
    
} // Dune namespace  

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVE_DISCONTINUOUSGALERKIN_HH
