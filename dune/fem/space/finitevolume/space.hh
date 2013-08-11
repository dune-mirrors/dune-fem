#ifndef DUNE_FEM_SPACE_FVSPACE_SPACE_HH
#define DUNE_FEM_SPACE_FVSPACE_SPACE_HH

#include <dune/fem/function/localfunction/average.hh>
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/version.hh>

#include <dune/fem/space/finitevolume/basisfunctionset.hh>
#include <dune/fem/space/finitevolume/declaration.hh>

/*
  @file
  @brief Finite Volume space
  @author Christoph Gersbacher
*/

namespace Dune
{

  namespace Fem
  {

    // FiniteVolumeSpaceTraits
    // -----------------------

    /* \class FiniteVolumeSpaceTraits
     *
     * \brief Finite volume space
     *
     * \tparam  FunctionSpace  Function space
     * \tparam  GridPart       Grid part
     * \tparam  codim          codimension
     * \tparam  Storage        Caching storage policy
     *
     * \note This shape function set has fixed polynomial order 0.
     */
    template< class FunctionSpace, class GridPart, int codim,
              template< class > class Storage >
    struct FiniteVolumeSpaceTraits
    {
      typedef FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = codim;
      
    private:
      static const int dimLocal = GridPartType::dimension; 
      static const int dimRange = FunctionSpaceType::dimRange;

      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

    public:
      typedef FiniteVolumeBasisFunctionSet< EntityType, typename FunctionSpaceType::RangeType > BasisFunctionSetType;

      static const int localBlockSize = dimRange;

      typedef CodimensionMapper< GridPartType, codimension > BlockMapperType;
      typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;

      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };
    };



    // FiniteVolumeSpace
    // -----------------
 
    template< class FunctionSpace, class GridPart, int codim = 0,
              template< class > class Storage = SimpleStorage >
    struct FiniteVolumeSpace
    : public DiscreteFunctionSpaceDefault< FiniteVolumeSpaceTraits< FunctionSpace, GridPart, codim, Storage > >
    {
      typedef FiniteVolumeSpaceTraits< FunctionSpace, GridPart, codim, Storage > Traits;

    private:
      typedef FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;
      static const CommunicationDirection defaultDirection =  ForwardCommunication;

    public:
      static const int polynomialOrder DUNE_DEPRECATED = 0;

      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::IteratorType IteratorType;
      typedef typename IteratorType::Entity EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::MapperType MapperType;
      typedef typename BaseType::BlockMapperType BlockMapperType;

      explicit FiniteVolumeSpace( GridPartType &gridPart,
                                  const InterfaceType commInterface = defaultInterface,
                                  const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        blockMapper_( gridPart ),
        mapper_( blockMapper_ )
      {}

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      DFSpaceIdentifier type () const
      {
        return FiniteVolumeSpace_id;
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet */
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity );
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous () const
      {
        return false;
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType &intersection ) const
      {
        return false;
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const
      {
        return 0;
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order ( const EntityType & ) const
      {
        return order();          
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::mapper */
      DUNE_VERSION_DEPRECATED(1,4,remove)
      MapperType &mapper () const
      {
        return mapper_;
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const
      {
        return blockMapper_;
      }

      ///////////////////////////
      // Non-interface methods //
      ///////////////////////////

      /** \brief local interpolation
       *
       *  \param[in]  localFunction  local function to interpolate
       *  \param[in]  dofs           local degrees of freedom of the interpolation
       */
      template< class LocalFunction, class LocalDofVector >
      void interpolate ( const LocalFunction &localFunction, LocalDofVector &dofs ) const
      {
        typename LocalFunction::RangeType value;
        LocalAverage< LocalFunction, GridPartType >::apply( localFunction, value );
        for( int i = 0; i < FunctionSpaceType::dimRange; ++i )
          dofs[ i ] = value[ i ];
      }

    private:
      mutable BlockMapperType blockMapper_;
      mutable MapperType mapper_; 
    };



    // DefaultLocalRestrictProlong for FiniteVolumeSpace
    // -------------------------------------------------

    template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
    struct DefaultLocalRestrictProlong< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
    : public ConstantLocalRestrictProlong< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
    {
      DefaultLocalRestrictProlong ( const FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > & )
      {}
    };



    // Dune::Fem::Capabilities for FiniteVolumeSpace
    // ---------------------------------------------

    namespace Capabilities
    {

      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct hasFixedPolynomialOrder< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct hasStaticPolynomialOrder< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
        static const int order = 0; 
      };


      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct isContinuous< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct isLocalized< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false; // there is no method 'shapeFunctionSet( const EntityType & )'
      };


      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct isParallel< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::isParallel< GridPart >::v;
      };


      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct isAdaptive< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct threadSafe< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct viewThreadSafe< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FVSPACE_SPACE_HH
