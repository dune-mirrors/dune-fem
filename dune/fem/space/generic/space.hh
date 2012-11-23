#ifndef DUNE_FEM_SPACE_GENERICDISCRETE_HH
#define DUNE_FEM_SPACE_GENERICDISCRETE_HH

#if HAVE_DUNE_LOCALFUNCTIONS

// C++ includes
#include <cassert>
#include <vector>

// dune-common includes
#include <dune/common/nullptr.hh>
#include <dune/common/static_assert.hh>

// dune-geometry types
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/localfunctions.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

// local includes
#include "declaration.hh"
#include "dofmappercode.hh"
#include "localinterpolation.hh"

/**
  @file
  @author Stefan Girke
  @brief  Provides space based on general local finite element.
*/



namespace Dune
{

  namespace Fem
  {

    // GenericDiscreteFunctionSpaceTraits
    // -----------------------------------------

    template< class GridPart, class LocalFiniteElement, int polOrder, template< class > class Storage >
    struct GenericDiscreteFunctionSpaceTraits
    {
      dune_static_assert( Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::v,
                          "GridPart has more than one geometry type." );

      typedef GenericDiscreteFunctionSpace< GridPart, LocalFiniteElement, polOrder, Storage > DiscreteFunctionSpaceType;

      typedef Dune::Fem::FunctionSpace< typename LocalFiniteElement::Traits::LocalBasisType::Traits::DomainFieldType,
                                   typename LocalFiniteElement::Traits::LocalBasisType::Traits::RangeFieldType,
                                   LocalFiniteElement::Traits::LocalBasisType::Traits::dimDomain,
                                   LocalFiniteElement::Traits::LocalBasisType::Traits::dimRange >
			  FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = 0;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      static const int dimLocal = GridPartType::dimension;

    public:
      typedef LocalFiniteElement LocalFiniteElementType;
      typedef typename LocalFiniteElementType::LocalBasisType LocalBasisType;
      typedef typename LocalFiniteElementType::LocalCoefficientsType LocalCoefficientsType;
      typedef typename LocalFiniteElementType::LocalInterpolationType LocalInterpolationType;

      typedef GenericBlockMapperFactory< GridPartType, LocalCoefficientsType > BlockMapperFactoryType;
      typedef typename BlockMapperFactoryType::BlockMapperType BlockMapperType;

      static const int localBlockSize = FunctionSpaceType::dimRange;
      typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;

      typedef LocalFunctionsShapeFunctionSet< LocalBasisType > LocalFunctionsShapeFunctionSetType;
			typedef LocalFiniteElement ShapeFunctionSetType;

      struct ShapeFunctionSetFactory
      {
        static ShapeFunctionSetType *createObject ( const GeometryType &type )
        {
          assert( type.isCube() );
          return new ShapeFunctionSetType( type, LocalFunctionsShapeFunctionSetType( LocalBasisType() ) );
        }
        
        static void deleteObject ( ShapeFunctionSetType *object ) { delete object; }
      };  
		
			template< class GeometryImp, class LocalFiniteElementKey >
			class ShapeFunctionSetSingletonKey
			{
				typedef ShapeFunctionSetSingletonKey< GeometryImp, LocalFiniteElementKey > ThisType;

			public:
				typedef GeometryImp GeometryType;

				ShapeFunctionSetSingletonKey( const GeometryType &geo, const LocalFiniteElementKey& key )
				: geo_( geo ),
					key_( key )
				{}

				const GeometryType& type() const { return geo_; }

				const LocalFiniteElementKey key() const { return key_; }
			
				bool operator== ( const ThisType &other ) const
				{
					return ( key_ == other.key_ && geo_ == other.geo_  );
				}
				
				bool operator!= ( const ThisType &other ) const
				{
					return !( *this == other );
				}

			private:
				const GeometryType& geo_;
				const LocalFiniteElementKey& key_;
			};

      typedef ShapeFunctionSetFactory ShapeFunctionSetFactoryType;

			typedef ShapeFunctionSetSingletonKey< GeometryType, int > ShapeFunctionSetSingletonKeyType;

      typedef DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };
    };



    // GenericDiscreteFunctionSpace
    // -----------------------------------

    template< class GridPart, class LocalFiniteElement, int polOrder, template< class > class Storage = CachingStorage >
    struct GenericDiscreteFunctionSpace
    : public DiscreteFunctionSpaceDefault< GenericDiscreteFunctionSpaceTraits< GridPart, LocalFiniteElement, polOrder, Storage > >
    {
      typedef GenericDiscreteFunctionSpace< GridPart, LocalFiniteElement, polOrder, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< GenericDiscreteFunctionSpaceTraits< GridPart, LocalFiniteElement, polOrder, Storage > > BaseType;

      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;
      static const CommunicationDirection defaultDirection =  ForwardCommunication;

    public:
      static const int polynomialOrder = polOrder;

      typedef typename BaseType::Traits Traits;

      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;
     
      typedef typename BaseType::Traits::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
      
      typedef typename BaseType::MapperType MapperType;
      typedef typename BaseType::BlockMapperType BlockMapperType;

    private:
      typedef typename Traits::ShapeFunctionSetSingletonKeyType ShapeFunctionSetSingletonKeyType;
      typedef SingletonList< ShapeFunctionSetSingletonKeyType, ShapeFunctionSetType, typename Traits::ShapeFunctionSetFactoryType > ShapeFunctionSetProviderType;

      typedef GenericBlockMapperSingletonKey< GridPartType > BlockMapperSingletonKeyType;
      typedef typename Traits::BlockMapperFactoryType BlockMapperFactoryType;
      typedef SingletonList< BlockMapperSingletonKeyType, BlockMapperType, BlockMapperFactoryType > BlockMapperProviderType;

    public:
      using BaseType::order;

      explicit GenericDiscreteFunctionSpace ( GridPartType &gridPart,
                                     const InterfaceType commInterface = defaultInterface,
                                     const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        shapeFunctionSet_( nullptr ),
        blockMapper_( nullptr ),
        mapper_( nullptr )
      {
        // create shape function set
        GeometryType type = GeometryType( Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPartType >::topologyId, GridPartType::dimension );
        shapeFunctionSet_ = &( ShapeFunctionSetProviderType::getObject( type, polOrder ) );

        // create block mapper
        BlockMapperSingletonKeyType key( gridPart );
        blockMapper_ = &( BlockMapperProviderType::getObject( key ) );

        // create mapper
        mapper_ = new MapperType( blockMapper() );
      }

      ~GenericDiscreteFunctionSpace ()
      {
        if( mapper_ )
          delete mapper_;
        mapper_ = nullptr;

        if( blockMapper_ )
          BlockMapperProviderType::removeObject( *blockMapper_ );
        blockMapper_ = nullptr;

        if( shapeFunctionSet_ )
          ShapeFunctionSetProviderType::removeObject( *shapeFunctionSet_ );
        shapeFunctionSet_ = nullptr;
      }
      
      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      DFSpaceIdentifier type () const { return GenericSpace_id; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet */
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity, shapeFunctionSet( entity ) );
      }

      /** \brief return shape function set for given entity
       *
       * \param[in]  entity  entity (of codim 0) for which shape function set 
       *                     is requested
       *
       * \returns  ShapeFunctionSetType  shape function set                     
       */
      ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity ) const
      {
        assert( shapeFunctionSet_);
        return ShapeFunctionSetType( shapeFunctionSet_ );
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous () const { return false; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType &intersection ) const { return false; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const { return polynomialOrder; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::multipleGeometryTypes */
      inline bool multipleGeometryTypes () const { return false; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const
      {
        assert( blockMapper_ );
        return *blockMapper_;
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::mapper */
      DUNE_VERSION_DEPRECATED(1,4,remove)
      MapperType &mapper () const
      {
        assert( mapper_ );
        return *mapper_;
      }      

    private:
      GenericDiscreteFunctionSpace ( const ThisType & );
      ThisType &operator= ( const ThisType & );

      ShapeFunctionSetType *shapeFunctionSet_;
      BlockMapperType *blockMapper_;
      MapperType *mapper_;
    };


  } // namespace Fem

} // end namespace Dune

#endif // #if HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef DUNE_FEM_SPACE_GENERICDISCRETE_HH
