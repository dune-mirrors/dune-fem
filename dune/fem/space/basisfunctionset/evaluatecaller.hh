#ifndef DUNE_FEM_EVALUATECALLER_HH
#define DUNE_FEM_EVALUATECALLER_HH

#include <cstdlib>
#include <iostream>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/fem/misc/threads/threadsafevalue.hh>

#ifdef USE_BASEFUNCTIONSET_CODEGEN
#ifndef DUNE_FEM_INCLUDE_AUTOGENERATEDCODE_FILENAME_SPEC
// default filename for autogenerated code is simply autogeneratedcode.hh
// this filename is overloaded by the codegeneration for the python modules
#include <autogeneratedcode.hh>
#endif

#ifdef DUNE_FEM_INCLUDE_AUTOGENERATEDCODE_FILENAME_SPEC
#define CODEGEN_INCLUDEMAXNUMS
// include max number definitions
#include DUNE_FEM_INCLUDE_AUTOGENERATEDCODE_FILENAME_SPEC
#undef CODEGEN_INCLUDEMAXNUMS
#endif
#endif

////////////////////////////////////////////
//
// pre-define these values for faster compilation
//
////////////////////////////////////////////
#ifndef MAX_NUMBER_OF_QUAD_POINTS
#define MAX_NUMBER_OF_QUAD_POINTS 20
#endif

#ifndef MAX_NUMBER_OF_BASE_FCT
#define MAX_NUMBER_OF_BASE_FCT 20
#endif

#ifndef MIN_NUMBER_OF_QUAD_POINTS
#define MIN_NUMBER_OF_QUAD_POINTS 1
#endif

#ifndef MIN_NUMBER_OF_BASE_FCT
#define MIN_NUMBER_OF_BASE_FCT 1
#endif

namespace Dune
{

  namespace Fem
  {

    // empty class for specialization of evaluation classes in basefunctionsets.hh
    class EmptyGeometry {};

    // forward declaration
    template <class Traits,
              int quadNop,
              int numBaseFct >
    class EvaluateCaller;

    template< class QuadratureImp,
              class FactorImp,
              class LocalDofVectorImp,
              class GeometryImp = EmptyGeometry >
    struct EvaluateCallerInterfaceTraits
    {
      typedef QuadratureImp     QuadratureType;
      typedef FactorImp         FactorType;
      typedef LocalDofVectorImp LocalDofVectorType;
      typedef GeometryImp       Geometry;
    };

    template <class Traits,
              class BaseFunctionSet,
              class RangeVectorImp>
    struct EvaluateCallerTraits
    {
      typedef Traits  BaseTraits;
      typedef typename Traits :: QuadratureType      QuadratureType ;
      typedef typename Traits :: FactorType          FactorType ;
      typedef typename Traits :: LocalDofVectorType  LocalDofVectorType ;
      typedef typename Traits :: Geometry            Geometry ;

      typedef BaseFunctionSet   BaseFunctionSetType;
      typedef RangeVectorImp    RangeVectorType;
    };


    //- base function evaluation interface
    template <class Traits>
    class EvaluateCallerInterface
    {
      typedef EvaluateCallerInterface< Traits >  ThisType;
    public:
      typedef std::unique_ptr< ThisType > StoragePointerType;
      typedef std::pair< bool, StoragePointerType > StorageItemType;

    protected:
      enum { maxNumBaseFunctions = MAX_NUMBER_OF_BASE_FCT };
      enum { maxQuadratures = 50 };
      enum { maxQuadNop = MAX_NUMBER_OF_QUAD_POINTS };

      class EvaluatorStorage
      {
      protected:
        //typedef std::pair< std::unique_ptr< ThisType >, bool > StorageItemType;
        std::vector< StorageItemType > storage_;
      public:
        EvaluatorStorage() :
          storage_( maxQuadratures )
        {
          for( auto& item : storage_ )
          {
            item.first = false ;
            item.second.reset();
          }

        }

        StorageItemType& operator [] ( const int i ) { return storage_[ i ]; }
        const StorageItemType& operator [] ( const int i ) const { return storage_[ i ]; }
      };


      EvaluateCallerInterface() {}

    public:
      typedef typename Traits :: QuadratureType      QuadratureType ;
      typedef typename Traits :: FactorType          FactorType ;
      typedef typename Traits :: LocalDofVectorType  LocalDofVectorType ;
      typedef typename Traits :: Geometry            Geometry ;

      virtual ~EvaluateCallerInterface() {}

      virtual void* storageAddress () const = 0;
      virtual size_t storageSize () const = 0;

      virtual void axpyRanges( const QuadratureType&,
                               const FactorType& ,
                               LocalDofVectorType & ) const = 0;

      virtual void evaluateRanges( const QuadratureType& quad,
                                   const LocalDofVectorType & dofs,
                                   FactorType& factors) const = 0;

      virtual void axpyJacobians( const QuadratureType&,
                                  const Geometry&,
                                  const FactorType& ,
                                  LocalDofVectorType & ) const = 0;

      virtual void evaluateJacobians( const QuadratureType&,
                                      const Geometry&,
                                      const LocalDofVectorType&,
                                      FactorType&) const = 0;

      template < class BaseFunctionSet, class Storage >
      static const StoragePointerType& storage(const BaseFunctionSet& baseSet,
                                               const Storage& dataCache,
                                               const QuadratureType& quad)
      {
        // assert that max numbers are big enough
        assert( baseSet.numDifferentBaseFunctions() <= maxNumBaseFunctions );
        assert( quad.nop()  <= maxQuadNop );
        assert( quad.id()   < maxQuadratures );

        // static vector holding all evaluator instances
        static ThreadSafeValue< EvaluatorStorage > evaluatorStorage;
        EvaluatorStorage& evaluators = *evaluatorStorage;

        // check if object already created
        const size_t quadId = quad.id();
        if( ! evaluators[ quadId ].first )
        {
          typedef EvaluateCallerTraits< Traits, BaseFunctionSet, Storage> NewTraits;
          auto& item = evaluators[ quadId ];
          // create appropriate evaluator
          item.second.reset(
              EvaluateCaller< NewTraits, maxQuadNop, maxNumBaseFunctions >
                 :: create( dataCache , quad.nop(), baseSet.numDifferentBaseFunctions() ) );
          // if pointer was checked, set flag to true, pointer may still be null
          item.first = true;
        }

        // make sure the storage is the same
        //assert( dataCache.size() == evaluators[ quadId ]->storageSize() );
        // return pointer to evaluator, if null then a fallback to default impl happens
        return evaluators[ quadId ].second;
      }
    };

    template <class Traits,
              int quadNop,
              int numBaseFct >
    class EvaluateRealImplementation
      : public EvaluateCallerInterface< typename Traits :: BaseTraits >
    {
    protected:
      typedef typename Traits :: BaseFunctionSetType BaseFunctionSetType;
      typedef typename Traits :: QuadratureType      QuadratureType ;
      typedef typename Traits :: FactorType          FactorType ;
      typedef typename Traits :: LocalDofVectorType  LocalDofVectorType ;
      typedef typename Traits :: Geometry            Geometry ;
      typedef typename Traits :: RangeVectorType     RangeVectorType ;
      typedef typename RangeVectorType :: value_type :: field_type FieldType;

      static const int dimRange = BaseFunctionSetType :: dimRange;
      typedef EvaluateRealImplementation< Traits, quadNop, numBaseFct > ThisType;
      typedef EvaluateCallerInterface< typename Traits :: BaseTraits >   BaseType;

      const RangeVectorType& rangeStorage_;
      std::vector< std::vector< FieldType > > rangeStorageTransposed_;
      std::vector< std::vector< FieldType > > rangeStorageFlat_;
      mutable std::vector< std::vector< std::vector< FieldType > > > rangeStorageTwisted_;

      template <class K>
      int getDim( const DenseVector< K >& vec) const
      {
        return vec.size();
      }

      template <class K>
      int getDim( const DenseMatrix< K >& mat) const
      {
        // dimRange == rows which is 1 for basis storage
        return getDim( mat[ 0 ] );
      }

      // initialize storage for ranges (i.e. scalar)
      void initRangeStorageTransposed( const std::integral_constant< bool, true > )
      {
        assert( rangeStorage_[ 0 ].size() == 1 );
        {
          const int quadPoints = rangeStorage_.size() / numBaseFct;
          const int faces = quadPoints / quadNop;
          rangeStorageTransposed_.resize( faces );
          for( int f=0; f<faces; ++f )
          {
            auto& rangeStorageTransposed = rangeStorageTransposed_[ f ];
            // rearrange such that we store for one basis functions all
            // evaluations for all quadrature points, i.e. numBaseFct * quadNop
            rangeStorageTransposed.resize( numBaseFct * quadNop );
            for( int i=0; i<numBaseFct; ++i )
            {
              const int idx  = i * quadNop;
              for( int j=0; j<quadNop; ++j )
              {
                int qp = f * quadNop + j ;
                assert( j*numBaseFct + i < int(rangeStorage_.size()) );
                rangeStorageTransposed[ idx + j ] = rangeStorage_[ qp*numBaseFct + i ][ 0 ];
              }
            }
          }
        }
      }

      // initialize storage for jacobians (i.e. vectors)
      void initRangeStorageTransposed( const std::integral_constant< bool, false > )
      {
        const int dim = rangeStorage_[ 0 ][ 0 ].size();
        {
          const int quadPoints = rangeStorage_.size() / numBaseFct;
          const int faces = quadPoints / quadNop;
          rangeStorageTransposed_.resize( faces );
          rangeStorageFlat_.resize( faces );
          for( int f=0; f<faces; ++f )
          {
            auto& rangeStorageTransposed = rangeStorageTransposed_[ f ];
            auto& rangeStorageFlat = rangeStorageFlat_[ f ];

            // rearrange such that we store for one basis functions all
            // evaluations for all quadrature points, i.e. numBaseFct * quadNop
            rangeStorageTransposed.resize( numBaseFct * quadNop * dim );
            rangeStorageFlat.resize( numBaseFct * quadNop * dim );
            for( int i=0; i<numBaseFct; ++i )
            {
              const int idx  = i * (quadNop * dim);
              for( int j=0; j<quadNop; ++j )
              {
                int qp = f * quadNop + j ;
                for( int d=0; d<dim; ++d )
                {
                  rangeStorageFlat[ j*numBaseFct*dim + (i * dim) + d ] = rangeStorage_[ qp*numBaseFct + i ][ 0 ][ d ];
                  rangeStorageTransposed[ idx + (j * dim) + d ] = rangeStorage_[ qp*numBaseFct + i ][ 0 ][ d ];
                }
              }
            }
          }
        }
      }

      template <class Quadrature>
      const std::vector< FieldType >&
      getTwistedStorage( const Quadrature& quad ) const
      {
        // for evaluation the range storage dimension should be 1 and therefore
        // rangeStorageTransposed should have been filled
        assert( ! rangeStorageTransposed_.empty() );

        // if we are in the ranges cases then basis can be stored transposed
        // quadrature points is the outer loop
        if( quad.twisted() )
        {
          auto& rangeStorageTwisted = rangeStorageTwisted_[ quad.twistId() ];
          if( rangeStorageTwisted.empty() )
          {
            // either 1 or dim of grid
            const int dim = getDim( rangeStorage_[ 0 ] );

            const int quadPoints = rangeStorage_.size() / numBaseFct;
            const int faces = quadPoints / quadNop;
            rangeStorageTwisted.resize( faces );
            for( int f=0; f<faces; ++f )
            {
              auto& rangeStorageFace = rangeStorageTwisted[ f ];
              const auto& rangeStorageTransposed = rangeStorageTransposed_[ f ];

              // rearrange such that we store for one basis functions all
              // evaluations for all quadrature points including the twisted mapping
              rangeStorageFace.resize( rangeStorageTransposed.size() );
              for( int i=0; i<numBaseFct; ++i )
              {
                const int idx = i * quadNop;
                for( int j=0; j<quadNop; ++j )
                {
                  const int qp = quad.localCachingPoint( j );
                  for( int d=0; d<dim; ++d )
                  {
                    rangeStorageFace[ idx + (j * dim) + d ] = rangeStorageTransposed[ idx + (qp * dim) + d ];
                  }
                }
              }
            }
          } // end if( rangeStorageTwisted.empty() )
          return rangeStorageTwisted[ quad.localFaceIndex() ];
        }
        else // no twist (i.e. twist = 0 and twistId == 5 (-4 is mapped to 0))
        {
          return rangeStorageTransposed_[ quad.localFaceIndex() ];
        }
      }
    public:
      // type of interface class
      typedef BaseType InterfaceType;

      EvaluateRealImplementation( const RangeVectorType& rangeStorage )
        : rangeStorage_( rangeStorage ), rangeStorageTwisted_( 8 ) // 8 different twists
      {
        initRangeStorageTransposed( std::integral_constant< bool,
                                    std::is_same< typename RangeVectorType::value_type,
                                                  Dune::FieldVector< double, 1 > > :: value > () );
      }

      virtual void* storageAddress() const { return (void *) &rangeStorage_ ; }
      virtual size_t storageSize() const { return rangeStorage_.size() ; }

      virtual void axpyRanges( const QuadratureType& quad,
                               const FactorType& rangeFactors,
                               LocalDofVectorType & dofs ) const
      {
        BaseFunctionSetType :: template AxpyRanges
          < BaseFunctionSetType, Geometry, dimRange, quadNop, numBaseFct > :: axpy
          ( quad, rangeStorage_, rangeFactors, dofs );
      }

      virtual void evaluateRanges( const QuadratureType& quad,
                                   const LocalDofVectorType & dofs,
                                   FactorType& rangeFactors) const
      {
        BaseFunctionSetType :: template EvaluateRanges
          < BaseFunctionSetType, Geometry, dimRange, quadNop, numBaseFct >
          :: eval ( quad, getTwistedStorage( quad ), dofs, rangeFactors );
      }

      virtual void axpyJacobians( const QuadratureType& quad,
                                  const Geometry& geometry,
                                  const FactorType& jacFactors,
                                  LocalDofVectorType& dofs) const
      {
        BaseFunctionSetType :: template AxpyJacobians
          < BaseFunctionSetType, Geometry, dimRange, quadNop, numBaseFct > :: axpy
          ( quad, geometry, rangeStorageFlat_[ quad.localFaceIndex() ], jacFactors, dofs );
      }

      virtual void evaluateJacobians( const QuadratureType& quad,
                                      const Geometry& geometry,
                                      const LocalDofVectorType& dofs,
                                      FactorType& jacFactors) const
      {
        BaseFunctionSetType :: template EvaluateJacobians
          < BaseFunctionSetType, Geometry, dimRange, quadNop, numBaseFct > :: eval
          ( quad, geometry, getTwistedStorage( quad ), dofs, jacFactors );
      }

      static InterfaceType* create( const RangeVectorType& rangeStorage )
      {
        return new ThisType( rangeStorage );
      }
    };

    // The default EvaluateImplementation is empty
    // to create this has to be specified and derived from EvaluateCallerDefault
    template <class Traits,
              int quadNop,
              int numBaseFct >
    class EvaluateImplementation
      : public EvaluateCallerInterface< typename Traits :: BaseTraits >
    {
    protected:
      typedef typename Traits :: BaseFunctionSetType BaseFunctionSetType;
      typedef typename Traits :: QuadratureType      QuadratureType ;
      typedef typename Traits :: FactorType          FactorType ;
      typedef typename Traits :: LocalDofVectorType  LocalDofVectorType ;
      typedef typename Traits :: Geometry            Geometry ;
      typedef typename Traits :: RangeVectorType     RangeVectorType ;

      typedef EvaluateImplementation< Traits, quadNop, numBaseFct > ThisType;

      typedef EvaluateCallerInterface< typename Traits :: BaseTraits >   BaseType;
      static const int dimRange = BaseFunctionSetType :: dimRange;
    public:
      // type of interface class
      typedef BaseType InterfaceType;

      EvaluateImplementation( const RangeVectorType& rangeStorage )
      {}

      virtual void axpyRanges( const QuadratureType& quad,
                               const FactorType& rangeFactors,
                               LocalDofVectorType & dofs ) const
      {
        std::cerr << "ERROR: EvaluateImplementation::axpyRanges not overloaded!" << std::endl;
        std::abort();
      }

      virtual void axpyJacobians( const QuadratureType& quad,
                                  const Geometry& geometry,
                                  const FactorType& jacFactors,
                                  LocalDofVectorType& dofs) const
      {
        std::cerr << "ERROR: EvaluateImplementation::axpyJacobians not overloaded!" << std::endl;
        std::abort();
      }

      virtual void evaluateRanges( const QuadratureType& quad,
                                   const LocalDofVectorType & dofs,
                                   FactorType& rangeFactors) const
      {
        std::cerr << "ERROR: EvaluateImplementation::evaluateRanges not overloaded!" << std::endl;
        std::abort();
      }

      virtual void evaluateJacobians( const QuadratureType& quad,
                                      const Geometry& geometry,
                                      const LocalDofVectorType& dofs,
                                      FactorType& jacFactors) const
      {
        std::cerr << "ERROR: EvaluateImplementation::evaluateJacobians not overloaded!" << std::endl;
        std::abort();
      }

      static InterfaceType* create( const RangeVectorType& )
      {
        std::cout << "Optimized EvaluateImplementation for < dimR="<<dimRange<< ", qp=" << quadNop << ", bases=" << numBaseFct << " > not created, falling back to default!" << std::endl;
        //DUNE_THROW(NotImplemented,"EvaluateImplementation for < " << quadNop << " , " << numBaseFct << " > not created!");
        //return (InterfaceType*) 0;
        return nullptr;
      }
    };

    template <class Traits,
              int quadNop,
              int numBaseFct >
    class EvaluateCaller
    {
    protected:
      typedef typename Traits :: RangeVectorType     RangeVectorType ;
      typedef EvaluateCallerInterface< typename Traits :: BaseTraits >  InterfaceType;
    public:
      static InterfaceType* createObj( const RangeVectorType& rangeStorage,
                                       const size_t numbase )
      {
        if( numBaseFct == numbase )
          return EvaluateImplementation< Traits, quadNop, numBaseFct > :: create( rangeStorage );
        else
          return EvaluateCaller< Traits, quadNop, numBaseFct - 1 > :: createObj( rangeStorage, numbase );
      }

      static InterfaceType* create( const RangeVectorType& rangeStorage,
                                    const size_t quadnop, const size_t numbase )
      {
        if( quadNop == quadnop )
          return EvaluateCaller< Traits, quadNop, numBaseFct > :: createObj( rangeStorage, numbase );
        else
          return EvaluateCaller< Traits, quadNop - 1, numBaseFct > :: create( rangeStorage, quadnop, numbase );
      }
    };

    template <class Traits,
              int numBaseFct >
    class EvaluateCaller< Traits, MIN_NUMBER_OF_QUAD_POINTS, numBaseFct >
    {
    protected:
      enum { quadNop = MIN_NUMBER_OF_QUAD_POINTS };
      typedef typename Traits :: RangeVectorType     RangeVectorType ;
      typedef EvaluateCallerInterface< typename Traits :: BaseTraits >  InterfaceType;
    public:
      static InterfaceType* createObj( const RangeVectorType& rangeStorage,
                                       const size_t numbase )
      {
        if( numBaseFct == numbase )
          return EvaluateImplementation< Traits, quadNop, numBaseFct > :: create( rangeStorage );
        else
          return EvaluateCaller< Traits, quadNop, numBaseFct - 1 > :: createObj( rangeStorage, numbase );
      }

      static InterfaceType* create( const RangeVectorType& rangeStorage,
                                    const size_t quadnop, const size_t numbase )
      {
        if( quadNop == quadnop )
          return EvaluateCaller< Traits, quadNop, numBaseFct > :: createObj( rangeStorage, numbase );
        else
        {
          std::cerr << "ERROR: EvaluateCaller< "<< quadNop << ", " << numBaseFct << " >::createObj: no working combination!" << std::endl;
          std::abort();
        }
      }
    };

    template <class Traits,
              int quadNop>
    class EvaluateCaller< Traits, quadNop, MIN_NUMBER_OF_BASE_FCT >
    {
    protected:
      enum { numBaseFct = MIN_NUMBER_OF_BASE_FCT };
      typedef typename Traits :: RangeVectorType     RangeVectorType ;
      typedef EvaluateCallerInterface< typename Traits :: BaseTraits >  InterfaceType;
    public:
      static InterfaceType* createObj( const RangeVectorType& rangeStorage,
                                       const size_t numbase )
      {
        if( numBaseFct == numbase )
          return EvaluateImplementation< Traits, quadNop, numBaseFct > :: create( rangeStorage );
        else
        {
          std::cerr << "ERROR: EvaluateCaller< "<< quadNop << ", " << numBaseFct << " >::createObj: no working combination!" << std::endl;
          std::abort();
        }
      }

      static InterfaceType* create( const RangeVectorType& rangeStorage,
                                    const size_t quadnop, const size_t numbase )
      {
        if( quadNop == quadnop )
          return EvaluateCaller< Traits, quadNop, numBaseFct > :: createObj( rangeStorage, numbase );
        else
        {
          return EvaluateCaller< Traits, quadNop - 1, numBaseFct > :: create( rangeStorage, quadnop, numbase );
        }
      }
    };

    template <class Traits>
    class EvaluateCaller< Traits, MIN_NUMBER_OF_QUAD_POINTS, MIN_NUMBER_OF_BASE_FCT>
    {
    protected:
      enum { quadNop = MIN_NUMBER_OF_QUAD_POINTS };
      enum { numBaseFct = MIN_NUMBER_OF_BASE_FCT };
      typedef typename Traits :: RangeVectorType     RangeVectorType ;
      typedef EvaluateCallerInterface< typename Traits :: BaseTraits >  InterfaceType;
    public:
      static InterfaceType* createObj( const RangeVectorType& rangeStorage,
                                       const size_t numbase )
      {
        if( numBaseFct == numbase )
          return EvaluateImplementation< Traits, quadNop, numBaseFct > :: create( rangeStorage );
        else
        {
          std::cerr << "ERROR: EvaluateCaller< "<< quadNop << ", " << numBaseFct << " >::createObj: no working combination!" << std::endl;
          std::abort();
        }
      }

      static InterfaceType* create( const RangeVectorType& rangeStorage,
                                    const size_t quadnop, const size_t numbase )
      {
        if( quadNop == quadnop )
          return EvaluateCaller< Traits, quadNop, numBaseFct > :: createObj( rangeStorage, numbase );
        else
        {
          std::cerr << "ERROR: EvaluateCaller< "<< quadNop << ", " << numBaseFct << " >::create: no working combination!" << std::endl;
          std::abort();
        }
      }
    };

#ifdef DUNE_FEM_INCLUDE_AUTOGENERATEDCODE_FILENAME_SPEC
#define CODEGEN_INCLUDEEVALCALLERS
// include specializations of EvaluateImplementation
#include DUNE_FEM_INCLUDE_AUTOGENERATEDCODE_FILENAME_SPEC
#undef CODEGEN_INCLUDEEVALCALLERS
#endif

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_EVALUATECALLER_HH
