#ifndef DUNE_FEM_STENCIL_HH
#define DUNE_FEM_STENCIL_HH

#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <type_traits>
#include <unordered_map>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/common/utility.hh>

namespace Dune
{
  namespace Fem
  {
    /** \class Stencil
     *  \brief default implementation for a general operator stencil
     *
     *  To assemble a matrix from an operator the method reserve has to be
     *  called on the linear operator class passing a stencil object as
     *  parameter. To setup a full stencil the method fill has to be
     *  called with each pair (en,nb) for which the localMatrix method is
     *  called during the assembly.
     *
     *  \tparam  DomainSpace  type of discrete function space for the domain
     *  \tparam  RangeSpace   type of discrete function space for the range
     *
     */
    template <class DomainSpace, class RangeSpace>
    class Stencil
    {
      // Domain = Row
      typedef typename DomainSpace::IteratorType        DomainIteratorType;
      typedef typename DomainSpace::BlockMapperType     DomainBlockMapper;

      // Range = Column
      typedef typename RangeSpace::IteratorType         RangeIteratorType;
      typedef typename RangeSpace::BlockMapperType      RangeBlockMapper;

    public:
      typedef typename DomainIteratorType::Entity        DomainEntityType;
      typedef typename RangeIteratorType::Entity         RangeEntityType;
      typedef typename DomainBlockMapper::GlobalKeyType  DomainGlobalKeyType;
      typedef typename RangeBlockMapper::GlobalKeyType   RangeGlobalKeyType;

      //! type for storing the stencil of one row
      typedef std::set< DomainGlobalKeyType >       LocalStencilType;

      //! type of std::vector for indexing
      typedef typename std::vector< std::size_t > :: size_type IndexType;

      static const bool indexIsSimple = std::is_convertible< RangeGlobalKeyType, IndexType >::value;
      //static const bool indexIsPOD = false ;//std::is_convertible< RangeGlobalKeyType, IndexType >::value;
      typedef typename std::conditional< indexIsSimple,
                            std::unordered_map< RangeGlobalKeyType, LocalStencilType >,
                            std::map< RangeGlobalKeyType, LocalStencilType > > :: type  GlobalStencilType;

    public:
      /** \brief Constructor
       *
       *  \param[in]  dSpace    domain space
       *  \param[in]  rSpace    range space
       *
       */
      Stencil(const DomainSpace &dSpace, const RangeSpace &rSpace)
        : domainSpace_(dSpace), rangeSpace_(rSpace)
        , domainBlockMapper_( dSpace.blockMapper() )
        , rangeBlockMapper_( rSpace.blockMapper() )
      {
      }

      const DomainSpace &domainSpace() const
      {
        return domainSpace_;
      }
      const RangeSpace &rangeSpace() const
      {
        return rangeSpace_;
      }

      /** \brief Create stencil entries for (dEntity,rEntity) pair
       *
       *  \param[in]  dEntity    domain entity
       *  \param[in]  rEntity    range entity
       *  \param[in]  fillGhost  setup stencil even for a ghost domain entity
       *
       */
      void fill ( const DomainEntityType &dEntity, const RangeEntityType &rEntity,
                  bool fillGhost=true ) const
      {
        if( (dEntity.partitionType() == GhostEntity) && !fillGhost )
          return;

        rangeBlockMapper_.mapEach( rEntity, [ this, &dEntity ] ( int localRow, auto globalRow ) {
            domainBlockMapper_.mapEach( dEntity, RowFillFunctor( globalStencil_[ globalRow ] ) );
          } );
      }

      /** \brief Return stencil for a given row of the matrix
       *
       *  \param[in]  key   key for matrix row
       *
       */
      const LocalStencilType &localStencil(const RangeGlobalKeyType &key) const
      {
        return globalStencil()[ key ];
      }

      /** \brief Return the full stencil
       */
      const GlobalStencilType &globalStencil() const
      {
        return globalStencil_;
      }

      /** \brief Return an upper bound for the maximum number of non-zero entries in all rows
       */
      int maxNonZerosEstimate() const
      {
        int maxNZ = 0;
        for( const auto& entry : globalStencil_ )
        {
          maxNZ = std::max( maxNZ, static_cast<int>( entry.second.size() ) );
        }
        return maxNZ;
      }

      int rows () const { return rangeBlockMapper_.size(); }
      int cols () const { return domainBlockMapper_.size(); }

      /** \brief clear previously computed entries such that a re-compute happens when used again */
      void update()
      {
        globalStencil_.clear();
        // compute stencil based on overloaded implementation
        setupStencil();
      }

      [[deprecated("Use Stencil::update instead()")]]
      void setup() {
        update();
      }

    protected:
      /** \brief method to setup stencil depending on entity set defined in derived class */
      virtual void setupStencil() const = 0;

    private:
      struct RowFillFunctor
      {
        explicit RowFillFunctor ( LocalStencilType &localStencil )
          : localStencil_( localStencil )
        {}

        void operator() ( const std::size_t, const DomainGlobalKeyType &domainGlobal ) const
        {
          localStencil_.insert( domainGlobal );
        }

      private:
        LocalStencilType &localStencil_;
      };

    protected:
      const DomainSpace &domainSpace_;
      const RangeSpace &rangeSpace_;
      const DomainBlockMapper &domainBlockMapper_;
      const RangeBlockMapper &rangeBlockMapper_;

      mutable GlobalStencilType globalStencil_;
    };

    /** \class SimpleStencil
     *  \brief a watered down stencil providing only the upper bound for the non-zero entries per row.
     *
     *  \tparam  DomainSpace  type of discrete function space for the domain
     *  \tparam  RangeSpace   type of discrete function space for the range
     *
     */
    template <class DomainSpace, class RangeSpace>
    class SimpleStencil
    {
      typedef Stencil<DomainSpace,RangeSpace> StencilType;
    public:
      typedef typename StencilType::DomainEntityType       DomainEntityType;
      typedef typename StencilType::RangeEntityType        RangeEntityType;
      typedef typename StencilType::DomainGlobalKeyType    DomainGlobalKeyType;
      typedef typename StencilType::RangeGlobalKeyType     RangeGlobalKeyType;
      typedef typename StencilType::LocalStencilType       LocalStencilType;
      typedef typename StencilType::GlobalStencilType      GlobalStencilType;

      SimpleStencil(int maxNZ)
      : maxNZ_(maxNZ)
      {}
      SimpleStencil()
      {
        maxNZ_ = 1;
      }
      int maxNonZerosEstimate() const
      {
        return maxNZ_;
      }
      const LocalStencilType &localStencil(const DomainGlobalKeyType &key) const
      {
         DUNE_THROW( Dune::NotImplemented, "SimpleStencil: exact stencil information is unavailable." );
         return localStencil_;
      }
      const GlobalStencilType &globalStencil() const
      {
        DUNE_THROW( Dune::NotImplemented, "SimpleStencil: global stencil is  unavailable." );
        return stencil_;
      }
    protected:
      int maxNZ_;
      GlobalStencilType stencil_;
      LocalStencilType localStencil_;
    };


    /** \class StencilWrapper
     *  \brief a simple wrapper class for sparsity patterns provide as vector< set< size_t > >
     *
     *  \tparam  DomainSpace  type of discrete function space for the domain
     *  \tparam  RangeSpace   type of discrete function space for the range
     *
     */
    template <class DomainSpace, class RangeSpace, class LocalStencil>
    class StencilWrapper
    {
      typedef Stencil<DomainSpace,RangeSpace> StencilType;
      typedef StencilWrapper< DomainSpace, RangeSpace, LocalStencil > ThisType;
    public:
      typedef typename StencilType::DomainEntityType       DomainEntityType;
      typedef typename StencilType::RangeEntityType        RangeEntityType;
      typedef typename StencilType::DomainGlobalKeyType    DomainGlobalKeyType;
      typedef typename StencilType::RangeGlobalKeyType     RangeGlobalKeyType;

      static_assert( Std::is_pod< DomainGlobalKeyType > :: value, "StencilWrapper only works for POD DomainGlobalKeyType");

      typedef LocalStencil                                 LocalStencilType;
      typedef std::vector< LocalStencilType >              GlobalStencilType;

      struct Iterator
      {
        typedef std::pair< DomainGlobalKeyType, const LocalStencilType& > ItemType;

        DomainGlobalKeyType index_;
        const GlobalStencilType& stencil_;


        Iterator( const DomainGlobalKeyType& init, const GlobalStencilType& stencil )
          : index_( init ), stencil_( stencil ) {}

        Iterator& operator++ ()
        {
          ++ index_ ;
          return *this;
        }

        ItemType operator*() const
        {
          assert( index_ < stencil_.size() );
          return ItemType( index_, stencil_[ index_ ] );
        }

        bool operator == ( const Iterator& other ) const
        {
          return index_ == other.index_;
        }

        bool operator != ( const Iterator& other ) const
        {
          return !this->operator==( other );
        }
      };

      StencilWrapper(const GlobalStencilType& stencil)
        : stencil_( stencil )
        , maxNZ_( computeMaxNZ() )
      {
      }

      int maxNonZerosEstimate() const
      {
        return maxNZ_;
      }

      const LocalStencilType &localStencil(const DomainGlobalKeyType &key) const
      {
         return stencil_[ key ];
      }

      const ThisType& globalStencil() const
      {
        return *this;
      }

      /** \brief Create stencil entries for (dEntity,rEntity) pair
       *
       *  \param[in]  dEntity    domain entity
       *  \param[in]  rEntity    range entity
       *  \param[in]  fillGhost  setup stencil even for a ghost domain entity
       *
       */
      void fill ( const DomainEntityType &dEntity, const RangeEntityType &rEntity,
                  bool fillGhost=true )
      {
      }

      Iterator begin() const { return Iterator(0, stencil_); }
      Iterator end() const { return Iterator(stencil_.size(), stencil_); }
      Iterator find( const DomainGlobalKeyType &key) const
      {
        assert( key < stencil_.size() );
        return Iterator( key, stencil_);
      }

    protected:
      int computeMaxNZ() const
      {
        int maxNZ = 0;
        for( const auto& row : stencil_ )
        {
          maxNZ = std::max( maxNZ, int(row.size()) );
        }
        return maxNZ;
      }

      const GlobalStencilType& stencil_;
      int maxNZ_;
    };


    /** \class DiagonalStencil
     *  \brief Stencil contaning the entries (en,en) for all entities in the space.
     *         Defailt for an operator over a Lagrange space or a DG mass operator.
     *
     *  \tparam  DomainSpace  type of discrete function space for the domain
     *  \tparam  RangeSpace   type of discrete function space for the range
     *
     */
    template <class DomainSpace, class RangeSpace, class Partition = Dune::Partitions::InteriorBorder>
    struct DiagonalStencil : public Stencil<DomainSpace,RangeSpace>
    {
      typedef Stencil<DomainSpace,RangeSpace> BaseType;
    public:
      typedef Partition                                 PartitionType;
      typedef typename BaseType::DomainEntityType       DomainEntityType;
      typedef typename BaseType::RangeEntityType        RangeEntityType;
      typedef typename BaseType::DomainGlobalKeyType    DomainGlobalKeyType;
      typedef typename BaseType::RangeGlobalKeyType     RangeGlobalKeyType;
      typedef typename BaseType::LocalStencilType       LocalStencilType;
      typedef typename BaseType::GlobalStencilType      GlobalStencilType;

      DiagonalStencil(const DomainSpace &dSpace, const RangeSpace &rSpace)
        : BaseType( dSpace, rSpace )
      {
        setupStencil();
      }

    protected:
      virtual void setupStencil () const override
      {
        const DomainSpace &dSpace = BaseType::domainSpace();
        for (const auto& entity : elements( dSpace.gridPart(), PartitionType{} ) )
          BaseType::fill(entity,entity);
      }
    };

    /** \class DiagonalAndNeighborStencil
     *  \brief Stencil contaning the entries (en,en) and (en,nb) for all entities en in the space
     *         and neighbors nb of en.
     *         Defailt for an operator over a DG space.
     *
     *  \tparam  DomainSpace  type of discrete function space for the domain
     *  \tparam  RangeSpace   type of discrete function space for the range
     *
     */
    template <class DomainSpace, class RangeSpace, class Partition = Dune::Partitions::InteriorBorder>
    struct DiagonalAndNeighborStencil : public Stencil<DomainSpace,RangeSpace>
    {
      typedef Stencil<DomainSpace,RangeSpace> BaseType;
    public:
      typedef Partition                                 PartitionType;
      typedef typename BaseType::DomainEntityType       DomainEntityType;
      typedef typename BaseType::RangeEntityType        RangeEntityType;
      typedef typename BaseType::DomainGlobalKeyType    DomainGlobalKeyType;
      typedef typename BaseType::RangeGlobalKeyType     RangeGlobalKeyType;
      typedef typename BaseType::LocalStencilType       LocalStencilType;
      typedef typename BaseType::GlobalStencilType      GlobalStencilType;

      DiagonalAndNeighborStencil(const DomainSpace &dSpace, const RangeSpace &rSpace,
                                 bool onlyNonContinuousNeighbors = false)
        : BaseType( dSpace, rSpace ),
          onlyNonContinuousNeighbors_(onlyNonContinuousNeighbors)
      {
        setupStencil();
      }

    protected:
      virtual void setupStencil() const override
      {
        const DomainSpace &dSpace = BaseType::domainSpace();
        const RangeSpace  &rSpace = BaseType::rangeSpace();
        for (const auto & entity: elements( dSpace.gridPart(), PartitionType{} ) )
        {
          BaseType::fill(entity,entity);
          for (const auto & intersection: intersections(dSpace.gridPart(), entity) )
          {
            if ( onlyNonContinuousNeighbors_
                && rSpace.continuous(intersection) && dSpace.continuous(intersection) )
              continue;
            if( intersection.neighbor() )
            {
              auto neighbor = intersection.outside();
              BaseType::fill(neighbor,entity);
            }
          }
        }
      }

      bool onlyNonContinuousNeighbors_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #if defined DUNE_FEM_STENCIL_HH

