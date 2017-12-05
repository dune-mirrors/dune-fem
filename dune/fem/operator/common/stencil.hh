#ifndef DUNE_FEM_STENCIL_HH
#define DUNE_FEM_STENCIL_HH

#include <algorithm>
#include <map>
#include <numeric>
#include <set>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/fem/misc/functor.hh>

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
      typedef typename DomainIteratorType::Entity       DomainEntityType;
      typedef typename RangeIteratorType::Entity        RangeEntityType;
      typedef typename DomainBlockMapper::GlobalKeyType DomainGlobalKeyType;
      typedef typename RangeBlockMapper::GlobalKeyType  RangeGlobalKeyType;

      //! type for storing the stencil of one row
      typedef std::set<RangeGlobalKeyType>                   LocalStencilType;
      //! type for storing the full stencil
      typedef std::map<DomainGlobalKeyType,LocalStencilType> GlobalStencilType;

    public:
      /** \brief Constructor
       *
       *  \param[in]  dSpace    domain space
       *  \param[in]  rSpace    range space
       *
       */
      Stencil(const DomainSpace &dSpace, const RangeSpace &rSpace)
        : domainBlockMapper_( dSpace.blockMapper() )
        , rangeBlockMapper_( rSpace.blockMapper() )
      {
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
      const LocalStencilType &localStencil(const DomainGlobalKeyType &key) const
      {
        return globalStencil_[key];
      }
      /** \brief Return the full stencil
       */
      const GlobalStencilType &globalStencil() const
      {
        return globalStencil_;
      }
      /** \brief Return an upper bound for the maximum number of non-zero entries in all row
       */
      int maxNonZerosEstimate() const
      {
        return std::accumulate( globalStencil().begin(), globalStencil().end(), 0,
         []( int ret, const auto& entry ){ return std::max( ret, static_cast<int>( entry.second.size() ) ); } );
      }

      int rows () const { return rangeBlockMapper_.size(); }
      int cols () const { return rangeBlockMapper_.size(); }

    private:
      struct RowFillFunctor
      {
        explicit RowFillFunctor ( LocalStencilType &localStencil )
          : localStencil_( localStencil )
        {}

        void operator() ( const std::size_t, const RangeGlobalKeyType &rangeGlobal ) const
        {
          localStencil_.insert( rangeGlobal );
        }

      private:
        LocalStencilType &localStencil_;
      };

      const DomainBlockMapper &domainBlockMapper_;
      const RangeBlockMapper &rangeBlockMapper_;
      GlobalStencilType globalStencil_;
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
    private:
      int maxNZ_;
      GlobalStencilType stencil_;
      LocalStencilType localStencil_;
    };

    /** \class DiagonalStencil
     *  \brief Stencil contaning the entries (en,en) for all entities in the space.
     *         Defailt for an operator over a Lagrange space or a DG mass operator.
     *
     *  \tparam  DomainSpace  type of discrete function space for the domain
     *  \tparam  RangeSpace   type of discrete function space for the range
     *
     */
    template <class DomainSpace, class RangeSpace>
    struct DiagonalStencil : public Stencil<DomainSpace,RangeSpace>
    {
      typedef Stencil<DomainSpace,RangeSpace> BaseType;
    public:
      typedef typename BaseType::DomainEntityType       DomainEntityType;
      typedef typename BaseType::RangeEntityType        RangeEntityType;
      typedef typename BaseType::DomainGlobalKeyType    DomainGlobalKeyType;
      typedef typename BaseType::RangeGlobalKeyType     RangeGlobalKeyType;
      typedef typename BaseType::LocalStencilType       LocalStencilType;
      typedef typename BaseType::GlobalStencilType      GlobalStencilType;

      DiagonalStencil(const DomainSpace &dSpace, const RangeSpace &rSpace)
        : BaseType( dSpace, rSpace )
      {
        for (const auto& entity : elements( dSpace.gridPart(), Dune::Partitions::interiorBorder ) )
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
    template <class DomainSpace, class RangeSpace>
    struct DiagonalAndNeighborStencil : public Stencil<DomainSpace,RangeSpace>
    {
      typedef Stencil<DomainSpace,RangeSpace> BaseType;
    public:
      typedef typename BaseType::DomainEntityType       DomainEntityType;
      typedef typename BaseType::RangeEntityType        RangeEntityType;
      typedef typename BaseType::DomainGlobalKeyType    DomainGlobalKeyType;
      typedef typename BaseType::RangeGlobalKeyType     RangeGlobalKeyType;
      typedef typename BaseType::LocalStencilType       LocalStencilType;
      typedef typename BaseType::GlobalStencilType      GlobalStencilType;

      DiagonalAndNeighborStencil(const DomainSpace &dSpace, const RangeSpace &rSpace,
                                 bool onlyNonContinuousNeighbors = false)
        : BaseType( dSpace, rSpace )
      {
        for (const auto & entity: elements( dSpace.gridPart(), Dune::Partitions::interiorBorder ) )
        {
          BaseType::fill(entity,entity);
          for (const auto & intersection: intersections(dSpace.gridPart(), entity) )
          {
            if ( onlyNonContinuousNeighbors
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
    };

  } // namespace Fem

} // namespace Dune

#endif // #if defined DUNE_FEM_STENCIL_HH

