// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_STENCIL_HH
#define DUNE_FEM_STENCIL_HH

#include <iostream>
#include <set>
#include <map>

namespace Dune 
{
  namespace Fem 
  {
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

      typedef std::set<RangeGlobalKeyType>                   LocalStencilType;
      typedef std::map<DomainGlobalKeyType,LocalStencilType> GlobalStencilType;

    public:
      Stencil(const DomainSpace &dSpace, const RangeSpace &rSpace)
        : domainBlockMapper_( dSpace.blockMapper() )
        , rangeBlockMapper_( rSpace.blockMapper() )
      {
      }

      //! create entries for element and neighbors
      void fill ( const DomainEntityType &dEntity, const RangeEntityType &rEntity )
      {
        domainBlockMapper_.mapEach(dEntity, 
                  MFunctor( rangeBlockMapper_, rEntity, FillFunctor(globalStencil_) ) );
      }

      const LocalStencilType &localStencil(const DomainGlobalKeyType &key) const
      { 
        return globalStencil_[key]; 
      }
      const GlobalStencilType &globalStencil() const
      { 
        return globalStencil_; 
      }
      int maxZerosEstimate() const
      {
        int ret = 0;
        typedef typename GlobalStencilType::const_iterator StencilIteratorType;
        const GlobalStencilType &glStencil = globalStencil();
        StencilIteratorType end = glStencil.end();
        for ( StencilIteratorType it = glStencil.begin(); it != end; ++it)
          ret = std::max(ret,(int)it->second.size());
        return ret;
      }

    private:

      struct FillFunctor
      {
        typedef DomainGlobalKeyType GlobalKey;
        FillFunctor(GlobalStencilType &stencil) 
        : stencil_(stencil),
          localStencil_(0)
        {}
        void set(const std::size_t, const DomainGlobalKeyType &domainGlobal)
        {
          localStencil_ = &(stencil_[ domainGlobal ]);
        }
        void operator() ( const std::size_t, const RangeGlobalKeyType &rangeGlobal)
        {
          localStencil_->insert( rangeGlobal );
        }
        private:
        GlobalStencilType &stencil_;
        LocalStencilType *localStencil_;
      };
      typedef typename Dune::Fem::MatrixFunctor<RangeBlockMapper,RangeEntityType,FillFunctor > MFunctor;

      const DomainBlockMapper &domainBlockMapper_;
      const RangeBlockMapper &rangeBlockMapper_;
      GlobalStencilType globalStencil_;
    };

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
        maxNZ_ = 0; // estimate optimal stencil....
      }
      int maxZerosEstimate() const
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

    template <class DomainSpace, class RangeSpace>
    struct DiagonalStencil : Stencil<DomainSpace,RangeSpace>
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
        typedef typename DomainSpace::IteratorType IteratorType;
        const IteratorType end = dSpace.end();
        for (IteratorType it = dSpace.begin(); it!=end; ++it)
        {
          const DomainEntityType &entity = *it;
          BaseType::fill(entity,entity);
        }
      }
    };

    template <class DomainSpace, class RangeSpace>
    struct DiagonalAndNeighborStencil : Stencil<DomainSpace,RangeSpace>
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
        typedef typename DomainSpace::IteratorType IteratorType;
        const IteratorType end = dSpace.end();
        for (IteratorType it = dSpace.begin(); it!=end; ++it)
        {
          const DomainEntityType &entity = *it;
          BaseType::fill(entity,entity);
          typedef typename DomainSpace::GridPartType GridPart;
          typedef typename GridPart :: IntersectionIteratorType IntersectionIteratorType;
          typedef typename IntersectionIteratorType :: Intersection IntersectionType;
          typedef typename GridPart :: template Codim<0> :: EntityPointerType EntityPointer;

          const IntersectionIteratorType endit = dSpace.gridPart().iend( entity );
          for( IntersectionIteratorType it = dSpace.gridPart().ibegin( entity );
               it != endit; ++it )
          {
            const IntersectionType& intersection = *it;
            if ( onlyNonContinuousNeighbors 
                && rSpace.continuous(intersection) && dSpace.continuous(intersection) )
              continue;
            if( intersection.neighbor() ) 
            {
              EntityPointer ep = intersection.outside(); 
              const DomainEntityType& neighbor = *ep;
              BaseType::fill(neighbor,entity);
            }
          }
        }
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #if defined DUNE_FEM_STENCIL_HH

