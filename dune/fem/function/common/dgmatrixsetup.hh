#ifndef DUNE_DGMATRIXSETUP_HH
#define DUNE_DGMATRIXSETUP_HH

#include <dune/common/timer.hh>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/misc/functor.hh>

#if 0 // HAVE_DUNE_ISTL
#include <dune/istl/operators.hh>
#include <dune/fem/operator/matrix/istlmatrix.hh>
#include <dune/fem/operator/matrix/preconditionerwrapper.hh>
#endif

namespace Dune 
{
  namespace Fem
  {
    ////////////////////////////////////////////////////////////
    //
    //  Setup of matrix structure 
    //
    ////////////////////////////////////////////////////////////
    /**  \brief Setup Matrix structure for DG operators by including
     * elements and it's neighbors. 
    */
    template <class SpaceImp >
    class ElementAndNeighbors
    {
    public:
      //! create entries for element and neighbors 
      template < class RowMapperType, class DiscreteFunctionType>
      static inline void setup(const SpaceImp& space,    
                               const RowMapperType& rowMapper,
                               const DiscreteFunctionType* )
      {}
    };
    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage>
    class ElementAndNeighbors<
      DiscontinuousGalerkinSpace<FunctionSpace,GridPart,polOrder,Storage> >
    {
      typedef DiscontinuousGalerkinSpace<FunctionSpace,GridPart,polOrder,Storage> SpaceImp;
    public:
      //! create entries for element and neighbors 
      template < class RowMapperType, class DiscreteFunctionType>
      static inline void setup(const SpaceImp& space,    
                               const RowMapperType& rowMapper,
                               const DiscreteFunctionType* )
      {
        typedef typename SpaceImp :: GridPartType GridPartType;
        const GridPartType& gridPart = space.gridPart();

        typedef Fem::ParallelScalarProduct<DiscreteFunctionType> ParallelScalarProductType;
        typedef typename ParallelScalarProductType :: BuildProxyType BuildProxyType;
        
        ParallelScalarProductType scp (space);

        std::auto_ptr<BuildProxyType> buildProxy = scp.buildProxy();

        // define used types 
        typedef typename GridPartType :: GridType GridType;
        typedef typename GridPartType :: template Codim<0> :: EntityType    EntityType;
        typedef typename GridPartType :: template Codim<0> :: IteratorType  IteratorType;

        // we need All_Partition here to insert overlap entities 
        // only for diagonal 
        IteratorType endit = gridPart.template end<0>(); 
        for(IteratorType it = gridPart.template begin<0>(); it != endit; ++it)
        {
          const EntityType & en = *it;
          // add all column entities to row  
          fill( gridPart,en,rowMapper,*buildProxy );
        }

        // insert size as last ghost 
        buildProxy->insert( rowMapper.size() );
      }

    protected:
      template <class GK, class SlaveDofProviderImp>
      struct FillSlaveFunctor
      {
        typedef SlaveDofProviderImp SlaveDofProvider;
        typedef GK GlobalKey;

        FillSlaveFunctor ( SlaveDofProvider &slaves ) 
        : slaves_( slaves ) 
        {}

        template < class ColGlobal >
        void operator() ( const std::size_t colLocal, const ColGlobal &colGlobal )
        {
          slaves_.insert( colGlobal );
        }

      private:
        SlaveDofProvider &slaves_;
      };

      //! create entries for element and neighbors 
      template <class GridPartImp,
                class EntityImp,
                class RowMapperImp,
                class ParallelScalarProductType>
      static inline void fill(const GridPartImp& gridPart,
                       const EntityImp& en,
                       const RowMapperImp& rowMapper,
                       ParallelScalarProductType& slaveDofs)
      {
#if HAVE_MPI 
        typedef FillSlaveFunctor< typename RowMapperImp::GlobalKeyType, ParallelScalarProductType > SlaveFillFunctorType;
        // if entity is not interior, insert into overlap entries 
        if(en.partitionType() != InteriorEntity && gridPart.indexSet().contains(en) )
          rowMapper.mapEach( en, SlaveFillFunctorType( slaveDofs ) );

        // insert neighbors 
        typedef typename GridPartImp::template Codim<0>::EntityPointerType EntityPointerType; 
        typedef typename GridPartImp:: IntersectionIteratorType IntersectionIteratorType;
        typedef typename IntersectionIteratorType :: Intersection IntersectionType;
        IntersectionIteratorType endnit = gridPart.iend(en);
        for(IntersectionIteratorType nit = gridPart.ibegin(en);
            nit != endnit; ++nit)
        {
          // get intersection 
          const IntersectionType& inter = *nit;

          if( inter.neighbor() )
          {
            // get neighbor 
            EntityPointerType ep = inter.outside();
            const EntityImp& nb = *ep;
            // check partition type 
            if( nb.partitionType() != InteriorEntity && gridPart.indexSet().contains(nb) )
            {
              rowMapper.mapEach( nb, SlaveFillFunctorType( slaveDofs ) );
            }
          }
        }
#endif
      }
    };

  } // end namespace Fem

} // end namespace Dune 
#endif
