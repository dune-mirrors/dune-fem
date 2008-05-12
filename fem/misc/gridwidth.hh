#ifndef DUNE_GRIDWIDTH_HH
#define DUNE_GRIDWIDTH_HH

//- Dune includes 
#include <dune/common/fvector.hh>
#include <dune/grid/common/capabilities.hh>

#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune
{

  template < class GridType >
  class GridWidthProvider;



  //! utility functions for calculating the maximum grid width 
  struct GridWidth 
  {
    template <class GridPartType, class EntityType> 
    static inline double maxEdgeWidth(const GridPartType& gridPart, const EntityType &en)
    {
      typedef typename GridPartType::IntersectionIteratorType
        IntersectionIteratorType;

      typedef typename EntityType :: Geometry Geometry;
      enum { dim = EntityType::dimension };

      Dune::FieldVector<double, dim-1> tmp1(0.5);
     
      const Geometry& geo = en.geometry();
      const double vol = geo.volume();
      const double fak = (dim == 3) ? 1 : 2;

      double minFace = 0;
      IntersectionIteratorType endit = gridPart.iend(en);
      for(IntersectionIteratorType it = gridPart.ibegin(en);
          it != endit; ++it)
      {
        double face = it->integrationOuterNormal(tmp1).two_norm();
        minFace = std::max(minFace,face);
      }
      return fak*vol/minFace;
    }

    template <class GridPartType> 
    static inline double calcGridWidth (const GridPartType & gridPart)
    {     
      double maxwidth = 0.0;
      typedef typename GridPartType :: GridType GridType;
      typedef typename GridPartType::template Codim<0> :: IteratorType IteratorType;
      
      // unstructured case 
      if( Capabilities::IsUnstructured<GridType>::v )
      {
        IteratorType endit = gridPart.template end<0> (); 
        for(IteratorType it = gridPart.template begin<0> (); 
            it != endit; ++it )
        {
          double w = maxEdgeWidth(gridPart,*it);
          if(w > maxwidth) maxwidth = w;
        }
      }
      else 
      {
        // here we only need to check one element 
        IteratorType it = gridPart.template begin<0> (); 
        if( it != gridPart.template end<0> () )
        {
          double w = maxEdgeWidth(gridPart,*it);
          if(w > maxwidth) maxwidth = w;
        }
      }
      return maxwidth;
    }   

    template <class EntityType> 
    static inline std::pair<double,double> maxLeafEdgeWidth(const EntityType &en)
    {
      std::pair<double,double> p(1e308, 0);

      typedef typename EntityType :: LeafIntersectionIterator IntersectionIteratorType;
      typedef typename EntityType :: Geometry Geometry;
      enum { dim = EntityType::dimension };

      Dune::FieldVector<double, dim-1> tmp1(0.5);
      
      const Geometry& geo = en.geometry();
      const double vol = (dim == 3 || geo.type().isCube() ) ? geo.volume() : (2.*geo.volume());
      
      IntersectionIteratorType endit = en.ileafend();
      for(IntersectionIteratorType it = en.ileafbegin();
          it != endit; ++it)
      {
        const double face = it->integrationOuterNormal(tmp1).two_norm();
        p.first  = std::min( p.first , face);
        p.second = std::max( p.second, face);
      }
      p.first  = vol/p.first;
      p.second = vol/p.second;
      return p;
    }

    template <class GridType> 
    static inline double calcLeafWidth (
        const GridType & grid, double& maxwidth, double& minwidth )
    {     
      minwidth = 1e308;
      maxwidth = 0;
      typedef typename GridType::template Codim<0> :: LeafIterator IteratorType; 
      
      // unstructured case 
      if( Capabilities::IsUnstructured<GridType>::v )
      {
        IteratorType endit = grid.template leafend<0> (); 
        for(IteratorType it = grid.template leafbegin<0> (); 
            it != endit; ++it )
        {
          std::pair< double, double> p = maxLeafEdgeWidth(*it);
          if(p.first > maxwidth) maxwidth = p.first;
          if(p.second < minwidth) minwidth = p.second;
        }
      }
      else 
      {
        // here we only need to check one element 
        IteratorType it = grid.template leafbegin<0> (); 
        if( it != grid.template leafend<0> () )
        {
          std::pair< double, double> p = maxLeafEdgeWidth(*it);
          if(p.first > maxwidth) maxwidth = p.first;
          if(p.second < minwidth) minwidth = p.second;
        }
      }
      return maxwidth;
    }   
  };



  //! utility functions for calculating the maximum grid width 
  template< class GridType >
  class GridWidthProvider 
  {
    typedef GridWidthProvider< GridType > ThisType;

  public:
    //! type of singleton provider 
    typedef SingletonList< const GridType *, ThisType > ProviderType;

  protected:
    typedef DofManager< GridType > DofManagerType;
    typedef DofManagerFactory< DofManagerType > DMFactoryType;

    const GridType &grid_;
    const DofManagerType &dm_;
    
    mutable double maxWidth_;
    mutable double minW_;
    mutable int sequence_;

    GridWidthProvider( const ThisType & );

  public:
    //! constructor taking grid part 
    inline explicit GridWidthProvider ( const GridType &grid )
    : grid_( grid ),
      dm_( DMFactoryType :: getDofManager( grid_ ) ),
      maxWidth_( -1.0 ),
      minW_( -1.0 ),
      sequence_( -1 )
    {}

    inline explicit GridWidthProvider ( const GridType *grid )
    : grid_( *grid ),
      dm_( DMFactoryType :: getDofManager( grid_ ) ),
      maxWidth_( -1.0 ),
      minW_( -1.0 ),
      sequence_( -1 )
    {}

    //! return characteristic grid width 
    inline double gridWidth () const
    {
      return maxGridWidth();
    }
    
    //! return characteristic grid width 
    inline double minGridWidth () const
    {
      calcWidths();
      return minW_;
    }
    
    //! return characteristic grid width 
    inline double maxGridWidth () const
    {
      calcWidths();
      return maxWidth_;
    }

  protected:  
    inline void calcWidths () const
    {
      if( dm_.sequence() != sequence_ )
      {
        GridWidth :: calcLeafWidth( grid_, this->maxWidth_, this->minW_ );
        assert( this->minW_ > 0 );
        assert( this->maxWidth_ > 0 );
        sequence_ = dm_.sequence();
      }
    }
  };

} // end namespace Dune

#endif

