#ifndef DUNE_FEM_ERRORDISPLAY_HH
#define DUNE_FEM_ERRORDISPLAY_HH

#include <dune/fem/function/common/discretefunctionadapter.hh>

namespace Dune
{

  template< class DiscreteFunction, class SolutionType, bool withTime = true >
  class DisplayErrorFunction;

  
  template< class DiscreteFunction, class SolutionType >
  class DisplayErrorFunction< DiscreteFunction, SolutionType, true >
  {
    typedef DisplayErrorFunction< DiscreteFunction, SolutionType > ThisType;

  public:
    typedef DiscreteFunction DiscreteFunctionType;

    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  protected:
    /** \cond */
    struct Error
    {
      typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
      
      typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType
        FunctionSpaceType;
      typedef typename FunctionSpaceType :: DomainType DomainType;
      typedef typename FunctionSpaceType :: RangeType RangeType;
      typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

      static const int dimDomain = DiscreteFunctionSpaceType :: dimDomain;
      static const int dimRange = DiscreteFunctionSpaceType :: dimRange;

    protected:
      typedef typename GridPartType :: template Codim< 0 > :: IteratorType :: Entity
        EntityType;
      typedef typename EntityType :: Geometry GeometryType;

    protected:
      LocalFunctionType lUh_;
      const SolutionType &initU0_;
      const GeometryType *geometry_;
      const double time_;
      bool initialized_;

    public:
      Error ( const DiscreteFunctionType &Uh,
              const SolutionType &solution,
              double time = 0 )
      : lUh_( Uh ), 
        initU0_( solution ),
        geometry_( 0 ),
        time_( time ),
        initialized_( false )
      {}

      template< class PointType >
      void evaluate ( const PointType &x, RangeType& ret) const
      {
        assert(initialized_);
        lUh_.evaluate( x, ret );
        DomainType global = geometry_->global( coordinate( x ) );
        RangeType phi;
        initU0_.evaluate( time_, global, phi );
        ret -= phi;
      }

      template< class PointType >
      void jacobian ( const PointType &x, JacobianRangeType &ret ) const
      {
        abort();
      }

      inline void init ( const EntityType &entity )
      {
        lUh_.init( entity );
        geometry_ = &( entity.geometry() );
        initialized_ = true;
      }
    };
    /** \endcond */

    typedef LocalFunctionAdapter< Error > ErrorFunctionType;
    
  protected:
    const GridPartType &gridPart_;
    Error error_;
    ErrorFunctionType errorFunction_;
    
  public:
    template< class GrapeDispType >
    DisplayErrorFunction ( GrapeDispType &disp,
                           const DiscreteFunctionType &Uh,
                           const SolutionType &solution,
                           const double time = 0 )
    : gridPart_( Uh.space().gridPart() ),
      error_( Uh, solution, time ),
      errorFunction_( "error", error_, gridPart_ )
    {
      disp.addData( errorFunction_, "error", time );
    }

  private:
    DisplayErrorFunction ( const ThisType & );
    ThisType &operator= ( const ThisType & );
  };



  template< class DiscreteFunction, class SolutionType >
  class DisplayErrorFunction< DiscreteFunction, SolutionType, false >
  {
    typedef DisplayErrorFunction< DiscreteFunction, SolutionType > ThisType;

  public:
    typedef DiscreteFunction DiscreteFunctionType;

    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  protected:
    /** \cond */
    struct Error
    {
      typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
      
      typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType
        FunctionSpaceType;
      typedef typename FunctionSpaceType :: DomainType DomainType;
      typedef typename FunctionSpaceType :: RangeType RangeType;
      typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

      static const int dimDomain = DiscreteFunctionSpaceType :: dimDomain;
      static const int dimRange = DiscreteFunctionSpaceType :: dimRange;

    protected:
      typedef typename GridPartType :: template Codim< 0 > :: IteratorType :: Entity
        EntityType;
      typedef typename EntityType :: Geometry GeometryType;

    protected:
      LocalFunctionType lUh_;
      const SolutionType &initU0_;
      const GeometryType *geometry_;
      bool initialized_;

    public:
      Error ( const DiscreteFunctionType &Uh,
              const SolutionType &solution )
      : lUh_( Uh ),
        initU0_( solution ),
        geometry_( 0 ),
        initialized_( false )
      {}

      template< class PointType >
      void evaluate ( const PointType &x, RangeType& ret) const
      {
        assert(initialized_);
        lUh_.evaluate( x, ret );
        DomainType global = geometry_->global( coordinate( x ) );
        RangeType phi;
        initU0_.evaluate( global, phi );
        ret -= phi;
      }

      template< class PointType >
      void jacobian ( const PointType &x, JacobianRangeType &ret ) const
      {
        abort();
      }

      inline void init ( const EntityType &entity )
      {
        lUh_.init( entity );
        geometry_ = &( entity.geometry() );
        initialized_ = true;
      }
    };
    /** \endcond */

    typedef LocalFunctionAdapter< Error > ErrorFunctionType;
    typedef DiscreteFunctionAdapter< SolutionType, GridPartType >
      GridSolutionType;
    
  protected:
    const GridPartType &gridPart_;
    GridSolutionType gridSolution_;
    Error error_;
    ErrorFunctionType errorFunction_;
    
  public:
    template< class GrapeDispType >
    DisplayErrorFunction ( GrapeDispType &disp,
                           const DiscreteFunctionType &Uh,
                           const SolutionType &solution )
    : gridPart_( Uh.space().gridPart() ),
      gridSolution_( "exact solution", solution, gridPart_ ),
      error_( Uh, solution ),
      errorFunction_( "error", error_, gridPart_ )
    {
      disp.addData( errorFunction_ );
      disp.addData( gridSolution_ );
    }

  private:
    DisplayErrorFunction ( const ThisType & );
    ThisType &operator= ( const ThisType & );
  };

}

#endif // DUNE_FEM_ERRORDISPLAY_HH

