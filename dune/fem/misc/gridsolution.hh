#ifndef DUNE_FEM_GRIDSOLUTION_HH
#define DUNE_FEM_GRIDSOLUTION_HH

#include <tuple>

#include <dune/common/exceptions.hh>
#include <dune/grid/utility/hierarchicsearch.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#endif

namespace Dune
{

  namespace Fem
  {

    // GridSolution
    //-------------

    /** \class GridSolution
     *  \brief creates a function with evaluate method from a check point
     *
     *  \tparam GridImp Grid type
     *  \tparam DiscreteFunctionImp Discrete function type
     */
    template <class GridImp, class DiscreteFunctionImp >
    class GridSolution
    {
    public:
      typedef GridImp  GridType;
      typedef DiscreteFunctionImp DiscreteFunctionType;
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: RangeType RangeType ;
      typedef typename DiscreteFunctionSpaceType :: DomainType DomainType ;
      typedef typename DiscreteFunctionSpaceType :: GridPartType   GridPartType ;
      typedef typename GridPartType :: IndexSetType IndexSetType;
      typedef CheckPointer< GridType >   CheckPointerType;

      typedef typename GridType :: template Codim<0> :: Entity        EntityType;

      typedef HierarchicSearch< GridType, IndexSetType > HierarchicSearchType;

      typedef std::tuple< DiscreteFunctionType* > IOTupleType;

    protected:
      GridType* grid_;
      GridPtr< GridType > gridPtr_;
      GridPartType gridPart_;
      DiscreteFunctionSpaceType space_;
      DiscreteFunctionType discreteFunction_;
      ConstLocalFunction< DiscreteFunctionType > lf_;
      IOTupleType data_;
      HierarchicSearchType hierarchicSearch_;

    public:
      GridType& grid() { assert( grid_ ); return *grid_ ; }
      const GridType& grid() const { assert( grid_ ); return *grid_ ; }

      //! Constructor
      explicit GridSolution(const std::string checkPointFile,
                            const int rank = -1 ) :
        grid_( CheckPointerType :: restoreGrid( checkPointFile, rank ) ),
        gridPtr_(),
        gridPart_( grid() ),
        space_( gridPart_ ),
        discreteFunction_("grid-sol", space_),
        lf_( discreteFunction_ ),
        data_( &discreteFunction_ ),
        hierarchicSearch_( grid(), gridPart_.indexSet() )
      {
        // store grid pointer
        gridPtr_ = grid_;
        DUNE_THROW(NotImplemented," GridSolution needs to be switched to DataWriter");
        // restore data from checkpoint
        CheckPointerType :: restoreData( grid(), data_, checkPointFile, rank );
      }

      /** \brief evaluates in a given space-time point
       *  \param[in] x Point in global coordinates
       *  \param[in] time Time
       *  \param[out] result The value of the discrete function in space-time point \f[ (x,time)\f]
       *
       *  \tparam PointType The point type
       */
      void evaluate(const DomainType& x, const double time, RangeType& result) const
      {
        evaluate(x, result );
      }

      /** \brief evaluates in a given space point
       *  \param[in] x Point in global coordinates
       *  \param[out] result The value of the discrete function in space point @a x
       *
       *  \tparam PointType The point type
       */
      void evaluate(const DomainType& x, RangeType& result) const
      {
        // search entity
        EntityType entity = hierarchicSearch_.findEntity( x );

        typedef typename EntityType :: Geometry Geometry;
        const Geometry& geo = entity.geometry();

        const DomainType local = geo.local( x );
#ifndef NDEBUG
        {
          // check that corners are within the reference element of the given type
          const auto &refElement = Dune::ReferenceElements< typename GridType::ctype, GridType::dimensionworld >::general( entity.type() );

          assert( refElement.checkInside( local ) );
        }
#endif

        // evaluate discrete function
        lf_.bind( entity );
        lf_.evaluate( local, result );
        lf_.unbind();
      }

      //! \brief writes a discrete function
      static void writeDiscreteFunction(const GridType& grid,
                                        const DiscreteFunctionType& discreteFunction,
                                        const double time,
                                        const int writeStep )
      {
        DUNE_THROW(NotImplemented," writeDiscreteFunction needs to be switched to DataWriter");
        //typedef std::tuple< const DiscreteFunctionType* > OutputTuple;
        //OutputTuple data( &discreteFunction );
        // false means don't backup persistent objects
        // CheckPointerType :: writeSingleCheckPoint( grid, data, time, false, writeStep );
      }

      const DiscreteFunctionType& discreteFunction() const { return discreteFunction_ ; }
    };

    ///////////////////////////////////////////////////////////////////
    //
    //
    //
    ///////////////////////////////////////////////////////////////////
    /** \class GridSolution
     *  \brief creates a function with evaluate method from a check point
     *
     *  \tparam GridImp Grid type
     *  \tparam DiscreteFunctionImp Discrete function type
     */
    template <class GridImp, class DiscreteFunctionImp >
    class GridSolutionVector
    {
    public:
      typedef GridImp  GridType;
      typedef DiscreteFunctionImp DiscreteFunctionType;

      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType     FunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: GridPartType          GridPartType;
      typedef typename FunctionSpaceType :: DomainType                    DomainType;
      typedef typename FunctionSpaceType :: DomainFieldType               DomainFieldType;
      typedef typename FunctionSpaceType :: RangeType                     RangeType;
      typedef typename FunctionSpaceType :: RangeFieldType                RangeFieldType;

      typedef GridSolution< GridType, DiscreteFunctionType > GridSolutionType ;

    protected:
      template <class DomainType, class Grid>
      struct CheckDomain
      {
        static bool isInside(const DomainType& x, const Grid& grid )
        {
          abort();
          typedef typename Grid :: LevelGridView MacroView ;
          typedef typename MacroView :: template Codim< 0 > :: Iterator Iterator ;
          typedef typename Iterator :: Entity Entity;
          const MacroView& macroView = grid.levelGridView( 0 );
          const Iterator end = macroView.template end< 0 > ();
          for( Iterator it = macroView.template begin< 0 > (); it != end; ++it )
          {
            const Entity& entity = * it ;
            // check that corners are within the reference element of the given type
            const auto &refElement = Dune::ReferenceElements< typename Grid::ctype, Grid::dimensionworld >::general( entity.type() );

            typedef typename Entity :: Geometry Geometry;
            const Geometry& geo = entity.geometry();

            if( refElement.checkInside( geo.local( geo.center() ) ) )
              return true ;
          }
          return false ;
        }
      };

#if HAVE_DUNE_SPGRID
      template <class DomainType, class ct, int dim, template< int > class Strategy , class Comm>
      struct CheckDomain< DomainType, SPGrid< ct, dim, Strategy, Comm > >
      {
        typedef SPGrid< ct, dim, Strategy, Comm >  Grid;
        static bool isInside(const DomainType& x, const Grid& grid )
        {
          return grid.domain().contains( x );
        }
      };
#endif

      const int numProcs_;
      std::vector< GridSolutionType* > solutions_;

      int numProcs(const std::string& checkPointFile) const
      {
        int numProc = MPIManager :: size();
        readParameter(checkPointFile, "NumberProcessors", numProc, Parameter::verbose () );
        return numProc;
      }
    public:
      //! Constructor
      explicit GridSolutionVector(const std::string checkPointFile) :
        numProcs_( numProcs( checkPointFile ) ),
        solutions_( numProcs_, (GridSolutionType *) 0 )
      {
        for(int p=0; p<numProcs_; ++p)
        {
          if( Parameter::verbose () )
            std::cout << "GridSolutionVector: Reading Grid " << p << " from checkpoint" << std::endl;
          solutions_[ p ] = new GridSolutionType( checkPointFile, p );
        }
      }

      ~GridSolutionVector()
      {
        for(int p=0; p<numProcs_; ++p)
        {
          delete solutions_[ p ];
          solutions_[ p ] = 0;
        }
      }

      /** \brief evaluates in a given space-time point
       *  \param[in] x Point in global coordinates
       *  \param[in] time Time
       *  \param[out] result The value of the discrete function in space-time point \f[ (x,time)\f]
       *
       *  \tparam PointType The point type
       */
      void evaluate(const DomainType& x, const double time, RangeType& result) const
      {
        evaluate(x, result );
      }

      /** \brief evaluates in a given space point
       *  \param[in] x Point in global coordinates
       *  \param[out] result The value of the discrete function in space point @a x
       *
       *  \tparam PointType The point type
       */
      void evaluate(const DomainType& x, RangeType& result) const
      {
        for(int p=0; p<numProcs_; ++p)
        {
          assert( solutions_[ p ] );
          const GridSolutionType& gridSolution = *(solutions_[ p ]);
          if( isInDomain( x, gridSolution.grid() ) )
          {
            //std::cout << "Found grid " << p << " for x = " << x << std::endl;
            gridSolution.evaluate( x, result );
            return ;
          }
        }

        std::cerr << "GridSolutionVector::evaluate: no grid found for point " << x << std::endl;
        assert( false );
        abort();
      }

      bool isInDomain( const DomainType& x, const GridType& grid ) const
      {
        return CheckDomain<DomainType, GridType> :: isInside( x, grid );
      }

      //! \brief writes a discrete function
      static void writeDiscreteFunction(const GridType& grid,
                                        const DiscreteFunctionType& discreteFunction,
                                        const double time,
                                        const int writeStep = 0 )
      {
        GridSolutionType :: writeDiscreteFunction( grid, discreteFunction, time, writeStep );
      }

      const DiscreteFunctionType& discreteFunction( const int rank ) const
      {
        assert( rank < numProcs_ );
        return solutions_[ rank ]->discreteFunction();
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDSOLUTION_HH
