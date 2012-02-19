#ifndef DUNE_FEM_GRIDPART_TEST_CHECKINDEXSET_HH
#define DUNE_FEM_GRIDPART_TEST_CHECKINDEXSET_HH

//- C++ includes
#include <algorithm>
#include <set>
#include <vector>

//- dune-common includes
#include <dune/common/fvector.hh>
#include <dune/common/forloop.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>

//- dune-geometry includes
#include <dune/geometry/referenceelements.hh>

//- dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/gridpart/test/failure.hh>


namespace Dune
{

  namespace Fem
  {
    // CheckIndexSet
    // -------------

    template< class GridPartType, class FailureHandler >
    struct CheckIndexSet
    {
      /** \brief dimension */
      static const int dimension = GridPartType::dimension;

      /** \brief single coordinate type */
      typedef typename GridPartType::ctype ctype;

      /** \brief type definitions */
      template< int codim >
      struct Codim
      : public GridPartType::template Codim< codim >
      { };

      /** \brief type of index set */ 
      typedef typename GridPartType::IndexSetType IndexSetType;

      /** \brie local failure type */
      struct SizeFailure;

    private:
      template< int codim >
      static void checkSize ( const IndexSetType &indexSet, 
                              const GridPartType &gridPart,
                              FailureHandler &failureHandler );

      static void checkGeomTypes ( const IndexSetType &indexSet, 
                                   const GridPartType &gridPart,
                                   FailureHandler &failureHandler );

      template< class GridPart, bool >
      struct If
      {
        template< int codim >
        static void count ( const IndexSetType &, const GridPartType &,
                            FailureHandler & )
        { }
      };

      template< class GridPart >
      struct If< GridPart, true >
      {
        template< int codim >
        static void count ( const IndexSetType &indexSet,
                            const GridPartType &gridPart,
                            FailureHandler &failureHandler )
        {
          checkSize< codim >( indexSet, gridPart, failureHandler );
        }
      };

      template< int codim >
      struct CheckSubEntityGeometry
      {
        template< class Entity >
        static void apply ( const Entity &entity, FailureHandler &failureHandler )
        {
          const bool hasEntity = Dune::Fem::GridPartCapabilities::hasEntity< GridPartType, codim >::v;
          If< GridPartType, hasEntity >::template check< codim, Entity >( entity, failureHandler );
        }
      };

      template< int codim >
      struct CheckSize
      {
        static void apply ( const IndexSetType &indexSet, const GridPartType &gridPart, 
                           FailureHandler &failureHandler )
        {
          const bool hasEntity = Dune::Fem::GridPartCapabilities::hasEntity< GridPartType, codim >::v;
          If< GridPartType, hasEntity >::template count< codim >( indexSet, gridPart, failureHandler );
        }
      };

    public:
      static void check ( const GridPartType &gridPart, FailureHandler &failureHandler )
      {
        const IndexSetType &indexSet = gridPart.indexSet();

        // check geometry types
        checkGeomTypes( indexSet, gridPart, failureHandler );

        // check size
        ForLoop< CheckSize, 0, dimension >::apply( indexSet, gridPart, failureHandler );
      }
    };



    // Implementation of CheckIndexSet
    // -------------------------------

    template< class GridPartType, class FailureHandler >
    template< int codim >
    inline void CheckIndexSet< GridPartType, FailureHandler >
      ::checkSize ( const IndexSetType &indexSet, 
                    const GridPartType &gridPart,
                    FailureHandler &failureHandler )
    {
      unsigned int count = 0;

      typedef typename Codim< codim >::IteratorType IteratorType;
      const IteratorType end = gridPart.template end< codim >();
      for( IteratorType it = gridPart.template begin< codim >(); it != end; ++it )
        ++count;

      if( (unsigned int) indexSet.size( codim ) != count )
      {
        SizeFailure failure( indexSet.size( codim ), count, codim );
        failureHandler( failure );
      }
    }


    template< class GridPartType, class FailureHandler >
    inline void CheckIndexSet< GridPartType, FailureHandler >
      ::checkGeomTypes ( const IndexSetType &indexSet, 
                         const GridPartType &gridPart,
                         FailureHandler &failureHandler )
    {
      // build vectors of geometry types
      std::vector< GeometryType > geomTypes[ dimension + 1 ];

      // find all geometry types for codimension 0
      if( GridPartCapabilities::hasSingleGeometryType< GridPartType >::v )
      {
        unsigned int topologyId 
          = GridPartCapabilities::hasSingleGeometryType< GridPartType >::topologyId;
        geomTypes[ 0 ].push_back( GeometryType::GeometryType( topologyId, dimension ) );
      }
      else
      {
        typedef typename Codim< 0 >::IteratorType IteratorType;
        const IteratorType end = gridPart.template end< 0 >();
        for( IteratorType it = gridPart.template begin< 0 >(); it != end; ++it )
          geomTypes[ 0 ].push_back( it->type() );
        std::sort( geomTypes[ 0 ].begin(), geomTypes[ 0 ].end() );
        std::vector< GeometryType > :: iterator it = std::unique( geomTypes[ 0 ].begin(), geomTypes[ 0 ].end() );
        geomTypes[ 0 ].erase( it, geomTypes[ 0 ].end() );
      }

      // build all other codimension from refernce elements
      const std::vector< GeometryType >::const_iterator end = geomTypes[ 0 ].end();
      for( std::vector< GeometryType >::const_iterator it = geomTypes[ 0 ].begin(); it !=  end; ++it )
      {
        const GenericReferenceElement< ctype, dimension > &referenceElement
          = GenericReferenceElements< ctype, dimension >::general( *it );
        for( int cd = 1; cd <= dimension; ++cd )
        {
          const int nSubEntities = referenceElement.size( cd );
          for( int i = 0; i < nSubEntities; ++i )
            geomTypes[ cd ].push_back( referenceElement.type( i, cd ) );
        }
      }

      // sort vectors and remove duplicates
      for( int cd = 1; cd <= dimension; ++cd )
      {
        std::sort( geomTypes[ cd ].begin(), geomTypes[ cd ].end() );
        std::vector< GeometryType > :: iterator it = std::unique( geomTypes[ cd ].begin(), geomTypes[ cd ].end() );
        geomTypes[ cd ].erase( it, geomTypes[ cd ].end() );
      }

      for( int cd = 0; cd <= dimension; ++cd )
      {
        if( !std::equal( geomTypes[ cd ].begin(), geomTypes[ cd ].end(), 
                         indexSet.geomTypes( cd ).begin() ) )
        {
          abort();
        }
      }
    }



    // Implementation of Failures used in CheckIndexSet
    // ------------------------------------------------

    template< class GridPartType, class FailureHandler >
    struct CheckIndexSet< GridPartType, FailureHandler >::SizeFailure
    : public Failure
    {
      /** \brief constructor */
      SizeFailure ( int size, int count, int codim )
      : size_( size ),
        count_( count ),
        codim_( codim )
      { }

      /** \brief write message to stream */
      virtual void writeTo ( std::ostream &out ) const
      {
        out <<  __FILE__
          << ":" << __LINE__ << ": Failure: "
          << "indexSet().size( " << codim_ << " ) = " << size_
          << ", but " << count_ << " entities found";
      }

    private:
      int size_, count_, codim_;
    };

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_TEST_CHECKINDEXSET_HH
