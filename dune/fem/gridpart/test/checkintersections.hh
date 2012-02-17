#ifndef DUNE_FEM_GRIDPART_TEST_CHECKINTERSECTIONS_HH
#define DUNE_FEM_GRIDPART_TEST_CHECKINTERSECTIONS_HH

//- dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>

//- dune-fem includes
#include <dune/fem/gridpart/test/failure.hh>


namespace Dune
{
  namespace Fem
  {
    template< class GridPartType, class FailureHandler >
    class CheckIntersections
    {
      struct AssignmentFailure
      : public Failure
      {
        virtual void writeTo ( std::ostream &out ) const
        {
          out <<  __FILE__
              << ":" << __LINE__ << ": Failure :"
              << "Assignment failure";
        }
      };

      template< class IntersectionIterator >
      static void checkIntersectionIteratorAssignment( const IntersectionIterator &it, 
                                                       FailureHandler &failureHandler );

    public:
      static void check ( const GridPartType &gridPart, FailureHandler &failureHandler )
      {
        typedef typename GridPartType::template Codim< 0 >::IteratorType IteratorType;
        typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

        const IteratorType end = gridPart.template end< 0>();
        for( IteratorType it = gridPart.template begin< 0 >(); it != end; ++it )
        {
          const EntityType &entity = *it;

          typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
          const IntersectionIteratorType iend = gridPart.iend( entity );
          for( IntersectionIteratorType iit = gridPart.ibegin( entity ); iit != iend; ++iit )
            checkIntersectionIteratorAssignment( iit, failureHandler );
        }
      }
    };



    // Implementation of CheckIntersections
    // ------------------------------------

    template< class GridPartType, class FailureHandler >
    template< class IntersectionIterator >
    inline void CheckIntersections< GridPartType, FailureHandler >
      ::checkIntersectionIteratorAssignment( const IntersectionIterator &it, FailureHandler &failureHandler )
    {
      // type of intersection iterator
      typedef IntersectionIterator IntersectionIteratorType;

      AssignmentFailure failure;
      // assert assignment
      IntersectionIteratorType it2 = it;
      if( it != it2 )
        failureHandler( failure );
      if( it->inside() != it2->inside() || it->outside() != it2->outside() )
        failureHandler( failure );

      // now increment second iterator
      ++it2;
      if( it == it2 )
        failureHandler( failure );
    }

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_TEST_CHECKINTERSECTIONS_HH
