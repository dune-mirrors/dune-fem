#ifndef DUNE_FEM_GRIDPART_TEST_FAILURE_HH
#define DUNE_FEM_GRIDPART_TEST_FAILURE_HH

#include <cassert>

//- dune-common includes
#include <dune/common/stdstreams.hh>

//- include this as very last
#include <dune/common/bartonnackmanifcheck.hh>


namespace Dune
{

  // Failure
  // -------

  struct Failure
  {
    /** \brief default constructor  */
    Failure () { }

  protected:
    friend std::ostream& operator<< ( std::ostream &, const Failure & );

    /** \brief print output to ostream */
    virtual void writeTo ( std::ostream &out ) const = 0;
  };

  std::ostream &operator<< ( std::ostream &s, const Failure &failure )
  {
    failure.writeTo( s );
    return s;
  }


  // FailureHandler
  // --------------

  template< class FailureHandlerImp >
  struct FailureHandler
  {
    /** \brief default constructor  */
    FailureHandler () { }

    /** \brief act on failure */
    template< class Failure >
    void operator() ( Failure &failure ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( ( asImp().operator()( failure ) ) );
      return asImp().operator()( failure );
    }

  private:
    // hide copy constructor
    FailureHandler ( const FailureHandler & );
    // hide assignment operator
    FailureHandler &operator= ( const FailureHandler & );

    //  Barton-Nackman trick
    FailureHandlerImp &asImp ()
    {
      return static_cast< FailureHandlerImp & >( *this );
    }

    //  Barton-Nackman trick
    const FailureHandlerImp &asImp () const
    {
      return static_cast< const FailureHandlerImp & >( *this );
    }
  };



  // DefaultFailureHandler
  // ---------------------

  struct DefaultFailureHandler
  : public FailureHandler< DefaultFailureHandler >
  {
    /** \brief act on failure */
    template< class Failure >
    void operator() ( Failure &failure ) const
    {
      std::cerr<< failure << std::endl;
      assert( false );
    }
  };

} // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_TEST_FAILURE_HH
