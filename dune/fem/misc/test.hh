// Test from Chuck Allison in C/C++ Users Journal, September 2000:
// "The Simplest Automated Unit Test Framework That Could Possibly Work"
// Numeric test feature courtesy to P. Frauenfelder (Concepts)

#ifndef TEST_HH
#define TEST_HH

#include <string>
#include <iostream>

using std::string;
using std::ostream;

namespace Dune {
/** @addtogroup Testing
  @{
  This is a basic testing framework to do black box testing in a systematic
  manner. The test framework consists of the two main classes Test and Suite.
  Test is the base class to derive user-defined tests. Suite is a container for
  test classes, working according to the Composite Pattern (GoF).
*/

/*! \file
  This file implements the class Test of the testing framework. Test is the
  base class for user-defined tests.
*/

// The following have underscores because they are macros
// (and it's impolite to usurp other users' functions!).
// For consistency, _succeed() also has an underscore.
#define _exec(str) doExec(*this, str,  __FILE__, __LINE__)
//! \brief Tests if \c cond is true (test passed) or not
//! This macro acts as the standard testing function, evaluating the result of 
//! boolean expression cond.
#define _test(cond) doTest(cond, #cond, __FILE__, __LINE__)

//! \brief Tests if \c arg1 and \c arg2 are identical up to a relative difference
//! Use this macro to test the equality of float types. A standard tolerance
//! of 10e-10 for the relative difference is allowed before the test fails.
#define _floatTest(arg1, arg2) \
  doFloatTest(arg1, arg2, #arg1, #arg2, __FILE__, __LINE__)

//! \brief Tests if \c arg1 and \c arg2 are identical up to a relative difference \c tol
//! Use this macro to test the equality of float types. A user-prescribed 
//! tolerance of \c tol for the relative difference is allowed before the test
//! fails.
#define _floatTestTol(arg1, arg2, tol) \
  doFloatTest(arg1, arg2, #arg1, #arg2, __FILE__, __LINE__, tol)

//! \brief Explicitly fails a test
//! Use this macro to explicitly signal a failure of the implemented algorithm
#define _fail(str) doFail(str, __FILE__, __LINE__)

/** Base class for tests. Writing a test case is done by deriving from this
    class.
    @author Chuck Allison, 2000.
    Adapted for float type comparisons by P. Frauenfelder
    @see Chuck Allison, <a href="http://www.cuj.com/documents/s=8035/cuj0009allison1/">The Simplest Automated Unit Test Framework That Could Possibly Work</a>, C/C++ Users Journal, September 2000.
*/
class Test
{
public:
  //! Constructor
  //! \param osptr Pointer to the output stream for reporting
  Test(ostream* osptr = &std::cout);
  //! Destructor
  virtual ~Test(){}
  //! \brief Runs the tests
  //! Must be overwritten in derived classes. run must explicitly call the 
  //! actual testing methods of a test class.
  virtual void run() = 0;
  
  //! Returns the number of passed tests
  long getNumPassed() const;
  //! Returns the number of failed tests
  long getNumFailed() const;
  //! Returns the output stream
  const ostream* getStream() const;
  //! Sets the output stream for the report
  //! \param osptr Pointer to the output stream
  void setStream(ostream* osptr);
  
  //! \brief Method signalling a successful test
  //! Use this method to explicitly signal a successful conclusion of a test.
  //! Consider using _test when you can express the fact using a boolean 
  //! condition.
  void _succeed();
  //! Display the report
  long report() const;
  //! Reset the test object's statistics
  virtual void reset();
  
protected:
  //! Internal class for compile-time checks
  template <bool cond>
  struct CompileAssertion {
    static inline void doExec(Test& tester, string str, 
                              const char* fname, long lineno);
  };

  //! Internal method for floating point comparisons
  void doFloatTest(double arg1, double arg2, const std::string& lbl1,
                   const std::string& lbl2, const char* fname, 
                   long lineno, double tol = 1e-10);
  //! Internal method for generic tests
  void doTest(bool cond, const string& lbl,
              const char* fname, long lineno);
  //! Internal method to signal a failure
  void doFail(const string& lbl,
              const char* fname, long lineno);

  //! Computer precision for doubles
  static const double eps;
  
private:
  ostream* m_osptr;
  long m_nPass;
  long m_nFail;
  
  // Disallowed:
  Test(const Test&);
  Test& operator=(const Test&);
};

/** @} end documentation group */

template<>
struct Test::CompileAssertion<true> {
  static inline void doExec(Test& tester, string str,
                            const char* fname, long lineno) {
    tester._succeed();
  }
};

template<>
struct Test::CompileAssertion<false> {
  static inline void doExec(Test& tester, string str,
                            const char* fname, long lineno) {
    tester.doFail(str, fname, lineno);
  }
};


inline
Test::Test(ostream* osptr)
{
    m_osptr = osptr;
    m_nPass = m_nFail = 0;
}

inline
long Test::getNumPassed() const
{
    return m_nPass;
}

inline
long Test::getNumFailed() const
{
    return m_nFail;
}

inline
const ostream* Test::getStream() const
{
    return m_osptr;
}

inline
void Test::setStream(ostream* osptr)
{
    m_osptr = osptr;
}

inline
void Test::_succeed()
{
    ++m_nPass;
}

inline
void Test::reset()
{
    m_nPass = m_nFail = 0;
}

} // end namespace Dune

#endif

