// Test suite from Chuck Allison in C/C++ Users Journal, September 2000:
// "The Simplest Automated Unit Test Framework That Could Possibly Work"

#ifndef SUITE_HH
#define SUITE_HH

#include "test.hh"   // includes <string>, <iostream>
#include <vector>
#include <stdexcept>
using std::string;
using std::ostream;
using std::vector;

namespace Dune {
/** @addtogroup Testing
    @{
*/

/*! \file
  This file implements the class Suite and TestSuiteError of the testing
  framework. Suites are used to gather together Tests.
*/

//! \brief Exception for test suites
class TestSuiteError : public std::logic_error
{
public:
    TestSuiteError(const string& s = "")
        : logic_error(s)
    {}
};

 
/** Suite of tests. Running a group of test cases is best be done by
    adding each test case to a test suite. This suite does then the
    run and report of the whole group.
    @author Chuck Allison, 2000.
    @see Chuck Allison, <a href="http://www.cuj.com/documents/s=8035/cuj0009allison1/">The Simplest Automated Unit Test Framework That Could Possibly Work</a>, C/C++ Users Journal, September 2000.
*/
class Suite
{
public:
  //! Constructor
  //! \param name Name of the test suite
  //! \param osptr Pointer to the output stream. Defaults to cout.
  Suite(const string& name, ostream* osptr = &std::cout);
  
  //! Returns name of the test suite
  string getName() const;
  //! Get number of passed tests of contained test objects
  long getNumPassed() const;
  //! Get number of failed tests of contained test objects
  long getNumFailed() const;
  //! Get output stream pointer
  const ostream* getStream() const;
  //! Set output stream for reporting
  void setStream(ostream* osptr);
  
  //! Add a test object to the this suite
  void addTest(Test* t) throw (TestSuiteError);
  //! Add a subsuite to this suite
  void addSuite(const Suite&) throw(TestSuiteError);
  //! Run all contained tests and subsuites
  void run();     // Calls Test::run() repeatedly
  //! Report test results
  long report() const;
  //! Delete all contained tests
  void free();    // deletes tests
  
private:
  string m_name;
  ostream* m_osptr;
  vector<Test*> m_tests;
  void reset();
  
  // Disallowed ops:
  Suite(const Suite&);
  Suite& operator=(const Suite&);
};

/** @} end documentation group */

inline
Suite::Suite(const string& name, ostream* osptr)
  : m_name(name)
{
  m_osptr = osptr;
}

inline
string Suite::getName() const
{
  return m_name;
}

inline
const ostream* Suite::getStream() const
{
  return m_osptr;
}

inline
void Suite::setStream(ostream* osptr)
{
  m_osptr = osptr;
}

} // end namespace Dune
#endif

