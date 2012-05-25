// Test from Chuck Allison in C/C++ Users Journal, September 2000:
// "The Simplest Automated Unit Test Framework That Could Possibly Work"
// Numeric test feature courtesy to P. Frauenfelder (Concepts)

#include "test.hh"
#include <cmath>
#include <iostream>
#include <typeinfo>     // Visual Studio requires /GR""
#include <limits>

#ifdef _MSC_VER
//Allow return-less mains:
#pragma warning(disable: 4541)
#endif

using namespace std;

namespace Dune {
const double Test::eps = std::numeric_limits<double>::epsilon();

void Test::doFloatTest(double arg1, double arg2, const std::string& lbl1,
                       const std::string& lbl2, const char* fname, 
                       long lineno, double relTol) {
  std::string lbl = lbl1 + std::string(" ~= ") + lbl2;
  double error;
  if (arg1 < 10.0*eps) {
    error = std::fabs(arg1-arg2);
  } else {
    error = std::fabs((arg1-arg2)/arg1);
  }
  if (error > relTol) {
    doFail(lbl, fname, lineno);
  } else {
    _succeed();
  }
}

void Test::doTest(bool cond, const std::string& lbl,
                  const char* fname, long lineno)
{
  if (!cond)
    doFail(lbl, fname, lineno);
  else
    _succeed();
}

void Test::doFail(const std::string& lbl,
                  const char* fname, long lineno)
{
  ++m_nFail;
  if (m_osptr)
    {
      *m_osptr << typeid(*this).name()
               << " failure: (" << lbl << "), "
               << fname
               << " (line " << lineno << ")\n";
    }
}

long Test::report() const
{
  if (m_osptr)
    {
      *m_osptr << "Test \"" 
               << typeid(*this).name() << "\":\n"
               << "\tPassed: " << m_nPass
               << "\tFailed: " << m_nFail
               << endl;
    }
  return m_nFail;
}

}
