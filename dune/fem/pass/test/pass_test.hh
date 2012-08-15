#ifndef DUNE_FEM_PASS_TEST_HH
#define DUNE_FEM_PASS_TEST_HH

#include <config.h>
#include <string>

#include <dune/fem/misc/test.hh>

namespace Dune 
{

  namespace Fem
  {
  
    class Pass_Test : public Test 
    {
    public:
      Pass_Test(std::string gridFile) :
        gridFile_(gridFile)
      {}

      virtual void run();

    private:
      void functorTest();
      void dummyTest();
      void dgTest();

    private:
      std::string gridFile_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_TEST_HH
