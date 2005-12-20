#ifndef DUNE_CACHEQUAD_TEST_HH
#define DUNE_CACHEQUAD_TEST_HH

#include "../../misc/test.hh"
#include "../cachequad.hh"

namespace Dune {

  class CachingQuadrature_Test : public Test 
  {
  public:
    CachingQuadrature_Test(std::string albertaGridFile,
                           std::string aluGridHexaFile,
                           std::string aluGridTetraFile) :
      albertaGridFile_(albertaGridFile),
      aluGridHexaFile_(aluGridHexaFile),
      aluGridTetraFile_(aluGridTetraFile)
    {}

    virtual void run();

    void codim0Test();
    void codim1AlbertaTest();
    void codim1ALUHexaTest();
    void codim1ALUTetraTest();
    void codim1SGridTest();

  private:
    std::string albertaGridFile_;
    std::string aluGridHexaFile_;
    std::string aluGridTetraFile_;
  };

} // end namespace Dune

#endif
