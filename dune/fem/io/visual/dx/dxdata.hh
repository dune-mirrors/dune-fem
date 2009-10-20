#define DUNE_IO_DXDATA_HH
#ifndef DUNE_IO_DXDATA_HH

//- System includes
#include <iostream>
#include <fstream>
#include <vector>

//- Dune includes
#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>

//- Local includes


namespace Dune {
  // Forward declaration
  template <int dim, typename CoordType>
  class VirtualRefinement;

  //! A helper struct to store the values for refinement
  struct RefinementData {
    //- Data members
    std::vector<int> startingCellIndex_;
    std::vector<int> startingVertexIndex_;
    
    int nFineCells_;
    int nFineVertices_;
    int cornersPerCell_;

    GeometryType geometryType_; // * should not be here
  };

  template <class Grid>
  void writeGridDX(Grid &grid, std::ostream& out);

  //! The DX output writer
  //! Maybe I am going to break this down into a data and mesh writer
  //! (or even more)
  template <class FunctionSpaceT, bool binary = true>
  class DXWriter {
  public:
    //- Typedefs
    typedef FunctionSpaceT FunctionSpaceType;

    //- Constructor and destructor
    //! Constructor taking the space on which the data lives
    DXWriter(const FunctionSpaceT& spc, 
             std::string outName,
             bool multiFile = false);

    //! Destructor
    ~DXWriter();

    //- Public methods
    //! Writes the output to a file
    template <class DiscreteFunctionT>
    void write(DiscreteFunctionT& fnc, std::string name);

  private:
    //- Local typedefs and enums
    typedef typename FunctionSpaceT::GridType GridType;
    typedef typename FunctionSpaceT::DomainType DomainType;
    typedef typename GridType::ctype CoordType;
    typedef typename GridType::template Codim<0>::LeafIterator LeafIterator;
    //typedef typename FunctionSpaceT::IteratorType Iterator;
    //- Local class
    template <int N>
    struct Int2Type {
      enum { value = N };
    };

    //- Local methods
    void init();
    void writeMesh();
    template <class DiscreteFunctionT>
    void writeData(DiscreteFunctionT& df, std::string name);
    template <class T>
    void writeScalar(T x, Int2Type<true>);
    template <class T>
    void writeScalar(T x, Int2Type<false>);

    //- Local data
    //! The space
    const FunctionSpaceT& spc_;
    //! The grid (obtained from spc_ and stored for convenience)
    const GridType& grid_;
    //! The output stream
    // * This will probably be one or more ofstream
    std::ofstream ofs_;
    //! Name of the output file (without extension)
    //! In case of multifile output this is used as the prefix for all output
    //! files
    std::string outName_;
    //! Flag triggering multifile output
    const bool multiFile_;
    //! The data used by the refinement
    RefinementData data_;
    //! The format string (binary or text)
    const std::string format_;
  };

} // end namespace Dune

#include "dxdata.cc"

#endif //DUNE_IO_DXDATA_HH
