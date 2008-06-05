#include "dxdata.hh"

//- System includes
#include <sstream>
#include <string>

//- Dune includes
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>
#include <dune/common/stdstreams.hh>
#include <dune/fem/space/common/geometryconversion.hh>

//- Local includes

namespace Dune {


  //!  Mapping from Dune::GeometryType to dx names
  //
  //!\param type The Dune representation of the geometry type, currently one of
  //!            Dune::triangle, Dune::tetrahedron, Dune::quadrilateral, Dune::hexahedron
  //
  //! Returns a reference to a static std::string.  Throws a Dune::NotImplemented for unsupported
  //! geometry types.
  const std::string &geometryName(GeometryType type, int dim)
  {
    static const std::string triangles = "triangles";
    static const std::string tetrahedrons = "tetrahedra";

    static const std::string quads = "quads";
    static const std::string cubes = "cubes";

    switch(type) {
    case triangle:      
      return triangles;
    case tetrahedron:
      return tetrahedrons;
    case simplex:
      return (dim == 2 ? triangles : tetrahedrons);

    case quadrilateral:
      return quads;
    case hexahedron: 
      return cubes;
    case cube:
      return (dim == 2 ? quads : cubes);

    default: DUNE_THROW(NotImplemented, "Dune::GeometryType " << type << " is not supported by dxwriter");
    }
  }

  // //////////////////////////////////////////////////////////////////////


  // * the new code follows here
  template <class FunctionSpaceT, bool binary>
  DXWriter<FunctionSpaceT, binary>::
  DXWriter(const FunctionSpaceT& spc, 
           std::string outName,
           bool multiFile) :
    spc_(spc),
    grid_(spc.grid()),
    outName_(outName),
    multiFile_(multiFile),
    format_(binary ? " ieee" : " text")
  {
    init();
    if (multiFile_) {
      ofs_.open(std::string(outName_ + "_mesh.dx").c_str(), 
                std::ios::out | std::ios::trunc);
    } else {
      ofs_.open(std::string(outName_ + ".dx").c_str(), 
                std::ios::out | std::ios::trunc);
    }
    writeMesh();
    if (multiFile_) {
      ofs_ << "end\n";
      ofs_.close();
    }
  }

  template <class FunctionSpaceT, bool binary>
  DXWriter<FunctionSpaceT, binary>::
  ~DXWriter() {
    try {
      if (ofs_.is_open()) {
        ofs_ << "end\n";
        ofs_.close();
      }
    } 
    // Avoid exceptions which leave the destructors
    catch (...) {}
  }

  template <class FunctionSpaceT, bool binary>
  template <class DiscreteFunctionT>
  void 
  DXWriter<FunctionSpaceT, binary>::
  write(DiscreteFunctionT& fnc, std::string name) {
    if (multiFile_) {
      assert(!ofs_.is_open());
      ofs_.open(std::string(outName_ + "_" + name + ".dx").c_str(),
                std::ios::out |std::ios::trunc);
    }
    writeData(fnc, name);
    if (multiFile_) {
      ofs_ << "end\n";
      ofs_.close();
    }
  }

  template <class FunctionSpaceT, bool binary>
  void DXWriter<FunctionSpaceT, binary>::init() {
    LeafIterator leafEnd = grid_.template leafend<0>();
    data_.cornersPerCell_ = 
      grid_.template leafbegin<0>()->geometry().corners();

    // Count leaves
    dverb << "Counting leaves" << std::endl;
    int nCells = 0;
    for(LeafIterator i = grid_.template leafbegin<0>();
        i != leafEnd; ++i)
      ++nCells;
 
    dverb << "Done counting leaves" << std::endl;

    // Record refinement levels
    dverb << "Preprocessing" << std::endl;

    // Record starting subindex for each cell
    data_.startingCellIndex_.reserve(nCells);
    data_.startingVertexIndex_.reserve(nCells);
    // Count fineCells
    data_.nFineCells_ = 0;
    // Count fineVertices
    data_.nFineVertices_ = 0;
    // Element type of the grid
    data_.geometryType_  = 
      grid_.template leafbegin<0>()->geometry().type();

    for(LeafIterator i = grid_.template leafbegin<0>(); i != leafEnd; ++i) {
      data_.startingVertexIndex_.push_back(data_.nFineVertices_);
      data_.startingCellIndex_.push_back(data_.nFineCells_);
      data_.nFineVertices_ += data_.cornersPerCell_;
      data_.nFineCells_ += 1;
    }
  }

  template <class FunctionSpaceT, bool binary>
  void DXWriter<FunctionSpaceT, binary>::writeMesh() {
    dverb << "Positions" << std::endl;

    const int dimension = GridType::dimension;

    if (binary) {
      ofs_ << "data mode lsb binary\n";
    }

    ofs_ << "object \"positions\" class array type float category real rank 1 "
        << "shape " 
        << GridType::dimensionworld 
        << " items " 
        << data_.nFineVertices_
        << format_
        << " data follows" 
        << std::endl;

    LeafIterator leafEnd = grid_.template leafend<0>();
    DomainType x; // temporary for later use
    int n = 0;
    for(LeafIterator i = grid_.template leafbegin<0>(); 
        i != leafEnd; ++n, ++i) {
 
      //     for(Iterator j=refinement_->vBegin(data_.refinementLevels_[n]);
      //  j != vertexEnd; ++j) {
      for (int j = 0; j < data_.cornersPerCell_; ++j) {
        x = i->geometry()[j];
        for (int k = 0; k < GridType::dimensionworld; ++k) {
          writeScalar(static_cast<float>(x[k]), Int2Type<binary>());
        }
        if (!binary) {
          ofs_ << "\n";
        }
      }
    }
    ofs_ << std::endl;
    dverb << "Positions done" << std::endl;

    // write connections
    dverb << "Connections" << std::endl;
    ofs_ << "\nobject \"connections\" class array type int category real rank 1 "
        << "shape " 
        << data_.cornersPerCell_ 
        << " items " 
        << data_.nFineCells_
        << format_
        << " data follows" 
        << std::endl;

    n = 0;
    for(LeafIterator i = grid_.template leafbegin<0>(); 
        i != leafEnd; ++n, ++i) {
      
      const int offset = data_.startingVertexIndex_[n];

      for (int k = 0; k < data_.cornersPerCell_; ++k) {
        writeScalar(static_cast<int>(offset + k), 
                    Int2Type<binary>());
      }
      if (!binary) {
        ofs_ << "\n";
      }
      
    }
    // Writing out attributes
    ofs_ << "\nattribute \"ref\" string \"positions\"\n";
    ofs_ << "attribute \"element type\" string \""
        << geometryName(data_.geometryType_, dimension) 
        << "\"\n"; 
    ofs_ << std::endl;
    dverb << "Connections done" << std::endl;
  }

  template <class FunctionSpaceT, bool binary>
  template <class DiscreteFunctionT>
  void 
  DXWriter<FunctionSpaceT, binary>::
  writeData(DiscreteFunctionT& df, std::string name) {
    //- Local typedefs
    typedef typename DiscreteFunctionT::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteFunctionT::RangeType RangeType;

    //- Actual code
    const int shape = RangeType::size;
    const int rank = (shape > 1) ? 1 : 0;
    std::string dataName = name + "_data";

    dverb << "Values" << std::endl;
    ofs_ << "object \"" 
        << dataName
        << "\" class array type float category real rank "
        << rank;
    if (rank > 0) {
      ofs_ << " shape "
          << shape;
    }
    ofs_ << " items "
        << data_.nFineVertices_
        << format_
        << " data follows"
        << std::endl;

    RangeType result(0.0);
    int n = 0;
    LeafIterator leafEnd = grid_.template leafend<0>();
    for(LeafIterator i = grid_.template leafbegin<0>(); 
        i != leafEnd; ++n, ++i) {
      LocalFunctionType lf = df.localFunction(*i);

      for (int j = 0; j < data_.cornersPerCell_; ++j) {
        lf.evaluateGlobal(*i, i->geometry()[j], result);
        //lf.evaluate(*i, i->geometry().global(j.coords()), result);
	
        for (int k =0; k < shape; ++k) {
          writeScalar(static_cast<float>(result[k]), 
                      Int2Type<binary>());
        }
        if (!binary) {
          ofs_ << "\n";
        }    
      }
    }

    // Temporary variables for mesh name
    std::string meshName = outName_ + "_mesh.dx";
    std::string dxMeshRef = "file \"" + meshName + "\",";
    ofs_ << "\nattribute \"dep\" string \"positions\"\n";
    ofs_ << std::endl;
    dverb << "Values done" << std::endl;    
    
    // write field defintion
    ofs_ << "object \""
         << name
         << "\" class field\n"
         << "\tcomponent \"positions\" value ";
    if (multiFile_) {
      ofs_ << dxMeshRef;
    }
    ofs_ << "\"positions\"\n"
         << "\tcomponent \"connections\" value ";
    if (multiFile_) {
      ofs_ << dxMeshRef;
    }
    ofs_ << "\"connections\"\n"
         << "\tcomponent \"data\" value \""
         << dataName
         << "\"\n";
    ofs_ << std::endl;
  }
  
  template <class FunctionSpaceT, bool binary>
  template <class T>
  inline void 
  DXWriter<FunctionSpaceT, binary>::
  writeScalar(T x, Int2Type<true> dummy) {
    ofs_.write(reinterpret_cast<char*>(&x), sizeof(T));
  }

  template <class FunctionSpaceT, bool binary>
  template <class T>
  inline void
  DXWriter<FunctionSpaceT, binary>::
  writeScalar(T x, Int2Type<false> dummy) {
    ofs_ << "\t" << x;
  }
} // namespace Dune
