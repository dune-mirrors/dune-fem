#include "dxdata.hh"

//- System includes
#include <sstream>
#include <string>

//- Dune includes
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>
#include <dune/common/stdstreams.hh>
#include <dune/grid/common/virtualrefinement.hh>

//#include <dune/grid/sgrid.hh>
//#include <dune/common/iteratorfacades.hh>

//- Local includes

namespace Dune {

  //
  //! Provide a dummy value for test output
  //
  //!\param pos Position to calculate the dummy value for
  template <class FieldVector>
  float dummyValue(const FieldVector &pos)
  {
    typedef typename FieldVector::ConstIterator Iterator;
    
    float result = 0;
    Iterator posEnd = pos.end();
    
    int sign = 1;
    for(Iterator i = pos.begin(); i != posEnd; ++i, sign *= -1)
      result += sign * *i * *i;
    
    return result;
  }
  
  //
  //!  Mapping from Dune::GeometryType to dx names
  //
  //!\param type The Dune representation of the geometry type, currently one of
  //!            Dune::triangle, Dune::tetrahedron, Dune::quadrilateral, Dune::hexahedron
  //
  //! Returns a reference to a static std::string.  Throws a Dune::NotImplemented for unsupported
  //! geometry types.
  const std::string &geometryName(GeometryType type)
  {
    static const std::string triangles = "triangles";
    static const std::string tetrahedrons = "tetrahedra";

    static const std::string quads = "quads";
    static const std::string cubes = "cubes";

    switch(type) {
    case triangle:      return triangles;
    case tetrahedron:   return tetrahedrons;

    case quadrilateral: return quads;
    case hexahedron:    return cubes;

    default: DUNE_THROW(NotImplemented, "Dune::GeometryType " << type << " is not supported by dxwriter");
    }
  }

  // //////////////////////////////////////////////////////////////////////

  int refinementLevel()
  {
    // * temp
    return 0;
  }

  template <class Scalar>
  std::vector<Scalar> operator+(const std::vector<Scalar> &V, const Scalar &S)
  {
    std::vector<Scalar> R;
    int size = V.size();
    R.reserve(size);
    for(int i = 0; i < size; ++i)
      R.push_back(V[i] + S);
    return R;
  }

  template<class Stream, class Value>
  Stream &operator<<(Stream &s, const std::vector<Value> &v)
  {
    int size = v.size();
    char *sep = "";
    for(int i = 0; i < size; ++i, sep = "\t")
      s << sep << v[i];
    return s;
  }

  template <class Grid>
  void writeGridDX(Grid &grid, std::ostream& out)
  {
    //- Type definitions
    typedef typename Grid::ctype CoordType;
    enum { dimension = Grid::dimension };

    typedef FieldVector<typename Grid::ctype, Grid::dimensionworld> PosType;
    typedef typename PosType::ConstIterator PosIterator;
    typedef std::vector<PosType> PosList;

    typedef std::vector<int> CellType;
    typedef std::vector<CellType> CellList;

    typedef float ValueType;
    typedef std::vector<ValueType> ValueList;

    typedef typename Grid::LeafIterator LeafIterator;

    typedef typename Grid::ctype CoordType;
    typedef VirtualRefinement<dimension, CoordType> Refinement;

    //- Initialisation
    LeafIterator leafEnd = grid.template leafend<0>();

    // Check if all elements are of the same type
    GeometryType geometryType = 
      grid.template leafbegin<0>()->geometry().type();
    for(LeafIterator i = grid.template leafbegin<0>(); i != leafEnd; ++i)
      if(geometryType != i->geometry().type())
	DUNE_THROW(NotImplemented, "dxwriter supports only uniformly typed grids");
    unsigned int cornersPerCell = 
      grid.template leafbegin<0>()->geometry().corners();

    // Count leaves
    std::cerr << "Counting leaves" << std::endl;
    unsigned int nCells = 0;
    for(LeafIterator i = grid.template leafbegin<0>(); i != leafEnd; ++i)
      ++nCells;
    CellList cellList;
    cellList.reserve(nCells);
    std::cerr << "Done counting leaves" << std::endl;

    // Record refinement levels
    std::cerr << "Preprocessing" << std::endl;
    std::vector<int> refinementLevels;
    refinementLevels.reserve(nCells);
    // Record starting subindex for each cell
    std::vector<int> startingCellIndex;
    std::vector<int> startingVertexIndex;
    startingCellIndex.reserve(nCells);
    startingVertexIndex.reserve(nCells);
    // Count fineCells
    unsigned int nFineCells = 0;
    // Count fineVertices
    unsigned int nFineVertices = 0;
    for(LeafIterator i = grid.template leafbegin<0>(); i != leafEnd; ++i) {
      Refinement &refinement = buildRefinement<dimension, CoordType>(geometryType, geometryType);
      int level = refinementLevel();
      refinementLevels.push_back(level);
      startingVertexIndex.push_back(nFineVertices);
      startingCellIndex.push_back(nFineCells);
      nFineVertices += refinement.nVertices(level);
      nFineCells += refinement.nElements(level);
    }
    std::cerr << "Preprocessing done" << std::endl;
    
    //- Write out mesh info
    // write positions
    std::cerr << "Positions" << std::endl;
    out << "object \"positions\" class array type float category real rank 1 "
        << "shape " << Grid::dimensionworld << " items " << nFineVertices
        << " text data follows" << std::endl;
    int n = 0;
    for(LeafIterator i = grid.template leafbegin<0>(); i != leafEnd; ++n, ++i) {
      typedef typename Refinement::VertexIterator Iterator;
      Refinement &refinement = buildRefinement<dimension, CoordType>(geometryType, geometryType);
      Iterator vertexEnd = refinement.vEnd(refinementLevels[n]);
      for(Iterator j = refinement.vBegin(refinementLevels[n]); j != vertexEnd; ++j)
	out << "\t" << i->geometry().global(j.coords()) << std::endl;
    }
    out << std::endl;
    std::cerr << "Positions done" << std::endl;

    // write connections
    std::cerr << "Connections" << std::endl;
    out << "object \"connections\" class array type int category real rank 1 "
        << "shape " << cornersPerCell << " items " << nFineCells
        << " text data follows" << std::endl;
    n = 0;
    for(LeafIterator i = grid.template leafbegin<0>(); i != leafEnd; ++n, ++i) {
      typedef typename Refinement::ElementIterator Iterator;
      std::cerr << "Building refinement" << std::endl;
      Refinement &refinement = buildRefinement<dimension, CoordType>(geometryType, geometryType);
      std::cerr << "Building elementEnd" << std::endl;
      Iterator elementEnd = refinement.eEnd(refinementLevels[n]);
      for(Iterator j = refinement.eBegin(refinementLevels[n]); j != elementEnd; ++j)
	out << "\t" << (j.vertexIndices()+startingVertexIndex[n]) << std::endl;
    }
    out << "attribute \"ref\" string \"positions\"" << std::endl;
    out << "attribute \"element type\" string \""<< geometryName(geometryType) << "\"" << std::endl;
    out << std::endl;
    std::cerr << "Connections done" << std::endl;
    
    //- Write data
    // write values
    std::cerr << "Values" << std::endl;
    out << "object \"data\" class array type float category real rank 0 "
        << "items " << nFineVertices
        << " text data follows" << std::endl;
    n = 0;
    for(LeafIterator i = grid.template leafbegin<0>(); i != leafEnd; ++n, ++i) {
      typedef typename Refinement::VertexIterator Iterator;
      Refinement &refinement = buildRefinement<dimension, CoordType>(geometryType, geometryType);
      Iterator vertexEnd = refinement.vEnd(refinementLevels[n]);
      for(Iterator j = refinement.vBegin(refinementLevels[n]); j != vertexEnd; ++j)
	out << "\t" << dummyValue(i->geometry().global(j.coords())) << std::endl;
    }
    out << "attribute \"dep\" string \"positions\"" << std::endl;
    out << std::endl;
    std::cerr << "Values done" << std::endl;
    
    // write container
    out << "object \"all\" class field" << std::endl
        << "\tcomponent \"positions\" value \"positions\"" << std::endl
        << "\tcomponent \"connections\" value \"connections\"" << std::endl
        << "\tcomponent \"data\" value \"data\"" << std::endl;
    out << std::endl;
    
    out << "end" << std::endl;
  }

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
    // * what to do with cellList -> comment out (not used)
    //CellList cellList;
    //cellList.reserve(nCells);
    dverb << "Done counting leaves" << std::endl;

    // Record refinement levels
    dverb << "Preprocessing" << std::endl;
    data_.refinementLevels_.reserve(nCells);
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
    // Get refinement object
    refinement_ = 
      &buildRefinement<GridType::dimension, CoordType>(data_.geometryType_, data_.geometryType_);

    for(LeafIterator i = grid_.template leafbegin<0>();
        i != leafEnd; ++i) {
      data_.level_ = refinementLevel();
      data_.refinementLevels_.push_back(data_.level_);
      data_.startingVertexIndex_.push_back(data_.nFineVertices_);
      data_.startingCellIndex_.push_back(data_.nFineCells_);
      data_.nFineVertices_ += refinement_->nVertices(data_.level_);
      data_.nFineCells_ += refinement_->nElements(data_.level_);
    }
  }

  template <class FunctionSpaceT, bool binary>
  void DXWriter<FunctionSpaceT, binary>::writeMesh() {
    dverb << "Positions" << std::endl;
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
      typedef typename Refinement::VertexIterator Iterator;
      Iterator vertexEnd = 
        refinement_->vEnd(data_.refinementLevels_[n]);
      for(Iterator j=refinement_->vBegin(data_.refinementLevels_[n]);
          j != vertexEnd; ++j) {
        x = i->geometry().global(j.coords());
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

    // temporary for later use
    std::vector<int> indices(data_.cornersPerCell_);
    n = 0;
    for(LeafIterator i = grid_.template leafbegin<0>(); 
        i != leafEnd; ++n, ++i) {
      typedef typename Refinement::ElementIterator Iterator;
      dverb << "Building elementEnd" << std::endl;
      
      const int offset = data_.startingVertexIndex_[n];
      Iterator elementEnd = 
        refinement_->eEnd(data_.refinementLevels_[n]);
      for(Iterator j=refinement_->eBegin(data_.refinementLevels_[n]);
          j != elementEnd; ++j) {
        indices = j.vertexIndices() + offset;
        for (int k = 0; k < data_.cornersPerCell_; ++k) {
          writeScalar(static_cast<int>(indices[k]), 
                      Int2Type<binary>());
        }
        if (!binary) {
          ofs_ << "\n";
        }
      }
    }
    // Writing out attributes
    ofs_ << "\nattribute \"ref\" string \"positions\"\n";
    ofs_ << "attribute \"element type\" string \""
        << geometryName(data_.geometryType_) 
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

    RangeType result;
    int n = 0;
    LeafIterator leafEnd = grid_.template leafend<0>();
    for(LeafIterator i = grid_.template leafbegin<0>(); 
        i != leafEnd; ++n, ++i) {
      typedef typename Refinement::VertexIterator Iterator;
      LocalFunctionType lf = df.localFunction(*i);

      Iterator vertexEnd = 
        refinement_->vEnd(data_.refinementLevels_[n]);
      for (Iterator j = 
             refinement_->vBegin(data_.refinementLevels_[n]);
           j != vertexEnd; ++j) {
        lf.evaluate(*i, i->geometry().global(j.coords()), result);
	
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
