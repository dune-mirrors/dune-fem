#ifndef DUNE_CACHEPOINTPROVIDER_HH
#define DUNE_CACHEPOINTPROVIDER_HH

namespace Dune {

  class PointMapper;

  template <class ct, int dim, int codim>
  class PointProvider {
    // * right?
    typedef CompileTimeChecker<false> Only_codim_1_implementation_exists;
  };

  // * add geo later
  template <class ct, int dim>
  class PointProvider<ct, dim, 1> {
  public:
    typedef FieldVector<ct, dim> PointType;

  public:
    // * return whole vector?
    static const PointMapper& getMapper(const QuadratureType& quad, 
                                        int faceIdx);
    
    static const PointType& getPoints(size_t id);

  private:
    typedef std::map<size_t, std::vector<PointMapper*> MapperStorageType;
    typedef std::map<size_t, std::vector<PointType> PointStorageType;
    typedef typename MapperStorageType::Iterator MapperIteratorType;
    typedef typename PointStorageType::Iterator PointIteratorType;

  private:
    // * define points anew
    static void addEntry(size_t id, GeometryType elemGeo,const Points& points);

  private:
    static MapperStorageType mappers_;
    static PointStorageType points_;
  }

  // * OLD CODE
  // * Find better name
  template <class ct, int dim, int codim>
  class CachePointCalculator 
  {
  private:
    typedef FieldMatrix<ct, dim-codim, dim> TransformationMatrixType;

  private:
    // this trick works as long as p_i+1 - p_0 = e_i on the subentities ref. element + the transformation stays linear
    // use mat_.umtv(xLocal, xRefElem) for the transformation
    void buildTransformationMatrix(int faceIndex) {
      // get global index of the 0 corner of face faceIndex
      int zeroLocal = refElem_.subEntity(faceIndex, codim, 0, dim-codim);
      for (int i = 0; i < dim-codim; ++i) {
        int indexLocal = refElem_.subEntity(faceIndex, codim, i, dim-codim);
        mat_[i] = refElem.position(indexLocal, dim) -
          refElem_.position(zeroLocal, dim);
      }
    }

  private:
    TransformationMatrixType mat_;

  public:
    
  };

}

#endif
