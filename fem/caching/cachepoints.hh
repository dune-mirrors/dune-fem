#ifndef DUNE_CACHEPOINTS_HH
#define DUNE_CACHEPOINTS_HH

namespace Dune {

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
