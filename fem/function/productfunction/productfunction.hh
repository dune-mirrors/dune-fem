#ifndef DUNE_PRODUCTFUNCTION_HH
#define DUNE_PRODUCTFUNCTION_HH

//- system includes 
#include <fstream>
#include <rpc/xdr.h>

//- Dune includes 
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/common/array.hh>
#include <dune/common/geometrytype.hh>

//#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/localfunction.hh>
#include <dune/fem/function/common/dofiterator.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune{
//  template <class DiscreteFunctionSpaceImp>    class ParaLocalFunctionAdapt;
 // template <class DofType, class DofArrayType>  class ParaDofIteratorAdapt;
  template <class DiscreteFunctionSpaceImp,class DiscreteFunctionSpace2Imp> class ProductDiscreteFunction;

  template <class DiscreteFunctionSpaceImp>
  struct ProductDiscreteFunctionTraits {
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
  
    typedef ProductDiscreteFunction<DiscreteFunctionSpaceImp,DiscreteFunctionSpaceImp> DiscreteFunctionType;
    // typedef ParaLocalFunctionAdapt<DiscreteFunctionType> LocalFunctionImp;
    // typedef LocalFunctionWrapper< DiscreteFunctionType > LocalFunctionType;

    typedef MutableArray<typename DiscreteFunctionSpaceImp::RangeFieldType> DofArrayType;
      
    typedef typename DofArrayType::DofIteratorType DofIteratorType;
    typedef typename DofArrayType::ConstDofIteratorType ConstDofIteratorType;
  };

//**********************************************************************
//
// --ProductDiscreteFunction 
//
//! This is one special implementation of a discrete function  \f$  u(x,y) = \sum _{i=1} ^{n_x}\sum _{j=1} ^{n_y} u_{ij}\varphi _i(x) \psi _j(y) \f$ for product spaces \f$ X \times Y \f$ using an
//! array for storing the dofs with \f$u_{ij}= u\left[j n_x +i\right]\f$.  
//!
//!
//**********************************************************************
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
class ProductDiscreteFunction 
//  : public DiscreteFunctionDefault<
//  ProductDiscreteFunctionTraits<DiscreteFunctionSpaceType> 
//>
{
public:
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef MutableArray< RangeFieldType > DofStorageType;

private:
//  typedef DiscreteFunctionDefault<
//    ProductDiscreteFunctionTraits<DiscreteFunctionSpaceType> 
//  > DiscreteFunctionDefaultType;
//  friend class DiscreteFunctionDefault< ProductDiscreteFunctionTraits<DiscreteFunctionSpaceType> >;
  
  enum { myId_ = 0};

public:
  typedef typename DiscreteFunctionSpaceType::GridType GridType;

  typedef DofManager<GridType> DofManagerType;
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

  typedef typename DiscreteFunctionSpaceType::Traits::MapperType MapperType;
  typedef typename DiscreteFunctionSpaceType::Traits::RangeFieldType DofType;
  typedef typename DofStorageType::DofIteratorType DofIteratorType;
  typedef typename DofStorageType::ConstDofIteratorType ConstDofIteratorType;

  typedef MemObjectInterface MemObjectInterfaceType;
    
  //! type of this class 
  typedef ProductDiscreteFunction <DiscreteFunctionSpaceType,DiscreteFunctionSpace2Type> DiscreteFunctionType;
  //typedef AdaptiveDiscreteFunction <DiscreteFunctionSpaceType> DiscreteFunctionType;
  typedef DiscreteFunctionType ThisType;
  //  LocalFunctionImp is the implementation 
  // typedef ParaLocalFunctionAdapt < DiscreteFunctionType > LocalFunctionImp;

  // LocalFunctionType is the exported lf type 
  // typedef LocalFunctionWrapper < DiscreteFunctionType> LocalFunctionType;
  
  //! type of class AdaptiveDiscreteFunction
  typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunction1Type;
  typedef typename DiscreteFunction1Type::DofStorageType DofStorage1Type;
  

  typedef DiscreteFunctionSpaceType FunctionSpaceType;
  typedef ProductDiscreteFunctionTraits<DiscreteFunctionSpaceType> Traits;

  /** \brief For ISTL-compatibility */
  typedef FieldVector<DofType,1> block_type;
  
  // typedef LocalFunctionStorage< DiscreteFunctionType > LocalFunctionStorageType; 

  //friend class ParaLocalFunctionAdapt< ThisType> ;
public:

  //! Constructor make Discrete Function for pair of spaces
  ProductDiscreteFunction(const DiscreteFunctionSpaceType& f, const DiscreteFunctionSpace2Type& f2) ;
  
  //! delete stack of free local functions belonging to this discrete
  //! function 
  ~ProductDiscreteFunction ();
       
  // ***********  Interface  *************************
  
  //! return local function \f$  u_{j_0}(x) = \sum _{i=1} ^{n_x}  u\left[j_0 n_x +i\right]\varphi _i(x) \f$ for global dof index \f$j_0\f$ of second space
  inline DiscreteFunction1Type
  localFunction(int dofIndex2) const; 
  
  //! return local function \f$  u_{j_0,en_2}(x) = \sum _{i=1} ^{n_x}  u\left[global(j_0,en_2) n_x +i\right]\varphi _i(x) \f$ for local dof index \f$j_0\f$ of second space and given entity2
  template <class Entity2Type>
  inline DiscreteFunction1Type
  localFunction(const Entity2Type &en2, int dofIndex2) const;
  
  //! local function \f$  u_{loc_2,en_2}(x) = \sum _{i=1} ^{n_x} \left[ \sum_{k=1} ^{\#dof_2}  u\left[global(k,en_2) n_x +i\right]\psi _k (loc_2) \right]\varphi _i(x) \f$ for given entity2 and local coordinate of second space and discrete function of first space
  template < class Entity2Type, class LocalCoord2Type>
  inline
  void localFunction(const Entity2Type &en2, const LocalCoord2Type &loc2,  DiscreteFunction1Type &discFunc) const;

  //! local function \f$  u_{qP_2,en_2}(x) = \sum _{i=1} ^{n_x} \left[ \sum_{k=1} ^{\#dof_2}  u\left[global(k,en_2) n_x +i\right]\psi _k (qP_2) \right]\varphi _i(x) \f$for given entity2, quadrature type and quadrature point number of second space and discrete function of first space
  template < class Entity2Type, class QuadratureType>
  inline
  void localFunction(const Entity2Type &en2, const QuadratureType &quad2, int pointNr, DiscreteFunction1Type &discFunc) const;

  //! points to the first dof of type cc
  inline
  DofIteratorType dbegin ( );
  
  //! points behind the last dof of type cc
  inline
  DofIteratorType dend   ( );

  //! const version of dof iterator  
  inline
  ConstDofIteratorType dbegin ( ) const;
  
  //! const version of dof iterator  
  inline
  ConstDofIteratorType dend   ( ) const;

  //! set all dofs to zero  
  inline
  void clear( );

  //! Add c*org to discrete function
  inline
  void addScaled (const ThisType & org, const RangeFieldType &c); 
      
  //! print all dofs 
  inline
  void print(std::ostream& s) const; 

  //! write data of discrete function to file filename|timestep 
  //! with xdr methods 
  inline
  bool write_xdr(std::string filename) const;

  //! write data of discrete function to file filename|timestep 
  //! with xdr methods 
  inline
  bool read_xdr(std::string filename);

  //! write function data to file filename|timestep in ascii Format
  inline
  bool write_ascii(std::string filename) const;

  //! read function data from file filename|timestep in ascii Format
  inline
  bool read_ascii(std::string filename);

  //! write function data in pgm fromat file
  inline
  bool write_pgm(std::string filename) const;

  //! read function data from pgm fromat file
  inline
  bool read_pgm(std::string filename); 

  //! return name of this discrete function
  std::string name () const { return name_; }

  //! return size fo this discrete function
  int size() const { return dofVec_.size(); }
  
  
  
  //! return pointer to internal array for use of BLAS routines 
  DofType * leakPointer () { return dofVec_.leakPointer();  };
  //! return pointer to internal array for use of BLAS routines 
  const DofType * leakPointer () const { return dofVec_.leakPointer(); };
  
  //! Get access to the related function space
  const DiscreteFunctionSpaceType& space() const { return functionSpace_; }
  //! Get access to the related function space2
  const DiscreteFunctionSpace2Type& space2() const { return functionSpace2_; }


private:  
  // ! return object pointer of type LocalFunctionImp 
  // LocalFunctionImp * newLocalFunctionObject () const;

  // name of this func
  std::string name_;

  // DofManager
  DofManager<GridType>& dm_;

  // MemObject that manages the memory for the dofs of this function
  std::pair<MemObjectInterface*, DofStorageType*> memPair_;
  
  //! array containing the dof of this function, see dofmanager.hh
  //! the array is stored within the mem object 
  DofStorageType & dofVec_;

  // one local function 
  // LocalFunctionType localFunc_;

  //! The related function space
  const DiscreteFunctionSpaceType & functionSpace_;
  //! The related function space2
  const DiscreteFunctionSpace2Type & functionSpace2_;
}; // end class ProductDiscreteFunction 



} // end namespace Dune

#include "productfunctionimp.cc"
#endif
