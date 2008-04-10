#ifndef DUNE_PRODUCTFUNCTION_HH
#define DUNE_PRODUCTFUNCTION_HH

//- system includes 
#include <fstream>
#include <rpc/xdr.h>

//- Dune includes 

#include <dune/common/geometrytype.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/dofiterator.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/io/file/xdrio.hh>

namespace Dune{
  template <class DiscreteFunctionSpaceImp,class DiscreteFunctionSpace2Imp> class ProductDiscreteFunction;

  template <class DiscreteFunctionSpaceImp>
  struct ProductDiscreteFunctionTraits {
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
  
    typedef ProductDiscreteFunction<DiscreteFunctionSpaceImp,DiscreteFunctionSpaceImp> DiscreteFunctionType;
    typedef MutableArray<typename DiscreteFunctionSpaceImp::RangeFieldType> DofArrayType;
      
    typedef typename DofArrayType::DofIteratorType DofIteratorType;
    typedef typename DofArrayType::ConstDofIteratorType ConstDofIteratorType;
  };

//**********************************************************************
//! @addtogroup ProductDFunction
// --ProductDiscreteFunction 
//
//! \brief This is one special implementation of a discrete function  \f$  u(x,y) = \sum _{i=1} ^{n_x}\sum _{j=1} ^{n_y} u_{ij}\varphi _i(x) \psi _j(y) \f$ for product spaces \f$ X \times Y \f$ using an
//! array for storing the dofs with \f$u_{ij}= u\left[j n_x +i\right]\f$.  
//!
//!
//**********************************************************************
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
class ProductDiscreteFunction 
{
public:
  //!Type of range field, (usually a float type)
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef MutableArray< RangeFieldType > DofStorageType;

private:
  enum { myId_ = 0};

public:
 //! Type of the underlying grid
  typedef typename DiscreteFunctionSpaceType::GridType GridType;

  typedef DofManager<GridType> DofManagerType;
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

  typedef typename DiscreteFunctionSpaceType::Traits::MapperType MapperType;
  typedef typename DiscreteFunctionSpaceType::Traits::RangeFieldType DofType;
//! Type of the dof iterator used in the discrete function implementation.
  typedef typename DofStorageType::DofIteratorType DofIteratorType;
 //! Type of the constant dof iterator used in the discrete function implementation
  typedef typename DofStorageType::ConstDofIteratorType ConstDofIteratorType;

  typedef MemObjectInterface MemObjectInterfaceType;
    
  //! Type of this class 
  typedef ProductDiscreteFunction <DiscreteFunctionSpaceType,DiscreteFunctionSpace2Type> DiscreteFunctionType;
  typedef DiscreteFunctionType ThisType;
  
  //! Type of class AdaptiveDiscreteFunction
  typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunction1Type;
  typedef typename DiscreteFunction1Type::DofStorageType DofStorage1Type;

  //! Type of the discrete function implementation
  typedef DiscreteFunctionSpaceType FunctionSpaceType;
  typedef ProductDiscreteFunctionTraits<DiscreteFunctionSpaceType> Traits;

  /** \brief For ISTL-compatibility */
  typedef FieldVector<DofType,1> block_type;
  
public:

 /** \brief Constructor storing discrete function spaces 
        \param[in] f discrete function space  
	\param[in] f2 discrete function space2
 */
  ProductDiscreteFunction(const DiscreteFunctionSpaceType& f, const DiscreteFunctionSpace2Type& f2) ;
  
  //! delete stack of free local functions belonging to this discrete
  //! function 
  virtual ~ProductDiscreteFunction ();
       
  // ***********  Interface  *************************
  
  /** \brief returns local function \f$  u_{j_0}(x) = \sum _{i=1} ^{n_x}  u\left[j_0 n_x +i\right]\varphi _i(x) \f$ for global dof index \f$j_0\f$ of second space
      \param[in] dofIndex2 global dof index of second space
      \returns \f$  u_{j_0}(x) = \sum _{i=1} ^{n_x}  u\left[j_0 n_x +i\right]\varphi _i(x) \f$
   */
  inline DiscreteFunction1Type
  localFunction(int dofIndex2) const; 
  
  /** \brief returns local function \f$  u_{j_0,en_2}(x) = \sum _{i=1} ^{n_x}  u\left[global(j_0,en_2) n_x +i\right]\varphi _i(x) \f$ for local dof index \f$j_0\f$ of second space and given entity2
      \param[in] en2 entity of second space
      \param[in] dofIndex2 global dof index of second space
      \returns \f$  u_{j_0,en_2}(x) = \sum _{i=1} ^{n_x}  u\left[global(j_0,en_2) n_x +i\right]\varphi _i(x) \f$
   */
  template <class Entity2Type>
  inline DiscreteFunction1Type
  localFunction(const Entity2Type &en2, int dofIndex2) const;
  
  /** \brief local function \f$  u_{loc_2,en_2}(x) = \sum _{i=1} ^{n_x} \left[ \sum_{k=1} ^{\#dof_2}  u\left[global(k,en_2) n_x +i\right]\psi _k (loc_2) \right]\varphi _i(x) \f$ for given entity2 and local coordinate of second space and discrete function of first space
      \param[in] en2 entity of second space
      \param[in] loc2 local coordinate of second space
      \param[out] discFunc \f$  = \sum _{i=1} ^{n_x} \left[ \sum_{k=1} ^{\#dof_2}  u\left[global(k,en_2) n_x +i\right]\psi _k (loc_2) \right]\varphi _i(x) \f$
   */
  template < class Entity2Type, class LocalCoord2Type>
  inline
  void localFunction(const Entity2Type &en2, const LocalCoord2Type &loc2,  DiscreteFunction1Type &discFunc) const;

  /** \brief local function \f$  u_{qP_2,en_2}(x) = \sum _{i=1} ^{n_x} \left[ \sum_{k=1} ^{\#dof_2}  u\left[global(k,en_2) n_x +i\right]\psi _k (qP_2) \right]\varphi _i(x) \f$ for given entity2, quadrature type and quadrature point number of second space and discrete function of first space
      \param[in] en2 entity of second space
      \param[in] quad2 quadrature of second space
      \param[in] pointNr number of quadrature point
      \param[out] discFunc \f$ = \sum _{i=1} ^{n_x} \left[ \sum_{k=1} ^{\#dof_2}  u\left[global(k,en_2) n_x +i\right]\psi _k (qP_2) \right]\varphi _i(x) \f$
   */
  template < class Entity2Type, class QuadratureType>
  inline
  void localFunction(const Entity2Type &en2, const QuadratureType &quad2, int pointNr, DiscreteFunction1Type &discFunc) const;

 
 /** \brief returns dof iterator pointing to the first degree of freedom of this discrete function 
     \return dof iterator pointing to first dof 
 */
  inline
  DofIteratorType dbegin ( );
  
 /** \brief returns dof iterator pointing behind the last degree of freedom of this discrete function 
     \return dof iterator pointing behind the last dof 
  */
  inline
  DofIteratorType dend   ( );

  /** \brief returns dof iterator pointing to the first degree of freedom of this discrete function 
      \return dof iterator pointing to first dof 
   */
   inline
  ConstDofIteratorType dbegin ( ) const;
  
  /** \brief returns dof iterator pointing behind the last degree of freedom of this discrete function 
      \return dof iterator pointing behind the last dof 
   */  
  inline
  ConstDofIteratorType dend   ( ) const;

   /** \brief set all degrees of freedom to zero
    */  
  inline
  void clear( );

   /** \brief axpy operation
       \param[in] g discrete function that is added 
       \param[in] c scalar value to scale 
    */  
  inline
  void addScaled (const ThisType & g, const RangeFieldType &c); 
      
   /** \brief print all degrees of freedom of this function to stream (for debugging purpose)
       \param[out] s std::ostream (e.g. std::cout)
    */
  inline
  void print(std::ostream& s) const; 

  /** \brief write discrete function to file with given filename using xdr encoding
      \param[in] filename name of file to which discrete function should be written using xdr 
      \return <b>true</b> if operation was successful 
    */  
  virtual bool write_xdr(const std::string filename) const;

  /** \brief read discrete function from file with given filename using xdr decoding
      \param[in] filename name of file from which discrete function should be read using xdr 
      \return <b>true</b> if operation was successful 
   */
  virtual bool read_xdr(const std::string filename);

  
  /** \brief write discrete function to file with given filename using ascii encoding
      \param[in] filename name of file to which discrete function should be written using ascii 
      \return <b>true</b> if operation was successful 
   */
  virtual bool write_ascii(const std::string filename) const;

  /** \brief read discrete function from file with given filename using ascii decoding
      \param[in] filename name of file from which discrete function should be read using ascii 
      \return <b>true</b> if operation was successful 
  */
  virtual bool read_ascii(const std::string filename);

  /** \brief write discrete function to file with given filename using pgm encoding
      \param[in] filename name of file to which discrete function should be written using pgm 
      \return <b>true</b> if operation was successful 
   */
  virtual bool DUNE_DEPRECATED write_pgm(const std::string filename) const;

  /** \brief read discrete function from file with given filename using pgm decoding
      \param[in] filename name of file from which discrete function should be read using pgm 
      \return <b>true</b> if operation was successful 
  */
  virtual bool DUNE_DEPRECATED read_pgm(const std::string filename); 

  /** \brief returns name of discrete function 
      \return string holding name of discrete function 
  */ 
  std::string name () const { return name_; }

  /** \brief returns total number of degrees of freedom, i.e. size of discrete function space 
      \return total number of dofs 
  */ 
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
 
  // name of this func
  std::string name_;

  // DofManager
  DofManager<GridType>& dm_;

  // MemObject that manages the memory for the dofs of this function
  std::pair<MemObjectInterface*, DofStorageType*> memPair_;
  
  //! array containing the dof of this function, see dofmanager.hh
  //! the array is stored within the mem object 
  DofStorageType & dofVec_;

  //! The related function space
  const DiscreteFunctionSpaceType & functionSpace_;
  //! The related function space2
  const DiscreteFunctionSpace2Type & functionSpace2_;
}; // end class ProductDiscreteFunction 



} // end namespace Dune

#include "productfunctionimp.cc"
#endif
