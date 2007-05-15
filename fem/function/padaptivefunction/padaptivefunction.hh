#ifndef DUNE_PADAPTIVEFUNCTION_HH
#define DUNE_PADAPTIVEFUNCTION_HH

//- system includes 
#include <fstream>
#include <rpc/xdr.h>

//- Dune includes 
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/common/array.hh>
#include <dune/common/geometrytype.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/localfunction.hh>
#include <dune/fem/function/common/dofiterator.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune{
  template <class DiscreteFunctionSpaceImp>    class ParaLocalFunctionAdapt;
  template <class DofType, class DofArrayType>  class ParaDofIteratorAdapt;
  template <class DiscreteFunctionSpaceImp,class DiscreteFunctionSpace2Imp> class PAdaptiveDiscreteFunction;

  template <class DiscreteFunctionSpaceImp>
  struct PAdaptiveDiscreteFunctionTraits {
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
  
    typedef PAdaptiveDiscreteFunction<DiscreteFunctionSpaceImp,DiscreteFunctionSpaceImp> DiscreteFunctionType;
    typedef ParaLocalFunctionAdapt<DiscreteFunctionType> LocalFunctionImp;
    typedef LocalFunctionWrapper< DiscreteFunctionType > LocalFunctionType;

    typedef typename DofArray<
      typename DiscreteFunctionSpaceImp::RangeFieldType
    >::DofIteratorType DofIteratorType;
    typedef typename DofArray<
      typename DiscreteFunctionSpaceImp::RangeFieldType
    >::ConstDofIteratorType ConstDofIteratorType;
  };

//**********************************************************************
//
//  --PAdaptiveDiscreteFunction 
//
//! this is one special implementation of a discrete function using an
//! array for storing the dofs.  
//!
//**********************************************************************
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
class PAdaptiveDiscreteFunction 
  : public DiscreteFunctionDefault<
  PAdaptiveDiscreteFunctionTraits<DiscreteFunctionSpaceType> 
>
{
public:
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef DofArray< RangeFieldType > DofStorageType;

private:
  typedef DiscreteFunctionDefault<
    PAdaptiveDiscreteFunctionTraits<DiscreteFunctionSpaceType> 
  > DiscreteFunctionDefaultType;
  friend class DiscreteFunctionDefault< PAdaptiveDiscreteFunctionTraits<DiscreteFunctionSpaceType> >;
  
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
  typedef PAdaptiveDiscreteFunction <DiscreteFunctionSpaceType,DiscreteFunctionSpace2Type> DiscreteFunctionType;
  //typedef AdaptiveDiscreteFunction <DiscreteFunctionSpaceType> DiscreteFunctionType;
  typedef DiscreteFunctionType ThisType;
  //! LocalFunctionImp is the implementation 
  typedef ParaLocalFunctionAdapt < DiscreteFunctionType > LocalFunctionImp;

  //! LocalFunctionType is the exported lf type 
  typedef LocalFunctionWrapper < DiscreteFunctionType> LocalFunctionType;
  
  typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunction1Type;
  typedef typename DiscreteFunction1Type::DofStorageType DofStorage1Type;
  

  typedef DiscreteFunctionSpaceType FunctionSpaceType;
  typedef PAdaptiveDiscreteFunctionTraits<DiscreteFunctionSpaceType> Traits;

  /** \brief For ISTL-compatibility */
  typedef FieldVector<DofType,1> block_type;
  
  typedef LocalFunctionStorage< DiscreteFunctionType > LocalFunctionStorageType; 

  friend class ParaLocalFunctionAdapt< ThisType> ;
public:

  //! Constructor make Discrete Function for pair of spaces
  PAdaptiveDiscreteFunction(const DiscreteFunctionSpaceType& f, const DiscreteFunctionSpace2Type& f2) ;
  
  //! delete stack of free local functions belonging to this discrete
  //! function 
  ~PAdaptiveDiscreteFunction ();
 
  DiscreteFunctionType & argument    () DUNE_DEPRECATED{ return *this; }
  const DiscreteFunctionType & argument () const DUNE_DEPRECATED{ return *this; }
  DiscreteFunctionType & destination () DUNE_DEPRECATED{ return *this; }
       
  // ***********  Interface  *************************
  //! return empty object of a local fucntion 
  //! old function, will be removed soon 
  
  //! return local function for global dof index dofIndex2 of space two
  inline DiscreteFunction1Type
  localFunction(int dofIndex2) const; 
  
  //! return local function for local dof index dofIndex2 of space two and given entity2
  template <class Entity2Type>
  inline DiscreteFunction1Type
  localFunction(const Entity2Type &en2, int dofIndex2) const;
  
  //! local function for given entity2 and local coordinate of second space and discrete function of first space
  template < class Entity2Type, class LocalCoord2Type>
  inline
  void localFunction(const Entity2Type &en2, const LocalCoord2Type &loc2,  DiscreteFunction1Type &discFunc) const;

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

  //! set all dof to value x 
  inline
  void set( RangeFieldType x )DUNE_DEPRECATED; 

  //! \todo Please do me!
  inline
  void addScaled (const PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type> & g,
      const RangeFieldType &scalar); 
      
      
      
      
/*
  //! \todo Please do me!
  template <class EntityType>
  inline
  void addScaledLocal (EntityType &en, 
      const PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type> & g,
      const RangeFieldType &scalar) DUNE_DEPRECATED; 
 
  //! add g to this on local entity
  template <class EntityType>
  inline
  void addLocal (EntityType &it, 
      const PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type> & g) DUNE_DEPRECATED; 
  
  //! add g to this on local entity
  template <class EntityType>
  inline
  void subtractLocal (EntityType &it, 
                      const PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type> & g)DUNE_DEPRECATED; 

  //! \todo Please do me!
  template <class EntityType>
  inline
  void setLocal (EntityType &it, const RangeFieldType &scalar) DUNE_DEPRECATED;
 */ 
 
 
 
 
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
  
  
   //! apply this*x = ret 
  void precondition(const DofType * x, DofType * ret) const
  {
    const int vecsize = size();
    for(int i=0; i<vecsize; ++i) 
    {
      ret[i] = x[i]*dofVec_[i];
    }
  }

  void multOEM (const DofType * x, DofType * ret) const 
  {
    precondition(x,ret);
  }

  //! return false as this is left precondition
  bool rightPrecondition () const { return false; }
  
  //! return pointer to internal array for use of BLAS routines 
  DofType * leakPointer () { return dofVec_.leakPointer();  };
  //! return pointer to internal array for use of BLAS routines 
  const DofType * leakPointer () const { return dofVec_.leakPointer(); };
  
  //! Get access to the related function space
  const DiscreteFunctionSpace2Type& space2() const { return functionSpace2_; }


private:  
  //! return object pointer of type LocalFunctionImp 
  LocalFunctionImp * newLocalFunctionObject () const;

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
  const DiscreteFunctionSpace2Type & functionSpace2_;
}; // end class PAdaptiveDiscreteFunction 


//**************************************************************************
//
//  --ParaLocalFunctionAdapt
//
//! Implementation of the local functions 
//
//**************************************************************************
template < class DiscreteFunctionType> 
class ParaLocalFunctionAdapt 
: public LocalFunctionDefault <typename DiscreteFunctionType::DiscreteFunctionSpaceType ,
  ParaLocalFunctionAdapt <DiscreteFunctionType>  >
{
public:
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef ParaLocalFunctionAdapt<DiscreteFunctionType> MyType;
  typedef DiscreteFunctionType DiscFuncType;
  friend class LocalFunctionWrapper< DiscreteFunctionType >; 

  enum { dimrange = DiscreteFunctionSpaceType::DimRange };
  enum { DimRange = DiscreteFunctionSpaceType::DimRange };

  friend class PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType>;
  typedef typename DiscreteFunctionSpaceType::Traits::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType::Traits::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::Traits::RangeType RangeType;
  typedef RangeFieldType DofType;
  typedef typename DiscreteFunctionSpaceType::Traits::JacobianRangeType JacobianRangeType;

  typedef typename DiscFuncType::DofStorageType DofStorageType;

  typedef typename DiscreteFunctionSpaceType::GridType GridType;
  typedef typename GridType::Traits::LocalIdSet LocalIdSetType;
  typedef typename LocalIdSetType::IdType IdType;
  //! Constructor 

public:
  //! Constructor taking DiscreteFunctionSpace and DofStorageType
  inline
  ParaLocalFunctionAdapt (const DiscreteFunctionType &df);

  //! Destructor 
  inline
  ~ParaLocalFunctionAdapt ();

  //! return name of local function 
  std::string name () const;

  //! access to dof number num, all dofs of the dof entity
  inline
  RangeFieldType & operator [] (int num);
  
  //! access to dof number num, all dofs of the dof entity
  inline
  const RangeFieldType & operator [] (int num) const;

  //! return number of degrees of freedom 
  inline
  int numberOfDofs () const;

  //! return number of degrees of freedom 
  inline
  int numDofs () const;

/*   #if OLDFEM
  //! sum over all local base functions 
  template <class EntityType> 
  inline
  void evaluate (EntityType &en, const DomainType & x, RangeType & ret) const ;
  
  template <class EntityType>
  inline
  void evaluateLocal(EntityType &en, const DomainType & x, RangeType & ret) const ;
  //! sum over all local base functions evaluated on given quadrature point
  template <class EntityType, class QuadratureType> 
  inline
  void evaluate (EntityType &en, QuadratureType &quad, int quadPoint , RangeType & ret) const;
   #endif
  //! sum over all local base functions evaluated on given quadrature point
  template <class EntityType, class QuadratureType> 
  inline
  void jacobian (EntityType &en, QuadratureType &quad, int quadPoint , JacobianRangeType & ret) const;
  #if OLDFEM
  template <class EntityType>
  inline
  void jacobianLocal(EntityType& en, const DomainType& x, JacobianRangeType& ret) const DUNE_DEPRECATED;
  #endif
  template <class EntityType>
  inline
  void jacobian(EntityType& en, const DomainType& x, JacobianRangeType& ret) const;
 */ 
/*******************************************************/
 /* //! evaluate the local function in a given local coordinate
  //! \param x: local coordinate
  //! \param ret: localfunction(x)
  inline
  void evaluate(const DomainType & x, RangeType & ret) const ;
  //! sum over all local base functions evaluated on given quadrature point
  template <class QuadratureType> 
  inline
  void evaluate (QuadratureType &quad, int quadPoint , RangeType & ret) const;
  */
  inline
  void assign(int numDof, const RangeType& dofs);
  
  inline
  const BaseFunctionSetType& baseFunctionSet() const;
  inline
  const BaseFunctionSetType& getBaseFunctionSet() const DUNE_DEPRECATED {return baseFunctionSet();}
protected: 
  //! update local function for given Entity  
  template <class EntityType > 
  inline
  void init ( const EntityType &en ) const;

  //! Forbidden! Would wreak havoc
  ParaLocalFunctionAdapt(const ParaLocalFunctionAdapt&);
  inline
  MyType& operator= (const MyType& other);

  //! needed once 
  mutable RangeType tmp_;
  mutable DomainType xtmp_;

  //! needed once 
  mutable JacobianRangeType tmpGrad_;

  //! diffVar for evaluate, is empty 
  const DiffVariable<0>::Type diffVar;

  //! number of all dofs 
  mutable int numOfDof_;

  //! DiscreteFunction we belong to
  const DiscreteFunctionType & df_; 
  
  //! the corresponding function space which provides the base function set
  const DiscreteFunctionSpaceType& fSpace_;
  
  //! Array holding pointers to the local dofs 
  mutable Array < RangeFieldType * > values_;

  //! dofVec from all levels of the discrete function 
  DofStorageType & dofVec_;

  //! the corresponding base function set 
  mutable const BaseFunctionSetType* baseSet_;
  
  //! is it initialised?
  mutable bool init_;

  //! actual geometry type 
  mutable GeometryType geoType_;

}; // end ParaLocalFunctionAdapt 

} // end namespace Dune

#include "padaptivefunctionimp.cc"
#endif
