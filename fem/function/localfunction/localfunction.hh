#ifndef DUNE_FEM_LOCALFUNCTION_HH
#define DUNE_FEM_LOCALFUNCTION_HH

#include <dune/fem/misc/engineconcept.hh>
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/misc/fieldmatrixhelper.hh>
#include <dune/fem/function/common/function.hh>

namespace Dune
{

  /** \addtogroup LocalFunction
   *
   * On every element from a discrete function the local funtion can be
   * accessed. With the local function one has access to the dof and on the
   * other hand to the base function set of this actual element. Therefore
   * this is called a local function.
   *
   * \remarks The interface for using a LocalFunction is defined by the class
   *          LocalFunction.
   *
   * \{
   */



  /** \class LocalFunction
   *  \brief interface for local functions
   *  
   *  Local functions are used to represend a discrete function on one entity.
   *  The LocalFunctionInterface defines the functionality that can be expected
   *  from such a local function.
   */
  template< class LocalFunctionTraits >
  class LocalFunction
  : public EngineWrapper< typename LocalFunctionTraits :: LocalFunctionImpType,
                          typename LocalFunctionTraits :: LocalFunctionUserType >
  {
  public:
    //! type of the traits
    typedef LocalFunctionTraits Traits;

    //! type of the discrete function space, the local function belongs to
    typedef typename Traits :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    //! type of the local function implementation (engine concept)
    typedef typename Traits :: LocalFunctionImpType LocalFunctionImpType;

    //! type of the user implementation (Barton-Nackman)
    typedef typename Traits :: LocalFunctionUserType LocalFunctionUserType;

  private:
    typedef EngineWrapper< LocalFunctionImpType, LocalFunctionUserType >
      BaseType;

  public:
    //! type of the local function (this type!)
    typedef LocalFunction< Traits > LocalFunctionType;

  private:
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

  public:
    //! type of the entity, the local function lives on
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;

  public:
    //! field type of the domain
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    //! field type of the range
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
    //! type of domain vectors, i.e., type of coordinates
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    //! type of range vectors, i.e., type of function values
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    //! type of Jacobian, i.e., type of evaluated Jacobian matrix
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;

    //! type of base function set  
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

  protected:
    using BaseType :: asImp;

  public:
    /** \brief access to local dofs (read-only)
     *
     *  \param[in]  num  local dof number 
     *  \return reference to dof 
     */
    inline const RangeFieldType &operator[] ( const int num ) const
    {
      return asImp()[ num ];
    }

    /** \brief access to local dofs (read-write)
     *
     *  \param[in]  num  local DoF number
     *  \return reference to DoF
     */
    inline RangeFieldType &operator[] ( const int num )
    {
      return asImp()[ num ];
    }

    /** \brief add another local function to this one
     *
     *  \note The local function to add may be any implementation of a
     *        LocalFunction.
     *
     *  \param[in]  lf  local function to add
     *
     *  \returns a reference to this local function (i.e., *this)
     */
    template< class T >
    inline LocalFunctionType &operator+= ( const LocalFunction< T > &lf )
    {
      asImp() += lf;
      return *this;
    }

    /** \brief assign all DoFs of this local function 
     *
     *  \param[in]  lf  local function to assign DoFs from 
     */
    template< class T >
    inline void assign ( const LocalFunction< T > &lf )
    {
      asImp().assign(lf);
    }

    /** \brief set all DoFs to zero */
    inline void clear ( )
    {
      asImp().clear();
    }

    /** \brief subtract another local function to this one
     *
     *  \note The local function to suctract may be any implementation of a
     *        LocalFunction.
     *
     *  \param[in]  lf  local function to subtract
     *
     *  \returns a reference to this local function (i.e., *this)
     */
    template< class T >
    inline LocalFunctionType &operator-= ( const LocalFunction< T > &lf )
    {
      asImp() -= lf;
      return *this;
    }

    /** \brief add a multiple of another local function to this one
     *
     *  \note The local function to add may be any implementation of a
     *        LocalFunction.
     *
     *  \param[in]  s   scalar factor to scale lf with
     *  \param[in]  lf  local function to add
     *
     *  \returns a reference to this local function (i.e., *this)
     */
    template< class T >
    inline LocalFunctionType &axpy ( const RangeFieldType s,
                                     const LocalFunction< T > &lf )
    {
      asImp().axpy( s, lf );
      return *this;
    }

    /** \brief axpy operation for local function
     *
     *  Denoting the DoFs of the local function by \f$u_i\f$ and the base
     *  functions by \f$\varphi_i\f$, this function performs the following
     *  operation:
     *  \f[
     *  u_i = u_i + factor \cdot \varphi_i( x )
     *  \f]
     *
     *  \param[in]  x       point to evaluate base functions in
     *  \param[in]  factor  axpy factor
     */
    template< class PointType >
    inline void axpy ( const PointType &x,
                       const RangeType &factor )
    {
      asImp().axpy( x, factor );
    }
  
    /** \brief axpy operation for local function
     *
     *  Denoting the DoFs of the local function by \f$u_i\f$ and the base
     *  functions by \f$\varphi_i\f$, this function performs the following
     *  operation:
     *  \f[
     *  u_i = u_i + factor \cdot \nabla\varphi_i( x )
     *  \f]
     *
     *  \param[in]  x       point to evaluate jacobian of base functions in
     *  \param[in]  factor  axpy factor
     */
    template< class PointType >
    inline void axpy ( const PointType &x,
                       const JacobianRangeType &factor)
    {
      asImp().axpy( x, factor );
    }

    /** \brief axpy operation for local function
     *
     *  Denoting the DoFs of the local function by \f$u_i\f$ and the base
     *  functions by \f$\varphi_i\f$, this function performs the following
     *  operation:
     *  \f[
     *  u_i = u_i + factor1 \cdot \varphi_i( x ) + factor2 \cdot \nabla\varphi_i( x )
     *  \f]
     *
     *  \param[in]  x        point to evaluate base functions in
     *  \param[in]  factor1  axpy factor for \f$\varphi( x )\f$
     *  \param[in]  factor2  axpy factor for \f$\nabla\varphi( x )\f$
     */
    template< class PointType >
    inline void axpy ( const PointType &x,
                       const RangeType &factor1,
                       const JacobianRangeType &factor2 )
    {
      asImp().axpy( x, factor1, factor2 );
    }

    /** \brief obtain the base function set for this local function
     *
     *  \returns reference to the base function set
     */
    inline const BaseFunctionSetType &baseFunctionSet () const 
    {
      return asImp().baseFunctionSet();
    }

    /** \brief obtain the entity, this local function lives on
     *
     *  \returns reference to the entity
     */
    inline const EntityType &entity () const
    {
      return asImp().entity();
    }
    
    /** \brief evaluate a partial deriviaive of the local function
     *
     *  \param[in]   diffVariable  vector describing the desired partial
     *                             derivative
     *  \param[in]   x             evaluation point in local coordinates 
     *  \param[out]  ret           value of the function in the given point
     */
    template< int diffOrder, class PointType >
    inline void evaluate ( const FieldVector< deriType, diffOrder > &diffVariable,
                           const PointType &x,
                           RangeType &ret ) const
    {
      asImp().evaluate( diffVariable, x, ret );
    }
   
    /** \brief evaluate the local function
     *
     *  \param[in]   x    evaluation point in local coordinates 
     *  \param[out]  ret  value of the function in the given point
     */
    template< class PointType >
    inline void evaluate ( const PointType &x,
                           RangeType &ret ) const
    {
      asImp().evaluate( x, ret );
    }

    inline void init ( const EntityType &entity )
    {
      asImp().init( entity );
    }

    /** \brief evaluate Jacobian of the local function
     *
     *  \note Though the Jacobian is evaluated on the reference element, the
     *        return value is the Jacobian with respect to the actual entity.
     *
     *  \param[in]   x    evaluation point in local coordinates
     *  \param[out]  ret  Jacobian of the function in the evaluation point
     */
    template< class PointType >
    inline void jacobian ( const PointType &x, 
                           JacobianRangeType &ret ) const
    {
      asImp().jacobian( x, ret );
    }

    /** \brief obtain the number of local DoFs
     *
     *  Obtain the number of local DoFs of this local function. The value is
     *  identical to the number of base functons on the entity.
     *  
     *  \returns number of local DoFs
     */
    inline int numDofs () const 
    {
      return asImp().numDofs();
    }
  };



  //- Forward declarations of Combined Space 
  template< class, int, DofStoragePolicy >
  class CombinedSpace;



  /** \class LocalFunctionDefault
   *  \brief default implementation of a LocalFunction
   */
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  class LocalFunctionDefault
  : public EngineDefault< LocalFunctionImp >
  {
  public:
    //! type of  discrete function space the local function belongs to
    typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

  private:
    typedef EngineDefault< LocalFunctionImp > BaseType;

  private:
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

  public:
    //! type of the entity, the local function lives on
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;

  public:
    //! field type of the domain
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    //! field type of the range
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
    //! type of domain vectors, i.e., type of coordinates
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    //! type of range vectors, i.e., type of function values
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    //! type of Jacobian, i.e., type of evaluated Jacobian matrix
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;

    //! type of base function set  
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

    //! dimension of the domain
    enum { dimDomain = DiscreteFunctionSpaceType :: dimDomain };
    //! dimension of the range
    enum { dimRange = DiscreteFunctionSpaceType :: dimRange };

  protected:
    //! type of entity's geometry
    typedef typename EntityType :: Geometry GeometryType;
    //! type of transposed of geometry's Jacobian Inverse
    typedef FieldMatrix
      < typename GridType :: ctype, GridType :: dimension, GridType :: dimension >
      GeometryJacobianInverseType;

  protected:
    using BaseType :: asImp;

  public:
    /** \copydoc Dune::LocalFunction::operator+=(const LocalFunction<T> &lf) */
    template< class T >
    inline void operator+= ( const LocalFunction< T > &lf );

    /** \copydoc Dune::LocalFunction::operator-=(const LocalFunction<T> &lf) */
    template< class T >
    inline void operator-= ( const LocalFunction< T > &lf );

    template< class T >
    inline void assign ( const LocalFunction< T > &lf );
    inline void clear ( );
    
    /** \copydoc Dune::LocalFunction::axpy(const RangeFieldType s,const LocalFunction<T> &lf) */
    template< class T >
    inline void axpy ( const RangeFieldType s,
                       const LocalFunction< T > &lf );

    /** \copydoc Dune::LocalFunction::evaluate(const FieldVector<deriType,diffOrder> &diffVariable,const PointType &x,RangeType &ret) const */
    template< int diffOrder, class PointType >
    inline void evaluate ( const FieldVector< deriType, diffOrder > &diffVariable,
                           const PointType &x,
                           RangeType &ret ) const;

    /** \copydoc Dune::LocalFunction::evaluate(const PointType &x,RangeType &ret) const */
    template< class PointType >
    inline void evaluate ( const PointType &x,
                           RangeType &ret ) const;

    /** \copydoc Dune::LocalFunction::jacobian(const PointType &x,JacobianRangeType &ret) const */
    template< class PointType >
    inline void jacobian ( const PointType &x,
                           JacobianRangeType &ret ) const;

    /** \copydoc Dune::LocalFunction::axpy(const PointType &x,const RangeType &factor) */
    template< class PointType >
    inline void axpy( const PointType &x,
                      const RangeType &factor );

    /** \copydoc Dune::LocalFunction::axpy(const PointType &x,const JacobianRangeType &factor) */
    template< class PointType >
    inline void axpy( const PointType &x,
                      const JacobianRangeType &factor );

    /** \copydoc Dune::LocalFunction::axpy(const PointType &x,const RangeType &factor1,const JacobianRangeType &factor2) */
    template< class PointType >
    inline void axpy( const PointType &x,
                      const RangeType &factor1,
                      const JacobianRangeType &factor2 );

  protected:
    inline void rightMultiply ( const JacobianRangeType &factor,
                                const DomainType &x,
                                JacobianRangeType &result ) const;
  };



  /** \copydoc Dune::LocalFunctionDefault
   *
   *  Specialised version for CombinedSpaces
   */
  template< class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy,
            class LocalFunctionImp >
  class LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpaceImp, N, policy >, LocalFunctionImp >
  : public EngineDefault< LocalFunctionImp >
  {
  public:
    //! type of  discrete function space the local function belongs to
    typedef CombinedSpace< ContainedFunctionSpaceImp, N, policy >
      DiscreteFunctionSpaceType;

  private:
    typedef EngineDefault< LocalFunctionImp > BaseType;

  private:
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

  public:
    //! type of the entity, the local function lives on
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;

  public:
    //! field type of the domain
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    //! field type of the range
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
    //! type of domain vectors, i.e., type of coordinates
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    //! type of range vectors, i.e., type of function values
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    //! type of Jacobian, i.e., type of evaluated Jacobian matrix
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;

    //! type of base function set  
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

    //! dimension of the domain
    enum { dimDomain = DiscreteFunctionSpaceType :: dimDomain };
    //! dimension of the range
    enum { dimRange = DiscreteFunctionSpaceType :: dimRange };

  protected:
    //! type of entity's geometry
    typedef typename EntityType :: Geometry GeometryType;
    //! type of transposed of geometry's Jacobian Inverse
    typedef FieldMatrix
      < typename GridType :: ctype, GridType :: dimension, GridType :: dimension >
      GeometryJacobianInverseType;

  protected:
    typedef typename DiscreteFunctionSpaceType :: ContainedRangeType
      ScalarRangeType;
    typedef typename DiscreteFunctionSpaceType :: ContainedJacobianRangeType
      ScalarJacobianRangeType;

  protected:
    using BaseType :: asImp;
      
  public:
    /** \copydoc Dune::LocalFunction::operator+=(const LocalFunction<T> &lf) */
    template< class T >
    inline void operator+= ( const LocalFunction< T > &lf );

    /** \copydoc Dune::LocalFunction::operator-=(const LocalFunction<T> &lf) */
    template< class T >
    inline void operator-= ( const LocalFunction< T > &lf );

    template< class T >
    inline void assign ( const LocalFunction< T > &lf );
    inline void clear ( );

    /** \copydoc Dune::LocalFunction::axpy(const RangeFieldType s,const LocalFunction<T> &lf) */
    template< class T >
    inline void axpy ( const RangeFieldType s,
                       const LocalFunction< T > &lf );

    /** \copydoc Dune::LocalFunction::axpy(const PointType &x,const RangeType &factor) */
    template< class PointType >
    inline void axpy ( const PointType &x,
                       const RangeType &factor );

    /** \copydoc Dune::LocalFunction::axpy(const PointType &x,const JacobianRangeType &factor) */
    template< class PointType >
    inline void axpy ( const PointType &x,
                       const JacobianRangeType &factor );

    /** \copydoc Dune::LocalFunction::axpy(const PointType &x,const RangeType &factor1,const JacobianRangeType &factor2) */
    template< class PointType >
    inline void axpy ( const PointType &x,
                       const RangeType &factor1,
                       const JacobianRangeType &factor2 );

    /** \copydoc Dune::LocalFunction::evaluate(const FieldVector<deriType,diffOrder> &diffVariable,const PointType &x,RangeType &ret) const */
    template< int diffOrder, class PointType >
    inline void evaluate ( const FieldVector< deriType, diffOrder > &diffVariable,
                           const PointType &x,
                           RangeType &ret ) const;

    /** \copydoc Dune::LocalFunction::evaluate(const PointType &x,RangeType &ret) const */
    template< class PointType >
    inline void evaluate( const PointType &x,
                          RangeType &ret ) const;

    /** \copydoc Dune::LocalFunction::jacobian(const PointType &x,JacobianRangeType &ret) const */
    template< class PointType >
    inline void jacobian ( const PointType &x,
                           JacobianRangeType &ret ) const;

    inline int numScalarDofs () const
    {
      const int numDofs = asImp().numDofs();
      assert( numDofs % N == 0 );
      return numDofs / N;
    }

  protected:
    inline void rightMultiply ( const JacobianRangeType &factor,
                                const DomainType &x,
                                JacobianRangeType &result ) const;
  };

  /** \} */

} // end namespace Dune 

#include "localfunction_inline.hh"

#endif
