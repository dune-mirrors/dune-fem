#ifndef DUNE_DISCRETEMODEL_HH
#define DUNE_DISCRETEMODEL_HH

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/fvector.hh>
#include "selection.hh"

namespace Dune {

  /**
   * @brief Interface class for problem definition in the LDG context.
   *
   * The problem interface prescribes the methods needed by the implementation
   * of local DG passes. Users need to derive from either this class or the
   * DiscreteModelDefault class. The methods provided by this class are used by
   * LocalDGPass (and possibly other classes) to solve for the space 
   * discretisation of the equation d_t u + div f(u) = s(u), where f
   * represents a flux function and s a source term. Even non-conservative
   * fluxes can be treated (for details see paper by Dedner et al).
   * \note The definition of f(u) differs from the usual definition for 
   * conservation laws in the leading sign.
   */
  template <class DiscreteModelTraits>
  class DiscreteModelInterface 
  {
  public:
    //! Traits class defined by the user
    typedef DiscreteModelTraits Traits;
    //! Implementation type for Barton-Nackman trick
    typedef typename Traits::DiscreteModelType DiscreteModelType;
    //! Function space (we take the discrete one even though the continuous one
    //! would suffice...)
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! Coordinate type (world coordinates)
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    //! Vector type of the function space's range
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    //! Vector type of the function space's range field type 
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
    //! Matrix type of the function space's jacobian
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    //! Type of GridPart 
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType; 
    //! Type of the grid
    typedef typename GridPartType::GridType GridType;
    //! Intersection iterator of the grid
    typedef typename GridPartType::IntersectionIteratorType IntersectionIterator;
    //! Element (codim 0 entity) of the grid
    typedef typename GridType::template Codim<0>::Entity EntityType;

    //! dimRange 
    enum { dimRange = DiscreteFunctionSpaceType :: DimRange };

    // mass factor type 
    template <class MatrixType> struct ExtractMatrix;

    template <template <class,int,int> class MatrixType,
              class ctype, int dimD, int dimR 
              > 
    struct ExtractMatrix<MatrixType<ctype,dimD,dimR> >
    {
      typedef MatrixType<RangeFieldType,dimRange,dimRange> Type;
    };

    //! type of mass factor 
    typedef typename ExtractMatrix<JacobianRangeType> :: Type MassFactorType ;

  public:
    //! Returns true if problem has a flux contribution.
    //! If you program this method to return true, make sure to implement
    //! numericalFlux and analyticalFlux as well.
    bool hasFlux() const { return asImp().hasFlux(); }

    //! Returns true if problem has a source term.
    //! If you program this method to return true, make sure to implement
    //! method source well.
    bool hasSource() const { return asImp().hasSource(); }

    //! Returns true if problem has a mass matrix factor disctinct from the identity
    bool hasMass() const { return asImp().hasMass(); }

    //! \brief Computes the numerical flux at a cell interface.
    //! Interface method for calculating the numerical flux at cell interfaces.
    //! The method has distinct return values for the contribution to the inner
    //! and outer element, since the framework allows for non-conservative
    //! terms. The non-conservative terms don't get evaluated over separate
    //! methods, but need to be implemented using the flux and source methods.
    //! \param it Iterator of the face in consideration.
    //! \param time Global time.
    //! \param x Local coordinate of the point where the numerical flux gets
    //! evaluated.
    //! \param uLeft Tuple of the states on the inner side of the intersection.
    //! \param uRight Tuple of the states on the outer side of the
    //! intersection.
    //! \param gLeft The numerical flux contribution to the inner element 
    //! (return value).
    //! \param gRight The numerical flux contribution to the outer element
    //! (return value).
    //! \return Speed of the fastest wave. This information is used to compute
    //! the maximum admissible timestep size.
    template <class ArgumentTuple, class FaceDomainType>
    double numericalFlux(IntersectionIterator& it,
                         const double time, 
                         const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { 
      CHECK_INTERFACE_IMPLEMENTATION( asImp().numericalFlux(it, time, x, 
                                   uLeft, uRight, gLeft, gRight) );
      return asImp().numericalFlux(it, time, x, 
                                   uLeft, uRight, gLeft, gRight); 
    }

    //! \brief Computes the flux at the boundary
    //! Special kind of numerical flux. The intersection iterator provides
    //! the necessary information to identify the type of boundary.
    //! \param it Iterator of the face in consideration.
    //! \param time Global time.
    //! \param x Local coordinate of the point where the numerical flux gets
    //! evaluated.
    //! \param uLeft Tuple of the states on the inner side of the intersection.
    //! \param gLeft The numerical flux contribution to the inner element 
    //! (return value).
    //! \return Speed of the fastest wave. This information is used to compute
    //! the maximum admissible timestep size.
    template <class ArgumentTuple, class FaceDomainType>
    double boundaryFlux(IntersectionIterator& it,
                        const double time, 
                        const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    { 
      CHECK_INTERFACE_IMPLEMENTATION( 
        asImp().boundaryFlux(it, time, x, uLeft, gLeft) );
      return asImp().boundaryFlux(it, time, x, uLeft, gLeft);
    }

    //! \brief Computes the analytical flux of the problem.
    //! Analytical flux of the problem.
    //! \param en Entity where the flux is to be evaluated
    //! \param time Global time.
    //! \param x Reference element coordinate where the flux is evaluated
    //! \param u Tuple of the states of the previous passes and of the global
    //! argument.
    //! \param f The analytical flux (return value)
    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        const double time, 
                        const DomainType& x,
                        const ArgumentTuple& u, 
                        JacobianRangeType& f) 
    { 
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
              asImp().analyticalFlux(en, time, x, u, f) ); 
    }

    //! \brief Implements the source term of the problem.
    //! Source term contribution. The source term can contain parts of the
    //! non-conservative contributions.
    //! \param en Entity where the flux is to be evaluated
    //! \param time Global time.
    //! \param x Reference element coordinate where the flux is evaluated
    //! \param u Tuple of the states of the previous passes and of the global
    //! argument.
    //! \param jac Tuple of the gradients of the precedings passes' results 
    //! (needed for the non-conservative contributions)
    //! \param s The source contribution (return value).
    template <class ArgumentTuple, class JacobianTuple>
    void source(EntityType& en, 
                const double time, 
                const DomainType& x,
                const ArgumentTuple& u, 
                const JacobianTuple& jac, 
                RangeType& s) 
    { 
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
          asImp().source(en, time, x, u, jac, s) ); 
    }

    //! \brief Implements the mass factor term of the problem.
    //! Mass term contribution. 
    //! \param en Entity where the mass is to be evaluated
    //! \param time Global time.
    //! \param x Reference element coordinate where the mass is evaluated
    //! \param u Tuple of the states of the previous passes and of the global
    //! \param m The mass contribution (return value).
    //! default implementation sets this factor to 1.0 
    template <class ArgumentTuple>
    void mass(const EntityType& en, 
              const double time, 
              const DomainType& x,
              const ArgumentTuple& u, 
              MassFactorType& m)
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
          asImp().mass(en, time, x, u, m) ); 
    }

    //! \brief Passes the active entity to the model.
    //! This can be used, to set local functions required as data function
    //! in the model.
    //! \param en active Entity 
    void setEntity(EntityType& en) 
    { 
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().setEntity(en) ); 
    }

    //! \brief Passes the active neigbor entity to the model.
    //! This can be used, to set local functions required as data functions 
    //! in the model.
    //! \param nb active neighbor Entity 
    void setNeighbor(EntityType& nb) 
    { 
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().setNeighbor(nb) ); 
    }

  protected:
    DiscreteModelType& asImp() { return static_cast<DiscreteModelType&>(*this); }
    const DiscreteModelType& asImp() const { return static_cast<const DiscreteModelType&>(*this); }
  };
  
  //! Default implementation of the DiscreteModelInterface where methods for 
  //! the fluxes and the source term do nothing, so that the user needn't
  //! implement them if not needed.
  //! N1, ..., N9 are passIds on which model depends
  template< class DiscreteModelTraits
            , int N1 = -1 
            , int N2 = -1 
            , int N3 = -1 
            , int N4 = -1 
            , int N5 = -1 
            , int N6 = -1 
            , int N7 = -1 
            , int N8 = -1 
            , int N9 = -1 
            >
  class DiscreteModelDefault : 
    public DiscreteModelInterface<DiscreteModelTraits> 
  {
    typedef DiscreteModelInterface<DiscreteModelTraits> BaseType;
  public:
    typedef DiscreteModelTraits Traits;
    typedef typename Traits::DiscreteModelType DiscreteModelType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType; 
    typedef typename GridPartType::GridType GridType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename GridPartType::IntersectionIteratorType IntersectionIterator;

    typedef typename BaseType :: MassFactorType MassFactorType;

    //! Selector for data tuple to use as arguments for all methods;
    //! this fixes the template type ArgumentTuple.
    //! If this discrete model is used for a pass n+1, i.e., follwoing
    //! passes p0,p1,..,pn then the return type of pass i (i=0,..,n)
    //! can be used by adding the interger number i in the Selector.
    //! Assume the following: \$ u_{n+1} = p_{n+1}(u_n,u_{n-1},..,u_1,u_0) \$
    //! where $u_0=u$ is the global argument of the combined passes.
    //! If \$ p_{n+1} \$ only depends on \$ u_0,u_2,u_n \$ then the
    //! following selector can be used: \c Selector<0,n-1,1>. Then
    //! ArgumentTuple is now filled with the values of these three
    //! functions and can be accessed by...
    //! Other way of filling the ArgumentTuple with corresponding pass
    //! results is when one uses passIds. In this case if \$ u_{n+1} \$
    //! depends on the passes with following passIds: firstPassId , passId2 
    //! , passId5 then the desired Selector is 
    //! Selector< firstPassId , passId2 , passId5 > ...
    //! If there's no SelectorType in user-implemented DiscreteModel then
    //! this Selector is used. Therefore it's good to pass passIds to this class
    //! and avoid writing SelectorType in user-implemented DiscreteModel.
    //! The point where a user specifies what's going to be in the Selector is
    //! in the template declaration of the DiscreteModel where one names
    //! passIds necessary for this DiscreteModel
    typedef Selector< N1 , N2 , N3 , N4 , N5 , N6 , N7 , N8 , N9 > SelectorType;
  public:
    /** \copydoc Dune::DiscreteModelInterface::hasFlux() const
     *
     *  The default implementation always returns false
     */
    inline bool hasFlux () const
    {
      return false;
    }
    
    /** \copydoc Dune::DiscreteModelInterface::hasSource() const
     *
     *  The default implementation always returns false
     */
    inline bool hasSource () const
    {
      return false;
    }
    
    /** \copydoc Dune::DiscreteModelInterface::hasMass() const
     *
     *  The default implementation always returns false
     */
    inline bool hasMass () const
    {
      return false;
    }
    
    //! Empty implementation that fails if problem claims to have a flux
    //! contribution.
    template <class ArgumentTuple, class FaceDomainType>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { 
      assert(!this->asImp().hasFlux()); 
      gLeft = 0.0;
      gRight = 0.0;
      return 0.0;
    }

    //! Empty implementation that fails if problem claims to have a flux
    //! contribution.
    template <class ArgumentTuple, class FaceDomainType>
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    {
      assert(!this->asImp().hasFlux());
      gLeft = 0.0;
      return 0.0;
    }

    //! Empty implementation that fails if problem claims to have a flux
    //! contribution.
    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    { 
      assert(!this->asImp().hasFlux()); 
      f = 0.0;
    }

    //! Empty implementation that fails if problem claims to have a source 
    //! term.
    template <class ArgumentTuple, class JacobianTuple>
    void source(EntityType& en, 
                double time, const DomainType& x,
                const ArgumentTuple& u, const JacobianTuple& jac, 
                RangeType& s)
    { 
      assert(!this->asImp().hasSource()); 
      s = 0.0;
    }

    //! empty implementation for mass factor 
    //! default implementation sets this factor to 1.0 
    template <class ArgumentTuple>
    void mass(const EntityType& en, 
              const double time, 
              const DomainType& x,
              const ArgumentTuple& u, 
              MassFactorType& m)
    {
      enum { rows = MassFactorType :: rows };
      enum { cols = MassFactorType :: cols };
      // default implementation sets mass factor to identity 
      for(int i=0; i<rows; ++i) 
      {
        // set diagonal to 1 
        m[i][i] = 1;
        // and all other values to 0
        for(int j=i+1; j<cols; ++j) 
          m[i][j] = m[j][i] = 0;
      }
    }

    //! Empty implementation 
    void setEntity(EntityType& en)
    { }

    //! Empty implementation 
    void setNeighbor(EntityType& nb)
    { }
  };

  //! Default implementation of the DiscreteModelInterface where methods for 
  //! the fluxes and the source term do nothing, so that the user needn't
  //! implement them if not needed.
  template <class DiscreteModelTraits
            , int N1 = -1 
            , int N2 = -1 
            , int N3 = -1 
            , int N4 = -1 
            , int N5 = -1 
            , int N6 = -1 
            , int N7 = -1 
            , int N8 = -1 
            , int N9 = -1 
            >
  class DiscreteModelDefaultWithInsideOutSide : 
    public DiscreteModelDefault<DiscreteModelTraits,N1,N2,N3,N4,N5,N6,N7,N8,N9> 
  {
  public:
    typedef DiscreteModelTraits Traits;
    typedef typename Traits::DiscreteModelType DiscreteModelType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType; 
    typedef typename GridPartType::GridType GridType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;

  public:
    //! \brief default constructor 
    DiscreteModelDefaultWithInsideOutSide() 
      : enVol_(-1.0) , nbVol_(-1.0) , en_(0) , nb_(0) 
    {}

    //! set entity and get volume  
    //! \brief method setting pointer of inside entity and getting volume 
    //! \param[in] en reference to inside entity 
    void setEntity(EntityType& en)
    { 
      en_ = &en;
      enVol_ = en.geometry().volume();
    }

    //! \brief method seting pointer of outside entity and getting volume 
    //! \param[in] nb reference to outside entity 
    void setNeighbor(EntityType& nb)
    { 
      nb_ = &nb;
      nbVol_ = nb.geometry().volume();
    }
    
    //! \brief method returning reference to inside entity 
    //! \return reference to inside entity 
    const EntityType & inside() const  
    {
      assert( en_ );
      return *en_;
    }
    
    //! \brief method returning reference to outside entity 
    //! \return reference to outside entity 
    const EntityType & outside() const 
    {
      assert( nb_ );
      return *nb_;
    }

    //! \brief return volume of entity
    double enVolume() const 
    { 
      assert(enVol_ > 0.0);
      return enVol_; 
    }

    //! \brief return volume of neighbor
    double nbVolume() const 
    {
      assert( nbVol_ > 0.0 );
      return nbVol_;
    }
  private:
    double enVol_;
    double nbVol_;

    EntityType* en_;
    EntityType* nb_;
  };

}  // end namespace Dune
#endif
