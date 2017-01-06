#ifndef DUNE_FEM_LOCALDG_DISCRETEMODEL_HH
#define DUNE_FEM_LOCALDG_DISCRETEMODEL_HH

#include <cassert>

#include <dune/common/typetraits.hh>
#include <dune/common/bartonnackmanifcheck.hh>

#include <dune/fem/pass/common/selector.hh>

namespace Dune
{
  namespace Fem
  {

    // DGDiscreteModelInterface
    // ------------------------

    /**
     * \brief Interface class for problem definition in the LDG context.
     *
     *  The problem interface prescribes the methods needed by the implementation
     *  of local DG passes. Users need to derive from either this class or the
     *  DGDiscreteModelDefault class. The methods provided by this class are used by
     *  LocalDGPass (and possibly other classes) to solve for the space
     *  discretization of the equation d_t u + div f(u) = s(u), where f
     *  represents a flux function and s a source term. Even non-conservative
     *  fluxes can be treated (for details see paper by Dedner et al).
     *
     *  \note The definition of f(u) differs from the usual definition for
     *        conservation laws in the leading sign.
     */
    template <class DGDiscreteModelTraits>
    struct DGDiscreteModelInterface
    {
      //! \brief Traits class defined by the user
      typedef DGDiscreteModelTraits Traits;
      //! \brief Implementation type for Barton-Nackman trick
      typedef typename Traits::DGDiscreteModelType DGDiscreteModelType;

      //! \brief Discrete function space
      typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      //! \brief Function space type
      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
      //! \brief Coordinate type (world coordinates)
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! \brief Vector type of the function space's range
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief Vector type of the function space's range field type
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      //! \brief Matrix type of the function space's jacobian
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

      //! \brief Type of GridPart
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      //! \brief Element (codim 0 entity) of the grid
      typedef typename DiscreteFunctionSpaceType::EntityType EntityType;
      //! \brief Intersection type
      typedef typename DiscreteFunctionSpaceType::IntersectionType IntersectionType;

      //! \brief Type of local coordinate
      typedef typename EntityType::Geometry::LocalCoordinate  LocalCoordinateType;

      //! \brief dimRange
      enum { dimRange = DiscreteFunctionSpaceType::dimRange };

      //! \brief Mass factor type
      template< class MatrixType > struct ExtractMatrix;

      template< template < class, int, int > class MatrixType,
                class ctype, int dimD, int dimR
              >
      struct ExtractMatrix<MatrixType< ctype, dimD, dimR> >
      {
        typedef MatrixType< RangeFieldType, dimRange, dimRange > Type;
      };

      //! \brief Type of mass factor
      typedef typename ExtractMatrix< JacobianRangeType >::Type MassFactorType;

    public:
      /** \brief Returns true if problem has a flux contribution.
       *         If you program this method to return true, make sure to implement
       *         numericalFlux and analyticalFlux as well.
       */
      bool hasFlux () const { return asImp().hasFlux(); }

      /** \brief Returns true if problem has a source term.
       *         If you program this method to return true, make sure to implement
       *         method source well.
       */
      bool hasSource () const { return asImp().hasSource(); }

      //! Returns true if problem has a mass matrix factor disctinct from the identity
      bool hasMass () const { return asImp().hasMass(); }

      /** \brief Computes the numerical flux at a cell interface.
       *
       *  Interface method for calculating the numerical flux at cell interfaces.
       *  The method has distinct return values for the contribution to the inner
       *  and outer element, since the framework allows for non-conservative
       *  terms. The non-conservative terms don't get evaluated over separate
       *  methods, but need to be implemented using the flux and source methods.
       *
       *  \param intersection The intersection in consideration.
       *  \param time Global time.
       *  \param x Local coordinate of the point where the numerical flux gets
       *           evaluated.
       *  \param uLeft Tuple of the states on the inner side of the intersection.
       *  \param uRight Tuple of the states on the outer side of the
       *                intersection.
       *  \param gLeft The numerical flux contribution to the inner element
       *               (return value).
       *  \param gRight The numerical flux contribution to the outer element
       *                (return value).
       *
       * \return Speed of the fastest wave. This information is used to compute
       * the maximum admissible timestep size.
       */
      template <class ArgumentTuple, class FaceDomainType>
      double numericalFlux ( const IntersectionType &intersection,
                             const double time,
                             const FaceDomainType &x,
                             const ArgumentTuple &uLeft,
                             const ArgumentTuple &uRight,
                             RangeType &gLeft,
                             RangeType &gRight)
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().numericalFlux( intersection, time, x, uLeft, uRight, gLeft, gRight) );
        return asImp().numericalFlux( intersection, time, x, uLeft, uRight, gLeft, gRight );
      }


      /** \brief Computes the flux at the boundary
       *         Special kind of numerical flux. The intersection iterator provides
       *         the necessary information to identify the type of boundary.
       *  \param intersection Intersection in consideration.
       *  \param time Global time.
       *  \param x Local coordinate of the point where the numerical flux gets
       *           evaluated.
       *  \param uLeft Tuple of the states on the inner side of the intersection.
       *  \param gLeft The numerical flux contribution to the inner element
       *               (return value).
       *
       *  \return Speed of the fastest wave. This information is used to compute
       *         the maximum admissible timestep size.
       */
      template <class ArgumentTuple, class FaceDomainType>
      double boundaryFlux ( const IntersectionType &intersection,
                            const double time,
                            const FaceDomainType &x,
                            const ArgumentTuple &uLeft,
                            RangeType &gLeft )
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().boundaryFlux( intersection, time, x, uLeft, gLeft ) );
        return asImp().boundaryFlux( intersection, time, x, uLeft, gLeft );
      }


      /** \brief Computes the analytical flux of the problem.
       *         Analytical flux of the problem.
       *
       *  \param en Entity where the flux is to be evaluated
       *  \param time Global time.
       *  \param x Reference element coordinate where the flux is evaluated
       *  \param u Tuple of the states of the previous passes and of the global
       *           argument.
       *  \param f The analytical flux (return value)
       */
      template <class ArgumentTuple>
      void analyticalFlux ( const EntityType& en,
                            const double time,
                            const LocalCoordinateType& x,
                            const ArgumentTuple& u,
                            JacobianRangeType& f )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
                asImp().analyticalFlux(en, time, x, u, f) );
      }

      /** \brief Implements the source term of the problem.
       *         Source term contribution. The source term can contain parts of the
       *         non-conservative contributions.
       *
       *  \param en Entity where the flux is to be evaluated
       *  \param time Global time.
       *  \param x Reference element coordinate where the flux is evaluated
       *  \param u Tuple of the states of the previous passes and of the global
       *           argument.
       *  \param jac Tuple of the gradients of the precedings passes' results
       *             (needed for the non-conservative contributions)
       *  \param s The source contribution (return value).
       */
      template <class ArgumentTuple, class JacobianTuple>
      double source ( const EntityType& en,
                      const double time,
                      const LocalCoordinateType& x,
                      const ArgumentTuple& u,
                      const JacobianTuple& jac,
                      RangeType& s )
      {
        CHECK_INTERFACE_IMPLEMENTATION(
            asImp().source(en, time, x, u, jac, s) );
        return asImp().source(en, time, x, u, jac, s);
      }

      /** \brief Implements the mass factor term of the problem.
       *         Mass term contribution.
       *
       *  \param en Entity where the mass is to be evaluated
       *  \param time Global time.
       *  \param x Reference element coordinate where the mass is evaluated
       *  \param u Tuple of the states of the previous passes and of the global
       *  \param m The mass contribution (return value). default implementation
       *           sets this factor to 1.0
       */
      template <class ArgumentTuple>
      void mass ( const EntityType& en,
                  const double time,
                  const LocalCoordinateType& x,
                  const ArgumentTuple& u,
                  MassFactorType& m )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
            asImp().mass(en, time, x, u, m) );
      }

      /** \brief Passes the active entity to the model.
       *         This can be used, to set local functions required as data function
       *         in the model.
       *
       *  \param en active Entity
       */
      void setEntity ( EntityType& en )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().setEntity(en) );
      }

      /** \brief Passes the active neigbor entity to the model.
       *         This can be used, to set local functions required as data functions
       *         in the model.
       *
       *  \param nb active neighbor Entity
       */
      void setNeighbor ( EntityType& nb )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().setNeighbor(nb) );
      }

    protected:
      DGDiscreteModelType& asImp() { return static_cast<DGDiscreteModelType&>(*this); }
      const DGDiscreteModelType& asImp() const { return static_cast<const DGDiscreteModelType&>(*this); }
    };



    // DGDiscreteModelDefault
    // ----------------------

    /** \brief Default implementation of the DGDiscreteModelInterface where methods for
     *         the fluxes and the source term do nothing, so that the user needn't
     *         implement them if not needed.
     *
     *  \note I are passIds on which model depends
     */
    template< class DGDiscreteModelTraits, int... I >
    class DGDiscreteModelDefault :
      public DGDiscreteModelInterface<DGDiscreteModelTraits>
    {
      typedef DGDiscreteModelInterface< DGDiscreteModelTraits > BaseType;

    public:
      /** Selector for data tuple to use as arguments for all methods;
       * this fixes the template type ArgumentTuple.
       * If this discrete model is used for a pass n+1, i.e., following
       * passes p0,p1,..,pn then the return type of pass i (i=0,..,n)
       * can be used by adding the integer number i in the Selector.
       * Assume the following: \$ u_{n+1} = p_{n+1}(u_n,u_{n-1},..,u_1,u_0) \$
       * where $u_0=u$ is the global argument of the combined passes.
       * If \$ p_{n+1} \$ only depends on \$ u_0,u_2,u_n \$ then the
       * following selector can be used: \c Selector<0,n-1,1>. Then
       * ArgumentTuple is now filled with the values of these three
       * functions and can be accessed by...
       * Other way of filling the ArgumentTuple with corresponding pass
       * results is when one uses passIds. In this case if \$ u_{n+1} \$
       * depends on the passes with following passIds: firstPassId , passId2
       * , passId5 then the desired Selector is
       * Selector< firstPassId , passId2 , passId5 > ...
       * If there's no SelectorType in user-implemented DGDiscreteModel then
       * this Selector is used. Therefore it's good to pass passIds to this class
       * and avoid writing SelectorType in user-implemented DGDiscreteModel.
       * The point where a user specifies what's going to be in the Selector is
       * in the template declaration of the DGDiscreteModel where one names
       * passIds necessary for this DGDiscreteModel
       */
      typedef typename Dune::Fem::VariadicSelector< I... >::Type Selector;

    public:
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;
      typedef typename BaseType::LocalCoordinateType LocalCoordinateType;
      typedef typename BaseType::RangeType RangeType;
      typedef typename BaseType::JacobianRangeType JacobianRangeType;
      typedef typename BaseType::MassFactorType MassFactorType;

      /** \copydoc Dune::DGDiscreteModelInterface::hasFlux() const
       *
       *  The default implementation always returns false
       */
      inline bool hasFlux () const { return false; }

      /** \copydoc Dune::DGDiscreteModelInterface::hasSource() const
       *
       *  The default implementation always returns false
       */
      inline bool hasSource () const { return false; }

      /** \copydoc Dune::DGDiscreteModelInterface::hasMass() const
       *
       *  The default implementation always returns false
       */
      inline bool hasMass () const { return false; }

      /** \brief Empty implementation that fails if problem claims to have a flux
       *          contribution.
       */
      template <class ArgumentTuple, class FaceDomainType>
      double numericalFlux ( const IntersectionType& it,
                             const double time,
                             const FaceDomainType& x,
                             const ArgumentTuple& uLeft,
                             const ArgumentTuple& uRight,
                             RangeType& gLeft,
                             RangeType& gRight )
      {
        assert(!this->asImp().hasFlux());
        gLeft = 0.0;
        gRight = 0.0;
        return 0.0;
      }

      /** \brief Empty implementation that fails if problem claims to have a flux
       *         contribution.
       */
      template <class ArgumentTuple, class FaceDomainType>
      double boundaryFlux ( const IntersectionType &intersection,
                            const double time,
                            const FaceDomainType &x,
                            const ArgumentTuple &uLeft,
                            RangeType &gLeft )
      {
        assert( !this->asImp().hasFlux() );
        gLeft = 0.0;
        return 0.0;
      }

      /** \brief Empty implementation that fails if problem claims to have a flux
       *         contribution.
       */
      template< class ArgumentTuple >
      void analyticalFlux ( const EntityType &entity,
                            const double time,
                            const LocalCoordinateType &x,
                            const ArgumentTuple &u,
                            JacobianRangeType &f )
      {
        assert( !this->asImp().hasFlux() );
        f = 0.0;
      }

      /** \brief Empty implementation that fails if problem claims to have a source
       *         term.
       */
      template <class ArgumentTuple, class JacobianTuple>
      double source ( const EntityType& en,
                      const double time,
                      const LocalCoordinateType& x,
                      const ArgumentTuple& u,
                      const JacobianTuple& jac,
                      RangeType& s )
      {
        assert(!this->asImp().hasSource());
        s = 0.0;
        return 0.0;
      }

      /** \brief empty implementation for mass factor
       *         default implementation sets this factor to 1.0
       */
      template <class ArgumentTuple>
      void mass (const EntityType& en,
                 const double time,
                 const LocalCoordinateType& x,
                 const ArgumentTuple& u,
                 MassFactorType& m )
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

      //! \brier empty implementation
      void setEntity( const EntityType& en )
      { }

      //! \brief empty implementation
      void setNeighbor( const EntityType& nb )
      { }
    };

    // DGDGAdaptiveDiscreteModel
    // ---------------------------------------

    class DGAdaptiveDiscreteModel
    {
    public:
      //! \brief default method for setting adaptation handle to discrete model
      template <class Adaptation>
      void setAdaptation( Adaptation&,
                          const double weight = 1.0 ) {}

      //! \brief default method for setting adaptation handle and domain filter to discrete model
      template <class Adaptation, class DomainFilter>
      void setAdaptation( Adaptation&,
                          const DomainFilter&,
                          const double weight = 1.0 ) {}

      //! remove pointer to adaptation handle
      void removeAdaptation() {}
    };




    // DGDiscreteModelDefaultWithInsideOutside
    // ---------------------------------------

    //! Default implementation of the DGDiscreteModelInterface where methods for
    //! the fluxes and the source term do nothing, so that the user needn't
    //! implement them if not needed.
    template <class DGDiscreteModelTraits, int... I >
    class DGDiscreteModelDefaultWithInsideOutside
    : public DGDiscreteModelDefault< DGDiscreteModelTraits, I... >,
      public DGAdaptiveDiscreteModel
    {
      typedef DGDiscreteModelDefault< DGDiscreteModelTraits, I... > BaseType;

    public:
      typedef typename BaseType::EntityType EntityType;

      DGDiscreteModelDefaultWithInsideOutside ()
      : enVol_(-1.0) , nbVol_(-1.0) , en_(0) , nb_(0)
      {}

      DGDiscreteModelDefaultWithInsideOutside ( const DGDiscreteModelDefaultWithInsideOutside& other )
      : enVol_(other.enVol_) , nbVol_(other.nbVol_),
        en_(other.en_) , nb_(other.nb_)
      {}

      /** \brief method setting pointer of inside entity and getting volume
       *
       *  \param[in] en reference to inside entity
       */
      void setEntity ( const EntityType& en )
      {
        en_ = &en;
        enVol_ = en.geometry().volume();
      }

      /** \brief method seting pointer of outside entity and getting volume
       *
       *  \param[in] nb reference to outside entity
       */
      void setNeighbor ( const EntityType& nb )
      {
        nb_ = &nb;
        nbVol_ = nb.geometry().volume();
      }

      /** \brief method returning reference to inside entity
       *
       *  \return reference to inside entity
       */
      const EntityType &inside () const
      {
        assert( en_ );
        return *en_;
      }

      /** \brief method returning reference to outside entity
       *
       *  \return reference to outside entity
       */
      const EntityType &outside () const
      {
        assert( nb_ );
        return *nb_;
      }

      //! \brief return volume of entity
      double enVolume () const
      {
        assert(enVol_ > 0.0);
        return enVol_;
      }

      //! \brief return volume of neighbor
      double nbVolume () const
      {
        assert( nbVol_ > 0.0 );
        return nbVol_;
      }

    private:
      double enVol_;
      double nbVol_;

      const EntityType* en_;
      const EntityType* nb_;
    };


    //Deprecation warning for old usage of DGDiscreteModelDefault (and derived class DGDiscreteModelDefaultWithInsideOutside )
    template <class T,int N1,int N2,int N3,int N4,int N5,int N6,int N7,int N8 >
    class DGDiscreteModelDefault< T,N1,N2,N3,N4,N5,N6,N7,N8,-1>
    { static_assert( AlwaysFalse<T>::value,"Deprecated: You can use a variadic number of ids now: Avoid trailing '-1'!" ); };
    template <class T,int N1,int N2,int N3,int N4,int N5,int N6,int N7 >
    class DGDiscreteModelDefault< T,N1,N2,N3,N4,N5,N6,N7,-1 >
    { static_assert( AlwaysFalse<T>::value,"Deprecated: You can use a variadic number of ids now: Avoid trailing '-1'!" ); };
    template <class T,int N1,int N2,int N3,int N4,int N5,int N6 >
    class DGDiscreteModelDefault< T,N1,N2,N3,N4,N5,N6,-1 >
    { static_assert( AlwaysFalse<T>::value,"Deprecated: You can use a variadic number of ids now: Avoid trailing '-1'!" ); };
    template <class T,int N1,int N2,int N3,int N4,int N5 >
    class DGDiscreteModelDefault< T,N1,N2,N3,N4,N5,-1 >
    { static_assert( AlwaysFalse<T>::value,"Deprecated: You can use a variadic number of ids now: Avoid trailing '-1'!" ); };
    template <class T,int N1,int N2,int N3,int N4 >
    class DGDiscreteModelDefault< T,N1,N2,N3,N4,-1 >
    { static_assert( AlwaysFalse<T>::value,"Deprecated: You can use a variadic number of ids now: Avoid trailing '-1'!" ); };
    template <class T,int N1,int N2,int N3 >
    class DGDiscreteModelDefault< T,N1,N2,N3,-1 >
    { static_assert( AlwaysFalse<T>::value,"Deprecated: You can use a variadic number of ids now: Avoid trailing '-1'!" ); };
    template <class T,int N1,int N2 >
    class DGDiscreteModelDefault< T,N1,N2,-1 >
    { static_assert( AlwaysFalse<T>::value,"Deprecated: You can use a variadic number of ids now: Avoid trailing '-1'!" ); };
    template <class T,int N1 >
    class DGDiscreteModelDefault< T,N1,-1 >
    { static_assert( AlwaysFalse<T>::value,"Deprecated: You can use a variadic number of ids now: Avoid trailing '-1'!" ); };
    template <class T >
    class DGDiscreteModelDefault< T,-1 >
    { static_assert( AlwaysFalse<T>::value,"Deprecated: You can use a variadic number of ids now: Avoid trailing '-1'!" ); };
  } // namespace Fem

}  // namespace Dune

#endif // #ifndef DUNE_FEM_LOCALDG_DISCRETEMODEL_HH
