#ifndef DUNE_FEM_PASS_LOCAL_HH
#define DUNE_FEM_PASS_LOCAL_HH

#if HAVE_DUNE_FEM_DG
#error "Outdated header, #include <dune/fem-dg/pass/pass.hh> instead!"
#endif


#include <memory>
#include <sstream>
#include <string>

#include <dune/fem/common/memory.hh>
#include <dune/fem/pass/common/pass.hh>

namespace Dune
{

  namespace Fem
  {

    // External forward declaration
    // ----------------------------

    template <class DiscreteModelImp, class PreviousPassImp , int passIdImp >
    class Pass;

    // LocalPass
    // ---------

    /** \brief Specialisation of Pass which provides a grid walk-through,
     *         but leaves open what needs to be done on each elements.
     *
     *  \tparam  DiscreteModelImp  discrete model
     *  \tparam  PreviousPassImp   previous pass
     *  \tparam  passIdImp         id for this pass
     */
    template< class DiscreteModelImp, class PreviousPassImp , int passIdImp >
    class LocalPass
    : public Pass< DiscreteModelImp , PreviousPassImp , passIdImp>
    {
    public:
      //! \brief type of the preceding pass
      typedef PreviousPassImp PreviousPassType;

      //! \brief base class
      typedef Pass< DiscreteModelImp , PreviousPassImp , passIdImp > BaseType;

      /** \brief The type of the argument (and destination) type of
       *         the overall operator
       */
      typedef typename BaseType::TotalArgumentType ArgumentType;

      //! \brief the discrete function representing the return value of this pass
      typedef typename DiscreteModelImp::Traits::DestinationType DestinationType;
      //! \brief the discrete function space belonging to destinationtype
      typedef typename DiscreteModelImp::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      //! \brief iterator over the space
      typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
      //! \brief the codim 0 entity
      typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

      /** \brief constructor
       *  \param pass Previous pass
       *  \param spc Space belonging to the discrete function of this pass.
       *  \param passName an identifier for this pass
       */
      LocalPass (PreviousPassImp &pass,
                 const DiscreteFunctionSpaceType &spc,
                 std::string passName = "LocalPass")
      : BaseType(pass),
        spc_( referenceToSharedPtr( spc ) ),
        passName_(passName),
        computeTime_(0.0),
        numberOfElements_( 0 ),
        passIsActive_( true )
      {}

      //! \brief destructor
      virtual ~LocalPass () {}

      //! \brief build up local memory
      virtual void allocateLocalMemory ()
      {
        if (!this->destination_)
        {
          std::ostringstream funcName;
          funcName << passName_ << "_" << this->passNumber();
          this->destination_ = new DestinationType(funcName.str(), space());

          // set mem handle for deleting destination_
          this->deleteHandler_ = &(BaseType::DeleteHandlerType::instance());
        }
      }

      //! \brief return reference to space
      const DiscreteFunctionSpaceType &space () const { return *spc_; }

      /** \brief return accumulated time needed by pass's operator () this method
       *         also resets the compute time to zero
       */
      virtual double computeTime () const
      {
        double ct = computeTime_;
        computeTime_ = 0.0;
        return ct;
      }

      /** \brief return number of elements visited during operator computation */
      virtual size_t numberOfElements() const { return numberOfElements_; }

      /** \brief return true if pass is active */
      bool active () const { return passIsActive_; }

      /** \brief set pass status to active */
      void enable () const { passIsActive_ = true ; }

      /** \brief set pass status to inactive */
      void disable() const { passIsActive_ = false ; }

    protected:
      //! Actions to be carried out before a global grid walkthrough.
      //! To be overridden in a derived class.
      virtual void prepare (const ArgumentType &arg, DestinationType &dest) const = 0;
      //! Actions to be carried out after a global grid walkthrough.
      //! To be overridden in a derived class.
      virtual void finalize (const ArgumentType &arg, DestinationType &dest) const = 0;
      //! Actions to be taken on every element. To be overridden in a derived
      //! class.
      virtual void applyLocal (const EntityType &en ) const = 0;

    public:
      //! The actual computations are performed as follows. First, prepare
      //! the grid walkthrough, then call applyLocal on each entity and then
      //! call finalize.
      void compute (const ArgumentType &arg, DestinationType &dest) const
      {
        // if pass was disable, don't do computation
        if( ! active() ) return ;

        // get stopwatch
        Dune::Timer timer;

        prepare(arg, dest);

        numberOfElements_ = 0 ;
        const IteratorType endit = space().end();
        for (IteratorType it = space().begin(); it != endit; ++it, ++ numberOfElements_ )
        {
          applyLocal(*it);
        }

        finalize(arg, dest);

        // accumulate time
        computeTime_ += timer.elapsed();
      }

    protected:
      std::shared_ptr< const DiscreteFunctionSpaceType > spc_;
      const std::string passName_;
      mutable double computeTime_;
      mutable size_t numberOfElements_;
      mutable bool passIsActive_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_LOCAL_HH
