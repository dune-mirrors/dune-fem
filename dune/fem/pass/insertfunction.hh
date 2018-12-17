#ifndef DUNE_FEM_PASS_INSERTFUNCTION_HH
#define DUNE_FEM_PASS_INSERTFUNCTION_HH

#if HAVE_DUNE_FEM_DG
#error "Outdated header, #include <dune/fem-dg/pass/insertfunction.hh> instead!"
#endif


#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/pass/common/pass.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class DiscreteFunction >
    struct InsertFunctionPassDiscreteModel;

    template< class DiscreteFunction, class PreviousPass, int passId = -1 >
    class InsertFunctionPass;



    // InsertFunctionPassDiscreteModelTraits
    // -------------------------------------

    //! Traits for InsertFunctionPass to create dummy discrete model
    template< class DiscreteFunction >
    struct InsertFunctionPassDiscreteModelTraits
    {
      typedef DiscreteFunction DiscreteFunctionType;

      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef DiscreteFunctionType DestinationType;
      typedef InsertFunctionPassDiscreteModel< DiscreteFunction > DiscreteModelType;
    };



    // InsertFunctionPassDiscreteModel
    // -------------------------------

    template< class DiscreteFunction >
    struct InsertFunctionPassDiscreteModel
    {
      typedef InsertFunctionPassDiscreteModelTraits< DiscreteFunction > Traits;
    };



    // InsertFunctionPass
    // ------------------

    /**
     * @brief Base class for specific pass implementations.
       InsertFunctionPass simply inserts a discrete function from outside of the pass tree
       into the current pass tree, for example when calculating the species
       transport the velocity function comes from a different pass but has to
       be inserted into the species pass.
     */
    template< class DiscreteFunction, class PreviousPass, int passId >
    class InsertFunctionPass
    : public Pass< InsertFunctionPassDiscreteModel< DiscreteFunction >, PreviousPass, passId >
    {
      typedef InsertFunctionPass< DiscreteFunction, PreviousPass, passId > ThisType;
      typedef Pass< InsertFunctionPassDiscreteModel< DiscreteFunction >, PreviousPass, passId > BaseType;

      static const bool hasLocalFunction = std::is_convertible< DiscreteFunction, HasLocalFunction >::value;
      static_assert( hasLocalFunction, "InsertFunctionPass can only insert grid functions." );

    protected:
      template <class DFType>
      struct LocalFunctionInitializer
      {
        template <class ArgType>
        static void init( const ArgType&, const double time, DiscreteFunction& ) {}
      };

      template <class LFType>
      struct LocalFunctionInitializer< LocalFunctionAdapter< LFType > >
      {
        template <class ArgType>
        static void init( const ArgType& arg, const double time, DiscreteFunction& dest )
        {
          // call initialize on LocalFunctionAdapter
          dest.initialize( arg, time );
        }
      };

    public:
      //! type of discrete model for this class
      typedef InsertFunctionPassDiscreteModel< DiscreteFunction > DiscreteModelType;

      //! type of traits for this class
      typedef typename DiscreteModelType::Traits Traits;

      //! Repetition of template arguments
      typedef PreviousPass PreviousPassType;

      typedef typename BaseType::TotalArgumentType ArgumentType;
      typedef typename BaseType::GlobalArgumentType GlobalArgumentType;

      //! discrete function representing the return value of this pass
      typedef typename Traits::DestinationType DestinationType;
      //! discrete function space belonging to DestinationType
      typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    public:
      /** constructor
       *
       *  \param[in]  destination   to be stored in this pass
       *  \param[in]  previousPass  previous pass in the tree
       */
      InsertFunctionPass ( const DestinationType &destination, PreviousPassType &previousPass )
      : BaseType( previousPass )
      {
        destination_ = const_cast< DestinationType * >( &destination );
      }

      /** constructor
       *
       *  \param[in]  destination   to be stored in this pass
       *  \param[in]  previousPass  previous pass in the tree
       *  \param[in]  space         discrete function space
       *
       *  \note The argument space is ignored; it just makes the constructor
       *        look like that from other passes.
       */
      InsertFunctionPass ( const DestinationType &destination,
                           PreviousPassType &previousPass,
                           const DiscreteFunctionSpaceType &space )
      : BaseType( previousPass )
      {
        destination_ = const_cast< DestinationType * >( &destination );
      }

      /** constructor
       *
       *  \param[in]  previousPass  previous pass in the tree
       */
      explicit InsertFunctionPass ( PreviousPassType &previousPass )
      : BaseType( previousPass )
      {
        assert( destination_ == 0 );
      }

      /** default constructor
       *  \note This constructor creates and instance of previous pass
       */
      InsertFunctionPass ( const DestinationType* destination, const std::shared_ptr< PreviousPassType >& previousPass )
      : BaseType( *previousPass ),
        prevPassPtr_( previousPass )
      {
        destination_ = const_cast< DestinationType * >( destination );
      }

      //! destructor
      ~InsertFunctionPass ()
      {
        destination_ = 0;
      }

      // empty method here
      void allocateLocalMemory ()
      {}

      //! return reference to space
      const DiscreteFunctionSpaceType &space () const
      {
        assert( destination_ != 0 );
        return destination_->space();
      }

      //! set internal destination pointer to dest
      void setDestination ( const DestinationType &destination )
      {
        destination_ = const_cast< DestinationType * >( &destination );
      }

      using BaseType::time;

    protected:
      // empty method here
      void compute ( const ArgumentType &arg, DestinationType &dest ) const
      {
        // in case DestinationType is a LocalFunctionAdapter
        // call initialize
        LocalFunctionInitializer< DestinationType >::init( arg, time(), dest );
      }

      std::shared_ptr< PreviousPassType > prevPassPtr_;
      using BaseType::destination_;
    }; // end class InsertFunctionPass

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_INSERTFUNCTION_HH
