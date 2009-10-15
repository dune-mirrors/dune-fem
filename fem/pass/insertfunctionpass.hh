#ifndef DUNE_INSERTFUNCTIONPASS_HH
#define DUNE_INSERTFUNCTIONPASS_HH

#include <dune/fem/pass/pass.hh>

namespace Dune
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

    static const bool hasLocalFunction = Conversion< DiscreteFunction, HasLocalFunction >::exists;
    dune_static_assert( hasLocalFunction, "InsertFunctionPass can only insert grid functions." );
    
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
    
  protected:
    // empty method here
    void compute ( const ArgumentType &arg, DestinationType &dest ) const
    {}

    using BaseType::destination_;
  }; // end class InsertFunctionPass


} // end namespace Dune

#endif // #ifndef DUNE_INSERTFUNCTIONPASS_HH
