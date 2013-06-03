#ifndef DUNE_FEM_PRECONDITIONING_HH
#define DUNE_FEM_PRECONDITIONING_HH

namespace OEMSolver 
{

  /** \brief Interface class for use of Preconditioners with the OEM solvers.
  */
  class PreconditionInterface 
  {
  public:
    //! type of this class 
    typedef PreconditionInterface ThisType; 

    //! return reference to precondition matrix 
    const ThisType & preconditionMatrix() const { return *this; }
    
    //! returns true, if preconditioning should be used 
    //! default is false 
    bool hasPreconditionMatrix() const { return false; }
  };

} // namespace OEMSolver 

#endif // #ifndef DUNE_FEM_PRECONDITIONING_HH
