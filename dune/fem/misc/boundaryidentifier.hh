#ifndef DUNE_BOUNDARYIDENTIFIER_HH 
#define DUNE_BOUNDARYIDENTIFIER_HH 

namespace Dune { 

/** \brief Simple class for boundary identifiers 
  use in numerical problems.
*/
class BoundaryIdentifier 
{
  public:
    //! \brief tpyes of boundary 
    enum BoundaryType { 
      //! Dirichlet boundary with non-zero boundary values 
      DirichletNonZero , 
      //! Dirichlet boundary with zero boundary values 
      DirichletZero , 
      //! Neumann boundary with non-zero boundary fluxes 
      NeumannNonZero, 
      //! Neuman boundary with zero boundary fluxes 
      NeumannZero , 
      //! mixed boundary type 
      Robin ,
      //! undefined boundary type 
      undefined };

    BoundaryType bnd_;  
  public:
    // create undefined bnd 
    BoundaryIdentifier() : bnd_(undefined) {}
    // create bnd with given id 
    BoundaryIdentifier(const BoundaryType bnd) : bnd_(bnd) {} 
    // copy bnd 
    BoundaryIdentifier(const  BoundaryIdentifier & org) : bnd_(org.bnd_) {} 

    //! returns true if bnd is either Dirichlet or DirichletZero
    bool isDirichletType() const 
    { 
      return (bnd_ == DirichletNonZero) || 
             (bnd_ == DirichletZero); 
    }
    
    //! returns true if bnd is DirichletZero
    bool isDirichletZero() const 
    { 
      return (bnd_ == DirichletZero);
    }

    //! returns true if bnd is DirichletNonZero
    bool isDirichletNonZero() const 
    { 
      return (bnd_ == DirichletNonZero);
    }

    //! returns true, if bnd is either NeumannNonZero or NeumannZero type
    bool isNeumannType() const 
    { 
      return (bnd_ == NeumannNonZero) || 
             (bnd_ == NeumannZero); 
    }
    
    //! returns true, if bnd is NeumannZero type
    bool isNeumannZero() const 
    { 
      return (bnd_ == NeumannZero); 
    }
   
    //! returns true, if bnd is Neumann type
    bool isNeumannNonZero() const 
    { 
      return (bnd_ == NeumannNonZero); 
    }
};

}// end namespace Dune 
#endif
