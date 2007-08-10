#ifndef DUNE_LOCALMATRIX_HH
#define DUNE_LOCALMATRIX_HH

//- system includes 
#include <vector> 

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune 
{ 

  /** @ingroup Matrix  
      @{ 
  **/

  /** \brief Interface for local matrix classes. */
  template <class LocalMatrixTraits> 
  class LocalMatrixInterface 
  {
  public:  
    //! type of traits class 
    typedef LocalMatrixTraits Traits; 
    
    //! type of local matrix implementation 
    typedef typename Traits :: LocalMatrixType LocalMatrixType; 

    //! type of range field 
    typedef typename Traits :: RangeFieldType RangeFieldType;

    //! 
    typedef typename LocalMatrixType:: block_type LittleBlockType;
    
  protected:  
    //! constructor 
    LocalMatrixInterface () {} 

  public:
    /*! \brief add value to matrix entry (row,col) where row and col are
        local row and local column 
        \param[in] localRow local row 
        \param[in] localCol local column 
        \param[in] value value to add 
    */
    void add(const int localRow, 
             const int localCol, 
             const RangeFieldType& value)
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().add(localRow,localCol,value));
    }

    /*! \brief set value of matrix entry (row,col) where row and col are
        local row and local column 
        \param[in] localRow local row 
        \param[in] localCol local column 
        \param[in] value value to set  
    */
    void set(const int localRow, 
             const int localCol, 
             const RangeFieldType& value)
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().set(localRow,localCol,value));
    }

    /*! \brief get value of matrix entry (row,col) where row and col are
        local row and local column 
        \param[in] localRow local row 
        \param[in] localCol local column 
        \return value of matrix entry   
    */
    const RangeFieldType get(const int localRow, 
                             const int localCol ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(
        asImp().get(localRow,localCol));
      return asImp().get(localRow,localCol);
    }

    /*! \brief set all entries of local matrix to zero */
    void clear ()
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().clear());
    }

    /*! \brief resort ordering in global matrix (if possible) */
    void resort() 
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().resort());
    }
  private:
    //! Barton-Nackman Trick
    LocalMatrixType& asImp() 
    {
      return static_cast<LocalMatrixType&> (*this);
    }
    
    //! Barton-Nackman Trick
    const LocalMatrixType& asImp() const 
    {
      return static_cast<const LocalMatrixType&> (*this);
    }
  };
  
  /** \brief Default implementation for local matrix classes. */
  template <class LocalMatrixTraits> 
  class LocalMatrixDefault : public LocalMatrixInterface<LocalMatrixTraits> 
  {
    protected:
      LocalMatrixDefault () {}
    public:
  };

///@} 
} // end namespace Dune 
#endif
