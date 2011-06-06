#ifndef DUNE_FEM_DIAGONALPRECONDITIONER_HH
#define DUNE_FEM_DIAGONALPRECONDTIIONER_HH

#include <dune/common/static_assert.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>



namespace Dune
{


// DiagonalPreconditioner
//   // -----------------------
//
//     /** \class DiagonalPreconditioner
//       *  \ingroup OEMSolver
//       *  \brief   Precondtioner, multiplies with inverse of the diagonal
//       *
//       *  \param  DFImp type of the disctete function
//       *  \param  MatrixImp type of the diagonal matrix
//       //                    */
//
template<class DFImp,class MatrixImp>
class DiagonalPreconditioner
  :public Operator< typename DFImp::RangeFieldType,typename DFImp::RangeFieldType,DFImp,DFImp>
{
  typedef DFImp DiscreteFunctionType;
  typedef MatrixImp MatrixType;
  
private:
  const MatrixType& matrix_;

public:
  DiagonalPreconditioner(const MatrixType &matrix)
    :matrix_(matrix)
  {}

 virtual  void operator()(const DFImp &u,DFImp &res) const
  {
    apply(u,res);
  }

//   void apply(const DFImp u,DFImp res)
//       {
//         typedef typename DFImp::DofIteratorType DofIteratorType;  
//         typedef typename DFImp::ConstDofIteratorType ConstDofIteratorType;  
     
// 	DofIteratorType ret_it = res.dbegin(); 
// 	ConstDofIteratorType f_it = u.dbegin();  
// 	for(ret_it; ret_it!=res.dend(); ++ret_it)
// 	  {
// 	    int i=ret_it.index();
// 	    ret_it=f_it;
// 	    ret_it/=matrix_(i,i);
// 	    ++f_it;    
// 	  }
//       } 

protected:
  void apply(const DFImp &u,DFImp &res) const
  {
    typedef typename DFImp::DofIteratorType DofIteratorType;  
    typedef typename DFImp::ConstDofIteratorType ConstDofIteratorType;  
    
    DofIteratorType ret_it = res.dbegin(); 
    ConstDofIteratorType f_it = u.dbegin();  
    for(int i=0; ret_it!=res.dend(); ++ret_it)
      {
	(*ret_it)=(*f_it);
        (*ret_it)/=matrix_(i,i);
        ++f_it;    
	++i;
      }
      }

};
 

}
#endif
