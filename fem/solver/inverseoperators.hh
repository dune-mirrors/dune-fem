#ifndef DUNE_INVERSE_OPERATORS_HH
#define DUNE_INVERSE_OPERATORS_HH

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>

namespace Dune {

  /** \brief Implementation of the CG algorithm */
  template <class OperatorType, class DiscreteFunctionType>
  struct CGAlgorithm 
  {
    /** solve Op(arg) - dest = 0 */      
    static void cg (const OperatorType & op, 
                    const DiscreteFunctionType& arg,
                    DiscreteFunctionType& dest ,
                    double epsilon , int maxIter , bool verbose ) 
    {
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType FunctionSpaceType;
      typedef typename FunctionSpaceType::RangeFieldType Field;
      typedef typename FunctionSpaceType :: GridType GridType; 

      typedef typename GridType :: Traits :: CollectiveCommunication
        CommunicatorType; 

      const CommunicatorType & comm = arg.space().grid().comm();

      int count = 0;
      Field spa=0, spn, q, quad;

      Field b = arg.scalarProductDofs( arg );
      b = comm.sum( b );
      const Field err = epsilon * b;

      DiscreteFunctionType r ( arg );
      DiscreteFunctionType h ( arg );
      DiscreteFunctionType p ( arg );
    
      op( dest, h );

      r.assign(h) ;
      r -= arg;

      p.assign(arg);
      p -= h;

      spn = r.scalarProductDofs( r );

      // global sum 
      spn = comm.sum( spn );
   
      while((spn > err ) && (count++ < maxIter)) 
      {
        // fall ab der zweiten iteration *************
    
        if(count > 1)
        { 
          const Field e = spn / spa;
          p *= e;
          p -= r;
        }

        // grund - iterations - schritt **************
        op( p, h );
    
        quad = p.scalarProductDofs( h );

        // global sum 
        quad = comm.sum( quad );
        
        q    = spn / quad;

        dest.addScaled( p, q );
        r.addScaled( h, q );

        spa = spn;
    
        // residuum neu berechnen *********************
    
        spn = r.scalarProductDofs( r ); 
        // global sum 
        spn = comm.sum( spn );
        
        if( verbose && (comm.rank() == 0))
          std::cerr << count << " cg-Iterationen  " << count << " Residuum:" << spn << "        \r";
      }
      if( verbose && (comm.rank() == 0))
        std::cerr << "\n";
    }
  };

  /** @ingroup OEMSolver
      @{
   **/
  
  /** \brief Inversion operator using CG algorithm, operator type is
      Mapping. 
   */
  template <class DiscreteFunctionType>
  class CGInverseOperator : public Operator<
    typename DiscreteFunctionType::DomainFieldType,
    typename DiscreteFunctionType::RangeFieldType,
    DiscreteFunctionType,DiscreteFunctionType> 
  {
    typedef Operator<
            typename DiscreteFunctionType::DomainFieldType,
            typename DiscreteFunctionType::RangeFieldType,
            DiscreteFunctionType,DiscreteFunctionType> BaseType;
    
    typedef Mapping<typename DiscreteFunctionType::DomainFieldType ,
                    typename DiscreteFunctionType::RangeFieldType ,
                    DiscreteFunctionType,DiscreteFunctionType> MappingType;

  public:
    /** \brief constructor of CGInverseOperator
      \param[in] op Mapping describing operator to invert 
      \param[in] redEps reduction epsilon 
      \param[in] absLimit absolut limit of residual  
      \param[in] maxIter maximal iteration steps 
      \param[in] verbose verbosity 
    */
    CGInverseOperator( const MappingType & op, 
                       double redEps,
                       double absLimit,
                       int maxIter,
                       int verbose ) 
      : op_(op), _redEps ( redEps ), epsilon_ ( absLimit*absLimit ) , 
        maxIter_ (maxIter ) , _verbose ( verbose ) {
    } 

    /** \brief solve the system 
        \param[in] arg right hand side 
        \param[out] dest solution 
    */
    virtual void operator()(const DomainType& arg, 
                            RangeType& dest ) const 
    {
      CGAlgorithm<MappingType,DomainType>::cg(op_,arg,dest,epsilon_,maxIter_,(_verbose >0));
    }

  private:
    //! reference to operator which should be inverted 
    const MappingType &op_;
  
    //! reduce error by this factor 
    double _redEps; 

    //! minial error to reach 
    typename DiscreteFunctionType::RangeFieldType epsilon_;

    //! number of maximal iterations
    int maxIter_;

    //! level of output 
    int _verbose ;
  };

  /** \brief Inversion operator using CG algorithm, operator type is
      a template parameter. 
   */
  template <class DiscreteFunctionType, class OperatorType>
  class CGInverseOp : public Operator<
    typename DiscreteFunctionType::DomainFieldType,
    typename DiscreteFunctionType::RangeFieldType,
    DiscreteFunctionType,DiscreteFunctionType> 
  {
  public:
    /** \brief constructor of CGInverseOperator
      \param[in] op Mapping describing operator to invert 
      \param[in] redEps reduction epsilon 
      \param[in] absLimit absolut limit of residual  
      \param[in] maxIter maximal iteration steps 
      \param[in] verbose verbosity 
    */
    CGInverseOp( OperatorType & op , double  redEps , double absLimit , int maxIter , int verbose ) : 
      op_(op),
      _redEps ( redEps ),
      epsilon_ ( absLimit*absLimit ) , 
      maxIter_ (maxIter ) ,
      _verbose ( verbose )  
    {} 

    /** \brief solve the system 
        \param[in] arg right hand side 
        \param[out] dest solution 
    */
    virtual void operator() (const DiscreteFunctionType& arg,
                             DiscreteFunctionType& dest ) const 
    {
      CGAlgorithm<OperatorType,DiscreteFunctionType>::cg(op_,arg,dest,epsilon_,maxIter_,(_verbose >0));
    }
  
  private:
    //! operator to invert 
    OperatorType &op_;
  
    //! relative reduction 
    double _redEps; 

    //! absolut reduction  
    typename DiscreteFunctionType::RangeFieldType epsilon_;

    //! number of maximal iterations
    int maxIter_;

    //! level of output 
    int _verbose ;
  };

///@}
} // end namespace Dune
#endif
