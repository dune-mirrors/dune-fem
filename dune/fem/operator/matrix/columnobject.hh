#ifndef DUNE_FEM_COLUMNOBJECT_HH
#define DUNE_FEM_COLUMNOBJECT_HH

namespace Dune
{

  namespace Fem
  {

    template< class LinearOperator >
    struct ColumnObject
    {
      typedef typename LinearOperator::DomainEntityType DomainEntityType;
      typedef typename LinearOperator::RangeEntityType RangeEntityType;

      typedef typename LinearOperator::DomainSpaceType DomainSpaceType;
      typedef typename LinearOperator::RangeSpaceType RangeSpaceType;

      typedef typename LinearOperator::LocalMatrixType LocalMatrixType;

      ColumnObject( const LinearOperator &linOp, const DomainEntityType &domainEntity )
      :
        linOp_( linOp ),
        domainEntity_( domainEntity )
      {}

      //! return local matrix 
      inline LocalMatrixType localMatrix( const RangeEntityType &rangeEntity ) const
      {
        return linOp_.localMatrix( domainEntity_, rangeEntity );
      }

      //! return domain space (i.e. space that builds the rows)
      const DomainSpaceType& domainSpace() const { return linOp_.domainSpace(); }

      //! return range space (i.e. space that builds the columns)
      const RangeSpaceType& rangeSpace() const { return linOp_.rangeSpace(); }

    private:
      const LinearOperator &linOp_;
      const DomainEntityType &domainEntity_;
    };

  } // namespace Fem
} // namespace Dune

#endif //#ifndef DUNE_FEM_COLUMNOBJECT_HH
