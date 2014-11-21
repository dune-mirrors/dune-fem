#ifndef DUNE_FEM_COLUMNOBJECT_HH
#define DUNE_FEM_COLUMNOBJECT_HH

namespace Dune
{

  namespace Fem
  {

    template< class LinearOperator >
    struct ColumnObject
    {
      typedef typename LinearOperator::ColumnEntityType ColumnEntityType;
      typedef typename LinearOperator::RowEntityType RowEntityType;

      typedef typename LinearOperator::DomainSpaceType DomainSpaceType;
      typedef typename LinearOperator::RangeSpaceType RangeSpaceType;

      typedef typename LinearOperator::LocalMatrixType LocalMatrixType;

      ColumnObject( const LinearOperator &linOp, const ColumnEntityType &colEntity )
      :
        linOp_( linOp ),
        colEntity_( colEntity )
      {}

      //! return local matrix
      inline LocalMatrixType localMatrix( const RowEntityType &rowEntity ) const
      {
        return linOp_.localMatrix( rowEntity, colEntity_ );
      }

      //! return domain space (i.e. space that builds the rows)
      const DomainSpaceType& domainSpace() const { return linOp_.domainSpace(); }

      //! return range space (i.e. space that builds the columns)
      const RangeSpaceType& rangeSpace() const { return linOp_.rangeSpace(); }

    private:
      const LinearOperator &linOp_;
      const ColumnEntityType &colEntity_;
    };

  } // namespace Fem
} // namespace Dune

#endif //#ifndef DUNE_FEM_COLUMNOBJECT_HH
