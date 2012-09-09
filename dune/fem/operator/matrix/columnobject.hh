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

      ColumnObject( const LinearOperator &linOp, const RowEntityType &rowEntity )
      :
        linOp_( linOp ),
        rowEntity_( rowEntity )
      {}

      //! return local matrix 
      inline LocalMatrixType localMatrix( const ColumnEntityType &colEntity ) const
      {
        return linOp_.localMatrix( rowEntity_, colEntity );
      }

      //! return domain space (i.e. space that builds the rows)
      const DomainSpaceType& domainSpace() const { return linOp_.domainSpace(); }

      //! return range space (i.e. space that builds the columns)
      const RangeSpaceType& rangeSpace() const { return linOp_.rangeSpace(); }

    private:
      const LinearOperator &linOp_;
      const RowEntityType &rowEntity_;
    };

  } // namespace Fem
} // namespace Dune

#endif //#ifndef DUNE_FEM_COLUMNOBJECT_HH
