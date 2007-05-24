#ifndef DUNE_LAGRANGESPACE_DOFHANDLER_HH
#define DUNE_LAGRANGESPACE_DOFHANDLER_HH

namespace Dune 
{

  /** \class LagrangeDofHandler
   *  \brief Provides access to the local coordinates of a Lagrange point
   *
   *  The LagrangeDofHandler provides quadrature-like access to the local 
   *  coordinates of a DoF. Moreover, it can provide access to the Lagrange
   *  points lying on a face of the entity.
   *
   *  \ņote This class is deprecated and only wraps the Lagrange point set
   *        found in the Lagrange discrete function space.
   */
  template< class DiscreteFunctionSpaceImp >
  class LagrangeDofHandler
  { 
  public:
    //! type of the discrete function space
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
  
    //! type of the grid
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;
  
    //! type of entities of codimension 0
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;

    //! type of Lagrange point set
    typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType
      LagrangePointSetType;

    //! field type of coordinates of Lagrange points
    typedef typename LagrangePointSetType :: FieldType FieldType;
    //! type of Lagrange points
    typedef typename LagrangePointSetType :: PointType PointType;
    //! dimension of Lagrange points
    enum { dimension = LagrangePointSetType :: dimension };

  private:
    typedef LagrangeDofHandler< DiscreteFunctionSpaceType > ThisType;

    typedef ReferenceElementContainer< FieldType, dimension >
      ReferenceElementContainerType;
    typedef typename ReferenceElementContainerType :: value_type
      ReferenceElementType;
 
  private:
    //! Lagrange point set implementing the actual dof handler
    const LagrangePointSetType &lagrangePointSet_;

    //! access to the reference elements
    ReferenceElementContainerType refElementContainer_;
    //! reference element corresponding to the entity
    const ReferenceElementType &refElement_;

  public:
    /** constructor: Initialization with discrete function space and entity
     *
     *  \param[in] space Lagrange function space
     *
     *  \param[in] entity entity the handler operates on
     */
    inline LagrangeDofHandler ( const DiscreteFunctionSpaceType& space,
                                const EntityType& entity )
    : lagrangePointSet_( space.lagrangePointSet( entity ) ),
      refElementContainer_(),
      refElement_( refElementContainer_( lagrangePointSet_.geometryType() ) )
    {
    }
  
    /** return the Lagrange point
     *
     *  \param[in] i number of the Lagrange point within the entity. The number
     *               corresponds directly to the index of the DoF within a
     *               local function.
     *
     *  \return Reference to the Lagrange Point
     */
    inline const PointType& point( int i ) const
    {
      return lagrangePointSet_[ i ];
    }

    /** determine number of lagrange points on a face
     * 
     *  This method can be used to determine base functions with support on
     *  the face.
     *
     *  \param[in] face local number of the face within the entity
     *
     *  \param[in] codim codimension of the face
     *
     *  \return number of Lagrange points on the specified face
     */
    inline int numDofsOnFace( int face, int codim ) const
    {
      int count = 0;
      for( int subCodim = codim; subCodim <= dimension; ++subCodim )
      {
        const int subFaces = refElement_.size( face, codim, subCodim );
        for( int i = 0; i < subFaces; ++i ) {
          const int subFace
            = refElement_.subEntity( face, codim, i, subCodim );
          count += lagrangePointSet_.numDofs( subCodim, subFace );
        }
      }
      return count;
    }
  
    /** determine local DoF number for a given DoF on a face
     *
     * By iterating over \$0 \le faceDof \le numDofsOnFace\$ all base functions
     * with support on the face can be accessed.
     * 
     * \param[in] face local number of the face within the entity
     *
     * \param{in] codim codimension of the face
     *
     * \param[in] faceDofNumber number of the DoF with respect to the face
     *
     * \return number of the DoF with respect to the entity
     */
    inline int entityDofNum( int face, int codim, int faceDofNumber ) const
    {
      for( int subCodim = codim; subCodim <= dimension; ++subCodim )
      {
        const int subFaces = refElement_.size( face, codim, subCodim );
        for( int i = 0; i < subFaces; ++i ) {
          const int subFace
            = refElement_.subEntity( face, codim, i, subCodim );
          const int count = lagrangePointSet_.numDofs( subCodim, subFace );
          if( faceDofNumber < count )
            return lagrangePointSet_.entityDofNumber
                     ( subCodim, subFace, faceDofNumber );
          faceDofNumber -= count;
        }
      }
      assert( false );
      return -1;
    }
  };

}

#endif
