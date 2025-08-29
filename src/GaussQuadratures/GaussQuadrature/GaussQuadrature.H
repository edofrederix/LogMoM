#ifndef GaussQuadrature_H
#define GaussQuadrature_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "List.H"
#include "scalar.H"
#include "label.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class GaussQuadrature Declaration
\*---------------------------------------------------------------------------*/

class GaussQuadrature
{
protected:

    // Protected data

        //- List of quadrature nodes
        List<scalar> x_;

        //- List of quadrature weights
        List<scalar> w_;


public:

    //- Runtime type information
    TypeName("GaussQuadrature");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            GaussQuadrature,
            dictionary,
            (),
            ()
        );


    // Constructor
    GaussQuadrature();


    //- Destructor
    virtual ~GaussQuadrature();


    // Selectors

        static autoPtr<GaussQuadrature> New
        (
            const word quadrature,
            const label numNodes
        );


    // Member Functions

        //- Nodes
        inline const List<scalar>& x() const { return x_; }

        //- Weights
        inline const List<scalar>& w() const { return w_; }

        //- Number of nodes
        inline label size() const { return x_.size(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
