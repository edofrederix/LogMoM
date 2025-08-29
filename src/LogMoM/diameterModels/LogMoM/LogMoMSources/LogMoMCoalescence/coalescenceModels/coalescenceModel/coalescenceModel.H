#ifndef coalescenceModel_H
#define coalescenceModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "LogMoM.H"
#include "runTimeSelectionTables.H"
#include "phaseCompressibleMomentumTransportModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class coalescenceModel Declaration
\*---------------------------------------------------------------------------*/

class coalescenceModel
{
protected:

    // Protected data

        //- Reference to the LogMoM diameter model
        const diameterModels::LogMoM& logmom_;

public:

    //- Runtime type information
    TypeName("coalescenceModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            coalescenceModel,
            dictionary,
            (
                const diameterModels::LogMoM& logmom,
                const dictionary& dict
            ),
            (logmom, dict)
        );


    // Constructors

        // Construct from ditionary
        coalescenceModel
        (
            const diameterModels::LogMoM& logmom,
            const dictionary& dict
        );


    //- Destructor
    virtual ~coalescenceModel();


    // Selectors

        static autoPtr<coalescenceModel> New
        (
            const diameterModels::LogMoM& logmom,
            const dictionary& dict
        );


    // Member Functions

        // Return a field of coalescence rates given two particle sizes
        virtual tmp<volScalarField> rate
        (
            const volScalarField& di,
            const volScalarField& dj
        ) const = 0;

        // Return the continuous phase turbulence model
        virtual const phaseCompressible::momentumTransportModel&
            continuousTurbulence() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
