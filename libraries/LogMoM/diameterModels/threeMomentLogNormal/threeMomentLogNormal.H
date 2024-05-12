#ifndef threeMomentLogNormal_H
#define threeMomentLogNormal_H

#include "diameterModel.H"
#include "coalescenceEfficiencyModel.H"
#include "coalescenceFrequencyModel.H"
#include "breakupModel.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class LogMoMInterfacialModel;

namespace functionObjects
{
    class threeMomentLogNormalDiameter;
    class threeMomentLogNormalSigma;
}

namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                           Class threeMomentLogNormal Declaration
\*---------------------------------------------------------------------------*/

class threeMomentLogNormal
:
    public diameterModel
{
    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- p and q diameter coefficients
        scalar p_;
        scalar q_;

        //- Number concentration (per cm^3)
        volScalarField N_;

        //- Interfacial area
        volScalarField A_;

        //- Limited and scaled number concentration (per cm^3)
        volScalarField lambda_;

        //- Limited and scaled interfacial area
        volScalarField kappa_;

        //- Maximum diameter used for stabilization
        dimensionedScalar dMax_;

        //- Minimum diameter used for stabilization
        dimensionedScalar dMin_;

        //- Number of correctors. Defaults to 2.
        label nCorr_;

        //- The Sauter-mean diameter of the phase
        volScalarField d_;

        //- The width of the size distribution
        volScalarField sigma_;

        //- Coalescence switch
        Switch coalescence_;

        //- Breakup switch
        Switch breakup_;

        //- Number concentration coalescence rate
        volScalarField NCoaRate_;

        //- Interfacial area concentration coalescence rate
        volScalarField ACoaRate_;

        //- Number concentration breakage rate
        volScalarField NBreakRate_;

        //- Interfacial area concentration breakage rate
        volScalarField ABreakRate_;

        //- Coalescence efficiency model
        autoPtr<coalescenceEfficiencyModel> coaEffModel_;

        //- Coalescence frequency model
        autoPtr<coalescenceFrequencyModel> coaFreqModel_;

        //- Breakup model
        autoPtr<breakupModel> breakupModel_;

        //- Fields to store the right-hand side implicit source term
        //  coefficients
        volScalarField NSource_;
        volScalarField ASource_;


protected:

    // Protected member functions

        //- Correct the scaled moments from the unscaled ones and apply
        //  realizability limiters.
        void correctLimitedScaledMoments();

        //- Compute the Sauter mean diameter
        tmp<volScalarField> dsm() const;

        //- Return the diameter based on moment p and q
        virtual tmp<volScalarField> d
        (
            const scalar p,
            const scalar q
        ) const;

        //- Width of the size distribution
        virtual tmp<volScalarField> sigma() const;

        //- Compute coalescence source term for moment gamma
        virtual void updateCoalescenceSources
        (
            const volScalarField& dcm,
            const volScalarField& s
        );

        //- Compute breakup source term for moment gamma using breakup models
        virtual void updateBreakupSources
        (
            const volScalarField& dcm,
            const volScalarField& s
        );

        //- Read the models
        virtual void readModels();

        //- Update the models
        virtual void updateSources();

public:

    //- Runtime type information
    TypeName("threeMomentLogNormal");


    // Constructors

        //- Construct from components
        threeMomentLogNormal
        (
            const dictionary& dict,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~threeMomentLogNormal();


    // Member Functions

        //- Return the local fields

        inline volScalarField& N()
        {
            return N_;
        }

        inline const volScalarField& N() const
        {
            return N_;
        }

        inline volScalarField& A()
        {
            return A_;
        }

        inline const volScalarField& A() const
        {
            return A_;
        }

        inline volScalarField& kappa()
        {
            return kappa_;
        }

        inline const volScalarField& kappa() const
        {
            return kappa_;
        }

        inline volScalarField& lambda()
        {
            return lambda_;
        }

        inline const volScalarField& lambda() const
        {
            return lambda_;
        }

        inline const volScalarField& NSource() const
        {
            return NSource_;
        }

        inline const volScalarField& ASource() const
        {
            return ASource_;
        }

        inline const dimensionedScalar& dMin() const
        {
            return dMin_;
        }

        inline const dimensionedScalar& dMax() const
        {
            return dMax_;
        }

        //- Evolve the moment transport equations
        virtual void correct();

        //- Get the diameter field
        virtual tmp<volScalarField> d() const;

        //- Interfacial area
        virtual tmp<volScalarField> a() const;

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);

        //- Return reference to the continuous phase
        const phaseModel& continuousPhase() const;

        //- Return the dispersed interface
        dispersedPhaseInterface interface() const;

        //- Return reference to turbulence model of the continuous phase
        const phaseCompressible::momentumTransportModel&
            continuousTurbulence() const;

        //- Return reference to turbulence model of the dispersed phase
        const phaseCompressible::momentumTransportModel&
            dispersedTurbulence() const;

        //- Befriend LogMoM interfacial model
        friend class Foam::LogMoMInterfacialModel;

        //- Befriend function objects
        friend class functionObjects::threeMomentLogNormalDiameter;
        friend class functionObjects::threeMomentLogNormalSigma;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
