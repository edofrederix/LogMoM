#ifndef LogMoM_H
#define LogMoM_H

#include "diameterModel.H"
#include "LogMoMSource.H"
#include "multivariateSelectionScheme.H"
#include "dispersedPhaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseCompressibleMomentumTransportModel;
class LogMoMInterfacialModel;

namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                            Class LogMoM Declaration
\*---------------------------------------------------------------------------*/

class LogMoM
:
    public regIOobject,
    public diameterModel
{
    // Private Data

        //- The representative diameter
        volScalarField d_;

        //- Representative diameter coefficients (defaults to p = 3 and q = 2,
        //  i.e., the Sauter-mean diameter)
        scalar p_;
        scalar q_;

        //- Continuous phase
        mutable const phaseModel* continuousPhasePtr_;

        // Main moment fields which are solved for

            //- Number concentration (per cm^3)
            volScalarField N_;

            //- Interfacial area
            volScalarField A_;

        // Limited moment fields shadowing the main ones

            //- Limited and scaled number concentration (per cm^3)
            volScalarField lambda_;

            //- Limited and scaled interfacial area
            volScalarField kappa_;

        // Log-normal size distribution parameters

            //- The count median diameter
            volScalarField dcm_;

            //- The width of the size distribution
            volScalarField sigma_;

        // Log-normal size distribution limits for stability

            //- Minimum/maximum diameter used for stabilization
            dimensionedScalar dMin_;
            dimensionedScalar dMax_;

            //- Minimum/maximum distribution width (defaults to 0 and 2)
            dimensionedScalar sigmaMin_;
            dimensionedScalar sigmaMax_;

        //- Number of correctors. Defaults to 2.
        label nCorr_;

        //- LogMoM sources
        PtrList<LogMoMSource> sources_;

        //- Moment order of integration. This value will be used by LogMoM
        //  interfacial closure models that perform integration. Defaults to 3,
        //  but can be (temporarily) changed to evaluate interfacial closure
        //  models for other moments.
        scalar gamma_;

        //- Multivariate convection scheme field table
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields_;

    // Private Member Functions

        // Correction Functions

            //- Correct limited moments from the normal ones
            void correctLimitedMoments();

            //- Correct the size distribution parameters (dcm and sigma) from
            //  the limited moments
            void correctDistribution();

            //- Correct the representative diameter
            void correctRepresentativeDiameter();

public:

    //- Cast a diameter model to LogMoM
    static const LogMoM* cast(const diameterModel& model)
    {
        return dynamic_cast<const LogMoM*>(&model);
    }

    //- Runtime type information
    TypeName("LogMoM");

    // Constructors

        //- Construct from dictionary and phase
        LogMoM
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~LogMoM();


    // Member Functions

        //- Number concentration (per cm^3)
        inline const volScalarField& N() const
        {
            return N_;
        }

        //- Interfacial area
        inline const volScalarField& A() const
        {
            return A_;
        }

        //- Limited and scaled number concentration (per cm^3)
        inline const volScalarField& lambda() const
        {
            return lambda_;
        }

        //- Limited and scaled interfacial area
        inline const volScalarField& kappa() const
        {
            return kappa_;
        }

        //- Count median diameter
        inline const volScalarField& dcm() const
        {
            return dcm_;
        }

        //- The width of the size distribution
        inline const volScalarField& sigma() const
        {
            return sigma_;
        }

        //- Return the moment order of integration
        inline scalar gamma() const
        {
            return gamma_;
        }

        //- Get the diameter field
        virtual tmp<volScalarField> d() const
        {
            return d_;
        }

        //- Get the surface area per unit volume field
        virtual tmp<volScalarField> Av() const
        {
            return A_;
        }

        //- Return continuous phase
        inline const phaseModel& continuousPhase() const
        {
            if (continuousPhasePtr_ == nullptr)
            {
                continuousPhasePtr_ =
                   &phase().mesh().lookupObject<phaseModel>
                    (
                        IOobject::groupName
                        (
                            "alpha",
                            diameterProperties().lookup("continuousPhase")
                        )
                    );
            }

            return *continuousPhasePtr_;
        }

        //- Return reference to momentumTransport model of the continuous phase
        const phaseCompressibleMomentumTransportModel&
            continuousTurbulence() const;

        //- Return the dispersed interface
        dispersedPhaseInterface interface() const;

        //- Return the diameter based on moment p and q
        virtual tmp<volScalarField> d
        (
            const scalar p,
            const scalar q
        ) const;

        //- Correct the model
        virtual void correct();

        //- Set the moment order of integration and return its previous value
        scalar setGamma(const scalar);

        //- Dummy write for regIOobject
        bool writeData(Ostream&) const;

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);

    // Befriend the LogMoM interfacial model class
    friend class ::Foam::LogMoMInterfacialModel;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
