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

        //- The Sauter-mean diameter
        volScalarField d_;

        //- Continuous phase
        mutable const phaseModel* continuousPhasePtr_;

        // Scaled moment fields

            //- Scaled number concentration (per cm^3)
            volScalarField lambda_;

            //- Scaled interfacial area
            volScalarField kappai_;

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

        //- Limit the moments
        void limitMoments();

        //- Correct the size distribution parameters (dcm and sigma) from the
        //  scaled moments
        void correctDistribution();

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

        //- Scaled number concentration (per cm^3)
        inline const volScalarField& lambda() const
        {
            return lambda_;
        }

        //- Scaled interfacial area
        inline const volScalarField& kappai() const
        {
            return kappai_;
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

        //- Get the Sauter-mean diameter field
        virtual tmp<volScalarField> d() const
        {
            return d_;
        }

        //- Get the surface area per unit volume field
        virtual tmp<volScalarField> Av() const
        {
            return phase()*kappai_;
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

        //- Compute the diameter based on moment p and q of the distribution
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
