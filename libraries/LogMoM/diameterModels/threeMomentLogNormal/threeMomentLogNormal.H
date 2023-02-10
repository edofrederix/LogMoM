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
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                           Class threeMomentLogNormal Declaration
\*---------------------------------------------------------------------------*/

class threeMomentLogNormal
:
    public diameterModel
{
public:

    // Public data types

        //- Scale transfer type enumerations
        enum class closingMomentType
        {
            interfacialArea,
            squaredVolume
        };

        //- Closing moment type names
        static const NamedEnum<closingMomentType, 2> closingMomentTypeNames_;

private:

    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- p and q diameter coefficients
        scalar p_;
        scalar q_;

        //- Additional moment used to close the log-normal size distribution
        closingMomentType closingMoment_;

        //- Scaled bubble number concentration (per cm^3, to avoid large numbers!)
        volScalarField lambda_;

        //- Scaled interfacial area concentration
        volScalarField kappa_;

        //- Scaled squared volume concentration
        volScalarField beta_;

        //- Maximum diameter used for stabilisation
        dimensionedScalar dMax_;

        //- Minimum diameter used for stabilisation
        dimensionedScalar dMin_;

        //- The Sauter-mean diameter of the phase
        volScalarField d_;

        //- The width of the size distribution
        volScalarField sigma_;

        //- Coalescence switch
        Switch coalescence_;

        //- Breakup switch
        Switch breakup_;

        //- Scaled number concentration coalescence rate
        volScalarField lambdaCoaRate_;

        //- Scaled interfacial area concentration coalescence rate
        volScalarField kappaCoaRate_;

        //- Scaled squared volume concentration coalescence rate
        volScalarField betaCoaRate_;

        //- Scaled number concentration breakage rate
        volScalarField lambdaBreakRate_;

        //- Scaled interfacial area concentration breakage rate
        volScalarField kappaBreakRate_;

        //- Scaled squared volume concentration breakage rate
        volScalarField betaBreakRate_;

        //- Coalescence efficiency model
        autoPtr<coalescenceEfficiencyModel> coaEffModel_;

        //- Coalescence frequency model
        autoPtr<coalescenceFrequencyModel> coaFreqModel_;

        //- Breakup model
        autoPtr<breakupModel> breakupModel_;


    // Private member functions

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
        virtual void updateModels();

        //- Limit moments so that mean diameters are in [dMin,dMax]. Option to
        //  prevent the boundary condition update, which is needed when using
        //  this function in the constructor
        virtual void limitMoments(bool = true);

public:

    //- Runtime type information
    TypeName("threeMomentLogNormal");


    // Constructors

        //- Construct from components
        threeMomentLogNormal
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~threeMomentLogNormal();


    // Member Functions

        //- Return the local fields
        volScalarField& kappa()
        {
            return kappa_;
        }

        const volScalarField& kappa() const
        {
            return kappa_;
        }

        volScalarField& lambda()
        {
            return lambda_;
        }

        const volScalarField& lambda() const
        {
            return lambda_;
        }

        volScalarField& beta()
        {
            return beta_;
        }

        const volScalarField& beta() const
        {
            return beta_;
        }

        //- Evolve the moment transport equations
        virtual void correct();

        //- Get the diameter field
        virtual tmp<volScalarField> d() const;

        virtual tmp<volScalarField> a() const;

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
