#include "LimitedPhaseChange.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BaseModel>
void Foam::fv::LimitedPhaseChange<BaseModel>::correctMDot(Switch info) const
{
    volScalarField::Internal& mDot =
        phase1_.mesh().lookupObjectRef<volScalarField::Internal>(mDotName_);

    const dimensionedScalar deltaT = phase1_.mesh().time().deltaT();

    // Limit the interfacial mass rate of change to a rate that is smaller than
    // what can be exchanged in a time step. Positive mDot is defined to be from
    // phase1 to phase2.

    mDot =
        max
        (
            min
            (
                mDot,
                (phase1_*phase1_.rho()/deltaT/mDotLimiter_)().internalField()
            ),
          - (phase2_*phase2_.rho()/deltaT/mDotLimiter_)().internalField()
        );

    if (info)
    {
        Info<< "correction:" << endl << incrIndent;

        Info<< indent << "min/mean/max mDot"
            << " = " << gMin(mDot.primitiveField())
            << '/' << gAverage(mDot.primitiveField())
            << '/' << gMax(mDot.primitiveField())
            << endl;

        Info<< decrIndent;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BaseModel>
Foam::fv::LimitedPhaseChange<BaseModel>::LimitedPhaseChange
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    BaseModel(name, modelType, mesh, dict),
    mDotLimiter_(dict.lookupOrDefault<scalar>("mDotLimiter", 4.0)),
    mDotName_(name + ":mDot"),
    fluid_(this->mesh().template lookupObject<phaseSystem>("phaseProperties")),
    phase1_(fluid_.phases()[this->phaseNames().first()]),
    phase2_(fluid_.phases()[this->phaseNames().second()])
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BaseModel>
void Foam::fv::LimitedPhaseChange<BaseModel>::correct()
{
    BaseModel::correct();
    this->correctMDot(true);
}

template<class BaseModel>
Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::LimitedPhaseChange<BaseModel>::mDot() const
{
    correctMDot();
    return BaseModel::mDot();
}

// ************************************************************************* //
