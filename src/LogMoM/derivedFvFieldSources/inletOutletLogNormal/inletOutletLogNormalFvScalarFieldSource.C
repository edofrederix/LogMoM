#include "DimensionedField.H"
#include "inletOutletLogNormalFvScalarFieldSource.H"
#include "GeometricField.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::inletOutletLogNormalFvScalarFieldSource::q() const
{
    const scalar pi(constant::mathematical::pi);

    const dimensionSet dims(internalField().dimensions());

    if
    (
        dims[dimensionSet::MASS] != 0.0
     || dims[dimensionSet::TIME] != 0.0
     || dims[dimensionSet::TEMPERATURE] != 0.0
     || dims[dimensionSet::MOLES] != 0.0
     || dims[dimensionSet::CURRENT] != 0.0
     || dims[dimensionSet::LUMINOUS_INTENSITY] != 0.0
    )
    {
        FatalErrorInFunction
            << "Dimensions of this moment field are incorrect"
            << abort(FatalError);
    }

    const scalar gamma(dims[dimensionSet::LENGTH]+3);

    if (gamma == 0.0)
    {
        return 1e-6;
    }
    else if (gamma == 2.0)
    {
        return pi;
    }
    else if (gamma == 3.0)
    {
        return pi/6.0;
    }
    else
    {
        return 1.0;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inletOutletLogNormalFvScalarFieldSource::
inletOutletLogNormalFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvScalarFieldSource(iF, dict),
    sigma_(readScalar(dict.lookup("sigma"))),
    dsm_(readScalar(dict.lookup("dsm")))
{}


Foam::inletOutletLogNormalFvScalarFieldSource::
inletOutletLogNormalFvScalarFieldSource
(
    const inletOutletLogNormalFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF),
    sigma_(field.sigma_),
    dsm_(field.dsm_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inletOutletLogNormalFvScalarFieldSource::
~inletOutletLogNormalFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::inletOutletLogNormalFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    const scalar pi(constant::mathematical::pi);
    const scalar gamma(internalField().dimensions()[dimensionSet::LENGTH]+3);

    return
        DimensionedField<scalar,volMesh>::New
        (
            model.name() + ":" + this->internalField().name() + "SourceValue",
            this->internalField().mesh(),
            dimensionedScalar
            (
                this->internalField().dimensions(),
                q()*6.0/pi
              * pow(dsm_, gamma-3.0)
              * exp((0.5*sqr(gamma)-2.5*gamma+3.0)*sqr(sigma_))
            )
        );
}


Foam::tmp<Foam::scalarField>
Foam::inletOutletLogNormalFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    const scalar pi(constant::mathematical::pi);
    const scalar gamma(internalField().dimensions()[dimensionSet::LENGTH]+3);

    return
        tmp<scalarField>
        (
            new scalarField
            (
                cells.size(),
                q()*6.0/pi
              * pow(dsm_, gamma-3.0)
              * exp((0.5*sqr(gamma)-2.5*gamma+3.0)*sqr(sigma_))
            )
        );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::inletOutletLogNormalFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return neg0(source);
}


Foam::tmp<Foam::scalarField>
Foam::inletOutletLogNormalFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return neg0(source);
}


void Foam::inletOutletLogNormalFvScalarFieldSource::write(Ostream& os) const
{
    fvScalarFieldSource::write(os);
    writeEntry(os, "sigma", sigma_);
    writeEntry(os, "dsm", dsm_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        inletOutletLogNormalFvScalarFieldSource
    );
}

// ************************************************************************* //
