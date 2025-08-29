#include "interfacialGrowthMomentFvScalarFieldSource.H"
#include "addToRunTimeSelectionTable.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * *  //

Foam::scalar Foam::interfacialGrowthMomentFvScalarFieldSource::gamma() const
{
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

    return dims[dimensionSet::LENGTH]+3;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacialGrowthMomentFvScalarFieldSource::
interfacialGrowthMomentFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvScalarFieldSource(iF, dict)
{}


Foam::interfacialGrowthMomentFvScalarFieldSource::
interfacialGrowthMomentFvScalarFieldSource
(
    const interfacialGrowthMomentFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfacialGrowthMomentFvScalarFieldSource::
~interfacialGrowthMomentFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::interfacialGrowthMomentFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    NotImplemented;
    return tmp<DimensionedField<scalar, volMesh>>(nullptr);
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::interfacialGrowthMomentFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    // When the source is negative, we have shrinkage. In that case, all moment
    // rates need to be proportional to each other, including the zeroth order
    // one.

    return neg(source);
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::interfacialGrowthMomentFvScalarFieldSource::sourceCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    const scalar gamma = this->gamma();

    if (gamma == 0.0)
    {
        // The zeroth-order moment has, in principle, no phase change source.
        // Only when it is negative, but that's handled by the internalCoeff
        // function as an implicit function.

        return source*0.0;
    }
    else if (gamma == 2.0)
    {
        return pos0(source)*this->internalField()/3.0;
    }
    else
    {
        FatalErrorInFunction
            << "Invalid moment order" << endl << abort(FatalError);

        return tmp<DimensionedField<scalar, volMesh>>(nullptr);
    }
}


void Foam::interfacialGrowthMomentFvScalarFieldSource::write
(
    Ostream& os
) const
{
    fvScalarFieldSource::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        interfacialGrowthMomentFvScalarFieldSource
    );
}

// ************************************************************************* //
