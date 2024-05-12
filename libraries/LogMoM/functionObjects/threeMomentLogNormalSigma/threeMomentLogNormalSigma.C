#include "threeMomentLogNormalSigma.H"
#include "addToRunTimeSelectionTable.H"
#include "threeMomentLogNormal.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(threeMomentLogNormalSigma, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        threeMomentLogNormalSigma,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::threeMomentLogNormalSigma::
threeMomentLogNormalSigma
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phase_
    (
        mesh_.lookupObject<phaseModel>
        (
            IOobject::groupName("alpha", dict.lookup("phase"))
        )
    ),
    sigma_
    (
        IOobject
        (
            IOobject::groupName("sigma", phase_.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0)
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::threeMomentLogNormalSigma::
~threeMomentLogNormalSigma()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::threeMomentLogNormalSigma::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::threeMomentLogNormalSigma::execute()
{
    const diameterModel& model = *phase_.dPtr();

    if (model.type() != diameterModels::threeMomentLogNormal::typeName)
    {
        FatalErrorInFunction
            << "Phase " << phase_.name() << " does not have "
            << diameterModels::threeMomentLogNormal::typeName
            << " as diameter model. Cannot calculate diameter." << endl
            << abort(FatalError);
    }
    else
    {
        const diameterModels::threeMomentLogNormal& LogMoM =
            static_cast<const diameterModels::threeMomentLogNormal&>(model);

        sigma_ = LogMoM.sigma();
    }

    return true;
}

bool Foam::functionObjects::threeMomentLogNormalSigma::write()
{
    return true;
}

// ************************************************************************* //
