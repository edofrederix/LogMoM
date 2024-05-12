#include "threeMomentLogNormalDiameter.H"
#include "addToRunTimeSelectionTable.H"
#include "threeMomentLogNormal.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(threeMomentLogNormalDiameter, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        threeMomentLogNormalDiameter,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::threeMomentLogNormalDiameter::
threeMomentLogNormalDiameter
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    p_(dict.lookupOrDefault<scalar>("p", 3.0)),
    q_(dict.lookupOrDefault<scalar>("q", 2.0)),
    phase_
    (
        mesh_.lookupObject<phaseModel>
        (
            IOobject::groupName("alpha", dict.lookup("phase"))
        )
    ),
    d_
    (
        IOobject
        (
            IOobject::groupName
            (
                "d("+Foam::name(p_)+","+Foam::name(q_)+")", phase_.name()
            ),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength, 0)
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::threeMomentLogNormalDiameter::
~threeMomentLogNormalDiameter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::threeMomentLogNormalDiameter::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::threeMomentLogNormalDiameter::execute()
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

        d_ = LogMoM.d(p_,q_);
    }

    return true;
}

bool Foam::functionObjects::threeMomentLogNormalDiameter::write()
{
    return true;
}

// ************************************************************************* //
