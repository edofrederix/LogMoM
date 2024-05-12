#include "populationBalanceForces.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "dragModel.H"
#include "virtualMassModel.H"
#include "liftModel.H"
#include "wallLubricationModel.H"
#include "turbulentDispersionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(populationBalanceForces, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        populationBalanceForces,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::populationBalanceForces::populationBalanceForces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    popBal_
    (
        obr_.lookupObject<Foam::diameterModels::populationBalanceModel>
        (
            dict.lookup("populationBalance")
        )
    ),
    gamma_(dict.lookupOrDefault<label>("gamma", 3))
{
    read(dict);

    const phaseSystem& fluid = popBal_.fluid();

    const word phaseName =
        popBal_.name() + (gamma_ == 3 ? "" : Foam::name(gamma_));

    // Check all phases that belong to this population but only add the force
    // field once

    forAllConstIter
    (
        HashTable<const diameterModels::velocityGroup*>,
        popBal_.velocityGroupPtrs(),
        iter
    )
    {
        const phaseModel& phase = iter()->phase();

        forAll(fluid.phases(), phasei)
        {
            const phaseModel& otherPhase = fluid.phases()[phasei];

            if (&otherPhase == &phase) continue;

            const phaseInterface interface(phase, otherPhase);

            if
            (
                !forceFields_.found(dragModel::typeName)
             && fluid.foundInterfacialModel<blendedDragModel>(interface)
            )
            {
                forceFields_.insert
                (
                    dragModel::typeName,
                    new volVectorField
                    (
                        IOobject
                        (
                            IOobject::groupName("dragForce", phaseName),
                            mesh_.time().timeName(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }

            if
            (
                !forceFields_.found(virtualMassModel::typeName)
             && fluid.foundInterfacialModel<blendedVirtualMassModel>(interface)
            )
            {
                forceFields_.insert
                (
                    virtualMassModel::typeName,
                    new volVectorField
                    (
                        IOobject
                        (
                            IOobject::groupName("virtualMassForce", phaseName),
                            mesh_.time().timeName(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }

            if
            (
                !forceFields_.found(liftModel::typeName)
             && fluid.foundInterfacialModel<blendedLiftModel>(interface)
            )
            {
                forceFields_.insert
                (
                    liftModel::typeName,
                    new volVectorField
                    (
                        IOobject
                        (
                            IOobject::groupName("liftForce", phaseName),
                            mesh_.time().timeName(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }

            if
            (
                !forceFields_.found(wallLubricationModel::typeName)
             && fluid.foundInterfacialModel
                <blendedWallLubricationModel>(interface)
            )
            {
                forceFields_.insert
                (
                    wallLubricationModel::typeName,
                    new volVectorField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "wallLubricationForce",
                                phaseName
                            ),
                            mesh_.time().timeName(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }

            if
            (
                !forceFields_.found(turbulentDispersionModel::typeName)
             && fluid.foundInterfacialModel
                <blendedTurbulentDispersionModel>(interface)
            )
            {
                forceFields_.insert
                (
                    turbulentDispersionModel::typeName,
                    new volVectorField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "turbulentDispersionForce",
                                phaseName
                            ),
                            mesh_.time().timeName(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::populationBalanceForces::~populationBalanceForces()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::populationBalanceForces::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    return true;
}

bool Foam::functionObjects::populationBalanceForces::execute()
{
    // Zero the force fields
    forAllConstIter
    (
        HashPtrTable<volVectorField>,
        forceFields_,
        forceFieldIter
    )
    {
        *forceFieldIter() = Zero;
    }

    const phaseSystem& fluid = popBal_.fluid();

    volScalarField M
    (
        volScalarField::New
        (
            "M",
            fluid.mesh(),
            dimensionedScalar(pow(dimLength, gamma_)/dimVolume, 0.0)
        )
    );

    PtrList<volScalarField> Ms(popBal_.velocityGroupPtrs().size());

    label i = 0;

    forAllConstIter
    (
        HashTable<const diameterModels::velocityGroup*>,
        popBal_.velocityGroupPtrs(),
        iter
    )
    {
        const phaseModel& phase = iter()->sizeGroups()[0].phase();

        const volScalarField d(iter()->d());

        Ms.set
        (
            i,
            new volScalarField(pow(d,gamma_-3))
        );

        M += max(phase, phase.residualAlpha()/Ms.size())*Ms[i];

        i++;
    }

    // Check all phases that belong to this population but only add the force
    // field once

    i = 0;

    forAllConstIter
    (
        HashTable<const diameterModels::velocityGroup*>,
        popBal_.velocityGroupPtrs(),
        iter
    )
    {
        const phaseModel& phase = iter()->phase();
        const volScalarField& Mi = Ms[i++];

        // Add the forces from all the interfaces which contain this phase
        forAll(fluid.phases(), phasei)
        {
            const phaseModel& otherPhase = fluid.phases()[phasei];

            if (&otherPhase == &phase) continue;

            const phaseInterface interface(phase, otherPhase);

            if (fluid.foundInterfacialModel<blendedDragModel>(interface))
            {
                *forceFields_[dragModel::typeName] +=
                    fluid.lookupInterfacialModel
                    <blendedDragModel>(interface).K()
                  * (otherPhase.U() - phase.U())
                  * popBal_.alphas()
                  * Mi/M;
            }

            if (fluid.foundInterfacialModel<blendedVirtualMassModel>(interface))
            {
                *forceFields_[virtualMassModel::typeName] +=
                    fluid.lookupInterfacialModel
                    <blendedVirtualMassModel>(interface).K()
                  * (otherPhase.DUDt() - phase.DUDt())
                  * popBal_.alphas()
                  * Mi/M;
            }

            if (fluid.foundInterfacialModel<blendedLiftModel>(interface))
            {
                *forceFields_[liftModel::typeName] +=
                    (&interface.phase1() == &phase ? -1 : +1)
                  * fluid.lookupInterfacialModel
                    <blendedLiftModel>(interface).F()
                  * popBal_.alphas()
                  * Mi/M;
            }

            if
            (
                fluid.foundInterfacialModel
                <blendedWallLubricationModel>(interface)
            )
            {
                *forceFields_[wallLubricationModel::typeName] +=
                    (&interface.phase1() == &phase ? -1 : +1)
                  * fluid.lookupInterfacialModel
                    <blendedWallLubricationModel>(interface).F()
                  * popBal_.alphas()
                  * Mi/M;
            }

            if
            (
                fluid.foundInterfacialModel
                <blendedTurbulentDispersionModel>(interface)
            )
            {
                *forceFields_[turbulentDispersionModel::typeName] +=
                    fluid.lookupInterfacialModel
                    <blendedTurbulentDispersionModel>(interface).D()
                  * fvc::grad
                    (
                        otherPhase
                      / max(phase + otherPhase, otherPhase.residualAlpha())
                    )
                  * popBal_.alphas()
                  * Mi/M;
            }
        }
    }

    return true;
}

bool Foam::functionObjects::populationBalanceForces::write()
{
    forAllConstIter
    (
        HashPtrTable<volVectorField>,
        forceFields_,
        forceFieldIter
    )
    {
        writeObject(forceFieldIter()->name());
    }

    return true;
}

// ************************************************************************* //
