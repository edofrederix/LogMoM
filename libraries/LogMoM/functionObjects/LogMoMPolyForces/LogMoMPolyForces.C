#include "LogMoMPolyForces.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

#include "polyDragModel.H"
#include "polyVirtualMassModel.H"
#include "polyLiftModel.H"
#include "polyWallLubricationModel.H"
#include "polyTurbulentDispersionModel.H"

#include "dragModel.H"
#include "virtualMassModel.H"
#include "liftModel.H"
#include "wallLubricationModel.H"
#include "turbulentDispersionModel.H"

#include "PolyPhaseModel.H"
#include "rhoReactionThermo.H"
#include "MovingPhaseModel.H"
#include "ThermoPhaseModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(LogMoMPolyForces, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        LogMoMPolyForces,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::LogMoMPolyForces::LogMoMPolyForces
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
    fluid_(mesh_.lookupObject<phaseSystem>(phaseSystem::propertiesName)),
    gamma_(dict.lookupOrDefault<label>("gamma", 3))
{
    read(dict);

    const polyPhaseModel* polyPhasePtr =
        dynamic_cast<const polyPhaseModel*>(&phase_);

    if (!polyPhasePtr)
    {
        WarningInFunction
            << "Phase " << phase_.name() << " is not a poly phase. "
            << "Calculating normal forces." << endl;
    }

    const word phaseName =
        phase_.name() + (gamma_ == 3 ? "" : Foam::name(gamma_));

    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& otherPhase = fluid_.phases()[phasei];

        if (&otherPhase == &phase_) continue;

        const phaseInterface interface(phase_, otherPhase);

        if
        (
            fluid_.foundInterfacialModel<blendedPolyDragModel>(interface)
         || fluid_.foundInterfacialModel<blendedDragModel>(interface)
        )
        {
            forceFields_.insert
            (
                polyPhasePtr
              ? polyDragModel::typeName
              : dragModel::typeName,
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
            fluid_.foundInterfacialModel
            <blendedPolyVirtualMassModel>(interface)
         || fluid_.foundInterfacialModel
            <blendedVirtualMassModel>(interface)
        )
        {
            forceFields_.insert
            (
                polyPhasePtr
              ? polyVirtualMassModel::typeName
              : virtualMassModel::typeName,
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "virtualMassForce",
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
            fluid_.foundInterfacialModel<blendedPolyLiftModel>(interface)
         || fluid_.foundInterfacialModel<blendedLiftModel>(interface)
        )
        {
            forceFields_.insert
            (
                polyPhasePtr
              ? polyLiftModel::typeName
              : liftModel::typeName,
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
            fluid_.foundInterfacialModel
            <blendedPolyWallLubricationModel>(interface)
         || fluid_.foundInterfacialModel
            <blendedWallLubricationModel>(interface)
        )
        {
            forceFields_.insert
            (
                polyPhasePtr
              ? polyWallLubricationModel::typeName
              : wallLubricationModel::typeName,
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
            fluid_.foundInterfacialModel
            <blendedPolyTurbulentDispersionModel>(interface)
         || fluid_.foundInterfacialModel
            <blendedTurbulentDispersionModel>(interface)
        )
        {
            forceFields_.insert
            (
                polyPhasePtr
              ? polyTurbulentDispersionModel::typeName
              : turbulentDispersionModel::typeName,
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

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::LogMoMPolyForces::~LogMoMPolyForces()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::LogMoMPolyForces::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    return true;
}

bool Foam::functionObjects::LogMoMPolyForces::execute()
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

    const polyPhaseModel* polyPhasePtr =
        dynamic_cast<const polyPhaseModel*>(&phase_);

    // Add the forces from all the interfaces which contain this phase
    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& otherPhase = fluid_.phases()[phasei];

        if (&otherPhase == &phase_) continue;

        const phaseInterface interface(phase_, otherPhase);

        if
        (
            polyPhasePtr
         && fluid_.foundInterfacialModel<blendedPolyDragModel>(interface)
        )
        {
            *forceFields_[polyDragModel::typeName] +=
                fluid_.lookupInterfacialModel<blendedPolyDragModel>
                (interface).K(gamma_)
              * (otherPhase.U() - polyPhasePtr->U(gamma_));
        }
        else if
        (
            fluid_.foundInterfacialModel<blendedDragModel>(interface)
        )
        {
            *forceFields_[dragModel::typeName] +=
                fluid_.lookupInterfacialModel<blendedDragModel>
                (interface).K()
              * (otherPhase.U() - phase_.U());
        }

        if
        (
            polyPhasePtr
         && fluid_.foundInterfacialModel<blendedPolyVirtualMassModel>(interface)
        )
        {
            *forceFields_[polyVirtualMassModel::typeName] +=
                fluid_.lookupInterfacialModel
                <blendedPolyVirtualMassModel>(interface).K(gamma_)
              * (otherPhase.DUDt() - polyPhasePtr->DUDt(gamma_));
        }
        else if
        (
            fluid_.foundInterfacialModel<blendedVirtualMassModel>(interface)
        )
        {
            *forceFields_[virtualMassModel::typeName] +=
                fluid_.lookupInterfacialModel
                <blendedVirtualMassModel>(interface).K()
              * (otherPhase.DUDt() - phase_.DUDt());
        }

        if
        (
            polyPhasePtr
         && fluid_.foundInterfacialModel<blendedPolyLiftModel>(interface)
        )
        {
            *forceFields_[polyLiftModel::typeName] +=
                (&interface.phase1() == &phase_ ? -1 : +1)
              * fluid_.lookupInterfacialModel
                <blendedPolyLiftModel>(interface).F(gamma_);
        }
        else if
        (
            fluid_.foundInterfacialModel<blendedLiftModel>(interface)
        )
        {
            *forceFields_[liftModel::typeName] +=
                (&interface.phase1() == &phase_ ? -1 : +1)
              * fluid_.lookupInterfacialModel
                <blendedLiftModel>(interface).F();
        }

        if
        (
            polyPhasePtr
         && fluid_.foundInterfacialModel
            <blendedPolyWallLubricationModel>(interface)
        )
        {
            *forceFields_[polyWallLubricationModel::typeName] +=
                (&interface.phase1() == &phase_ ? -1 : +1)
               *fluid_.lookupInterfacialModel
                <blendedPolyWallLubricationModel>(interface).F(gamma_);
        }
        else if
        (
            fluid_.foundInterfacialModel<blendedWallLubricationModel>(interface)
        )
        {
            *forceFields_[wallLubricationModel::typeName] +=
                (&interface.phase1() == &phase_ ? -1 : +1)
               *fluid_.lookupInterfacialModel
                <blendedWallLubricationModel>(interface).F();
        }

        if
        (
            polyPhasePtr
         && fluid_.foundInterfacialModel
            <blendedPolyTurbulentDispersionModel>(interface)
        )
        {
            if (gamma_ == 3)
            {
                *forceFields_[polyTurbulentDispersionModel::typeName] +=
                    fluid_.lookupInterfacialModel
                    <blendedPolyTurbulentDispersionModel>(interface).D(gamma_)
                  * fvc::grad
                    (
                        otherPhase
                      / max(phase_ + otherPhase, otherPhase.residualAlpha())
                    );
            }
            else
            {
                *forceFields_[polyTurbulentDispersionModel::typeName] +=
                    fluid_.lookupInterfacialModel
                    <blendedPolyTurbulentDispersionModel>(interface).D(gamma_)
                  * fvc::grad(polyPhasePtr->M(gamma_))
                  / polyPhasePtr->moment(gamma_);
            }
        }
        else if
        (
            fluid_.foundInterfacialModel
            <blendedTurbulentDispersionModel>(interface)
        )
        {
            *forceFields_[turbulentDispersionModel::typeName] +=
                fluid_.lookupInterfacialModel
                <blendedTurbulentDispersionModel>(interface).D()
              * fvc::grad
                (
                    otherPhase
                  / max(phase_ + otherPhase, otherPhase.residualAlpha())
                );
        }
    }

    return true;
}

bool Foam::functionObjects::LogMoMPolyForces::write()
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
