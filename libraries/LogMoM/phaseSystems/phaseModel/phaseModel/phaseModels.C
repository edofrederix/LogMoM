#include "addToRunTimeSelectionTable.H"

#include "rhoThermo.H"
#include "rhoReactionThermo.H"

#include "combustionModel.H"

#include "phaseModel.H"
#include "ThermoPhaseModel.H"
#include "IsothermalPhaseModel.H"
#include "AnisothermalPhaseModel.H"
#include "PurePhaseModel.H"
#include "MultiComponentPhaseModel.H"
#include "InertPhaseModel.H"
#include "ReactingPhaseModel.H"
#include "MovingPhaseModel.H"
#include "PolyPhaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        AnisothermalPhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    PolyPhaseModel
                    <
                        MovingPhaseModel
                        <
                            ThermoPhaseModel<phaseModel, rhoThermo>
                        >
                    >
                >
            >
        >
        purePolyPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        purePolyPhaseModel,
        phaseSystem,
        purePolyPhaseModel
    );

    typedef
        IsothermalPhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    PolyPhaseModel
                    <
                        MovingPhaseModel
                        <
                            ThermoPhaseModel<phaseModel, rhoThermo>
                        >
                    >
                >
            >
        >
        pureIsothermalPolyPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        pureIsothermalPolyPhaseModel,
        phaseSystem,
        pureIsothermalPolyPhaseModel
    );

    typedef
        AnisothermalPhaseModel
        <
            MultiComponentPhaseModel
            <
                InertPhaseModel
                <
                    PolyPhaseModel
                    <
                        MovingPhaseModel
                        <
                            ThermoPhaseModel<phaseModel, rhoReactionThermo>
                        >
                    >
                >
            >
        >
        multiComponentPolyPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        multiComponentPolyPhaseModel,
        phaseSystem,
        multiComponentPolyPhaseModel
    );

    typedef
        IsothermalPhaseModel
        <
            MultiComponentPhaseModel
            <
                InertPhaseModel
                <
                    PolyPhaseModel
                    <
                        MovingPhaseModel
                        <
                            ThermoPhaseModel<phaseModel, rhoReactionThermo>
                        >
                    >
                >
            >
        >
        multiComponentIsothermalPolyPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        multiComponentIsothermalPolyPhaseModel,
        phaseSystem,
        multiComponentIsothermalPolyPhaseModel
    );

    typedef
        AnisothermalPhaseModel
        <
            MultiComponentPhaseModel
            <
                ReactingPhaseModel
                <
                    PolyPhaseModel
                    <
                        MovingPhaseModel
                        <
                            ThermoPhaseModel<phaseModel, rhoReactionThermo>
                        >
                    >
                >
            >
        >
        reactingPolyPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        reactingPolyPhaseModel,
        phaseSystem,
        reactingPolyPhaseModel
    );
}

// ************************************************************************* //
