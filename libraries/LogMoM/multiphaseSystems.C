#include "phaseSystem.H"
#include "MomentumTransferPhaseSystem.H"
#include "PolyPhaseSystem.H"
#include "OneResistanceHeatTransferPhaseSystem.H"
#include "TwoResistanceHeatTransferPhaseSystem.H"
#include "PhaseTransferPhaseSystem.H"
#include "InterfaceCompositionPhaseChangePhaseSystem.H"
#include "ThermalPhaseChangePhaseSystem.H"
#include "PopulationBalancePhaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        PhaseTransferPhaseSystem
        <
            OneResistanceHeatTransferPhaseSystem
            <
                PolyPhaseSystem
                <
                    MomentumTransferPhaseSystem<phaseSystem>
                >
            >
        >
        basicPolyPhaseSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        basicPolyPhaseSystem,
        dictionary,
        basicPolyPhaseSystem
    );

    typedef
        InterfaceCompositionPhaseChangePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                TwoResistanceHeatTransferPhaseSystem
                <
                    PolyPhaseSystem
                    <
                        MomentumTransferPhaseSystem<phaseSystem>
                    >
                >
            >
        >
        interfaceCompositionPhaseChangePolyPhaseSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        interfaceCompositionPhaseChangePolyPhaseSystem,
        dictionary,
        interfaceCompositionPhaseChangePolyPhaseSystem
    );

    typedef
        ThermalPhaseChangePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                TwoResistanceHeatTransferPhaseSystem
                <
                    PolyPhaseSystem
                    <
                        MomentumTransferPhaseSystem<phaseSystem>
                    >
                >
            >
        >
        thermalPhaseChangePolyPhaseSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        thermalPhaseChangePolyPhaseSystem,
        dictionary,
        thermalPhaseChangePolyPhaseSystem
    );
}


// ************************************************************************* //
