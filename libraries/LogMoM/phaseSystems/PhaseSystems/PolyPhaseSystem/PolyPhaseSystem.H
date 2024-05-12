#ifndef PolyPhaseSystem_H
#define PolyPhaseSystem_H

#include "phaseSystem.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blendedPolyDragModel;
class blendedPolyVirtualMassModel;
class blendedPolyLiftModel;
class blendedPolyWallLubricationModel;
class blendedPolyTurbulentDispersionModel;

template<class> class PolyPhaseModel;
template<class> class MovingPhaseModel;
template<class, class> class ThermoPhaseModel;

/*---------------------------------------------------------------------------*\
                 Class PolyPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class PolyPhaseSystem
:
    public BasePhaseSystem
{
private:

    // Private typedefs

        // Poly phase model type
        typedef PolyPhaseModel
        <
            MovingPhaseModel
            <
                ThermoPhaseModel<phaseModel, rhoThermo>
            >
        > polyPhaseModel;

        typedef HashPtrTable
        <
            volScalarField,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > KdTable;

        typedef HashPtrTable
        <
            volScalarField,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > VmTable;

        typedef HashTable
        <
            autoPtr<blendedPolyDragModel>,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > dragModelTable;

        typedef HashTable
        <
            autoPtr<blendedPolyVirtualMassModel>,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > virtualMassModelTable;

        typedef HashTable
        <
            autoPtr<blendedPolyLiftModel>,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > liftModelTable;

        typedef HashTable
        <
            autoPtr<blendedPolyWallLubricationModel>,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > wallLubricationModelTable;

        typedef HashTable
        <
            autoPtr<blendedPolyTurbulentDispersionModel>,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > turbulentDispersionModelTable;


    // Private Data

        //- Poly phase models
        phaseSystem::phaseModelPartialList polyPhaseModels_;

        //- Drag coefficients
        KdTable Kd0s_;
        KdTable Kd2s_;

        //- Virtual mass coefficients
        VmTable Vm0s_;
        VmTable Vm2s_;

        // Sub models

            //- Drag models
            dragModelTable dragModels_;

            //- Virtual mass models
            virtualMassModelTable virtualMassModels_;

            //- Lift models
            liftModelTable liftModels_;

            //- Wall lubrication models
            wallLubricationModelTable wallLubricationModels_;

            //- Turbulent dispersion models
            turbulentDispersionModelTable turbulentDispersionModels_;


protected:

    //- Dynamically cast phase model to poly. Returns nullptr if not possible
    inline polyPhaseModel* polyCast(phaseModel& phase) const
    {
        return dynamic_cast<polyPhaseModel*>(&phase);
    }

    //- Dynamically cast const phase model to poly. Returns nullptr if not
    //  possible
    inline const polyPhaseModel* polyCast(const phaseModel& phase) const
    {
        return dynamic_cast<const polyPhaseModel*>(&phase);
    }

    //- Check if the given phase is a poly phase
    inline bool isPoly(const phaseModel& phase) const
    {
        return polyCast(phase);
    }

    //- Check if the given phase interface contains a poly phase
    inline bool isPoly(const phaseInterface& interface) const
    {
        return isPoly(interface.phase1()) || isPoly(interface.phase2());
    }

    //- Check if the provided gamma value is compatible. Gamma can only be zero
    //  or two.
    inline void checkGamma(const label gamma) const
    {
        if (gamma != 0 && gamma != 2)
            FatalErrorInFunction
                << "Invalid gamma" << endl << abort(FatalError);
    }

    //- Return the dimensions of the given moment
    inline dimensionSet dim(const label gamma) const
    {
        checkGamma(gamma);
        return gamma == 0 ? dimless/dimVolume : dimArea/dimVolume;
    }

public:

    // Constructors

        //- Construct from fvMesh
        PolyPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~PolyPhaseSystem();

    // Member Functions

        //- Use explicitly and then overload with moment variant
        using BasePhaseSystem::momentumTransfer;

        //- Return the momentum transfer matrices for the cell-based algorithm.
        //  This includes implicit and explicit forces that add into the cell
        //  UEqn in the normal way.
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransfer
        (
            const label gamma
        );

        //- Return the explicit force fluxes for the cell-based algorithm, that
        //  do not depend on phase mass/volume fluxes, and can therefore be
        //  evaluated outside the corrector loop. This includes things like
        //  lift, turbulent dispersion, and wall lubrication.
        virtual PtrList<surfaceScalarField> phiFs
        (
            const label gamma,
            const PtrList<volScalarField>& rAUs
        );

        //- Return the explicit drag force fluxes for the cell-based algorithm.
        //  These depend on phase mass/volume fluxes, and must therefore be
        //  evaluated inside the corrector loop.
        virtual PtrList<surfaceScalarField> phiKdPhis
        (
            const label gamma,
            const PtrList<volScalarField>& rAUs
        ) const;

        //- Return the explicit part of the drag force for the cell-based
        //  algorithm. This is the cell-equivalent of phiKdPhis. These depend on
        //  phase velocities, and must therefore be evaluated inside the
        //  corrector loop.
        virtual PtrList<volVectorField> KdUByAs
        (
            const label gamma,
            const PtrList<volScalarField>& rAUs
        ) const;

        //- Return the flux corrections for the cell-based algorithm. These
        //  depend on phase mass/volume fluxes, and must therefore be evaluated
        //  inside the corrector loop.
        virtual PtrList<surfaceScalarField> ddtCorrByAs
        (
            const label gamma,
            const PtrList<volScalarField>& rAUs
        ) const;

        //- Solve for the phase fractions
        virtual void solve
        (
            const PtrList<volScalarField>& rAUs,
            const PtrList<surfaceScalarField>& rAUfs
        );

        //- Correct the kinematics
        virtual void correctKinematics();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PolyPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
