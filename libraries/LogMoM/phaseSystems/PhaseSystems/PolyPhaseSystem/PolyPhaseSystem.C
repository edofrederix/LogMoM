#include "PolyPhaseSystem.H"

#include "polyDragModel.H"
#include "polyVirtualMassModel.H"
#include "polyLiftModel.H"
#include "polyWallLubricationModel.H"
#include "polyTurbulentDispersionModel.H"

#include "HashPtrTable.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcAverage.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcMeshPhi.H"
#include "fvMatrix.H"

#include "fvModels.H"
#include "fvConstraints.H"

#include "pimpleControl.H"

#include "convectionScheme.H"

#include "rhoThermo.H"
#include "rhoReactionThermo.H"
#include "threeMomentLogNormal.H"
#include "PolyPhaseModel.H"
#include "MovingPhaseModel.H"
#include "ThermoPhaseModel.H"

#include "constrainHbyA.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PolyPhaseSystem<BasePhaseSystem>::PolyPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    // Create poly phase model list

    label nPolyPhases = 0;

    List<polyPhaseModel*> polyPhases;

    forAll(this->movingPhaseModels_, phasei)
    {
        polyPhaseModel* phase =
            this->polyCast(this->movingPhaseModels_[phasei]);

        if (phase)
        {
            nPolyPhases++;
            polyPhases.append(phase);
        }
    }

    polyPhaseModels_.resize(nPolyPhases);

    forAll(polyPhaseModels_, phasei)
    {
        polyPhaseModels_.set(phasei, polyPhases[phasei]);

        Info<< "Added poly phase " << polyPhaseModels_[phasei].name() << endl;
    }

    polyPhases.clear();

    // Generate interfacial models. We need to specify the interfacial
    // dictionary explicitly.

    this->generateInterfacialModels
    (
        this->template interfacialDict<dictionary>
        (
            this->template modelName<dragModel>()
        ),
        dragModels_
    );

    this->generateInterfacialModels
    (
        this->template interfacialDict<dictionary>
        (
            this->template modelName<virtualMassModel>()
        ),
        virtualMassModels_
    );

    this->generateInterfacialModels
    (
        this->template interfacialDict<dictionary>
        (
            this->template modelName<liftModel>()
        ),
        liftModels_
    );

    this->generateInterfacialModels
    (
        this->template interfacialDict<dictionary>
        (
            this->template modelName<wallLubricationModel>()
        ),
        wallLubricationModels_
    );

    this->generateInterfacialModels
    (
        this->template interfacialDict<dictionary>
        (
            this->template modelName<turbulentDispersionModel>()
        ),
        turbulentDispersionModels_
    );

    forAllConstIter
    (
        dragModelTable,
        dragModels_,
        dragModelIter
    )
    if (isPoly(dragModelIter()->interface()))
    {
        const phaseInterface& interface = dragModelIter()->interface();

        Kd0s_.insert
        (
            dragModelIter.key(),
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Kd0", interface.name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(dragModel::dimK, 0)
            )
        );

        Kd2s_.insert
        (
            dragModelIter.key(),
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Kd2", interface.name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(dragModel::dimK, 0)
            )
        );
    }

    forAllConstIter
    (
        virtualMassModelTable,
        virtualMassModels_,
        virtualMassModelIter
    )
    if (isPoly(virtualMassModelIter()->interface()))
    {
        const phaseInterface& interface = virtualMassModelIter()->interface();

        Vm0s_.insert
        (
            interface,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Vm0", interface.name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(virtualMassModel::dimK, 0)
            )
        );

        Vm2s_.insert
        (
            interface,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Vm2", interface.name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(virtualMassModel::dimK, 0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PolyPhaseSystem<BasePhaseSystem>::
~PolyPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PolyPhaseSystem<BasePhaseSystem>::momentumTransfer(const label gamma)
{
    this->checkGamma(gamma);

    KdTable& Kds = gamma == 0.0 ? Kd0s_ : Kd2s_;
    VmTable& Vms = gamma == 0.0 ? Vm0s_ : Vm2s_;

    // Create a momentum transfer matrix for each poly phase
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr
    (
        new phaseSystem::momentumTransferTable()
    );

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    forAll(polyPhaseModels_, phasei)
    {
        polyPhaseModel& phase =
            *this->polyCast(polyPhaseModels_[phasei]);

        eqns.insert
        (
            phase.name(),
            new fvVectorMatrix
            (
                phase.U(gamma),
                phase.MRef(gamma).dimensions()*dimMass*dimVelocity/dimTime
            )
        );
    }

    // Update the drag coefficients
    forAllConstIter
    (
        dragModelTable,
        dragModels_,
        dragModelIter
    )
    if (isPoly(dragModelIter()->interface()))
    {
        *Kds[dragModelIter.key()] = dragModelIter()->K(gamma);
    }

    // Add the implicit part of the drag force
    forAllConstIter(KdTable, Kds, KdIter)
    {
        const volScalarField& K(*KdIter());
        const phaseInterface interface(*this, KdIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            if (!iter().stationary() && isPoly(iter()))
            {
                const polyPhaseModel& phase = *polyCast(iter());

                fvVectorMatrix& eqn = *eqns[phase.name()];

                eqn -= fvm::Sp(phase.moment(gamma)*K, eqn.psi());
            }
        }
    }

    // Update the virtual mass coefficients
    forAllConstIter
    (
        virtualMassModelTable,
        virtualMassModels_,
        virtualMassModelIter
    )
    if (isPoly(virtualMassModelIter()->interface()))
    {
        *Vms[virtualMassModelIter.key()] = virtualMassModelIter()->K(gamma);
    }

    // Add the virtual mass force
    forAllConstIter(VmTable, Vms, VmIter)
    {
        const volScalarField& Vm(*VmIter());
        const phaseInterface interface(*this, VmIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            if (!iter().stationary() && isPoly(iter()))
            {
                const polyPhaseModel& phase = *polyCast(iter());
                const phaseModel& otherPhase = iter.otherPhase();

                fvVectorMatrix& eqn = *eqns[phase.name()];

                const volVectorField& U = phase.U(gamma);
                const surfaceScalarField& phi = phase.phi(gamma);
                const tmp<surfaceScalarField> taphi(fvc::absolute(phi, U));
                const surfaceScalarField& aphi(taphi());

                eqn -=
                    phase.moment(gamma)*Vm
                  * (
                        fvm::ddt(U)
                      + fvm::div(aphi, U) - fvm::Sp(fvc::div(aphi), U)
                      - otherPhase.DUDt()
                    )
                  + this->MRF_.DDt(phase.moment(gamma)*Vm, U - otherPhase.U());
            }
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::PolyPhaseSystem<BasePhaseSystem>::phiFs
(
    const label gamma,
    const PtrList<volScalarField>& rAUs
)
{
    this->checkGamma(gamma);

    PtrList<surfaceScalarField> phiFs(this->phaseModels_.size());

    // Add the lift force
    forAllConstIter
    (
        liftModelTable,
        liftModels_,
        liftModelIter
    )
    if (isPoly(liftModelIter()->interface()))
    {
        const phaseInterface& interface = liftModelIter()->interface();

        const volVectorField F(liftModelIter()->F(gamma));

        if (isPoly(interface.phase1()))
        {
            const polyPhaseModel& phase = *polyCast(interface.phase1());

            addField
            (
                phase,
                "phiF",
                fvc::flux(rAUs[phase.index()]*F*phase.moment(gamma)),
                phiFs
            );
        }

        if (isPoly(interface.phase2()))
        {
            const polyPhaseModel& phase = *polyCast(interface.phase2());

            addField
            (
                phase,
                "phiF",
              - fvc::flux(rAUs[phase.index()]*F*phase.moment(gamma)),
                phiFs
            );
        }
    }

    // Add the wall lubrication force
    forAllConstIter
    (
        wallLubricationModelTable,
        wallLubricationModels_,
        wallLubricationModelIter
    )
    if (isPoly(wallLubricationModelIter()->interface()))
    {
        const phaseInterface& interface =
            wallLubricationModelIter()->interface();

        const volVectorField F(wallLubricationModelIter()->F(gamma));

        if (isPoly(interface.phase1()))
        {
            const polyPhaseModel& phase = *polyCast(interface.phase1());

            addField
            (
                phase,
                "phiF",
                fvc::flux(rAUs[phase.index()]*F*phase.moment(gamma)),
                phiFs
            );
        }

        if (isPoly(interface.phase2()))
        {
            const polyPhaseModel& phase = *polyCast(interface.phase2());

            addField
            (
                phase,
                "phiF",
              - fvc::flux(rAUs[phase.index()]*F*phase.moment(gamma)),
                phiFs
            );
        }
    }

    // Add the phase pressure
    forAll(this->polyPhaseModels_, phasei)
    {
        polyPhaseModel& phase = *polyCast(polyPhaseModels_[phasei]);

        const surfaceScalarField pPrimeByAf
        (
            fvc::interpolate(rAUs[phase.index()]*phase.pPrime())
        );

        const surfaceScalarField snGradM
        (
            fvc::snGrad(phase.MRef(gamma))*this->mesh_.magSf()
        );

        addField(phase, "phiF", pPrimeByAf*snGradM, phiFs);
    }

    // Add the turbulent dispersion force
    forAllConstIter
    (
        turbulentDispersionModelTable,
        turbulentDispersionModels_,
        turbulentDispersionModelIter
    )
    if (isPoly(turbulentDispersionModelIter()->interface()))
    {
        const phaseInterface& interface =
            turbulentDispersionModelIter()->interface();

        const volScalarField D(turbulentDispersionModelIter()->D(gamma));

        if (isPoly(interface.phase1()))
        {
            const polyPhaseModel& phase = *polyCast(interface.phase1());

            const surfaceScalarField DByA1f
            (
                fvc::interpolate(rAUs[phase.index()]*D)
            );

            const surfaceScalarField snGradM
            (
                fvc::snGrad(phase.M(gamma))*this->mesh_.magSf()
            );

            addField(phase, "phiF", DByA1f*snGradM, phiFs);
        }

        if (isPoly(interface.phase2()))
        {
            const polyPhaseModel& phase = *polyCast(interface.phase2());

            const surfaceScalarField DByA1f
            (
                fvc::interpolate(rAUs[phase.index()]*D)
            );

            const surfaceScalarField snGradM
            (
                fvc::snGrad(phase.M(gamma))*this->mesh_.magSf()
            );

            addField(phase, "phiF", DByA1f*snGradM, phiFs);
        }
    }

    if (this->fillFields_)
    {
        this->fillFields
        (
            "phiF",
            dim(gamma)*dimForce/dimDensity/dimVelocity,
            phiFs
        );
    }

    return phiFs;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::PolyPhaseSystem<BasePhaseSystem>::phiKdPhis
(
    const label gamma,
    const PtrList<volScalarField>& rAUs
) const
{
    this->checkGamma(gamma);

    const KdTable& Kds = gamma == 0.0 ? Kd0s_ : Kd2s_;

    PtrList<surfaceScalarField> phiKdPhis(this->phaseModels_.size());

    // Add the explicit part of the drag force
    forAllConstIter(KdTable, Kds, KdIter)
    {
        const volScalarField& K(*KdIter());
        const phaseInterface interface(*this, KdIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            if (isPoly(iter()))
            {
                const polyPhaseModel& phase = *polyCast(iter());
                const phaseModel& otherPhase = iter.otherPhase();

                addField
                (
                    phase,
                    "phiKdPhi",
                  - fvc::interpolate(rAUs[phase.index()]*K*phase.moment(gamma))
                  * fvc::absolute
                    (
                        this->MRF().absolute(otherPhase.phi()),
                        otherPhase.U()
                    ),
                    phiKdPhis
                );
            }
        }
    }

    if (this->fillFields_)
    {
        this->fillFields
        (
            "phiKdPhi",
            dim(gamma)*dimForce/dimDensity/dimVelocity,
            phiKdPhis
        );
    }

    return phiKdPhis;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volVectorField>
Foam::PolyPhaseSystem<BasePhaseSystem>::KdUByAs
(
    const label gamma,
    const PtrList<volScalarField>& rAUs
) const
{
    this->checkGamma(gamma);

    const KdTable& Kds = gamma == 0.0 ? Kd0s_ : Kd2s_;

    PtrList<volVectorField> KdUByAs(this->phaseModels_.size());

    // Add the explicit part of the drag force
    forAllConstIter(KdTable, Kds, KdIter)
    {
        const volScalarField& K(*KdIter());
        const phaseInterface interface(*this, KdIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            if (isPoly(iter()))
            {
                const polyPhaseModel& phase = *polyCast(iter());
                const phaseModel& otherPhase = iter.otherPhase();

                addField
                (
                    phase,
                    "KdUByA",
                  - rAUs[phase.index()]*K*phase.moment(gamma)*otherPhase.U(),
                    KdUByAs
                );
            }
        }
    }

    if (this->fillFields_)
    {
        this->fillFields("KdUByA", dim(gamma)*dimVelocity, KdUByAs);
    }

    return KdUByAs;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::PolyPhaseSystem<BasePhaseSystem>::ddtCorrByAs
(
    const label gamma,
    const PtrList<volScalarField>& rAUs
) const
{
    this->checkGamma(gamma);

    PtrList<surfaceScalarField> ddtCorrByAs(this->phaseModels_.size());

    // Construct phi differences
    PtrList<surfaceScalarField> phiCorrs(this->phaseModels_.size());

    forAll(this->polyPhaseModels_, phasei)
    {
        const polyPhaseModel& phase =
            *this->polyCast(polyPhaseModels_[phasei]);

        phiCorrs.set
        (
            phase.index(),
            this->MRF().zeroFilter
            (
                (
                    phase.Uf(gamma).valid()
                  ? (this->mesh_.Sf() & phase.Uf(gamma)().oldTime())()
                  : phase.phi(gamma)().oldTime()
                )
              - fvc::flux(phase.U(gamma)().oldTime())
            )()
        );
    }

    // Add correction
    forAll(this->polyPhaseModels_, phasei)
    {
        const polyPhaseModel& phase = *polyCast(polyPhaseModels_[phasei]);

        const volScalarField& alpha = phase;

        // Apply ddtPhiCorr filter in pure(ish) phases
        surfaceScalarField alphafBar
        (
            fvc::interpolate(fvc::average(fvc::interpolate(alpha)))
        );

        tmp<surfaceScalarField> phiCorrCoeff = pos0(alphafBar - 0.99);

        surfaceScalarField::Boundary& phiCorrCoeffBf =
            phiCorrCoeff.ref().boundaryFieldRef();

        forAll(this->mesh_.boundary(), patchi)
        {
            // Set ddtPhiCorr to 0 on non-coupled boundaries
            if
            (
                !this->mesh_.boundary()[patchi].coupled()
             || isA<cyclicAMIFvPatch>(this->mesh_.boundary()[patchi])
            )
            {
                phiCorrCoeffBf[patchi] = 0;
            }
        }

        addField
        (
            phase,
            "ddtCorrByA",
          - phiCorrCoeff*phiCorrs[phase.index()]*fvc::interpolate
            (
                byDt
                (
                    alpha.oldTime()
                  * phase.moment(gamma)().oldTime()
                  * phase.rho()().oldTime()
                  * rAUs[phase.index()]
                )
            ),
            ddtCorrByAs
        );
    }

    return ddtCorrByAs;
}

template<class BasePhaseSystem>
void Foam::PolyPhaseSystem<BasePhaseSystem>::solve
(
    const PtrList<volScalarField>& rAUs,
    const PtrList<surfaceScalarField>& rAUfs
)
{
    BasePhaseSystem::solve(rAUs, rAUfs);

    // At this point alpha is updated. Update the poly-celerity fluxes.

    forAll(polyPhaseModels_, phasei)
    {
        polyPhaseModel& phase =
            *this->polyCast(polyPhaseModels_[phasei]);

        const volScalarField& alpha = phase;

        const surfaceScalarField& phi0 = phase.phi0Ref();
        const surfaceScalarField& phi2 = phase.phi2Ref();

        const word divSchemeName
        (
            "div("+phase.phiRef().name()+","+alpha.name()+")"
        );

        tmp<fv::convectionScheme<scalar>> convPhi0
        (
            fv::convectionScheme<scalar>::New
            (
                this->mesh_,
                phi0,
                this->mesh_.schemes().div(divSchemeName)
            )
        );

        tmp<fv::convectionScheme<scalar>> convPhi2
        (
            fv::convectionScheme<scalar>::New
            (
                this->mesh_,
                phi2,
                this->mesh_.schemes().div(divSchemeName)
            )
        );

        phase.alphaPhi0Ref() = convPhi0->flux(phi0, alpha);
        phase.alphaPhi2Ref() = convPhi2->flux(phi2, alpha);

        phase.alphaRhoPhi0Ref() =
            fvc::interpolate(phase.rho())*phase.alphaPhi0Ref();

        phase.alphaRhoPhi2Ref() =
            fvc::interpolate(phase.rho())*phase.alphaPhi2Ref();
    }
}

template<class BasePhaseSystem>
void Foam::PolyPhaseSystem<BasePhaseSystem>::correctKinematics()
{
    BasePhaseSystem::correctKinematics();

    // Update velocities

    PtrList<fvVectorMatrix> U0Eqns(polyPhaseModels_.size());
    PtrList<fvVectorMatrix> U2Eqns(polyPhaseModels_.size());

    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    {
        autoPtr<phaseSystem::momentumTransferTable>
            momentumTransfer0Ptr(this->momentumTransfer(0));

        autoPtr<phaseSystem::momentumTransferTable>
            momentumTransfer2Ptr(this->momentumTransfer(2));

        phaseSystem::momentumTransferTable&
            momentumTransfer0(momentumTransfer0Ptr());

        phaseSystem::momentumTransferTable&
            momentumTransfer2(momentumTransfer2Ptr());

        forAll(polyPhaseModels_, phasei)
        {
            polyPhaseModel& phase =
                *this->polyCast(polyPhaseModels_[phasei]);

            const volScalarField& rho = phase.rho();

            volScalarField& N = phase.NRef();
            volScalarField& A = phase.ARef();

            volVectorField& U0 = phase.U0Ref();
            volVectorField& U2 = phase.U2Ref();

            U0Eqns.set
            (
                phasei,
                new fvVectorMatrix
                (
                    phase.U0Eqn()
                ==
                  * momentumTransfer0[phase.name()]
                  + fvModels.source(N, rho, U0)
                )
            );

            U2Eqns.set
            (
                phasei,
                new fvVectorMatrix
                (
                    phase.U2Eqn()
                ==
                  * momentumTransfer2[phase.name()]
                  + fvModels.source(A, rho, U2)
                )
            );

            U0Eqns[phasei].relax();
            U2Eqns[phasei].relax();

            fvConstraints.constrain(U0Eqns[phasei]);
            fvConstraints.constrain(U2Eqns[phasei]);

            U0.correctBoundaryConditions();
            U2.correctBoundaryConditions();

            fvConstraints.constrain(U0);
            fvConstraints.constrain(U2);
        }
    }

    // Diagonal coefficients

    PtrList<volScalarField> rAUs0(polyPhaseModels_.size());
    PtrList<volScalarField> rAUs2(polyPhaseModels_.size());

    forAll(polyPhaseModels_, phasei)
    {
        polyPhaseModel& phase =
            *this->polyCast(polyPhaseModels_[phasei]);

        const volScalarField& alpha = phase;

        rAUs0.set
        (
            phasei,
            new volScalarField
            (
                IOobject::groupName("rAU0", phase.name()),
                1.0
              / (
                    U0Eqns[phasei].A()
                  + byDt
                    (
                        max(phase.residualAlpha() - alpha, scalar(0))
                      * phase.lambda()
                      * phase.rho()
                    )
                )
            )
        );

        rAUs2.set
        (
            phasei,
            new volScalarField
            (
                IOobject::groupName("rAU2", phase.name()),
                1.0
              / (
                    U2Eqns[phasei].A()
                  + byDt
                    (
                        max(phase.residualAlpha() - alpha, scalar(0))
                      * phase.kappa()
                      * phase.rho()
                    )
                )
            )
        );
    }

    // Phase diagonal coefficients

    PtrList<surfaceScalarField> alpharAUfs0(polyPhaseModels_.size());
    PtrList<surfaceScalarField> alpharAUfs2(polyPhaseModels_.size());

    forAll(polyPhaseModels_, phasei)
    {
        polyPhaseModel& phase =
            *this->polyCast(polyPhaseModels_[phasei]);

        const volScalarField& alpha = phase;

        alpharAUfs0.set
        (
            phasei,
            fvc::interpolate
            (
                max(alpha, phase.residualAlpha())
              * phase.lambda()
              * rAUs0[phasei]
            ).ptr()
        );

        alpharAUfs2.set
        (
            phasei,
            fvc::interpolate
            (
                max(alpha, phase.residualAlpha())
              * phase.kappa()
              * rAUs2[phasei]
            ).ptr()
        );
    }

    // Explicit force fluxes
    PtrList<surfaceScalarField> phiFs0(this->phiFs(0, rAUs0));
    PtrList<surfaceScalarField> phiFs2(this->phiFs(2, rAUs2));

    volScalarField rho("rho", this->rho());

    // Combined buoyancy and force fluxes

    const uniformDimensionedVectorField& g =
        this->mesh_.template lookupObject<uniformDimensionedVectorField>("g");

    const surfaceScalarField& ghf =
        this->mesh_.template lookupObject<surfaceScalarField>("ghf");

    PtrList<surfaceScalarField> phigFs0(polyPhaseModels_.size());
    PtrList<surfaceScalarField> phigFs2(polyPhaseModels_.size());

    {
        const surfaceScalarField ghSnGradRho
        (
            "ghSnGradRho",
            ghf*fvc::snGrad(rho)*this->mesh_.magSf()
        );

        forAll(polyPhaseModels_, phasei)
        {
            polyPhaseModel& phase =
                *this->polyCast(polyPhaseModels_[phasei]);

            phigFs0.set
            (
                phasei,
                (
                    alpharAUfs0[phasei]
                  * (
                        ghSnGradRho
                      - fvc::interpolate(phase.rho()-rho)*(g & this->mesh_.Sf())
                      - this->surfaceTension(phase)*this->mesh_.magSf()
                    )
                ).ptr()
            );

            phigFs2.set
            (
                phasei,
                (
                    alpharAUfs2[phasei]
                  * (
                        ghSnGradRho
                      - fvc::interpolate(phase.rho()-rho)*(g & this->mesh_.Sf())
                      - this->surfaceTension(phase)*this->mesh_.magSf()
                    )
                ).ptr()
            );

            phigFs0[phasei] += phiFs0[phasei];
            phigFs2[phasei] += phiFs2[phasei];
        }
    }

    const volScalarField& p_rgh =
        this->mesh_.template lookupObject<volScalarField>("p_rgh");

    // Correct velocity

    const pimpleControl& pimple =
        this->mesh().template lookupObject<pimpleControl>("solutionControl");

    for (label corr = 0; corr < pimple.nCorrPiso(); corr++)
    {
        // Predicted velocities and fluxes for each phase

        PtrList<volVectorField> HbyAs0(polyPhaseModels_.size());
        PtrList<volVectorField> HbyAs2(polyPhaseModels_.size());

        PtrList<surfaceScalarField> phiHbyAs0(polyPhaseModels_.size());
        PtrList<surfaceScalarField> phiHbyAs2(polyPhaseModels_.size());

        {
            // Correction force fluxes
            PtrList<surfaceScalarField> ddtCorrByAs0
            (
                this->ddtCorrByAs(0,rAUs0)
            );

            PtrList<surfaceScalarField> ddtCorrByAs2
            (
                this->ddtCorrByAs(2,rAUs2)
            );

            forAll(polyPhaseModels_, phasei)
            {
                polyPhaseModel& phase =
                    *this->polyCast(polyPhaseModels_[phasei]);

                const volScalarField& alpha = phase;

                HbyAs0.set
                (
                    phasei,
                    constrainHbyA
                    (
                        rAUs0[phasei]
                      * (
                            U0Eqns[phasei].H()
                          + byDt
                            (
                                max(phase.residualAlpha() - alpha, scalar(0))
                              * phase.lambda()
                              * phase.rho()
                            )
                          * phase.U0Ref().oldTime()
                        ),
                        phase.U0Ref(),
                        p_rgh
                    )
                );

                HbyAs2.set
                (
                    phasei,
                    constrainHbyA
                    (
                        rAUs2[phasei]
                      * (
                            U2Eqns[phasei].H()
                          + byDt
                            (
                                max(phase.residualAlpha() - alpha, scalar(0))
                              * phase.kappa()
                              * phase.rho()
                            )
                          * phase.U2Ref().oldTime()
                        ),
                        phase.U2Ref(),
                        p_rgh
                    )
                );

                phiHbyAs0.set
                (
                    phasei,
                    new surfaceScalarField
                    (
                        IOobject::groupName("phiHbyA0", phase.name()),
                        fvc::flux(HbyAs0[phasei])
                      - phigFs0[phasei]
                      - ddtCorrByAs0[phasei]
                    )
                );

                phiHbyAs2.set
                (
                    phasei,
                    new surfaceScalarField
                    (
                        IOobject::groupName("phiHbyA2", phase.name()),
                        fvc::flux(HbyAs2[phasei])
                      - phigFs2[phasei]
                      - ddtCorrByAs2[phasei]
                    )
                );
            }
        }

        // Add explicit drag forces and fluxes

        PtrList<volVectorField> KdUByAs0(this->KdUByAs(0, rAUs0));
        PtrList<volVectorField> KdUByAs2(this->KdUByAs(2, rAUs2));

        PtrList<surfaceScalarField> phiKdPhis0(this->phiKdPhis(0, rAUs0));
        PtrList<surfaceScalarField> phiKdPhis2(this->phiKdPhis(2, rAUs2));

        forAll(polyPhaseModels_, phasei)
        {
            const phaseModel& phase = polyPhaseModels_[phasei];

            HbyAs0[phasei] -= KdUByAs0[phase.index()];
            HbyAs2[phasei] -= KdUByAs2[phase.index()];

            phiHbyAs0[phasei] -= phiKdPhis0[phase.index()];
            phiHbyAs2[phasei] -= phiKdPhis2[phase.index()];
        }

        const surfaceScalarField mSfGradp
        (
            "mSfGradp",
          - fvc::snGrad(p_rgh)*this->mesh_.magSf()
        );

        forAll(polyPhaseModels_, phasei)
        {
            polyPhaseModel& phase =
                *this->polyCast(polyPhaseModels_[phasei]);

            phase.phi0Ref() = phiHbyAs0[phasei] + alpharAUfs0[phasei]*mSfGradp;
            phase.phi2Ref() = phiHbyAs2[phasei] + alpharAUfs2[phasei]*mSfGradp;
        }

        forAll(polyPhaseModels_, phasei)
        {
            polyPhaseModel& phase =
                *this->polyCast(polyPhaseModels_[phasei]);

            phase.U0Ref() =
                HbyAs0[phasei]
              + fvc::reconstruct
                (
                    alpharAUfs0[phasei]*mSfGradp
                  - phigFs0[phasei]
                );

            phase.U2Ref() =
                HbyAs2[phasei]
              + fvc::reconstruct
                (
                    alpharAUfs2[phasei]*mSfGradp
                  - phigFs2[phasei]
                );
        }

        forAll(polyPhaseModels_, phasei)
        {
            polyPhaseModel& phase =
                *this->polyCast(polyPhaseModels_[phasei]);

            this->MRF_.makeRelative(phase.phi0Ref());
            this->MRF_.makeRelative(phase.phi2Ref());

            fvc::makeRelative(phase.phi0Ref(), phase.U0());
            fvc::makeRelative(phase.phi2Ref(), phase.U2());
        }

        forAll(polyPhaseModels_, phasei)
        {
            polyPhaseModel& phase =
                *this->polyCast(polyPhaseModels_[phasei]);

            volVectorField& U0 = phase.U0Ref();
            volVectorField& U2 = phase.U2Ref();

            surfaceScalarField& phi0 = phase.phi0Ref();
            surfaceScalarField& phi2 = phase.phi2Ref();

            U0.correctBoundaryConditions();
            U2.correctBoundaryConditions();

            // Correct boundaries

            volVectorField::Boundary& U0bf = U0.boundaryFieldRef();
            volVectorField::Boundary& U2bf = U2.boundaryFieldRef();

            surfaceScalarField::Boundary& phi0bf = phi0.boundaryFieldRef();
            surfaceScalarField::Boundary& phi2bf = phi2.boundaryFieldRef();

            forAll(U0bf, patchi)
            {
                if (U0bf[patchi].fixesValue())
                {
                    phi0bf[patchi] =
                        U0bf[patchi]
                      & this->mesh().Sf().boundaryField()[patchi];
                }
            }

            forAll(U2bf, patchi)
            {
                if (U2bf[patchi].fixesValue())
                {
                    phi2bf[patchi] =
                        U2bf[patchi]
                      & this->mesh().Sf().boundaryField()[patchi];
                }
            }

            phase.correctU0f();
            phase.correctU2f();

            fvConstraints.constrain(U0);
            fvConstraints.constrain(U2);
        }
    }
}


template<class BasePhaseSystem>
bool Foam::PolyPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Read models ...

        return readOK;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
