#include "PolyPhaseModel.H"
#include "phaseSystem.H"

#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcGrad.H"
#include "fvcMeshPhi.H"

#undef NoRepository
#include "PolyPhaseSystem.H"
#include "MomentumTransferPhaseSystem.H"
#define NoRepository

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::PolyPhaseModel<BasePhaseModel>::phi
(
    const volVectorField& U,
    const word phiName
) const
{
    typeIOobject<surfaceScalarField> phiHeader
    (
        phiName,
        U.mesh().time().timeName(),
        U.mesh(),
        IOobject::NO_READ
    );

    if (phiHeader.headerOk())
    {
        Info<< "Reading face flux field " << phiName << endl;

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    U.mesh().time().timeName(),
                    U.mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                U.mesh()
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U.boundaryField(), patchi)
        {
            if (!U.boundaryField()[patchi].assignable())
            {
                phiTypes[patchi] = fixedValueFvPatchScalarField::typeName;
            }
        }

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    U.mesh().time().timeName(),
                    U.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::flux(U),
                phiTypes
            )
        );
    }
}

template<class BasePhaseModel>
bool Foam::PolyPhaseModel<BasePhaseModel>::isPoly(const phaseSystem& fluid)
const
{
    return isA<polyPhaseSystem>(fluid);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::PolyPhaseModel<BasePhaseModel>::PolyPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
    U0_
    (
        IOobject
        (
            IOobject::groupName("U0", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    U2_
    (
        IOobject
        (
            IOobject::groupName("U2", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    phi0_(phi(U0_, IOobject::groupName("phi0", this->name()))),
    phi2_(phi(U2_, IOobject::groupName("phi2", this->name()))),
    alphaPhi0_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi0", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), 0)
    ),
    alphaPhi2_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi2", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), 0)
    ),
    alphaRhoPhi0_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi0", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(1, 0, -1, 0, 0), 0)
    ),
    alphaRhoPhi2_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi2", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(1, 0, -1, 0, 0), 0)
    ),
    U0f_(nullptr),
    U2f_(nullptr)
{
    phi0_.writeOpt() = IOobject::AUTO_WRITE;
    phi2_.writeOpt() = IOobject::AUTO_WRITE;

    if (fluid.mesh().dynamic())
    {
        U0f_ = new surfaceVectorField
        (
            IOobject
            (
                IOobject::groupName("U0f", this->name()),
                fluid.mesh().time().timeName(),
                fluid.mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U0_)
        );

        U2f_ = new surfaceVectorField
        (
            IOobject
            (
                IOobject::groupName("U2f", this->name()),
                fluid.mesh().time().timeName(),
                fluid.mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U2_)
        );
    }

    if (this->dPtr()->type() != diameterModels::threeMomentLogNormal::typeName)
    {
        FatalErrorInFunction
            << "PolyPhaseModel selected but the diameter model is not of type "
            << diameterModels::threeMomentLogNormal::typeName
            << endl << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::PolyPhaseModel<BasePhaseModel>::~PolyPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::diameterModels::threeMomentLogNormal&
Foam::PolyPhaseModel<BasePhaseModel>::LogMoM()
{
    return
        const_cast<diameterModels::threeMomentLogNormal&>
        (
            static_cast<const diameterModels::threeMomentLogNormal&>
            (
                *this->dPtr()
            )
        );
}

template<class BasePhaseModel>
const Foam::diameterModels::threeMomentLogNormal&
Foam::PolyPhaseModel<BasePhaseModel>::LogMoM() const
{
    return
        static_cast<const diameterModels::threeMomentLogNormal&>
        (
            *this->dPtr()
        );
}

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::PolyPhaseModel<BasePhaseModel>::N() const
{
    return this->LogMoM().N();
}

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::PolyPhaseModel<BasePhaseModel>::A() const
{
    return this->LogMoM().A();
}

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::PolyPhaseModel<BasePhaseModel>::M(const label gamma) const
{
    if (gamma == 0)
    {
        return N();
    }
    else if (gamma == 2)
    {
        return A();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return N();
    }
}

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::PolyPhaseModel<BasePhaseModel>::lambda() const
{
    return this->LogMoM().lambda();
}

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::PolyPhaseModel<BasePhaseModel>::kappa() const
{
    return this->LogMoM().kappa();
}

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::PolyPhaseModel<BasePhaseModel>::moment(const label gamma) const
{
    if (gamma == 0)
    {
        return lambda();
    }
    else if (gamma == 2)
    {
        return kappa();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return lambda();
    }
}

template<class BasePhaseModel>
Foam::volScalarField&
Foam::PolyPhaseModel<BasePhaseModel>::NRef()
{
    return this->LogMoM().N();
}

template<class BasePhaseModel>
Foam::volScalarField&
Foam::PolyPhaseModel<BasePhaseModel>::ARef()
{
    return this->LogMoM().A();
}

template<class BasePhaseModel>
Foam::volScalarField&
Foam::PolyPhaseModel<BasePhaseModel>::MRef(const label gamma)
{
    if (gamma == 0)
    {
        return NRef();
    }
    else if (gamma == 2)
    {
        return ARef();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return NRef();
    }
}

template<class BasePhaseModel>
Foam::volScalarField&
Foam::PolyPhaseModel<BasePhaseModel>::lambdaRef()
{
    return this->LogMoM().lambda();
}

template<class BasePhaseModel>
Foam::volScalarField&
Foam::PolyPhaseModel<BasePhaseModel>::kappaRef()
{
    return this->LogMoM().kappa();
}

template<class BasePhaseModel>
Foam::volScalarField&
Foam::PolyPhaseModel<BasePhaseModel>::momentRef(const label gamma)
{
    if (gamma == 0)
    {
        return lambdaRef();
    }
    else if (gamma == 2)
    {
        return kappaRef();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return lambdaRef();
    }
}

template<class BasePhaseModel>
void Foam::PolyPhaseModel<BasePhaseModel>::correct()
{
    BasePhaseModel::correct();

    if (isPoly(this->fluid()))
    {
        this->fluid().MRF().correctBoundaryVelocity(U0_);
        this->fluid().MRF().correctBoundaryVelocity(U2_);
    }
    else
    {
        // If the system is not poly, set the moment velocities and fluxes to
        // the third moment ones.

        WarningInFunction
            << "Phase " << this->name() << " is a poly phase "
            << "but the phase system is not poly. Consider making either the "
            << "phase system poly, or the phase non-poly."
            << endl;

        U0_ = this->U_;
        U2_ = this->U_;

        phi0_ = this->phi_;
        phi2_ = this->phi_;

        if (this->fluid().mesh().dynamic())
        {
            U0f_.ref() = this->Uf_();
            U2f_.ref() = this->Uf_();
        }
    }
}

template<class BasePhaseModel>
void Foam::PolyPhaseModel<BasePhaseModel>::correctU0f()
{
    const fvMesh& mesh = this->fluid().mesh();

    if (mesh.dynamic())
    {
        U0f_.ref() = fvc::interpolate(U0_);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        U0f_.ref() +=
          n*(
                this->fluid().MRF().absolute(fvc::absolute(phi0_, U0_))
              / mesh.magSf()
              - (n & U0f_())
            );

        surfaceVectorField::Boundary& UfBf = U0f_.ref().boundaryFieldRef();
        const volVectorField::Boundary& UBf = U0_.boundaryField();

        forAll(mesh.boundary(), patchi)
        {
            // Remove the flux correction on AMI patches to compensate for
            // AMI non-conservation error
            if (isA<cyclicAMIFvPatch>(mesh.boundary()[patchi]))
            {
                UfBf[patchi] = UBf[patchi];
            }
        }
    }
}

template<class BasePhaseModel>
void Foam::PolyPhaseModel<BasePhaseModel>::correctU2f()
{
    const fvMesh& mesh = this->fluid().mesh();

    if (mesh.dynamic())
    {
        U2f_.ref() = fvc::interpolate(U2_);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        U2f_.ref() +=
          n*(
                this->fluid().MRF().absolute(fvc::absolute(phi2_, U2_))
              / mesh.magSf()
              - (n & U2f_())
            );

        surfaceVectorField::Boundary& UfBf = U2f_.ref().boundaryFieldRef();
        const volVectorField::Boundary& UBf = U2_.boundaryField();

        forAll(mesh.boundary(), patchi)
        {
            // Remove the flux correction on AMI patches to compensate for
            // AMI non-conservation error
            if (isA<cyclicAMIFvPatch>(mesh.boundary()[patchi]))
            {
                UfBf[patchi] = UBf[patchi];
            }
        }
    }
}

template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::PolyPhaseModel<BasePhaseModel>::U0Eqn()
{
    const volScalarField& rho = this->thermo().rho();
    const volScalarField& alpha = *this;

    const volScalarField nuEff(this->LogMoM().dispersedTurbulence().nuEff());

    const volScalarField& N = this->NRef();

    // Use non-limited scaled moment, for stability

    const volScalarField lambda
    (
        this->lambdaRef().name(),
        N/max(alpha, this->residualAlpha())
    );

    tmp<fv::convectionScheme<scalar>> conv
    (
        fv::convectionScheme<scalar>::New
        (
            this->mesh(),
            alphaRhoPhi0_,
            this->mesh().schemes().div
            (
                "div("+alphaRhoPhi0_.name()+","+lambda.name()+")"
            )
        )
    );

    const surfaceScalarField alphaLambdaRhoPhi0
    (
        "alphaLambdaRhoPhi0",
        conv->flux(alphaRhoPhi0_,lambda)
    );

    return
    (
        fvm::ddt(N, rho, U0_)
      + fvm::div(alphaLambdaRhoPhi0, U0_)
      + this->fluid().MRF().DDt(N*rho, U0_)
      + fvc::div((N*rho*nuEff)*dev2(T(fvc::grad(U0_))))
      + fvm::laplacian(N*rho*nuEff, U0_)
    );
}

template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::PolyPhaseModel<BasePhaseModel>::U2Eqn()
{
    const volScalarField& rho = this->thermo().rho();
    const volScalarField& alpha = *this;

    const volScalarField nuEff(this->LogMoM().dispersedTurbulence().nuEff());

    const volScalarField& A = this->ARef();

    // Use non-limited scaled moment, for stability

    const volScalarField kappa
    (
        this->kappaRef().name(),
        A/max(alpha, this->residualAlpha())
    );

    tmp<fv::convectionScheme<scalar>> conv
    (
        fv::convectionScheme<scalar>::New
        (
            this->mesh(),
            alphaRhoPhi2_,
            this->mesh().schemes().div
            (
                "div("+alphaRhoPhi2_.name()+","+kappa.name()+")"
            )
        )
    );

    const surfaceScalarField alphaKappaRhoPhi2
    (
        "alphaKappaRhoPhi2",
        conv->flux(alphaRhoPhi2_,kappa)
    );

    return
    (
        fvm::ddt(A, rho, U2_)
      + fvm::div(alphaKappaRhoPhi2, U2_)
      + this->fluid().MRF().DDt(A*rho, U2_)
      + fvc::div((A*rho*nuEff)*dev2(T(fvc::grad(U2_))))
      + fvm::laplacian(A*rho*nuEff, U2_)
    );
}

template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::PolyPhaseModel<BasePhaseModel>::UEqn(const label gamma)
{
    if (gamma == 0)
    {
        return U0Eqn();
    }
    else if (gamma == 2)
    {
        return U2Eqn();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return U0Eqn();
    }
}

template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::PolyPhaseModel<BasePhaseModel>::U(const label gamma) const
{
    if (gamma == 0)
    {
        return U0();
    }
    else if (gamma == 2)
    {
        return U2();
    }
    else if (gamma == 3)
    {
        return BasePhaseModel::U();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return U0();
    }
}

template<class BasePhaseModel>
Foam::volVectorField&
Foam::PolyPhaseModel<BasePhaseModel>::URef(const label gamma)
{
    if (gamma == 0)
    {
        return U0Ref();
    }
    else if (gamma == 2)
    {
        return U2Ref();
    }
    else if (gamma == 3)
    {
        return BasePhaseModel::URef();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return U0Ref();
    }
}

template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::PolyPhaseModel<BasePhaseModel>::phi(const label gamma) const
{
    if (gamma == 0)
    {
        return phi0();
    }
    else if (gamma == 2)
    {
        return phi2();
    }
    else if (gamma == 3)
    {
        return BasePhaseModel::phi();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return phi0();
    }
}

template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::PolyPhaseModel<BasePhaseModel>::phiRef(const label gamma)
{
    if (gamma == 0)
    {
        return phi0Ref();
    }
    else if (gamma == 2)
    {
        return phi2Ref();
    }
    else if (gamma == 3)
    {
        return BasePhaseModel::phiRef();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return phi0Ref();
    }
}

template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::PolyPhaseModel<BasePhaseModel>::alphaPhi(const label gamma) const
{
    if (gamma == 0)
    {
        return alphaPhi0();
    }
    else if (gamma == 2)
    {
        return alphaPhi2();
    }
    else if (gamma == 3)
    {
        return BasePhaseModel::alphaPhi();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return alphaPhi0();
    }
}

template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::PolyPhaseModel<BasePhaseModel>::alphaPhiRef(const label gamma)
{
    if (gamma == 0)
    {
        return alphaPhi0Ref();
    }
    else if (gamma == 2)
    {
        return alphaPhi2Ref();
    }
    else if (gamma == 3)
    {
        return BasePhaseModel::alphaPhiRef();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return alphaPhi0Ref();
    }
}

template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::PolyPhaseModel<BasePhaseModel>::alphaRhoPhi(const label gamma) const
{
    if (gamma == 0)
    {
        return alphaRhoPhi0();
    }
    else if (gamma == 2)
    {
        return alphaRhoPhi2();
    }
    else if (gamma == 3)
    {
        return BasePhaseModel::alphaRhoPhi();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return alphaRhoPhi0();
    }
}

template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::PolyPhaseModel<BasePhaseModel>::alphaRhoPhiRef(const label gamma)
{
    if (gamma == 0)
    {
        return alphaRhoPhi0Ref();
    }
    else if (gamma == 2)
    {
        return alphaRhoPhi2Ref();
    }
    else if (gamma == 3)
    {
        return BasePhaseModel::alphaRhoPhiRef();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return alphaRhoPhi0Ref();
    }
}

template<class BasePhaseModel>
Foam::tmp<Foam::surfaceVectorField>
Foam::PolyPhaseModel<BasePhaseModel>::U0f() const
{
    return
        U0f_.valid()
      ? tmp<surfaceVectorField>(U0f_())
      : tmp<surfaceVectorField>();
}

template<class BasePhaseModel>
Foam::tmp<Foam::surfaceVectorField>
Foam::PolyPhaseModel<BasePhaseModel>::U2f() const
{
    return
        U2f_.valid()
      ? tmp<surfaceVectorField>(U2f_())
      : tmp<surfaceVectorField>();
}

template<class BasePhaseModel>
Foam::tmp<Foam::surfaceVectorField>
Foam::PolyPhaseModel<BasePhaseModel>::Uf(const label gamma) const
{
    if (gamma == 0)
    {
        return U0f();
    }
    else if (gamma == 2)
    {
        return U2f();
    }
    else if (gamma == 3)
    {
        return BasePhaseModel::Uf();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return U0f();
    }
}

template<class BasePhaseModel>
Foam::surfaceVectorField&
Foam::PolyPhaseModel<BasePhaseModel>::U0fRef()
{
    if (U0f_.valid())
    {
        return U0f_.ref();
    }
    else
    {
        FatalErrorInFunction
            << "U0f has not been allocated."
            << exit(FatalError);

        return const_cast<surfaceVectorField&>(surfaceVectorField::null());
    }
}

template<class BasePhaseModel>
Foam::surfaceVectorField&
Foam::PolyPhaseModel<BasePhaseModel>::U2fRef()
{
    if (U2f_.valid())
    {
        return U2f_.ref();
    }
    else
    {
        FatalErrorInFunction
            << "U0f has not been allocated."
            << exit(FatalError);

        return const_cast<surfaceVectorField&>(surfaceVectorField::null());
    }
}

template<class BasePhaseModel>
Foam::surfaceVectorField&
Foam::PolyPhaseModel<BasePhaseModel>::UfRef(const label gamma)
{
    if (gamma == 0)
    {
        return U0fRef();
    }
    else if (gamma == 2)
    {
        return U2fRef();
    }
    else if (gamma == 3)
    {
        return BasePhaseModel::UfRef();
    }
    else
    {
        FatalErrorInFunction
            << "Invalid gamma" << endl << abort(FatalError);

        return U0fRef();
    }
}

template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::PolyPhaseModel<BasePhaseModel>::DUDt(const label gamma) const
{
    const surfaceScalarField& phi = this->phi(gamma);
    const volVectorField& U = this->U(gamma);

    const tmp<surfaceScalarField> taphi(fvc::absolute(phi, U));
    const surfaceScalarField& aphi(taphi());

    return fvc::ddt(U) + fvc::div(aphi, U) - fvc::div(aphi)*U;
}

template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::PolyPhaseModel<BasePhaseModel>::DUDtf(const label gamma) const
{
    const surfaceScalarField& phi = this->phi(gamma);

    return  byDt(phi - phi.oldTime());
}

// ************************************************************************* //
