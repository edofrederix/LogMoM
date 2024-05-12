using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class Model>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::LogMoMInterfacialModel::evaluate
(
    const label gamma,
    tmp<GeometricField<Type, fvPatchField, volMesh>> (Model::*f)() const,
    const Model& model
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    scalar p = 0;
    scalar q = 0;

    forAll(p_, i)
    {
        p += p_[i]*pow(gamma,i);
    }

    forAll(q_, i)
    {
        q += q_[i]*pow(gamma,i);
    }

    // Cache the diameter and velocity fields so we can manipulate and then
    // restore them. Manipulating the velocity field is required because it may
    // be used implicitly by the provided f() function.

    volScalarField& dRef = diameterModelPtr_->d_;
    volVectorField& URef =
        const_cast<phaseModel&>(diameterModelPtr_->phase()).URef();

    tmp<volScalarField> tdOrig(new volScalarField(dRef));
    tmp<volVectorField> tUOrig(new volVectorField(URef));

    if
    (
        integrationMode_ == integrationMode::quadrature
     || integrationMode_ == integrationMode::Riemann
    )
    {
        FieldType F
        (
            FieldType::New
            (
                "F",
                interface_.mesh(),
                dimensioned<Type>(dimless, Zero)
            )
        );

        volScalarField M
        (
            volScalarField::New
            (
                "M",
                interface_.mesh(),
                pow(dimLength,p)
            )
        );

        F = Zero;
        M = Zero;

        const volScalarField sigma
        (
            diameterModelPtr_->sigma()
        );

        const dimensionedScalar& dMin = diameterModelPtr_->dMin();
        const dimensionedScalar& dMax = diameterModelPtr_->dMax();

        const volScalarField dcm
        (
            min
            (
                max
                (
                    diameterModelPtr_->dsm()
                  * Foam::exp(-2.5*Foam::sqr(sigma)),
                    dMin
                ),
                dMax
            )
        );

        tmp<volScalarField> td0, td2, td3;
        tmp<volVectorField> tU0, tU2, tU3;

        if (blendVelocity_)
        {
            td3 = diameterModelPtr_->d(3,0);
            td2 = min(diameterModelPtr_->d(2,0), td3());
            td0 = min(diameterModelPtr_->d(0,0), td2());

            tU0 = this->U(0);
            tU2 = this->U(2);
            tU3 = this->U(3);
        }
        else if (gamma != 3)
        {
            URef = this->U(gamma);
        }

        // Integrate d^p*f()*n(d) using Guass-Hermite quadrature. Limit the
        // evaluation diameters for stability.

        if (integrationMode_ == integrationMode::quadrature)
        {
            for (label i = 0; i < quadrature_->size(); i++)
            {
                const scalar wi(quadrature_->w()[i]);
                const scalar x(quadrature_->x()[i]);

                dRef = exp(x*sqrt(2.0)*sigma)*dcm;

                if (blendVelocity_)
                {
                    URef =
                        this->UBlend(dRef,td0(),td2(),td3(),tU0(),tU2(),tU3());
                }

                const FieldType fi((model.*f)());

                if (i == 0)
                {
                    F.dimensions().reset(fi.dimensions()*pow(dimLength,p));
                }

                F += wi*fi*pow(dRef,p);
                M += wi*pow(dRef,p);
            }
        }
        else
        {
            const dimensionedScalar dd = (dMax-dMin)/scalar(Nr_);

            for (label i = 0; i < Nr_; i++)
            {
                const dimensionedScalar dl = dMin + dd*i;
                const dimensionedScalar du = dMin + dd*(i+1);
                const dimensionedScalar d = (dl+du)/2.0;

                dRef = d;

                if (blendVelocity_)
                {
                    URef =
                        this->UBlend(dRef,td0(),td2(),td3(),tU0(),tU2(),tU3());
                }

                const FieldType fi((model.*f)());

                if (i == 0)
                {
                    F.dimensions().reset(fi.dimensions()*pow(dimLength,p));
                }

                const volScalarField Ni
                (
                    Foam::erf(Foam::log(du/dcm)/(sqrt(2.0)*sigma))
                  - Foam::erf(Foam::log(dl/dcm)/(sqrt(2.0)*sigma))
                );

                F += Ni*fi*pow(d,p);
                M += Ni*pow(d,p);
            }
        }

        if (blendVelocity_)
        {
            td0.clear();
            td2.clear();
            td3.clear();

            tU0.clear();
            tU2.clear();
            tU3.clear();
        }

        // Restore the diameter and velocity from cache

        dRef = tdOrig();
        tdOrig.clear();

        if (gamma != 3)
        {
            URef = tUOrig();
        }

        tUOrig.clear();

        return F/max(M, dimensionedScalar(M.dimensions(), SMALL));
    }
    else
    {
        // Evaluate function at the specified diameter

        dRef = diameterModelPtr_->d(p,q);

        if (gamma != 3)
        {
            URef = this->U(gamma);
        }

        const tmp<FieldType> tF(new FieldType((model.*f)()));

        // Restore the diameter and velocity from cache

        dRef = tdOrig();
        tdOrig.clear();

        if (gamma != 3)
        {
            URef = tUOrig();
        }

        tUOrig.clear();

        return tF;
    }
}

// ************************************************************************* //
