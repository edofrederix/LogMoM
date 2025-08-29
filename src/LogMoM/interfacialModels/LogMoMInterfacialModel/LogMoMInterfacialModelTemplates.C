using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, class Model>
Foam::tmp<Foam::GeometricField<Type, GeoMesh>>
Foam::LogMoMInterfacialModel::evaluate
(
    const scalar beta,
    const scalar gamma,
    tmp<GeometricField<Type, GeoMesh>> (Model::*f)() const,
    const Model& model
) const
{
    typedef GeometricField<Type, GeoMesh> FieldType;

    const diameterModels::LogMoM& logmom = this->logmom();

    // Cache the diameter so we can manipulate and then restore it

    volScalarField& dRef =
        const_cast<diameterModels::LogMoM&>(logmom).d_;

    if (integrationMode_ == integrationMode::none)
    {
        return (model.*f)();
    }
    else if (integrationMode_ == integrationMode::quadrature)
    {
        tmp<volScalarField> tdOrig(new volScalarField(dRef));

        const volScalarField& sigma = logmom.sigma();
        const volScalarField& dcm = logmom.dcm();

        FieldType F
        (
            FieldType::New
            (
                "F",
                dispersedInterface_.mesh(),
                dimensioned<Type>(dimless, Zero)
            )
        );

        const volScalarField db(dcm*exp(beta*sqr(sigma)));

        // Integrate using centered Gauss-Hermite quadrature

        for (label i = 0; i < quadrature_->size(); i++)
        {
            const scalar wi(quadrature_->w()[i]);
            const scalar x(quadrature_->x()[i]);

            dRef = db*exp(x*sqrt(2.0)*sigma);

            const FieldType fi((model.*f)());

            if (i == 0)
            {
                F.dimensions().reset(fi.dimensions());
            }

            F += wi*fi;
        }

        F /= sqrt(pi);

        // Restore the diameter

        dRef = tdOrig();
        tdOrig.clear();

        // Return

        if (beta == gamma)
        {
            return tmp<FieldType>(new FieldType(F));
        }
        else
        {
            return
                F
              * pow(dcm, beta - gamma)
              * exp(0.5*sqr(sigma)*(beta*beta - gamma*gamma));
        }
    }
    else if (integrationMode_ == integrationMode::diameter)
    {
        // Evaluate function at the specified diameter

        tmp<volScalarField> tdOrig(new volScalarField(dRef));

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

        dRef = logmom.d(p,q);

        const FieldType F(pow(dRef, beta - gamma)*(model.*f)());

        // Restore the diameter

        dRef = tdOrig();
        tdOrig.clear();

        // Return

        return tmp<FieldType>(new FieldType(F));
    }
    else
    {
        FatalErrorInFunction
            << "Invalid integration mode" << endl << abort(FatalError);

        return tmp<FieldType>();
    }
}

template<class Type, class GeoMesh, class Model>
Foam::tmp<Foam::GeometricField<Type, GeoMesh>>
Foam::LogMoMInterfacialModel::evaluate
(
    const scalar gamma,
    tmp<GeometricField<Type, GeoMesh>> (Model::*f)() const,
    const Model& model
) const
{
    return evaluate(gamma, gamma, f, model);
}

// ************************************************************************* //
