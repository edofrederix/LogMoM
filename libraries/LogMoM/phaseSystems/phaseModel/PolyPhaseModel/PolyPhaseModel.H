#ifndef PolyPhaseModel_H
#define PolyPhaseModel_H

#include "phaseModel.H"
#include "threeMomentLogNormal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class PolyPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class> class PolyPhaseSystem;
template<class> class MomentumTransferPhaseSystem;
class phaseSystem;

template<class BasePhaseModel>
class PolyPhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Moment velocity fields
        volVectorField U0_;
        volVectorField U2_;

        //- Moment velocity fluxes
        surfaceScalarField phi0_;
        surfaceScalarField phi2_;

        //- Moment velocity volumetric fluxes
        surfaceScalarField alphaPhi0_;
        surfaceScalarField alphaPhi2_;

        //- Moment velocity mass fluxes
        surfaceScalarField alphaRhoPhi0_;
        surfaceScalarField alphaRhoPhi2_;

        //- Moment face velocity fields
        tmp<surfaceVectorField> U0f_;
        tmp<surfaceVectorField> U2f_;


private:

    // Private types

        typedef PolyPhaseSystem
        <
            MomentumTransferPhaseSystem<phaseSystem>
        > polyPhaseSystem;

    // Private static member functions

        //- Calculate and return the flux field
        tmp<surfaceScalarField> phi
        (
            const volVectorField& U,
            const word phiName
        ) const;

        //- Check if the phase system is poly
        bool isPoly(const phaseSystem& fluid) const;

public:

    // Constructors

        PolyPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~PolyPhaseModel();


    // Member Functions

        //- Correct the phase properties other than the thermo and turbulence
        virtual void correct();

        //- Correct the moment face velocities for moving meshes
        virtual void correctU0f();
        virtual void correctU2f();

        //- Access to the LogMoM diameter model
        virtual diameterModels::threeMomentLogNormal& LogMoM();

        //- Const access to the LogMoM diameter model
        virtual const diameterModels::threeMomentLogNormal& LogMoM() const;

        //- Return the number concentration
        virtual tmp<volScalarField> N() const;

        //- Return the interfacial area
        virtual tmp<volScalarField> A() const;

        //- Return the requested moment
        virtual tmp<volScalarField> M(const label) const;

        //- Return the scaled number concentration
        virtual tmp<volScalarField> lambda() const;

        //- Return the scaled interfacial area
        virtual tmp<volScalarField> kappa() const;

        //- Return the requested scaled moment
        virtual tmp<volScalarField> moment(const label) const;

        //- Access to the number concentration
        virtual volScalarField& NRef();

        //- Access to the interfacial area
        virtual volScalarField& ARef();

        //- Access to the requested moment
        virtual volScalarField& MRef(const label);

        //- Access to the scaled number concentration
        virtual volScalarField& lambdaRef();

        //- Access to the scaled interfacial area
        virtual volScalarField& kappaRef();

        //- Access to the requested scaled moment
        virtual volScalarField& momentRef(const label);

        //- Return reference to the continuous phase
        const phaseModel& continuousPhase() const
        {
            return this->LogMoM().continuousPhase();
        }

        // Momentum

            using BasePhaseModel::UEqn;
            using BasePhaseModel::U;
            using BasePhaseModel::URef;
            using BasePhaseModel::phiRef;
            using BasePhaseModel::alphaPhi;
            using BasePhaseModel::alphaPhiRef;
            using BasePhaseModel::alphaRhoPhi;
            using BasePhaseModel::alphaRhoPhiRef;
            using BasePhaseModel::Uf;
            using BasePhaseModel::UfRef;
            using BasePhaseModel::DUDt;
            using BasePhaseModel::DUDtf;


            //- Return the moment momentum equations
            virtual tmp<fvVectorMatrix> U0Eqn();
            virtual tmp<fvVectorMatrix> U2Eqn();
            virtual tmp<fvVectorMatrix> UEqn(const label);

            //- Return the moment velocities
            virtual tmp<volVectorField> U0() const
            {
                return U0_;
            }
            virtual tmp<volVectorField> U2() const
            {
                return U2_;
            }
            virtual tmp<volVectorField> U(const label) const;

            //- Access the moment velocities
            virtual volVectorField& U0Ref()
            {
                return U0_;
            }
            virtual volVectorField& U2Ref()
            {
                return U2_;
            }
            virtual volVectorField& URef(const label);


            //- Return the moment fluxes
            virtual tmp<surfaceScalarField> phi0() const
            {
                return phi0_;
            }
            virtual tmp<surfaceScalarField> phi2() const
            {
                return phi2_;
            }
            virtual tmp<surfaceScalarField> phi(const label) const;

            //- Access the moment fluxes
            virtual surfaceScalarField& phi0Ref()
            {
                return phi0_;
            }
            virtual surfaceScalarField& phi2Ref()
            {
                return phi2_;
            }
            virtual surfaceScalarField& phiRef(const label);

            //- Return the moment volumetric fluxes
            virtual tmp<surfaceScalarField> alphaPhi0() const
            {
                return alphaPhi0_;
            }
            virtual tmp<surfaceScalarField> alphaPhi2() const
            {
                return alphaPhi2_;
            }
            virtual tmp<surfaceScalarField> alphaPhi(const label) const;

            //- Access the moment volumetric fluxes
            virtual surfaceScalarField& alphaPhi0Ref()
            {
                return alphaPhi0_;
            }
            virtual surfaceScalarField& alphaPhi2Ref()
            {
                return alphaPhi2_;
            }
            virtual surfaceScalarField& alphaPhiRef(const label);

            //- Return the moment mass fluxes
            virtual tmp<surfaceScalarField> alphaRhoPhi0() const
            {
                return alphaRhoPhi0_;
            }
            virtual tmp<surfaceScalarField> alphaRhoPhi2() const
            {
                return alphaRhoPhi2_;
            }
            virtual tmp<surfaceScalarField> alphaRhoPhi(const label) const;

            //- Access the moment mass fluxes
            virtual surfaceScalarField& alphaRhoPhi0Ref()
            {
                return alphaRhoPhi0_;
            }
            virtual surfaceScalarField& alphaRhoPhi2Ref()
            {
                return alphaRhoPhi2_;
            }
            virtual surfaceScalarField& alphaRhoPhiRef(const label);

            //- Return the moment face velocity
            //  Required for moving mesh cases
            virtual tmp<surfaceVectorField> U0f() const;
            virtual tmp<surfaceVectorField> U2f() const;
            virtual tmp<surfaceVectorField> Uf(const label) const;

            //- Access the moment face velocity
            //  Required for moving mesh cases
            virtual surfaceVectorField& U0fRef();
            virtual surfaceVectorField& U2fRef();
            virtual surfaceVectorField& UfRef(const label);

            //- Return the substantive acceleration
            virtual tmp<volVectorField> DUDt(const label) const;

            //- Return the substantive acceleration on the faces
            virtual tmp<surfaceScalarField> DUDtf(const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PolyPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
