#ifndef LogMoMSource_H
#define LogMoMSource_H

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

class LogMoM;

/*---------------------------------------------------------------------------*\
                         Class LogMoMSource Declaration
\*---------------------------------------------------------------------------*/

class LogMoMSource
{

protected:

    // Protected data

        //- Reference to the LogMoM this source applies to
        const LogMoM& logmom_;


public:

    //- Runtime type information
    TypeName("LogMoMSource");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LogMoMSource,
            dictionary,
            (
                const LogMoM& logmom,
                const dictionary& dict
            ),
            (logmom, dict)
        );


    //- Class used for the read-construction of
    //  PtrLists of LogMoM sources
    class iNew
    {
        const LogMoM& logmom_;

    public:

        iNew(const LogMoM& logmom)
        :
            logmom_(logmom)
        {}

        autoPtr<LogMoMSource> operator()(Istream& is) const
        {
            word type(is);
            dictionary dict(is);
            return LogMoMSource::New(type, logmom_, dict);
        }
    };


    // Constructors

        LogMoMSource(const LogMoM& logmom)
        :
            logmom_(logmom)
        {}

        autoPtr<LogMoMSource> clone() const
        {
            NotImplemented;
            return autoPtr<LogMoMSource>(nullptr);
        }


    // Selectors

        static autoPtr<LogMoMSource> New
        (
            const word& type,
            const LogMoM& logmom,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LogMoMSource()
    {}


    // Member Functions

        // Moment source term
        virtual tmp<fvScalarMatrix> R(volScalarField&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
