#include "LimitedPhaseChange.H"
#include "heatTransferLimitedPhaseChange.H"
#include "massDiffusionLimitedPhaseChange.H"
#include "homogeneousCondensation.H"
#include "homogeneousLiquidPhaseSeparation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeLimitedPhaseChangeModel(Model)                                     \
                                                                               \
    typedef Foam::fv::LimitedPhaseChange<Model> Model##LimitedPhaseChange;     \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Model##LimitedPhaseChange,                                             \
        (                                                                      \
            word(Model##LimitedPhaseChange::typeName_()) + "<"                 \
          + Model::typeName + ">"                                              \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        fvModel,                                                               \
        Model##LimitedPhaseChange,                                             \
        dictionary                                                             \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

makeLimitedPhaseChangeModel(heatTransferLimitedPhaseChange)
makeLimitedPhaseChangeModel(massDiffusionLimitedPhaseChange)
makeLimitedPhaseChangeModel(homogeneousCondensation)
makeLimitedPhaseChangeModel(homogeneousLiquidPhaseSeparation)

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
