/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "capillaryTransportModelBase.H"
#include "fvMesh.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstream.H"

namespace Foam
{
    defineTypeNameAndDebug(capillaryTransportModelBase, 0);
    defineRunTimeSelectionTable(capillaryTransportModelBase, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillaryTransportModelBase::capillaryTransportModelBase
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    mesh_(mesh),
    alpha1_(alpha1),
    alpha2_(alpha2),
    coeffDict_(dict),
    pCapillary_
    (
        IOobject
        (
            "pCapillary",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("pCapillary", dimPressure, 0.0)
    ),
    enabled_(dict.getOrDefault<Switch>("enabled", true)),
    wettingPhase_(dict.getOrDefault<label>("wettingPhase", 1)),
    maxPc_
    (
        "maxPc", dimPressure,
        dict.getOrDefault<scalar>("maxPc", 1e5)
    )
    
{
    if (enabled_)
    {
        Info<< "Capillary transport model base enabled." << endl;
        Info<< "    Wetting Phase: " << wettingPhase_ << endl;
    }
    else
    {
        Info<< "Capillary transport model base disabled." << endl;
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::capillaryTransportModelBase> Foam::capillaryTransportModelBase::New
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting capillary transport model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "capillaryTransportModelBase",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<capillaryTransportModelBase>(ctorPtr(mesh, alpha1, alpha2, dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::capillaryTransportModelBase::update()
{
    if (!enabled_) return;

    calculateCapillaryPressure();
}

Foam::tmp<Foam::volVectorField>
Foam::capillaryTransportModelBase::capillaryMomentumSource() const
{
    return tmp<volVectorField>
    (
        new volVectorField
        (
            IOobject
            (
                "capillaryPressureGradient",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::grad(pCapillary_)
        )
    );
}

Foam::tmp<Foam::volVectorField>
Foam::capillaryTransportModelBase::capillaryPressureGradient() const
{
    return this->capillaryMomentumSource();
}

bool Foam::capillaryTransportModelBase::writeData(Ostream& os) const
{
    return true;
}

bool Foam::capillaryTransportModelBase::read(const dictionary& dict)
{
    coeffDict_ = dict;
    enabled_ = dict.getOrDefault<Switch>("enabled", true);
    wettingPhase_ = dict.getOrDefault<label>("wettingPhase", wettingPhase_);
    return true;
}

Foam::Ostream& operator<<(Foam::Ostream& os, const Foam::capillaryTransportModelBase& model)
{
    model.writeData(os);
    return os;
}