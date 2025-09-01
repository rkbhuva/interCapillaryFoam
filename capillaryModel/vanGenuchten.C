/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "vanGenuchten.H"
#include "fvMesh.H"
#include "fvc.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstream.H"

namespace Foam
{
    // Define the run-time selection entry
    defineTypeNameAndDebug(vanGenuchten, 0);
    addToRunTimeSelectionTable(capillaryTransportModelBase, vanGenuchten, dictionary);
}

// capillary pressure calculation - vanGenuchten
void Foam::vanGenuchten::calculateCapillaryPressure()
{
    forAll(pCapillary_, cellI)
    {
        scalar alphaW = (wettingPhase_ == 1) ? alpha1_[cellI] : alpha2_[cellI];

        scalar Srw = residualSaturationWetting_.value();
        scalar m = max(m_.value(), SMALL);
        scalar p0 = pc0_.value();

        scalar Se = (alphaW - Srw) / (1.0 - Srw + SMALL);
        Se = max(1e-4, min(1.0 - 1e-4, Se));

        scalar maxPcLimit = maxPc_.value();

        scalar Pc_calculated;
        Pc_calculated = p0 * pow(((pow(Se, -1.0/m)) - 1.0), 1-m);

        pCapillary_[cellI] = max(min(Pc_calculated, maxPcLimit), 0.0);
    }
    pCapillary_.correctBoundaryConditions();
}

Foam::vanGenuchten::vanGenuchten
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    capillaryTransportModelBase(mesh, alpha1, alpha2, dict), // constructor
    residualSaturationWetting_
    (
        "residualSaturationWetting", dimless, dict.getOrDefault<scalar>("residualSaturationWetting", 0.05)
    ),
    m_
    (
        "m", dimless, dict.getOrDefault<scalar>("m", 0.5)
    ),
    pc0_
    (
        "pc0", dimPressure, dict.getOrDefault<scalar>("pc0", 1.0)
    )
{
    if (enabled_)
    {
        Info<< "Capillary transport model enabled (vanGenuchten):" << endl;
        Info<< "    Residual Saturation (Wetting): " << residualSaturationWetting_.value() << endl;
        Info<< "    m : " << m_.value() << endl;
        Info<< "    Entry (Displacement) Pressure: " << pc0_.value() << " Pa" << endl;
    }

    calculateCapillaryPressure(); // Initial calculation
}

bool Foam::vanGenuchten::writeData(Ostream& os) const
{
    capillaryTransportModelBase::writeData(os);
    return os.good();
}

bool Foam::vanGenuchten::read(const dictionary& dict)
{
    capillaryTransportModelBase::read(dict);
    residualSaturationWetting_.value() = dict.getOrDefault<scalar>("residualSaturationWetting", residualSaturationWetting_.value());
    m_.value() = dict.getOrDefault<scalar>("m", m_.value());
    pc0_.value() = dict.getOrDefault<scalar>("pc0", pc0_.value());
    return true;
}

Foam::Ostream& operator<<(Foam::Ostream& os, const Foam::vanGenuchten& model)
{
    model.writeData(os);
    return os;
}