/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "brooksCorey.H"
#include "fvMesh.H"
#include "fvc.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstream.H"

namespace Foam
{
    // Define the run-time selection entry
    defineTypeNameAndDebug(brooksCorey, 0);
    addToRunTimeSelectionTable(capillaryTransportModelBase, brooksCorey, dictionary);
}

// capillary pressure calculation - Brooks Corey
void Foam::brooksCorey::calculateCapillaryPressure()
{
    forAll(pCapillary_, cellI)
    {
        scalar alphaW = (wettingPhase_ == 1) ? alpha1_[cellI] : alpha2_[cellI];

        scalar Srw = residualSaturationWetting_.value();
        //scalar Srnw = residualSaturationNonWetting_.value();
        scalar lambda = max(lambda_.value(), SMALL);
        scalar Pd = entryPressure_.value();

        scalar Se = (alphaW - Srw) / (1.0 - Srw + SMALL);
        Se = max(1e-4, min(1.0 - 1e-4, Se));

        scalar maxPcLimit = maxPc_.value();

        scalar Pc_calculated;
        Pc_calculated = Pd * pow(Se, -1.0/lambda);

       /* if (Se > 1.0 - 1e-9)
        {
            Pc_calculated = 0;
        }
        else if (Se < 1e-9)
        {
            Pc_calculated = Pd;
        }
        else
        {
            Pc_calculated = Pd * pow(Se, -1.0/lambda);
        } */

        pCapillary_[cellI] = max(min(Pc_calculated, maxPcLimit), 0);
    }
    pCapillary_.correctBoundaryConditions();
}

Foam::brooksCorey::brooksCorey
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
    //residualSaturationNonWetting_
    //(
    //    "residualSaturationNonWetting", dimless, dict.getOrDefault<scalar>("residualSaturationNonWetting", 0.1)
    //),
    lambda_
    (
        "lambda", dimless, dict.getOrDefault<scalar>("lambda", 2.0)
    ),
    entryPressure_
    (
        "entryPressure", dimPressure, dict.getOrDefault<scalar>("entryPressure", 1.0)
    )
{
    if (enabled_)
    {
        Info<< "Capillary transport model enabled (Brooks-Corey):" << endl;
        Info<< "    Residual Saturation (Wetting): " << residualSaturationWetting_.value() << endl;
       // Info<< "    Residual Saturation (Non-Wetting): " << residualSaturationNonWetting_.value() << endl;
        Info<< "    Lambda (Pore Size Distribution Index): " << lambda_.value() << endl;
        Info<< "    Entry (Displacement) Pressure: " << entryPressure_.value() << " Pa" << endl;
    }

    calculateCapillaryPressure(); // Initial calculation
}

bool Foam::brooksCorey::writeData(Ostream& os) const
{
    capillaryTransportModelBase::writeData(os);
    return os.good();
}

bool Foam::brooksCorey::read(const dictionary& dict)
{
    capillaryTransportModelBase::read(dict);
    residualSaturationWetting_.value() = dict.getOrDefault<scalar>("residualSaturationWetting", residualSaturationWetting_.value());
    //residualSaturationNonWetting_.value() = dict.getOrDefault<scalar>("residualSaturationNonWetting", residualSaturationNonWetting_.value());
    lambda_.value() = dict.getOrDefault<scalar>("lambda", lambda_.value());
    entryPressure_.value() = dict.getOrDefault<scalar>("entryPressure", entryPressure_.value());
    return true;
}

Foam::Ostream& operator<<(Foam::Ostream& os, const Foam::brooksCorey& model)
{
    model.writeData(os);
    return os;
}