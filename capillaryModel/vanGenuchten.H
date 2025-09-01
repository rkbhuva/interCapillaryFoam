/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

\*---------------------------------------------------------------------------*/

#ifndef vanGenuchten_H
#define vanGenuchten_H

#include "capillaryTransportModelBase.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class VanGenuchten Declaration
\*---------------------------------------------------------------------------*/

class vanGenuchten
:
    public capillaryTransportModelBase
{
    // Private Data (specific to VanGenuchten)
    dimensionedScalar residualSaturationWetting_;
    dimensionedScalar m_;
    dimensionedScalar pc0_;

    // Private Member Functions
    void calculateCapillaryPressure() override;

    //- No copy construct
    vanGenuchten(const vanGenuchten&) = delete;

    //- No copy assignment
    void operator=(const vanGenuchten&) = delete;

public:

    TypeName("vanGenuchten");

    // Constructor
    vanGenuchten
    (
        const fvMesh& mesh,
        const volScalarField& alpha1,
        const volScalarField& alpha2,
        const dictionary& dict
    );

    virtual ~vanGenuchten() = default;

    // Accessors specific to vanGenuchten
    const dimensionedScalar& residualSaturationWetting() const { return residualSaturationWetting_; }
    const dimensionedScalar& m() const { return m_; }
    const dimensionedScalar& pc0() const { return pc0_; }

    // I-O
    bool writeData(Ostream& os) const override;
    bool read(const dictionary& dict) override;
};

// Ostream operator (for Info<< model)
Foam::Ostream& operator<<(Foam::Ostream&, const vanGenuchten&);
}

#endif