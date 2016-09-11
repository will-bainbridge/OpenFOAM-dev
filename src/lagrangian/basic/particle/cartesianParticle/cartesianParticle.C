/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "cartesianParticle.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cartesianParticle, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cartesianParticle::cartesianParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    particleBase(mesh, celli, tetFacei, tetPti),
    position_(position)
{}


Foam::cartesianParticle::cartesianParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    bool doCellFacePt
)
:
    particleBase(mesh, celli),
    position_(position)
{
    if (doCellFacePt)
    {
        initCellFacePt();
    }
}


Foam::cartesianParticle::cartesianParticle(const cartesianParticle& p)
:
    particleBase(p),
    position_(p.position_)
{}


Foam::cartesianParticle::cartesianParticle
(
    const cartesianParticle& p,
    const polyMesh& mesh
)
:
    particleBase(p, mesh),
    position_(p.position_)
{}


// ************************************************************************* //
