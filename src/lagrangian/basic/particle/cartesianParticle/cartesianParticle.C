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

Foam::label Foam::cartesianParticle::particleCount_ = 0;

const Foam::scalar Foam::cartesianParticle::trackingCorrectionTol = 1e-5;

const Foam::scalar Foam::cartesianParticle::lambdaDistanceToleranceCoeff =
    1e3*SMALL;

const Foam::scalar Foam::cartesianParticle::minStepFractionTol = 1e5*SMALL;

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
    mesh_(mesh),
    position_(position),
    celli_(celli),
    facei_(-1),
    stepFraction_(0.0),
    tetFacei_(tetFacei),
    tetPti_(tetPti),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleID())
{}


Foam::cartesianParticle::cartesianParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    bool doCellFacePt
)
:
    mesh_(mesh),
    position_(position),
    celli_(celli),
    facei_(-1),
    stepFraction_(0.0),
    tetFacei_(-1),
    tetPti_(-1),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleID())
{
    if (doCellFacePt)
    {
        initCellFacePt();
    }
}


Foam::cartesianParticle::cartesianParticle(const cartesianParticle& p)
:
    mesh_(p.mesh_),
    position_(p.position_),
    celli_(p.celli_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    tetFacei_(p.tetFacei_),
    tetPti_(p.tetPti_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


Foam::cartesianParticle::cartesianParticle
(
    const cartesianParticle& p,
    const polyMesh& mesh
)
:
    mesh_(mesh),
    position_(p.position_),
    celli_(p.celli_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    tetFacei_(p.tetFacei_),
    tetPti_(p.tetPti_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cartesianParticle::transformProperties(const tensor&)
{}


void Foam::cartesianParticle::transformProperties(const vector&)
{}


Foam::scalar Foam::cartesianParticle::wallImpactDistance(const vector&) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

bool Foam::operator==(const cartesianParticle& pA, const cartesianParticle& pB)
{
    return (pA.origProc() == pB.origProc() && pA.origId() == pB.origId());
}


bool Foam::operator!=(const cartesianParticle& pA, const cartesianParticle& pB)
{
    return !(pA == pB);
}


// ************************************************************************* //
