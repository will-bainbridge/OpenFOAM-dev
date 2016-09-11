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

#include "particleBase.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::particleBase::particleCount_ = 0;

const Foam::scalar Foam::particleBase::trackingCorrectionTol = 1e-5;

const Foam::scalar Foam::particleBase::lambdaDistanceToleranceCoeff =
    1e3*SMALL;

const Foam::scalar Foam::particleBase::minStepFractionTol = 1e5*SMALL;

namespace Foam
{
    defineTypeNameAndDebug(particleBase, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleBase::particleBase
(
    const polyMesh& mesh,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    mesh_(mesh),
    celli_(celli),
    facei_(-1),
    stepFraction_(0.0),
    tetFacei_(tetFacei),
    tetPti_(tetPti),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleID())
{}


Foam::particleBase::particleBase
(
    const polyMesh& mesh,
    const label celli
)
:
    mesh_(mesh),
    celli_(celli),
    facei_(-1),
    stepFraction_(0.0),
    tetFacei_(-1),
    tetPti_(-1),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleID())
{}


Foam::particleBase::particleBase(const particleBase& p)
:
    mesh_(p.mesh_),
    celli_(p.celli_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    tetFacei_(p.tetFacei_),
    tetPti_(p.tetPti_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


Foam::particleBase::particleBase
(
    const particleBase& p,
    const polyMesh& mesh
)
:
    mesh_(mesh),
    celli_(p.celli_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    tetFacei_(p.tetFacei_),
    tetPti_(p.tetPti_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::particleBase::transformProperties(const tensor&)
{}


void Foam::particleBase::transformProperties(const vector&)
{}


Foam::scalar Foam::particleBase::wallImpactDistance(const vector&) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

bool Foam::operator==(const particleBase& pA, const particleBase& pB)
{
    return (pA.origProc() == pB.origProc() && pA.origId() == pB.origId());
}


bool Foam::operator!=(const particleBase& pA, const particleBase& pB)
{
    return !(pA == pB);
}


// ************************************************************************* //
