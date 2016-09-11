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

#include "barycentricParticle.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::barycentricParticle::particleCount_ = 0;

// <-- !!! These three tolerances do nothing in the barycentric method, and
// should be deleted. They turn up in derived classes, and (for some reason) the
// sampling library, though, so it's more work to delete them than I'm willing
// to do right now.

const Foam::scalar Foam::barycentricParticle::trackingCorrectionTol = 1e-5;

const Foam::scalar Foam::barycentricParticle::lambdaDistanceToleranceCoeff =
    1e3*SMALL;

const Foam::scalar Foam::barycentricParticle::minStepFractionTol = 1e5*SMALL;

namespace Foam
{
    defineTypeNameAndDebug(barycentricParticle, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::barycentricParticle::barycentricParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    mesh_(mesh),
    coordinates_(),
    celli_(celli),
    facei_(-1),
    stepFraction_(0.0),
    tetFacei_(tetFacei),
    tetPti_(tetPti),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleID())
{
    this->position(position);
}


Foam::barycentricParticle::barycentricParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    bool doCellFacePt
)
:
    mesh_(mesh),
    coordinates_(),
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
    this->position(position);
}


Foam::barycentricParticle::barycentricParticle(const barycentricParticle& p)
:
    mesh_(p.mesh_),
    coordinates_(p.coordinates_),
    celli_(p.celli_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    tetFacei_(p.tetFacei_),
    tetPti_(p.tetPti_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


Foam::barycentricParticle::barycentricParticle
(
    const barycentricParticle& p,
    const polyMesh& mesh
)
:
    mesh_(mesh),
    coordinates_(p.coordinates_),
    celli_(p.celli_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    tetFacei_(p.tetFacei_),
    tetPti_(p.tetPti_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::barycentricParticle::tetFaceIndices
(
    label& baseI,
    label& vertex1I,
    label& vertex2I
) const
{
    if (!hasCellFacePt())
    {
        FatalErrorInFunction
            << "Tet face indices were requested for a particle before the "
            << "tetrahedron was initialised."
            << exit(FatalError);
    }

    const Foam::face& f = mesh_.faces()[tetFacei_];

    baseI = max(0, mesh_.tetBasePtIs()[tetFacei_]);

    vertex1I = (baseI + tetPti_) % f.size();

    vertex2I = f.fcIndex(vertex1I);

    if (mesh_.faceOwner()[tetFacei_] != celli_)
    {
        Swap(vertex1I, vertex2I);
    }
}


void Foam::barycentricParticle::tetMeshIndices
(
    label& basei,
    label& vertex1i,
    label& vertex2i
) const
{
    if (!hasCellFacePt())
    {
        FatalErrorInFunction
            << "Tet mesh indices were requested for a particle before the "
            << "tetrahedron was initialised."
            << exit(FatalError);
    }

    const Foam::face& f = mesh_.faces()[tetFacei_];

    tetFaceIndices(basei, vertex1i, vertex2i);

    basei = f[basei];
    vertex1i = f[vertex1i];
    vertex2i = f[vertex2i];
}


void Foam::barycentricParticle::tetGeometry
(
    vector& centre,
    vector& base,
    vector& vertex1,
    vector& vertex2
) const
{
    if (!hasCellFacePt())
    {
        centre = vector::zero;
        base = vector(1, 0, 0);
        vertex1 = vector(0, 1, 0);
        vertex2 = vector(0, 0, 1);
    }
    else
    {
        label basei, vertex1i, vertex2i;
        tetMeshIndices(basei, vertex1i, vertex2i);

        centre = mesh_.cellCentres()[celli_];
        base = mesh_.points()[basei];
        vertex1 = mesh_.points()[vertex1i];
        vertex2 = mesh_.points()[vertex2i];
    }
}


void Foam::barycentricParticle::tetTransform
(
    vector& centre,
    tensor& A
) const
{
    vector base, vertex1, vertex2;
    tetGeometry(centre, base, vertex1, vertex2);

    A = tensor
    (
        base - centre,
        vertex1 - centre,
        vertex2 - centre
    ).T();
}


void Foam::barycentricParticle::tetReverseTransform
(
    vector& centre,
    scalar& detA,
    tensor& T
) const
{
    tensor A;
    tetTransform(centre, A);

    // <-- This transpose happens twice. It could be optimised out, but the
    //     compiler might be removing it anyway...
    A = A.T();

    detA = A.x() & (A.y() ^ A.z());

    T = tensor
    (
        A.y() ^ A.z(),
        A.z() ^ A.x(),
        A.x() ^ A.y()
    );
}


void Foam::barycentricParticle::movingTetGeometry
(
    Pair<vector>& centre,
    Pair<vector>& base,
    Pair<vector>& vertex1,
    Pair<vector>& vertex2
) const
{
    NotImplemented;
}


void Foam::barycentricParticle::movingTetTransform
(
    Pair<vector>& centre,
    Pair<tensor>& A
) const
{
    NotImplemented;
}


void Foam::barycentricParticle::movingTetReverseTransform
(
    Pair<vector>& centre,
    FixedList<scalar, 4>& detA,
    FixedList<tensor, 3>& A
) const
{
    NotImplemented;
}


void Foam::barycentricParticle::reflect()
{
    Swap(coordinates_.c(), coordinates_.d());
}


void Foam::barycentricParticle::rotate(const bool reverse)
{
    if (!reverse)
    {
        scalar temp = coordinates_.b();
        coordinates_.b() = coordinates_.c();
        coordinates_.c() = coordinates_.d();
        coordinates_.d() = temp;
    }
    else
    {
        scalar temp = coordinates_.d();
        coordinates_.d() = coordinates_.c();
        coordinates_.c() = coordinates_.b();
        coordinates_.b() = temp;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::barycentricParticle::transformProperties(const tensor&)
{}


void Foam::barycentricParticle::transformProperties(const vector&)
{}


Foam::scalar Foam::barycentricParticle::wallImpactDistance(const vector&) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

bool Foam::operator==
(
    const barycentricParticle& pA,
    const barycentricParticle& pB
)
{
    return (pA.origProc() == pB.origProc() && pA.origId() == pB.origId());
}


bool Foam::operator!=
(
    const barycentricParticle& pA,
    const barycentricParticle& pB
)
{
    return !(pA == pB);
}


// ************************************************************************* //
