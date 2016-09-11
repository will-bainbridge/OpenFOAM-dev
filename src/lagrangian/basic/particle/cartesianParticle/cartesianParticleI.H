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

#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::cartesianParticle::findTris
(
    const vector& position,
    DynamicList<label>& faceList,
    const tetPointRef& tet,
    const FixedList<vector, 4>& tetAreas,
    const FixedList<label, 4>& tetPlaneBasePtIs,
    const scalar tol
) const
{
    faceList.clear();

    const point Ct = tet.centre();

    for (label i = 0; i < 4; i++)
    {
        scalar lambda = tetLambda
        (
            Ct,
            position,
            i,
            tetAreas[i],
            tetPlaneBasePtIs[i],
            celli_,
            tetFacei_,
            tetPti_,
            tol
        );

        if ((lambda > 0.0) && (lambda < 1.0))
        {
            faceList.append(i);
        }
    }
}


inline Foam::scalar Foam::cartesianParticle::tetLambda
(
    const vector& from,
    const vector& to,
    const label triI,
    const vector& n,
    const label tetPlaneBasePtI,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const scalar tol
) const
{
    const pointField& pPts = mesh_.points();

    if (mesh_.moving())
    {
        return movingTetLambda
        (
            from,
            to,
            triI,
            n,
            tetPlaneBasePtI,
            celli,
            tetFacei,
            tetPti,
            tol
        );
    }

    const point& base = pPts[tetPlaneBasePtI];

    scalar lambdaNumerator = (base - from) & n;
    scalar lambdaDenominator = (to - from) & n;

    // n carries the area of the tet faces, so the dot product with a
    // delta-length has the units of volume.  Comparing the component of each
    // delta-length in the direction of n times the face area to a fraction of
    // the cell volume.

    if (mag(lambdaDenominator) < tol)
    {
        if (mag(lambdaNumerator) < tol)
        {
            // Track starts on the face, and is potentially
            // parallel to it.  +-tol/+-tol is not a good
            // comparison, return 0.0, in anticipation of tet
            // centre correction.
            return 0.0;
        }
        else
        {
            if (mag((to - from)) < tol/mag(n))
            {
                // 'Zero' length track (compared to the tolerance, which is
                // based on the cell volume, divided by the tet face area), not
                // along the face, face cannot be crossed.
                return GREAT;
            }
            else
            {
                // Trajectory is non-zero and parallel to face
                lambdaDenominator = sign(lambdaDenominator)*SMALL;
            }
        }
    }

    return lambdaNumerator/lambdaDenominator;
}


inline Foam::scalar Foam::cartesianParticle::movingTetLambda
(
    const vector& from,
    const vector& to,
    const label triI,
    const vector& n,
    const label tetPlaneBasePtI,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const scalar tol
) const
{
    const pointField& pPts = mesh_.points();
    const pointField& oldPPts = mesh_.oldPoints();

    // Base point of plane at end of motion
    const point& b = pPts[tetPlaneBasePtI];

    // n: Normal of plane at end of motion

    // Base point of plane at start of timestep
    const point& b00 = oldPPts[tetPlaneBasePtI];

    // Base point of plane at start of tracking portion (cast forward by
    // stepFraction)
    point b0 = b00 + stepFraction_*(b - b00);

    // Normal of plane at start of tracking portion
    vector n0 = Zero;

    {
        tetIndices tetIs(celli, tetFacei, tetPti, mesh_);

        // Tet at timestep start
        tetPointRef tet00 = tetIs.oldTet(mesh_);

        // Tet at timestep end
        tetPointRef tet = tetIs.tet(mesh_);

        point tet0PtA = tet00.a() + stepFraction_*(tet.a() - tet00.a());
        point tet0PtB = tet00.b() + stepFraction_*(tet.b() - tet00.b());
        point tet0PtC = tet00.c() + stepFraction_*(tet.c() - tet00.c());
        point tet0PtD = tet00.d() + stepFraction_*(tet.d() - tet00.d());

        // Tracking portion start tet (cast forward by stepFraction)
        tetPointRef tet0(tet0PtA, tet0PtB, tet0PtC, tet0PtD);

        switch (triI)
        {
            case 0:
            {
                n0 = tet0.Sa();
                break;
            }
            case 1:
            {
                n0 = tet0.Sb();
                break;
            }
            case 2:
            {
                n0 = tet0.Sc();
                break;
            }
            case 3:
            {
                n0 = tet0.Sd();
                break;
            }
            default:
            {
                break;
            }
        }
    }

    if (mag(n0) < SMALL)
    {
        // If the old normal is zero (for example in layer addition)
        // then use the current normal;
        n0 = n;
    }

    scalar lambdaNumerator = 0;
    scalar lambdaDenominator = 0;

    vector dP = to - from;
    vector dN = n - n0;
    vector dB = b - b0;
    vector dS = from - b0;

    if (mag(dN) > SMALL)
    {
        scalar a = (dP - dB) & dN;
        scalar b = ((dP - dB) & n0) + (dS & dN);
        scalar c = dS & n0;

        if (mag(a) > SMALL)
        {

            // Solve quadratic for lambda
            scalar discriminant = sqr(b) - 4.0*a*c;

            if (discriminant < 0)
            {
                // Imaginary roots only - face not crossed
                return GREAT;
            }
            else
            {
                scalar q = -0.5*(b + sign(b)*Foam::sqrt(discriminant));

                if (mag(q) < VSMALL)
                {
                    // If q is zero, then l1 = q/a is the required
                    // value of lambda, and is zero.
                    return 0.0;
                }

                scalar l1 = q/a;
                scalar l2 = c/q;

                // There will be two roots, a big one and a little
                // one, choose the little one.

                if (mag(l1) < mag(l2))
                {
                    return l1;
                }
                else
                {
                    return l2;
                }
            }
        }
        {
            // When a is zero, solve the first order polynomial
            lambdaNumerator = -c;
            lambdaDenominator = b;
        }
    }
    else
    {
        // When n = n0 is zero, there is no plane rotation, solve the
        // first order polynomial
        lambdaNumerator = -(dS & n0);
        lambdaDenominator = ((dP - dB) & n0);
    }

    if (mag(lambdaDenominator) < tol)
    {
        if (mag(lambdaNumerator) < tol)
        {
            // Track starts on the face, and is potentially
            // parallel to it.  +-tol)/+-tol is not a good
            // comparison, return 0.0, in anticipation of tet
            // centre correction.
            return 0.0;
        }
        else
        {
            if (mag((to - from)) < tol/mag(n))
            {
                // Zero length track, not along the face, face
                // cannot be crossed.
                return GREAT;
            }
            else
            {
                // Trajectory is non-zero and parallel to face
                lambdaDenominator = sign(lambdaDenominator)*SMALL;
            }
        }
    }

    return lambdaNumerator/lambdaDenominator;
}



inline void Foam::cartesianParticle::tetNeighbour(label triI)
{
    const labelList& pOwner = mesh_.faceOwner();
    const faceList& pFaces = mesh_.faces();

    bool own = (pOwner[tetFacei_] == celli_);

    const Foam::face& f = pFaces[tetFacei_];

    label tetBasePtI = mesh_.tetBasePtIs()[tetFacei_];

    if (tetBasePtI == -1)
    {
        FatalErrorInFunction
            << "No base point for face " << tetFacei_ << ", " << f
            << ", produces a valid tet decomposition."
            << abort(FatalError);
    }

    label facePtI = (tetPti_ + tetBasePtI) % f.size();
    label otherFacePtI = f.fcIndex(facePtI);

    switch (triI)
    {
        case 0:
        {
            // Crossing this triangle changes tet to that in the
            // neighbour cell over tetFacei

            // Modification of celli_ will happen by other indexing,
            // tetFacei_ and tetPti don't change.

            break;
        }
        case 1:
        {
            crossEdgeConnectedFace
            (
                celli_,
                tetFacei_,
                tetPti_,
                Foam::edge(f[facePtI], f[otherFacePtI])
            );

            break;
        }
        case 2:
        {
            if (own)
            {
                if (tetPti_ < f.size() - 2)
                {
                    tetPti_ = f.fcIndex(tetPti_);
                }
                else
                {
                    crossEdgeConnectedFace
                    (
                        celli_,
                        tetFacei_,
                        tetPti_,
                        Foam::edge(f[tetBasePtI], f[otherFacePtI])
                    );
                }
            }
            else
            {
                if (tetPti_ > 1)
                {
                    tetPti_ = f.rcIndex(tetPti_);
                }
                else
                {
                    crossEdgeConnectedFace
                    (
                        celli_,
                        tetFacei_,
                        tetPti_,
                        Foam::edge(f[tetBasePtI], f[facePtI])
                    );
                }
            }

            break;
        }
        case 3:
        {
            if (own)
            {
                if (tetPti_ > 1)
                {
                    tetPti_ = f.rcIndex(tetPti_);
                }
                else
                {
                    crossEdgeConnectedFace
                    (
                        celli_,
                        tetFacei_,
                        tetPti_,
                        Foam::edge(f[tetBasePtI], f[facePtI])
                    );
                }
            }
            else
            {
                if (tetPti_ < f.size() - 2)
                {
                    tetPti_ = f.fcIndex(tetPti_);
                }
                else
                {
                    crossEdgeConnectedFace
                    (
                        celli_,
                        tetFacei_,
                        tetPti_,
                        Foam::edge(f[tetBasePtI], f[otherFacePtI])
                    );
                }
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Tet tri face index error, can only be 0..3, supplied "
                << triI << abort(FatalError);

            break;
        }
    }
}


inline void Foam::cartesianParticle::crossEdgeConnectedFace
(
    const label& celli,
    label& tetFacei,
    label& tetPti,
    const edge& e
)
{
    const faceList& pFaces = mesh_.faces();
    const cellList& pCells = mesh_.cells();

    const Foam::face& f = pFaces[tetFacei];

    const Foam::cell& thisCell = pCells[celli];

    forAll(thisCell, cFI)
    {
        // Loop over all other faces of this cell and
        // find the one that shares this edge

        label fI = thisCell[cFI];

        if (tetFacei == fI)
        {
            continue;
        }

        const Foam::face& otherFace = pFaces[fI];

        label edDir = otherFace.edgeDirection(e);

        if (edDir == 0)
        {
            continue;
        }
        else if (f == pFaces[fI])
        {
            // This is a necessary condition if using duplicate baffles
            // (so coincident faces). We need to make sure we don't cross into
            // the face with the same vertices since we might enter a tracking
            // loop where it never exits. This test should be cheap
            // for most meshes so can be left in for 'normal' meshes.
            continue;
        }
        else
        {
            //Found edge on other face
            tetFacei = fI;

            label eIndex = -1;

            if (edDir == 1)
            {
                // Edge is in the forward circulation of this face, so
                // work with the start point of the edge
                eIndex = findIndex(otherFace, e.start());
            }
            else
            {
                // edDir == -1, so the edge is in the reverse
                // circulation of this face, so work with the end
                // point of the edge
                eIndex = findIndex(otherFace, e.end());
            }

            label tetBasePtI = mesh_.tetBasePtIs()[fI];

            if (tetBasePtI == -1)
            {
                FatalErrorInFunction
                    << "No base point for face " << fI << ", " << f
                    << ", produces a decomposition that has a minimum "
                    << "volume greater than tolerance."
                    << abort(FatalError);
            }

            // Find eIndex relative to the base point on new face
            eIndex -= tetBasePtI;

            if (neg(eIndex))
            {
                eIndex = (eIndex + otherFace.size()) % otherFace.size();
            }

            if (eIndex == 0)
            {
                // The point is the base point, so this is first tet
                // in the face circulation
                tetPti = 1;
            }
            else if (eIndex == otherFace.size() - 1)
            {
                // The point is the last before the base point, so
                // this is the last tet in the face circulation
                tetPti = otherFace.size() - 2;
            }
            else
            {
                tetPti = eIndex;
            }

            break;
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::vector& Foam::cartesianParticle::position() const
{
    return position_;
}


inline void Foam::cartesianParticle::position(const vector& position)
{
    position_ = position;
}


inline void Foam::cartesianParticle::initCellFacePt()
{
    particleBase::initCellFacePt(position_);
}


// ************************************************************************* //
