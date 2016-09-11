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

#include "IOPosition.H"

#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "processorPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wallPolyPatch.H"
#include "wedgePolyPatch.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackData>
Foam::label Foam::barycentricParticle::track
(
    const vector& endPosition,
    TrackData& td
)
{
    facei_ = -1;

    // Tracks to endPosition or stop on boundary
    while (!onBoundary() && stepFraction_ < 1.0 - SMALL)
    {
        stepFraction_ += trackToFace(endPosition, td)*(1.0 - stepFraction_);
    }

    return facei_;
}


template<class TrackData>
Foam::scalar Foam::barycentricParticle::trackToFace
(
    const vector& endPosition,
    TrackData& td
)
{
    scalar trackFraction = 0.0;

    td.tetTriI = -1;

    facei_ = -1;

    // Tracking loop
    do
    {
        trackFraction += trackToTri(endPosition, td)*(1 - trackFraction);

        if (td.tetTriI == -1)
        {
            // The track is complete
            return trackFraction;
        }
        else if (td.tetTriI == 0)
        {
            // Set the face so that the loop will exit for face/patch processing
            facei_ = tetFacei_;
        }
        else
        {
            // Move to the next tet and continue the track
            changeTet(td);
        }

    } while(facei_ == -1);

    // Face/patch processing
    typedef typename TrackData::cloudType::particleType particleType;
    particleType& p = static_cast<particleType&>(*this);
    p.hitFace(td);
    if (internalFace(facei_))
    {
        changeCell(td);
    }
    else
    {
        label origFacei = facei_;
        label patchi = patch(facei_);

        // No action is taken for tetPti_ for tetFacei_ here. These are handled
        // by the patch interaction call or later during processor transfer.

        // <-- !!! The proper hitWallFaces method hasn't been implemented here
        // yet. This is a poor-man's place-holder, which assumes a small
        // particle that hits the wall within the current tet. To be honest, the
        // "proper" method is a bit suspect as it only considers tets in the
        // current cell. If the particle can hit the wall at a location within
        // another tet then it can certainly hit a location outside the cell...
        const tetIndices faceHitTetIs =
            polyMeshTetDecomposition::triangleTetIndices
            (
                mesh_,
                tetFacei_,
                celli_,
                tetPti_
            );

        if
        (
            !p.hitPatch
            (
                mesh_.boundaryMesh()[patchi],
                td,
                patchi,
                trackFraction,
                faceHitTetIs
            )
        )
        {
            // Did patch interaction model switch patches?
            if (facei_ != origFacei)
            {
                patchi = patch(facei_);
            }

            const polyPatch& patch = mesh_.boundaryMesh()[patchi];

            if (isA<wedgePolyPatch>(patch))
            {
                p.hitWedgePatch
                (
                    static_cast<const wedgePolyPatch&>(patch), td
                );
            }
            else if (isA<symmetryPlanePolyPatch>(patch))
            {
                p.hitSymmetryPlanePatch
                (
                    static_cast<const symmetryPlanePolyPatch&>(patch), td
                );
            }
            else if (isA<symmetryPolyPatch>(patch))
            {
                p.hitSymmetryPatch
                (
                    static_cast<const symmetryPolyPatch&>(patch), td
                );
            }
            else if (isA<cyclicPolyPatch>(patch))
            {
                p.hitCyclicPatch
                (
                    static_cast<const cyclicPolyPatch&>(patch), td
                );
            }
            else if (isA<cyclicAMIPolyPatch>(patch))
            {
                p.hitCyclicAMIPatch
                (
                    static_cast<const cyclicAMIPolyPatch&>(patch),
                    td,
                    endPosition - position()
                );
            }
            else if (isA<processorPolyPatch>(patch))
            {
                p.hitProcessorPatch
                (
                    static_cast<const processorPolyPatch&>(patch), td
                );
            }
            else if (isA<wallPolyPatch>(patch))
            {
                p.hitWallPatch
                (
                    static_cast<const wallPolyPatch&>(patch), td, faceHitTetIs
                );
            }
            else
            {
                p.hitPatch(patch, td);
            }
        }
    }

    return trackFraction;
}


template<class TrackData>
Foam::scalar Foam::barycentricParticle::trackToTri
(
    const vector& endPosition,
    TrackData& td
)
{
    if (mesh_.moving())
    {
        NotImplemented; // <-- !!!
    }

    const vector x0 = position();
    // !!! <-- We really need the input as a displacement, not an end position
    const vector x1 = endPosition - x0;
    const barycentric y0 = coordinates_;

    if (debug)
    {
        Info<< "Tracking from " << x0 << " to " << x0 + x1 << endl;
    }

    // Get the tet geometry
    vector centre;
    scalar detA;
    tensor T;
    tetReverseTransform(centre, detA, T);

    if (debug)
    {
        vector o, b, v1, v2;
        tetGeometry(o, b, v1, v2);
        Info<< "Tet points o=" << o << ", b=" << b
            << ", v1=" << v1 << ", v2=" << v2 << endl
            << "Tet determinant = " << detA << endl
            << "Start local coordinates = " << y0 << endl;
    }

    // Calculate the local tracking displacement
    const barycentric Tx1 = toBarycentric(T, x1);

    if (debug)
    {
        Info<< "Local displacement = " << Tx1 << "/" << detA << endl;
    }

    // Calculate the hit fraction
    label iH = -1;
    scalar muH = detA == 0 ? VGREAT : mag(1/detA);
    for (label i = 0; i < 4; ++ i)
    {
        if (Tx1[i] < 0)
        {
            scalar mu = - y0[i]/Tx1[i];

            if (debug)
            {
                Info<< "Hit on tet face " << i << " at local coordinate "
                    << y0 + mu*Tx1 << ", " << mu*detA*100 << "\% of the "
                    << "way along the track" << endl;
            }

            if (0 <= mu && mu < muH)
            {
                iH = i;
                muH = mu;
            }
        }
    }

    // Set the new coordinates
    barycentric yH = y0 + muH*Tx1;

    // Remove tolerance issues in the event of a hit
    if (iH != -1)
    {
        yH.replace(iH, 0);
    }

    // Set the new position and hit index
    coordinates_ = yH;
    td.tetTriI = iH;

    if (debug)
    {
        if (iH != -1)
        {
            Info<< "Track hit tet face " << iH << " first" << endl;
        }
        else
        {
            Info<< "Track hit no tet faces" << endl;
        }
        Info<< "End local coordinates = " << yH << endl
            << "End global coordinates = " << position() << endl
            << "Tracking displacement = " << position() - x0 << endl
            << muH*detA*100 << "\% of the track completed" << endl;
    }

    // Return the proportion of the track that has been completed
    return muH*detA;
}


template<class TrackData>
void Foam::barycentricParticle::changeTet(const TrackData& td)
{
    const bool isOwner = mesh_.faceOwner()[tetFacei_] == celli_;

    const label firstTetPtI = 1;
    const label lastTetPtI = mesh_.faces()[tetFacei_].size() - 2;

    if (td.tetTriI == 1)
    {
        changeFace(td);
    }
    else if (td.tetTriI == 2)
    {
        if (isOwner)
        {
            if (tetPti_ == lastTetPtI)
            {
                changeFace(td);
            }
            else
            {
                reflect();
                tetPti_ += 1;
            }
        }
        else
        {
            if (tetPti_ == firstTetPtI)
            {
                changeFace(td);
            }
            else
            {
                reflect();
                tetPti_ -= 1;
            }
        }
    }
    else if (td.tetTriI == 3)
    {
        if (isOwner)
        {
            if (tetPti_ == firstTetPtI)
            {
                changeFace(td);
            }
            else
            {
                reflect();
                tetPti_ -= 1;
            }
        }
        else
        {
            if (tetPti_ == lastTetPtI)
            {
                changeFace(td);
            }
            else
            {
                reflect();
                tetPti_ += 1;
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Changing tet without changing cell should only happen when the "
            << "track is on triangle 1, 2 or 3."
            << exit(FatalError);
    }
}


template<class TrackData>
void Foam::barycentricParticle::changeFace(const TrackData& td)
{
    // Get the tet topology
    label basei, vertex1i, vertex2i;
    tetMeshIndices(basei, vertex1i, vertex2i);

    // Get the shared edge and the pre-rotation
    edge sharedEdge;
    if (td.tetTriI == 1)
    {
        sharedEdge = edge(vertex1i, vertex2i);
    }
    else if (td.tetTriI == 2)
    {
        sharedEdge = edge(vertex2i, basei);
    }
    else if (td.tetTriI == 3)
    {
        sharedEdge = edge(basei, vertex1i);
    }
    else
    {
        FatalErrorInFunction
            << "Changing face without changing cell should only happen when the"
            << " track is on triangle 1, 2 or 3."
            << exit(FatalError);
    }

    // Find the face in the same cell that shares the edge, and the
    // corresponding tetrahedra point
    tetPti_ = -1;
    const label nFaces = mesh_.cells()[celli_].size();
    for (label cellFaceI = 0; tetPti_ == -1 && cellFaceI < nFaces; ++ cellFaceI)
    {
        const label newFaceI = mesh_.cells()[celli_][cellFaceI];

        // Exclude the current face
        if (tetFacei_ == newFaceI)
        {
            continue;
        }

        const class face& newFace = mesh_.faces()[newFaceI];

        // Loop over the edges, looking for the shared one
        label edgeComp = 0;
        label edgeI = 0;
        for (; edgeComp == 0 && edgeI < newFace.size(); ++ edgeI)
        {
            edgeComp = edge::compare(sharedEdge, newFace.faceEdge(edgeI));
        }

        // If the face does not contain the edge, then move on to the next face
        if (edgeComp == 0)
        {
            continue;
        }

        // Correct the edge index based on whether the face is owned or
        // neighbours the current cell, and whether the comparison was in order
        // or not
        const bool isOwner = mesh_.faceOwner()[newFaceI] == celli_;
        if (isOwner && edgeComp == 1)
        {
            edgeI = newFace.prevLabel(edgeI);
        }
        else if (!isOwner && edgeComp == -1)
        {
            edgeI = newFace.nextLabel(edgeI);
        }

        // Set the new face and tet point. The loop will now exit as the tet
        // point has a value.
        const label newBaseI = max(0, mesh_.tetBasePtIs()[newFaceI]);
        tetFacei_ = newFaceI;
        tetPti_ = (edgeI + newFace.size() - newBaseI) % newFace.size();
    }

    if (tetPti_ == -1)
    {
        FatalErrorInFunction
            << "The search for an edge-connected face and tet-point failed."
            << exit(FatalError);
    }

    // Pre-rotation puts the shared edge opposite the base of the tetrahedron
    if (sharedEdge.otherVertex(vertex1i) == -1)
    {
        rotate(false);
    }
    else if (sharedEdge.otherVertex(vertex2i) == -1)
    {
        rotate(true);
    }

    // Update the new tet topology
    tetMeshIndices(basei, vertex1i, vertex2i);

    // Reflect to account for the change of triangle orientation on the new face
    reflect();

    // Post rotation puts the shared edge back in the correct location
    if (sharedEdge.otherVertex(vertex1i) == -1)
    {
        rotate(true);
    }
    else if (sharedEdge.otherVertex(vertex2i) == -1)
    {
        rotate(false);
    }
}


template<class TrackData>
void Foam::barycentricParticle::changeCell(const TrackData& td)
{
    // Set the cell to be the one on the other side of the face
    const label ownerCellI = mesh_.faceOwner()[tetFacei_];
    const bool isOwner = celli_ == ownerCellI;
    celli_ = isOwner ? mesh_.faceNeighbour()[tetFacei_] : ownerCellI;

    // Reflect to account for the change of triangle orientation in the new cell
    reflect();
}


// ************************************************************************* //
