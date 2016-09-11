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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::particleBase::getNewParticleID() const
{
    label id = particleCount_++;

    if (id == labelMax)
    {
        WarningInFunction
            << "Particle counter has overflowed. This might cause problems"
            << " when reconstructing particle tracks." << endl;
    }
    return id;
}


inline const Foam::polyMesh& Foam::particleBase::mesh() const
{
    return mesh_;
}


inline Foam::label Foam::particleBase::cell() const
{
    return celli_;
}


inline Foam::label& Foam::particleBase::cell()
{
    return celli_;
}


inline Foam::label Foam::particleBase::tetFace() const
{
    return tetFacei_;
}


inline Foam::label& Foam::particleBase::tetFace()
{
    return tetFacei_;
}


inline Foam::label Foam::particleBase::tetPt() const
{
    return tetPti_;
}


inline Foam::label& Foam::particleBase::tetPt()
{
    return tetPti_;
}


inline Foam::tetIndices Foam::particleBase::currentTetIndices() const
{
    return tetIndices(celli_, tetFacei_, tetPti_, mesh_);
}


inline Foam::tetPointRef Foam::particleBase::currentTet() const
{
    return currentTetIndices().tet(mesh_);
}


inline Foam::vector Foam::particleBase::normal() const
{
    return currentTetIndices().faceTri(mesh_).normal();
}


inline Foam::vector Foam::particleBase::oldNormal() const
{
    return currentTetIndices().oldFaceTri(mesh_).normal();
}


inline Foam::label Foam::particleBase::face() const
{
    return facei_;
}


inline Foam::label& Foam::particleBase::face()
{
    return facei_;
}


inline void Foam::particleBase::initCellFacePt(Foam::vector& position)
{
    if (celli_ == -1)
    {
        mesh_.findCellFacePt
        (
            position,
            celli_,
            tetFacei_,
            tetPti_
        );

        if (debug && celli_ == -1)
        {
            WarningInFunction
                << "cell, tetFace and tetPt search failure for position "
                << position
                << endl;
        }
    }
    else
    {
        mesh_.findTetFacePt(celli_, position, tetFacei_, tetPti_);

        if (tetFacei_ == -1 || tetPti_ == -1)
        {
            label oldCelli = celli_;

            mesh_.findCellFacePt
            (
                position,
                celli_,
                tetFacei_,
                tetPti_
            );

            if (celli_ == -1 || tetFacei_ == -1 || tetPti_ == -1)
            {
                // The particle has entered this function with a cell number,
                // but hasn't been able to find a cell to occupy.

                if (!mesh_.pointInCellBB(position, oldCelli, 0.1))
                {
                    // If the position is not inside the (slightly extended)
                    // bound-box of the cell that it thought it should be in,
                    // then this is considered an error.

                    if (debug)
                    {
                        WarningInFunction
                            << "position " << position
                            << " not in specified cell " << oldCelli
                            << " with centre " << mesh_.cellCentres()[oldCelli]
                            << endl;
                    }

                    celli_ = -1;
                    tetFacei_ = -1;
                    tetPti_ = -1;

                    return;
                }

                // The position is in the (slightly extended) bound-box of the
                // cell.  This situation may arise because the face
                // decomposition of the cell is not the same as when the
                // particle acquired the cell index.  For example, it has been
                // read into a mesh that has made a different face base-point
                // decision for a boundary face and now this particle is in a
                // position that is not in the mesh.  Gradually move the
                // particle towards the centre of the cell that it thought that
                // it was in.

                celli_ = oldCelli;

                point newPosition = position;

                const point& cC = mesh_.cellCentres()[celli_];

                label trap(1.0/trackingCorrectionTol + 1);

                label iterNo = 0;

                do
                {
                    newPosition += trackingCorrectionTol*(cC - position);

                    mesh_.findTetFacePt
                    (
                        celli_,
                        newPosition,
                        tetFacei_,
                        tetPti_
                    );

                    iterNo++;

                } while (tetFacei_ < 0  && iterNo <= trap);

                if (tetFacei_ == -1)
                {
                    FatalErrorInFunction
                        << "cell, tetFace and tetPt search failure at position "
                        << position << abort(FatalError);
                }

                if (debug)
                {
                    WarningInFunction
                        << "Particle moved from " << position
                        << " to " << newPosition
                        << " in cell " << celli_
                        << " tetFace " << tetFacei_
                        << " tetPt " << tetPti_ << nl
                        << "    (A fraction of "
                        << 1.0 - mag(cC - newPosition)/mag(cC - position)
                        << " of the distance to the cell centre)"
                        << " because a decomposition tetFace and tetPt "
                        << "could not be found."
                        << endl;
                }

                position = newPosition;
            }

            if (debug && celli_ != oldCelli)
            {
                WarningInFunction
                    << "Particle at position " << position
                    << " searched for a cell, tetFace and tetPt." << nl
                    << "    Found"
                    << " cell " << celli_
                    << " tetFace " << tetFacei_
                    << " tetPt " << tetPti_ << nl
                    << "    This is a different cell to that which was supplied"
                    << " (" << oldCelli << ")." << nl
                    << endl;
            }
        }
    }
}


inline bool Foam::particleBase::onBoundary() const
{
    return facei_ != -1 && facei_ >= mesh_.nInternalFaces();
}


inline Foam::scalar& Foam::particleBase::stepFraction()
{
    return stepFraction_;
}


inline Foam::scalar Foam::particleBase::stepFraction() const
{
    return stepFraction_;
}


inline Foam::label Foam::particleBase::origProc() const
{
    return origProc_;
}


inline Foam::label& Foam::particleBase::origProc()
{
    return origProc_;
}


inline Foam::label Foam::particleBase::origId() const
{
    return origId_;
}


inline Foam::label& Foam::particleBase::origId()
{
    return origId_;
}


inline bool Foam::particleBase::softImpact() const
{
    return false;
}


inline Foam::scalar Foam::particleBase::currentTime() const
{
    return
        mesh_.time().value()
      + stepFraction_*mesh_.time().deltaTValue();
}


inline bool Foam::particleBase::internalFace(const label facei) const
{
    return mesh_.isInternalFace(facei);
}


bool Foam::particleBase::boundaryFace(const label facei) const
{
    return !internalFace(facei);
}


inline Foam::label Foam::particleBase::patch(const label facei) const
{
    return mesh_.boundaryMesh().whichPatch(facei);
}


inline Foam::label Foam::particleBase::patchFace
(
    const label patchi,
    const label facei
) const
{
    return mesh_.boundaryMesh()[patchi].whichFace(facei);
}


inline Foam::label Foam::particleBase::faceInterpolation() const
{
    return facei_;
}


// ************************************************************************* //
