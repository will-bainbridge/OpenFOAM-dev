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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::barycentric& Foam::barycentricParticle::coordinates() const
{
    return coordinates_;
}


inline Foam::barycentric& Foam::barycentricParticle::coordinates()
{
    return coordinates_;
}


inline Foam::vector Foam::barycentricParticle::position() const
{
    vector centre;
    tensor A;
    tetTransform(centre, A);

    return fromBarycentric(A, centre, coordinates());
}


inline void Foam::barycentricParticle::position(const vector& position)
{
    vector centre;
    scalar detA;
    tensor T;
    tetReverseTransform(centre, detA, T);

    coordinates() = toBarycentric(T/detA, centre, position);
}


inline void Foam::barycentricParticle::initCellFacePt()
{
    vector pos = position();

    particleBase::initCellFacePt(pos);

    position(pos);
}

inline bool Foam::barycentricParticle::hasCellFacePt() const
{
    // This is a bit of a mess. The particle class is very fundamentally built
    // around the concept of staged construction. After construction, the
    // tetrahedra may be uninitialised. Therefore we have to consider the stored
    // coordinates as cartesian up until the point at which the tet is fully
    // specified. This function provides that logic.
    return celli_ != -1 && tetFacei_ != -1 && tetPti_ != -1;
}


// ************************************************************************* //
