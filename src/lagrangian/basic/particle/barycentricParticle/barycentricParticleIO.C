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
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::barycentricParticle::propertyList_ =
    Foam::barycentricParticle::propertyList();


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::barycentricParticle::barycentricParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particleBase(mesh, is, readFields),
    coordinates_()
{
    vector pos;
    if (readFields)
    {
        particleBase::readFields(is, pos);
    }
    else
    {
        particleBase::readPosition(is, pos);
    }
    position(pos);
}


void Foam::barycentricParticle::writePosition(Ostream& os) const
{
    particleBase::writePosition(os, position());
}


Foam::Ostream& Foam::operator<<(Ostream& os, const barycentricParticle& p)
{
    p.writeFields(os, p.position());

    return os;
}


// ************************************************************************* //
