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
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::particleBase::propertyList_ =
    Foam::particleBase::propertyList();

const std::size_t Foam::particleBase::sizeofFields_
(
    sizeof(particleBase) - offsetof(particleBase, celli_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleBase::particleBase
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    mesh_(mesh),
    celli_(-1),
    facei_(-1),
    stepFraction_(0.0),
    tetFacei_(-1),
    tetPti_(-1),
    origProc_(Pstream::myProcNo()),
    origId_(-1)
{}


void Foam::particleBase::writePosition
(
    Ostream& os,
    const Foam::vector& position
) const
{
    if (os.format() == IOstream::ASCII)
    {
        os  << position << token::SPACE << celli_;
    }
    else
    {
        os.write(reinterpret_cast<const char*>(&position), sizeof(vector));
        os.write(reinterpret_cast<const char*>(&celli_), sizeof(label));
    }

    // Check state of Ostream
    os.check("particleBase::writePosition(Ostream&, const vector&) const");
}


void Foam::particleBase::writeFields
(
    Ostream& os,
    const Foam::vector& position
) const
{
    if (os.format() == IOstream::ASCII)
    {
        os  << position
            << token::SPACE << celli_
            << token::SPACE << facei_
            << token::SPACE << stepFraction_
            << token::SPACE << tetFacei_
            << token::SPACE << tetPti_
            << token::SPACE << origProc_
            << token::SPACE << origId_;
    }
    else
    {
        os.write(reinterpret_cast<const char*>(&position), sizeof(vector));
        os.write(reinterpret_cast<const char*>(&celli_), sizeofFields_);
    }

    // Check state of Ostream
    os.check("particleBase::writeFields(Ostream&, const vector&) const");
}


void Foam::particleBase::readPosition
(
    Istream& is,
    Foam::vector& position
)
{
    if (is.format() == IOstream::ASCII)
    {
        is  >> position >> celli_;
    }
    else
    {
        is.read(reinterpret_cast<char*>(&position), sizeof(vector));
        is.read(reinterpret_cast<char*>(&celli_), sizeof(label));
    }

    // Check state of Istream
    is.check("particleBase::readPosition(Ostream&, vector& position) const");
}


void Foam::particleBase::readFields
(
    Istream& is,
    Foam::vector& position
)
{
    if (is.format() == IOstream::ASCII)
    {
        is  >> position
            >> celli_
            >> facei_
            >> stepFraction_
            >> tetFacei_
            >> tetPti_
            >> origProc_
            >> origId_;
    }
    else
    {
        is.read(reinterpret_cast<char*>(&position), sizeof(vector));
        is.read(reinterpret_cast<char*>(&celli_), sizeofFields_);
    }

    // Check state of Istream
    is.check("particleBase::readFields(Ostream&, vector& position) const");
}


// ************************************************************************* //
