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

Class
    Foam::cartesianParticle

Description
    Base Cartesian particle class

\*---------------------------------------------------------------------------*/

#ifndef cartesianParticle_H
#define cartesianParticle_H

#include "particleBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class cartesianParticle;

// Forward declaration of friend functions and operators

Ostream& operator<<
(
    Ostream&,
    const cartesianParticle&
);

/*---------------------------------------------------------------------------*\
                          Class Particle Declaration
\*---------------------------------------------------------------------------*/

class cartesianParticle
:
    public particleBase
{
protected:

    // Protected data

        //- Position of particle
        vector position_;


    // Private Member Functions

        //- Find the tet tri faces between position and tet centre
        void findTris
        (
            const vector& position,
            DynamicList<label>& faceList,
            const tetPointRef& tet,
            const FixedList<vector, 4>& tetAreas,
            const FixedList<label, 4>& tetPlaneBasePtIs,
            const scalar tol
        ) const;

        //- Find the lambda value for the line to-from across the
        //  given tri face, where p = from + lambda*(to - from)
        inline scalar tetLambda
        (
            const vector& from,
            const vector& to,
            const label triI,
            const vector& tetArea,
            const label tetPlaneBasePtI,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const scalar tol
        ) const;

        //- Find the lambda value for a moving tri face
        inline scalar movingTetLambda
        (
            const vector& from,
            const vector& to,
            const label triI,
            const vector& tetArea,
            const label tetPlaneBasePtI,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const scalar tol
        ) const;

        //- Modify the tet owner data by crossing triI
        inline void tetNeighbour(label triI);

        //- Cross the from the given face across the given edge of the
        //  given cell to find the resulting face and tetPti
        inline void crossEdgeConnectedFace
        (
            const label& celli,
            label& tetFacei,
            label& tetPti,
            const edge& e
        );

        //- Hit wall faces in the current cell if the
        //- wallImpactDistance is non-zero.  They may not be in
        //- Different tets to the current.
        template<class CloudType>
        inline void hitWallFaces
        (
            const CloudType& td,
            const vector& from,
            const vector& to,
            scalar& lambdaMin,
            tetIndices& closestTetIs
        );


public:

    // Static data members

        //- Runtime type information
        TypeName("cartesianParticle");

        //- String representation of properties
        DefinePropertyList
        (
            string("(Px Py Pz)") + particleBase::propertyList()
        );


    // Constructors

        //- Construct from components
        cartesianParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from components, tetFacei_ and tetPti_ are not
        //  supplied so they will be deduced by a search
        cartesianParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            bool doCellFacePt = true
        );

        //- Construct from Istream
        cartesianParticle
        (
            const polyMesh& mesh,
            Istream&,
            bool readFields = true
        );

        //- Construct as a copy
        cartesianParticle(const cartesianParticle& p);

        //- Construct as a copy with refernce to a new mesh
        cartesianParticle(const cartesianParticle& p, const polyMesh& mesh);


    //- Destructor
    virtual ~cartesianParticle()
    {}


    // Member Functions

        // Access

            //- Position as a noun. Return the cartesian position.
            inline const vector& position() const;

            //- Position as a verb. Set the current cartesian position.
            void position(const vector& x);


        // Check

            //- Check the stored cell value (setting if necessary) and
            //  initialise the tetFace and tetPt values
            inline void initCellFacePt();


        // Track

            //- Track particle to end of trajectory
            //  or until it hits the boundary.
            //  On entry 'stepFraction()' should be set to the fraction of the
            //  time-step at which the tracking starts and on exit it contains
            //  the fraction of the time-step completed.
            //  Returns the boundary face index if the track stops at the
            //  boundary, -1 otherwise.
            template<class TrackData>
            label track(const vector& endPosition, TrackData& td);

            //- Track particle to a given position and returns 1.0 if the
            //  trajectory is completed without hitting a face otherwise
            //  stops at the face and returns the fraction of the trajectory
            //  completed.
            //  on entry 'stepFraction()' should be set to the fraction of the
            //  time-step at which the tracking starts.
            template<class TrackData>
            scalar trackToFace(const vector& endPosition, TrackData& td);


        // I-O

            //- Write the particle position and cell
            void writePosition(Ostream&) const;


    // Friend Operators

        friend Ostream& operator<<(Ostream&, const cartesianParticle&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cartesianParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cartesianParticleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
