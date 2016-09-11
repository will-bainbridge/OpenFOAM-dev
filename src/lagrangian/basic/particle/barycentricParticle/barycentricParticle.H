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
    Foam::barycentricParticle

Description
    Base Barycentric particle class

\*---------------------------------------------------------------------------*/

#ifndef barycentricParticle_H
#define barycentricParticle_H

#include "particleBase.H"
#include "barycentric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class barycentricParticle;

// Forward declaration of friend functions and operators

Ostream& operator<<
(
    Ostream&,
    const barycentricParticle&
);

/*---------------------------------------------------------------------------*\
                          Class Particle Declaration
\*---------------------------------------------------------------------------*/

class barycentricParticle
:
    public particleBase
{
public:

    template<class CloudType>
    class TrackingData
    :
        public particleBase::TrackingData<CloudType>
    {
    public:

        // Public data

            //  the track has stopped inside the tet.
            label tetTriI;


        // Constructors

            // Inherit base constructors
            using particleBase::TrackingData<CloudType>::TrackingData;
    };


protected:

    // Protected data

        //- Barycentric coordinates of the particle within the current tet
        barycentric coordinates_;


private:

    // Private Member Functions

        // Tetrahedra functions

            //- Get indices into the current face for the face-bound vertices of
            //  the current tet.
            void tetFaceIndices
            (
                label& baseI,
                label& vertex1I,
                label& vertex2I
            ) const;

            //- Get indices into the mesh points for the face-bound vertices of
            //  the current tet.
            void tetMeshIndices
            (
                label& baseI,
                label& vertex1I,
                label& vertex2I
            ) const;

            //- Get the vertices of the current tet
            void tetGeometry
            (
                vector& centre,
                vector& base,
                vector& vertex1,
                vector& vertex2
            ) const;

            //- Get the transformation associated with the current tet. This
            //  will convert a barycentric position within the tet to a
            //  cartesian position in the global coordinate system. The
            //  conversion is x = centre + A&y, where x is the cartesian
            //  position, y is the barycentric position, centre is the origin,
            //  and A is the transformation tensor.
            void tetTransform
            (
                vector& centre,
                tensor& A
            ) const;

            //- Get the reverse transform associated with the current tet. The
            //  conversion is detA*y = T&(x - centre). The variables x, y and
            //  centre have the same meaning as for the forward transform. T is
            //  the inverse of the forward transform tensor, A, multiplied by
            //  its determinant, detA. This separation allows the barycentric
            //  tracking algorithm to function on inverted or degenerate
            //  tetrahedra.
            void tetReverseTransform
            (
                vector& centre,
                scalar& detA,
                tensor& T
            ) const;

            //- Get the vertices of the current moving tet. Two values are
            //  returned for each vertex. The first is a constant, and the
            //  second is a linear coefficient of the track fraction.
            void movingTetGeometry
            (
                Pair<vector>& centre,
                Pair<vector>& base,
                Pair<vector>& vertex1,
                Pair<vector>& vertex2
            ) const;

            //- Get the transformation associated with the current, moving, tet.
            //  This is of the same form as for the static case. As with the
            //  moving geometry, a linear function of the tracking fraction is
            //  returned for each component.
            void movingTetTransform
            (
                Pair<vector>& centre,
                Pair<tensor>& A
            ) const;

            //- Get the reverse transformation associated with the current,
            //  moving, tet. This is of the same form as for the static case. As
            //  with the moving geometry, a function of the tracking fraction is
            //  returned for each component. The functions are higher order than
            //  for the forward transform; the determinant is cubic, and the
            //  tensor is quadratic.
            void movingTetReverseTransform
            (
                Pair<vector>& centre,
                FixedList<scalar, 4>& detA,
                FixedList<tensor, 3>& T
            ) const;


        // Transformations

            //- Reflection transform. Corrects the coordinates when the particle
            //  moves between two tets which share a base vertex, but for which
            //  the other two non cell-centre vertices are reversed. All hits
            //  which retain the same face behave this way, as do face hits.
            void reflect();

            //- Rotation transform. Corrects the coordinates when the particle
            //  moves between two tets with different base vertices, but are
            //  otherwise similarly oriented. Hits which change the face within
            //  the cell make use of both this and the reflect transform.
            void rotate(const bool direction);


        // Track

            //- Change tet within a cell. Called after a triangle is hit.
            template<class TrackData>
            void changeTet(const TrackData& td);

            //- Change tet face within a cell. Called by changeTet.
            template<class TrackData>
            void changeFace(const TrackData& td);

            //- Change cell. Called when the particle hits an internal face.
            template<class TrackData>
            void changeCell(const TrackData& td);


public:

    // Static data members

        //- Runtime type information
        TypeName("barycentricParticle");

        //- String representation of properties
        DefinePropertyList
        (
            string("(Ca Cb Cc Cd)") + particleBase::propertyList()
        );


    // Constructors

        //- Construct from components
        barycentricParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from components, tetFacei_ and tetPti_ are not
        //  supplied so they will be deduced by a search
        barycentricParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            bool doCellFacePt = true
        );

        //- Construct from Istream
        barycentricParticle
        (
            const polyMesh& mesh,
            Istream&,
            bool readFields = true
        );

        //- Construct as a copy
        barycentricParticle(const barycentricParticle& p);

        //- Construct as a copy with refernce to a new mesh
        barycentricParticle(const barycentricParticle& p, const polyMesh& mesh);


    //- Destructor
    virtual ~barycentricParticle()
    {}


    // Member Functions

        // Access

            //- Return current particle barycentric coordinates
            inline const barycentric& coordinates() const;

            //- Return current particle barycentric coordinates
            inline barycentric& coordinates();

            //- Position as a noun. Calculate and return the cartesian position
            //  from the current barycentric coordinates.
            vector position() const;

            //- Position as a verb. Set the current barycentric coordinates from
            //  the given cartesian position.
            void position(const vector& x);


        // Check

            //- Check the stored cell value (setting if necessary) and
            //  initialise the tetFace and tetPt values
            inline void initCellFacePt();

            //- Have the cell, face, tet points been initialised?
            inline bool hasCellFacePt() const;


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

            //- Track to the next intersection with a triangle of the current
            //  tet, or to the end position, whichever is the nearer.
            template<class TrackData>
            scalar trackToTri(const vector& endPosition, TrackData& td);


    // I-O

        //- Write the particle position and cell
        void writePosition(Ostream&) const;


    // Friend Operators

        friend Ostream& operator<<(Ostream&, const barycentricParticle&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "barycentricParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "barycentricParticleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
