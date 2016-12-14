/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::Elmer

Description
    Coupling with FEM open source multiphysical simulation software Elmer

SourceFiles
    Elmer.C

Authors: Juris Vencels
Email:   juris.vencels@gmail.com

Original Date: 12.12.2016.

\*---------------------------------------------------------------------------*/

#ifndef Elmer_H
#define Elmer_H

#include "fvCFD.H"
#include "PstreamGlobals.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class Elmer Declaration
\*---------------------------------------------------------------------------*/

class Elmer
{

public:

    // Public data

        //- Construct and initialize
        Elmer(const fvMesh& meshOF);

        //- Get interpolated scalar field at cell centres from Elmer
        void getScalar(volScalarField& field);

        //- Get interpolated vector field at cell centres from Elmer
        void getVector(volVectorField& field);

        //- Data structure for each Elmer MPI rank
        class ElmerProc_t
        {
            public:
		        int *foundCellsIndx;
		        MPI_Request *reqSend, *reqRecv;
		        int globalRank;
		        int nFoundCells;
				double *recvValues;
        };

private:

    // Private data

        //- Number of cells for this OpenFOAM MPI rank
        const fvMesh& mesh_;

        //- Number of cells for this OpenFOAM MPI rank
        int nCells;

        //- My rank's cell centres
        double *cellCentres_x, *cellCentres_y, *cellCentres_z;

        //- Total global MPI ranks
        int totGlobalRanks;

        //- My rank at MPI_COMM_WORLD
        int myGlobalRank;

        //- Total OpenFOAM ranks
        int totLocalRanks;

        //- My rank at MPI_OF_WORLD
        int myLocalRank;

        //- Total ELmer ranks
        int totElmerRanks;
    
        //- Elmer ranks start in MPI_COMM_WORLD
        int ElmerRanksStart;

        //- OpenFOAM ranks start in MPI_COMM_WORLD
        int OFRanksStart;

        //- Data structure for each Elmer MPI rank
        ElmerProc_t *ELp;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //