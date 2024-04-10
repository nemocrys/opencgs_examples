/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 AUTHOR,AFFILIATION
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

Application
    azimuthalAverage

Description

Usage
   
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fvCFD.H"
#include "regionProperties.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    argList::addNote
    (
        "Note"
    );

    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << endl << "Reading azimuthalAverageDict" << endl;

    IOdictionary azimuthalAverageDict
    (
        IOobject
        (
            "azimuthalAverageDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    // 3D velocity field for azimuthal averaging; default is "UMean"
    const word U3DName_ = azimuthalAverageDict.getOrDefault<word>("U3DName", "UMean");
    // these variables depend on case setup in runAzimuthalAverage
    const word U2DName_ = azimuthalAverageDict.getOrDefault<word>("U2DName", "U");
    const word U3Dtime_ = azimuthalAverageDict.getOrDefault<word>("U3Dtime", "0");
    const word case3DName_ = azimuthalAverageDict.getOrDefault<word>("case3DName", "case3D");



    // read 2D and 3D meshes
    regionProperties rp(runTime);
    // const wordList mesh2Dnames(rp["case2D"]);
    const wordList mesh3Dnames(rp[case3DName_]);

    // Foam::autoPtr<Foam::fvMesh> mesh2DPtr(nullptr);
    Foam::autoPtr<Foam::fvMesh> mesh3DPtr(nullptr);

    Foam::fvMesh& mesh2D = mesh; // mesh2D is internalMesh in base case


    // mesh2DPtr.reset
    // (
    //     new Foam::fvMesh
    //     (
    //         Foam::IOobject
    //         (
    //             mesh2Dnames[0],
    //             runTime.timeName(),
    //             runTime,
    //             Foam::IOobject::MUST_READ
    //         ),
    //         false
    //     )
    // );

    // Foam::fvMesh& mesh2D = mesh2DPtr();


    Info << "2D mesh with elements: " << mesh2D.C().size() << endl;
    Info << "List of all available patches in 2D mesh:" << endl
         << (mesh2D.boundary().size()) << endl
        << "(" << endl;
    forAll (mesh2D.boundary(),patchI)
    {
        Info << mesh2D.boundary()[patchI].name() << endl;
    }    
    Info << ")" << endl;


    mesh3DPtr.reset
    (
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                mesh3Dnames[0],
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            ),
            false
        )
    );
    Foam::fvMesh& mesh3D = mesh3DPtr();


    Info << "3D mesh with elements: " << mesh3D.C().size() << endl;
    Info << "List of all available patches in 3D mesh:" << endl
         << (mesh3D.boundary().size()) << endl
         << "(" << endl;
    forAll (mesh3D.boundary(),patchI)
    {
        Info << mesh3D.boundary()[patchI].name() << endl;
    }    
    Info << ")" << endl;




    // The base fields required
    Info<< "Reading field U3D\n" << endl;
    volVectorField U3D
    (
        IOobject
        (
            U3DName_,
            U3Dtime_,
            mesh3D,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh3D
    );


    Info<< "Reading field U2D\n" << endl;
    volVectorField U2D
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh2D,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh2D
    );


    Info << "Starting interpolation... " << endl;

    label nNotFound = 0;

    // loop over all 2D cells
    forAll(U2D, cell2D)
    {

        Info << "2D Cell =  " << cell2D;

        vector p2D = mesh2D.C()[cell2D];
        label nCells = 0;

        U2D[cell2D] = vector(0, 0, 0);


        // loop over all 3D cells to find minimum distance

        scalar mindist = 1000;

        forAll(U3D, cell3D) // loop over all 3D cells to find minimum distance
        {
            vector p3D = mesh3D.C()[cell3D];

            // convert 3D mesh coordinates to cylindrical
            vector p3Dproj;
            p3Dproj[0] = Foam::sqrt( p3D[0]*p3D[0] + p3D[2]*p3D[2] ); // y is rotational axis
            p3Dproj[1] = p3D[1];
            p3Dproj[2] = 0;

            scalar dist = mag(p2D-p3Dproj);
            if (dist<mindist) mindist = dist;
        }


        // loop over all 3D cells for interpolation

        forAll(U3D, cell3D) 
        {
            vector p3D = mesh3D.C()[cell3D];
            scalar d3D = mesh3D.V()[cell3D];
            d3D = Foam::pow(d3D, 1.0/3.0); // approximate 3D cell size


            // convert 3D mesh coordinates to cylindrical
            vector p3Dproj;
            p3Dproj[0] = Foam::sqrt( p3D[0]*p3D[0] + p3D[2]*p3D[2] ); // y is rotational axis
            p3Dproj[1] = p3D[1];
            p3Dproj[2] = 0;


            // simplified test if 2D cell center is inside 3D cell
            scalar dist = mag(p2D-p3Dproj);
            if (dist < mindist + 0.5*d3D) 
                {

                // U projection to the ry plane
                if ( p3Dproj[0]<1e-12 ) continue;
                vector U3Dproj;
                U3Dproj[0] = ( U3D[cell3D][0]*p3D[0] + U3D[cell3D][2]*p3D[2] ) / p3Dproj[0]; // U_r = vec(U_xz)*vec(r_xz)/mag(r_xz)
                U3Dproj[1] = U3D[cell3D][1]; 
                U3Dproj[2] = 0;

                nCells += 1;
                U2D[cell2D] = U2D[cell2D] + U3Dproj;

                }
        }

        if (nCells>1) 
            {
            U2D[cell2D] /= nCells;
            }
        else nNotFound++;

        Info << " -> 3D cells = " << nCells << " U2D = " << U2D[cell2D] << endl;
    }

    // Set boundary fields to zero

    forAll(U2D.boundaryFieldRef(), patchesi)
        {
        forAll(U2D.boundaryFieldRef()[patchesi], facesi)
            {
            U2D.boundaryFieldRef()[patchesi][facesi] = vector(0, 0, 0); 
            }
        }

    Info << endl << nNotFound << " cells were not found on 3D grid" << endl;

    runTime.write();
    U2D.write();

    Info << endl << "End\n" << endl;
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
