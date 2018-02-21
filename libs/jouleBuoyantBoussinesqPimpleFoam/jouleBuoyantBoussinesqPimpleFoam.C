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

Application
    jouleBuoyantBoussinesqPimpleFoam = buoyantBoussinesqPimpleFoam + TEqn.H + Elmer
    \f[

Description
    Transient solver for buoyant, turbulent flow of incompressible fluids.

    Uses the Boussinesq approximation:
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]
    
    -------------------------------------------------------------------------------
Ogirinal buoyantBoussinesqPimpleFoam solver is part of OpenFOAM
//Thermal solver taken from Qingming Liu:
//http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2011/QingmingLiu/Project_QingmingLIU-final.pdf
Write the actual source
Modified by: Anil Kunwar
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "Elmer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    // Send fields to Elmer
    Elmer sending(mesh,1); // 1=send, -1=receive
    sending.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    //elcond = alpha1 * elcond_ref;
    elcond = 1.0 * elcond_ref;
    sending.sendScalar(elcond); //send elcond to Elmer

      // Receive fields from Elmer
    Elmer receiving(mesh,-1); // 1=send, -1=receive
    receiving.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    //receiving.recvVector(JxB_recv);
  receiving.recvScalar(qjoule_recv);

    while (runTime.run())
    {
        qjoule = qjoule_recv*1.0;

        #include "readTimeControls.H"

        if (LTS)
        {
            //#include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            //#include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;
   // while (runTime.run())
   // {
        //#include "readTimeControls.H"
        //#include "CourantNo.H"
        //#include "setDeltaT.H"

       // runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        //Check whether we need to update electromagnetic stuff with Elmer
        double maxRelDiff = (max(mag(alpha_old - alpha1f))).value();

        bool doElmer = false;
        if(maxRelDiff>0.5) {
            doElmer = true;
        }

        if(doElmer || !runTime.run()) {
            //alpha_old = alpha1f;
            double commTime = MPI_Wtime();

            // Send fields to Elmer
            sending.sendStatus(runTime.run());
            //elcond = alpha1f * elcond_ref;
           elcond = 1.0 * elcond_ref;
            sending.sendScalar(elcond);

            Info<< "OpenFOAM2Elmer = " << MPI_Wtime()-commTime << " s" << nl << endl;
            commTime = MPI_Wtime();

            // Receive fields form Elmer
            receiving.sendStatus(runTime.run());
            //receiving.recvVector(JxB_recv);
            receiving.recvScalar(qjoule_recv);

            Info<< "Elmer2OpenFOAM = " << MPI_Wtime()-commTime << " s" << nl << endl;
        }
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
