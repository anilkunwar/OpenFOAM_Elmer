//http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2010/highLevelProgramming.
//http://caefn.com/openfoam/solvers-
//http://www.foamacademy.com/wp-content/uploads/2016/11/thermo_training_handout_public.pdf
//https://www.cfd-online.com/Forums/openfoam-pre-processing/179786-buoyantboussinesqpimplefoam-heat-transfer.html
// Try with simulationType laminar; for convergence reasons and things.

//{
//    alphat = turbulence->nut()/Prt;
//    alphat.correctBoundaryConditions();

//    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

//   fvScalarMatrix TEqn
//    (
//        fvm::ddt(T)
//      + fvm::div(phi, T)
//     (- fvm::laplacian(alphaEff, T)
     //==
       // radiation->ST(rhoCpRef, T)
      //+ fvOptions(T)
        // Added for Joule Heating Source
//+ sigma*(fvc::grad(Vel) ) & fvc::grad(Vel))  //Vel is electric voltage
     //   + fvm::Sp(qjoule,T) //Implicit source term case is as such
//       + ((qjoule*alphaEff)/kappa))
//    );

//    TEqn.relax();

    //fvOptions.constrain(TEqn);

 //   TEqn.solve();

    //radiation->correct();

    //fvOptions.correct(T);

 //   rhok = 1.0 - beta*(T - TRef);
//}
{
    //volScalarField kappaf = kappa/cp ; // we need to define a volScalarField at the RHS of equation as both kappa and cp are dimensioned scalars

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - (
            fvm::laplacian(kappa/cp, T)
          + qjoule*(rho/cp)
        )
    );

    TEqn.relax();
    TEqn.solve();
    rhok = 1.0 - beta*(T - TRef);
}
