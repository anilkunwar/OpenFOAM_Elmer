//http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2010/highLevelProgramming.
//http://caefn.com/openfoam/solvers-
//http://www.foamacademy.com/wp-content/uploads/2016/11/thermo_training_handout_public.pdf

{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
      + fvOptions(T)
        // Added for Joule Heating Source
//+ sigma*(fvc::grad(Vel) ) & fvc::grad(Vel))  //Vel is electric voltage
     //   + fvm::Sp(qjoule,T) //Implicit source term case is as such
        + qjoule*alphaEff/kappa
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    radiation->correct();

    fvOptions.correct(T);

    rhok = 1.0 - beta*(T - TRef);
}
