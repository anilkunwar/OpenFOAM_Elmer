Info<< "Reading thermophysical properties\n" << endl;

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field qjoule\n" << endl;
volScalarField qjoule
(
    IOobject
    (
        "qjoule",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing field qjoule_recv\n" << endl;
volScalarField qjoule_recv
(
    IOobject
    (
        "qjoule_recv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    qjoule
);


Info<< "Initializing field alpha_old\n" << endl;
volScalarField alpha_old
(
    IOobject
    (
        "alpha_old",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
	dimensionedScalar("zero", dimless, 0)
);


Info<< "Initializing field alpha1f\n" << endl;
volScalarField alpha1f
(
    IOobject
    (
        "alpha1f",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
	dimensionedScalar("zero", dimless, 0)
);


//Info<< "Reading field kappa\n" << endl;
//volScalarField kappa
//(
    //IOobject
    //(
      //  "kappa",
      //  runTime.timeName(),
       // mesh,
       // IOobject::MUST_READ,
      //  IOobject::AUTO_WRITE
    //),
   // mesh
//);


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
//Add the field electric potential and electric conductivity for joule heatin
Info<< "Initializing field Vel\n" << endl;
volScalarField Vel
(
    IOobject
    (
        "Vel",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing field elcond\n" << endl;
volScalarField elcond
(
    IOobject
    (
        "elcond",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

#include "createPhi.H"

#include "readTransportProperties.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    1.0 - beta*(T - TRef)
);

// kinematic turbulent thermal thermal conductivity m2/s
Info<< "Reading field alphat\n" << endl;
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rhok*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p_rgh.name());

#include "createMRF.H"
#include "createIncompressibleRadiationModel.H"

//can define in readTransportProperties.H file as well
// http://openfoamwiki.net/index.php/How_to_add_temperature_to_icoFoam (A very important note for adding extra Physics)
// Reference conductivity
//dimensionedScalar elcond_ref("elcond_ref", dimless, mixture);
 dimensionedScalar elcond_ref("elcond_ref", dimless, laminarTransport);
//dimensionedScalar elcond_ref("elcond_ref", dimless, 0);
//dimensionedScalar elcond_ref
//(
//    transportProperties.lookup("elcond_ref")
//);

//dimensionedScalar kappa
//(
 //    transportProperties.lookup("kappa")
//);
dimensionedScalar kappa("kappa", dimless, laminarTransport); //please change later to dimKappa
//dimensionedScalar kappa("kappa", dimless, 0); //please change later to dimKappa
