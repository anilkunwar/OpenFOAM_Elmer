// Need to check more for the exact syntax
{
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            solve
            (
                fvm::laplacian(sigma, Vel)
            );
        }
    Info<< "Min/max V:" << min(Vel) << ' '
        << max(Vel) << endl;

}
