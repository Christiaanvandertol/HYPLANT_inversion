%pkg install "https://github.com/gnu-octave/octave-netcdf/releases/download/v1.0.19/netcdf-1.0.19.tar.gz"


FLOXdata = 'output';

FLEXdatafolder = [{'c:\Users\tol\projects\disc\l2pp\validation\Scenarios\L2PP\L2PP-TC-V1.0-A-L2PP-GlobalPerformancesE2ES-TDS6\Karpov\submoduleproducts\'}, ...
    {'c:\Users\tol\projects\disc\l2pp\validation\Scenarios\L2PP\L2PP-TC-V1.0-A-L2PP-GlobalPerformancesE2ES-TDS6\Karpov\submoduleproducts\'}];

 Var = matchup(FLEXdatafolder,FLOXdata);

