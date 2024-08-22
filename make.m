path0 = cd;

cd 'C:\Users\lynab\OneDrive\Bureau\UAE INTERNSHIP\gtw\lib';
mex cellss.cpp;
mex oness.cpp;
mex zeross.cpp;
cd(path0);

cd 'C:\Users\lynab\OneDrive\Bureau\UAE INTERNSHIP\gtw\src\ali\dtw';
mex dtwFord.cpp;
mex dtwBack.cpp;
mex dtwFordAsy.cpp;
cd(path0);

cd 'C:\Users\lynab\OneDrive\Bureau\UAE INTERNSHIP\gtw\src\ali\help';
mex rowBd.cpp;
cd(path0);

cd 'C:\Users\lynab\OneDrive\Bureau\UAE INTERNSHIP\gtw\src\ali\imw';
mex timewarp.cpp;
cd(path0);
