make clean

#!/bin/bash
# Before making, you need to copy all healix headers in the current folder.
# I don't understand why this should be needed, but it is needed.
cp /global/homes/e/eschaan/local/Healpix_3.31/include_f90/* ./

make
