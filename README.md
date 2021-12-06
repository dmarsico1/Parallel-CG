# Parallel-CG
Parallel conjugate gradient solvers in two and three dimensions.

In two the dimensions, mpi_cart_poisson solves Poisson's equation in parallel on a square domain using cartesian topologies where the domain is decomposed into rectangles.

In three dimensions, mpi_cart_poisson_3d solves the same equation in parallel on a cube, again using cartesian topologies.  Here, howeverm, the domain is decomposed into rectangular prisms.

To compile the code:

mpifort -o executable mpi_cart_poisson_3d.f95 -I/(netCDF include directory) -L/(path to netCDF lib directory) -lnetcdf -lnetcdff
