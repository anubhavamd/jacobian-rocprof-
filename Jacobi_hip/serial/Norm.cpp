//**************************************************************************
//* Copyright (c) 2019, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "Jacobi.hpp"


void NormKernel(const int N,
                const dfloat dx, const dfloat dy,
                const dfloat*__restrict__ U,
                      dfloat*__restrict__ norm2) {

  for (int id=0; id < N; id++) {
    *norm2 += U[id] * U[id] * dx * dy;
  }
}

dfloat Norm(grid_t& grid, mesh_t& mesh, dfloat *U) {

  dfloat h_norm=0;
  NormKernel(mesh.N, mesh.dx, mesh.dy, U, &h_norm);

  dfloat norm;
  MPI_Allreduce(&h_norm, &norm, 1, MPI_DFLOAT, MPI_SUM, grid.comm);

  return sqrt(norm)*mesh.invNtotal;
}
