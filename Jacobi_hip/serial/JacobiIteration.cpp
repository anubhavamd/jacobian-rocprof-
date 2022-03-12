//**************************************************************************
//* Copyright (c) 2019, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "Jacobi.hpp"

//Jacobi iterative method
// U = U + D^{-1}*(RHS - AU)
void JacobiIterationKernel(const int N,
                           const dfloat dx,
                           const dfloat dy,
                           const dfloat *__restrict__ RHS,
                           const dfloat *__restrict__ AU,
                                 dfloat *__restrict__ RES,
                                 dfloat *__restrict__ U) {

  for (int id=0;id<N;id++) {
    const dfloat r_res = RHS[id] - AU[id];

    RES[id] = r_res;

    U[id] += r_res/(2.0/(dx*dx) + 2.0/(dy*dy));
  }
}

void JacobiIteration(grid_t& grid, mesh_t& mesh,
                     dfloat* h_RHS,
                     dfloat* h_AU,
                     dfloat* h_RES,
                     dfloat* h_U) {

  JacobiIterationKernel(mesh.N,
                     mesh.dx, mesh.dy,
                     h_RHS, h_AU, h_RES, h_U);
}
