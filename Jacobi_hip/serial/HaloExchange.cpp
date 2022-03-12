//**************************************************************************
//* Copyright (c) 2019, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "Jacobi.hpp"

void HaloExchangeStart(grid_t& grid, mesh_t& mesh, dfloat* h_U) {

  //copy left and right sides to a contiguous send buffer
  if (grid.Neighbor[SIDE_LEFT]>-1) {
    dfloat *sbuf = h_U;
    dfloat *dbuf = mesh.sendBuffer;

    for (int j=0;j<mesh.Ny;j++) {
      dbuf[j] = sbuf[j*mesh.Nx];
    }
  }

  if (grid.Neighbor[SIDE_RIGHT]>-1) {
    dfloat *sbuf = h_U+mesh.Nx-1;
    dfloat *dbuf = mesh.sendBuffer+mesh.Ny;

    for (int j=0;j<mesh.Ny;j++) {
      dbuf[j] = sbuf[j*mesh.Nx];
    }
  }

  //post recvs & sends
  if (grid.Neighbor[SIDE_DOWN]>-1) {
    MPI_Irecv(mesh.haloBuffer, mesh.Nx, MPI_DFLOAT,
              grid.Neighbor[SIDE_DOWN], 0, grid.comm, mesh.requests+0);
    MPI_Isend(h_U, mesh.Nx, MPI_DFLOAT,
              grid.Neighbor[SIDE_DOWN], 0, grid.comm, mesh.requests+1);
  }

  if (grid.Neighbor[SIDE_UP]>-1) {
    MPI_Irecv(mesh.haloBuffer+mesh.Nx, mesh.Nx, MPI_DFLOAT,
              grid.Neighbor[SIDE_UP], 0, grid.comm, mesh.requests+2);
    MPI_Isend(h_U+(mesh.Ny-1)*mesh.Nx, mesh.Nx, MPI_DFLOAT,
              grid.Neighbor[SIDE_UP], 0, grid.comm, mesh.requests+3);
  }

  if (grid.Neighbor[SIDE_LEFT]>-1) {
    MPI_Irecv(mesh.haloBuffer+2*mesh.Nx, mesh.Ny, MPI_DFLOAT,
              grid.Neighbor[SIDE_LEFT], 0, grid.comm, mesh.requests+4);
    MPI_Isend(mesh.sendBuffer, mesh.Ny, MPI_DFLOAT,
              grid.Neighbor[SIDE_LEFT], 0, grid.comm, mesh.requests+5);
  }

  if (grid.Neighbor[SIDE_RIGHT]>-1) {
    MPI_Irecv(mesh.haloBuffer+2*mesh.Nx+mesh.Ny, mesh.Ny, MPI_DFLOAT,
              grid.Neighbor[SIDE_RIGHT], 0, grid.comm, mesh.requests+6);
    MPI_Isend(mesh.sendBuffer+mesh.Ny, mesh.Ny, MPI_DFLOAT,
              grid.Neighbor[SIDE_RIGHT], 0, grid.comm, mesh.requests+7);
  }
}

void HaloExchangeFinish(grid_t& grid, mesh_t& mesh, dfloat* h_U) {
  // Wait for all sent messages to have left and received messages to have arrived
  MPI_Waitall(8, mesh.requests, mesh.status);
}
