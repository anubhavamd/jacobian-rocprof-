//**************************************************************************
//* Copyright (c) 2019, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "Jacobi.hpp"

// AU_i,j = (-U_i+1,j + 2U_i,j - U_i-1,j)/dx^2 +
//          (-U_i,j+1 + 2U_i,j - U_i,j-1)/dy^2
void LocalLaplacianKernel(const int localNx,
                          const int localNy,
                          const int stride,
                          const dfloat dx,
                          const dfloat dy,
                          const dfloat *__restrict__ U,
                                dfloat *__restrict__ AU) {

  for (int j=0;j<localNy;j++) {
    for (int i=0;i<localNx;i++) {

      const int id = (i+1) + (j+1)*stride;

      const int id_l = id - 1;
      const int id_r = id + 1;
      const int id_d = id - stride;
      const int id_u = id + stride;

      AU[id] = (-U[id_l] + 2*U[id] - U[id_r])/(dx*dx) +
               (-U[id_d] + 2*U[id] - U[id_u])/(dy*dy);
    }
  }
}

void LocalLaplacian(grid_t& grid, mesh_t& mesh,
                    dfloat* h_U,
                    dfloat* h_AU) {

  int localNx = mesh.Nx-2;
  int localNy = mesh.Ny-2;

  LocalLaplacianKernel(localNx, localNy, mesh.Nx,
                     mesh.dx, mesh.dy,
                     h_U, h_AU);
}


// AU_i,j = (-U_i+1,j + 2U_i,j - U_i-1,j)/dx^2 +
//          (-U_i,j+1 + 2U_i,j - U_i,j-1)/dy^2
void HaloLaplacianKernel(const int Nx,
                         const int Ny,
                         const int stride,
                         const dfloat dx,
                         const dfloat dy,
                         const dfloat *__restrict__ U,
                         const dfloat *__restrict__ haloBuffer,
                               dfloat *__restrict__ AU) {

  for (int id=0;id<2*Nx+2*Ny-2;id++) {

    //get the (i,j) coordinates of node based on how large id is
    int i, j;
    if (id < Nx) { //bottom
      i = id;
      j = 0;
    } else if (id<2*Nx) { //top
      i = id - Nx;
      j = Ny-1;
    } else if (id < 2*Nx + Ny-1) { //left
      i = 0;
      j = id - 2*Nx + 1;
    } else { //right
      i = Nx-1;
      j = id - 2*Nx - Ny + 2;
    }

    const int iid = i+j*stride;

    //check if node is one the boundary and use haloBuffer's data if so
    const dfloat U_d = (j==0)    ?  haloBuffer[i+0]       : U[iid - stride];
    const dfloat U_u = (j==Ny-1) ?  haloBuffer[i+Nx]      : U[iid + stride];
    const dfloat U_l = (i==0)    ?  haloBuffer[j+2*Nx]    : U[iid - 1];
    const dfloat U_r = (i==Nx-1) ?  haloBuffer[j+2*Nx+Ny] : U[iid + 1];

    AU[iid] = (-U_l + 2*U[iid] - U_r)/(dx*dx) +
              (-U_d + 2*U[iid] - U_u)/(dy*dy);
  }
}

void HaloLaplacian(grid_t& grid, mesh_t& mesh,
                    dfloat* h_U,
                    dfloat* h_AU) {

  HaloLaplacianKernel(mesh.Nx, mesh.Ny, mesh.Nx,
                     mesh.dx, mesh.dy,
                     h_U, mesh.haloBuffer, h_AU);
}