//**************************************************************************
//* Copyright (c) 2019, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include <unistd.h>
#include <string>
#include "Jacobi.hpp"

Jacobi_t::Jacobi_t(grid_t& grid_, mesh_t& mesh_):
  grid(grid_),
  mesh(mesh_)
{
  ApplyTopology();

  CreateMesh();

  InitializeData();
}

Jacobi_t::~Jacobi_t() {
  SafeHostFree(mesh.x);
  SafeHostFree(mesh.y);

  SafeHostFree(h_U);
  SafeHostFree(h_AU);
  SafeHostFree(h_RHS);

  SafeHostFree(mesh.sendBuffer);
  SafeHostFree(mesh.haloBuffer);
}

// ====================
// Topology application
// ====================

/**
 * @brief Generates the 2D topology and establishes the neighbor relationships between MPI processes
 */
void Jacobi_t::ApplyTopology() {

  // The number of MPI processes must fill the topology
  int size = grid.Nrow * grid.Ncol;
  if (size != grid.size) {
    OneErrPrintf((!grid.rank),
        "Error: The number of MPI processes (%d) doesn't match "
        "the topology size (%d).\n", size, grid.size);
    MPI_Abort(grid.comm, STATUS_ERR);
  }

  grid.mycol = grid.rank % grid.Ncol;
  grid.myrow = grid.rank / grid.Ncol;

  grid.Neighbor[SIDE_LEFT ] = (grid.mycol==0)           ? -1 : grid.rank-1;
  grid.Neighbor[SIDE_RIGHT] = (grid.mycol==grid.Ncol-1) ? -1 : grid.rank+1;

  grid.Neighbor[SIDE_DOWN ] = (grid.myrow==0)           ? -1 : grid.rank-grid.Ncol;
  grid.Neighbor[SIDE_UP   ] = (grid.myrow==grid.Nrow-1) ? -1 : grid.rank+grid.Ncol;
}

// ====================
// Mesh Creation
// ====================

/**
 * @brief Generates the 2D mesh
 */
void Jacobi_t::CreateMesh() {

  mesh.N = mesh.Nx * mesh.Ny;

  mesh.Nhalo = 2*mesh.Nx + 2*mesh.Ny;

  //domain dimensions
  mesh.Lx = (X_MAX) - (X_MIN);
  mesh.Ly = (Y_MAX) - (Y_MIN);

  //mesh spacing
  mesh.dx = mesh.Lx/(mesh.Nx*grid.Ncol+1);
  mesh.dy = mesh.Ly/(mesh.Ny*grid.Nrow+1);

  dfloat DX = mesh.dx*(mesh.Nx-1);
  dfloat DY = mesh.dy*(mesh.Ny-1);

  uint64 Ntotal = grid.size*mesh.N;
  mesh.invNtotal = 1.0/Ntotal;

  //coordinates (including boundary points)
  mesh.x = (dfloat *) malloc((mesh.Nx+2)*sizeof(dfloat));
  mesh.y = (dfloat *) malloc((mesh.Ny+2)*sizeof(dfloat));

  for (int i=0;i<mesh.Nx;i++)
    mesh.x[i] = mesh.dx + grid.mycol*(DX+mesh.dx) + i*mesh.dx;

  for (int j=0;j<mesh.Ny;j++)
    mesh.y[j] = mesh.dy + grid.myrow*(DY+mesh.dy) + j*mesh.dy;

  //halo exchange storage
  mesh.sendBuffer = (dfloat *) malloc(mesh.Nhalo*sizeof(dfloat));
  mesh.haloBuffer = (dfloat *) calloc(mesh.Nhalo,sizeof(dfloat));

  for (int i=0;i<2*NSIDES;i++) mesh.requests[i] = MPI_REQUEST_NULL;

  //number of entries for each side in send/recv buffers
  mesh.sideLength[SIDE_DOWN ] = mesh.Nx;
  mesh.sideLength[SIDE_UP   ] = mesh.Nx;
  mesh.sideLength[SIDE_LEFT ] = mesh.Ny;
  mesh.sideLength[SIDE_RIGHT] = mesh.Ny;

  //offets of halo side data in send/recv buffers
  mesh.sideOffset[SIDE_DOWN ] = 0;
  mesh.sideOffset[SIDE_UP   ] = mesh.Nx;
  mesh.sideOffset[SIDE_LEFT ] = 2*mesh.Nx;
  mesh.sideOffset[SIDE_RIGHT] = 2*mesh.Nx + mesh.Ny;
}

/**
 * @brief This allocates and initializes all the relevant data buffers before the Jacobi run
 *
 */
void Jacobi_t::InitializeData() {

  //host buffers
  h_U   = (dfloat*) malloc(mesh.N*sizeof(dfloat));
  h_AU  = (dfloat*) malloc(mesh.N*sizeof(dfloat));
  h_RHS = (dfloat*) malloc(mesh.N*sizeof(dfloat));
  h_RES = (dfloat*) malloc(mesh.N*sizeof(dfloat));

  for (int j=0;j<mesh.Ny;j++) {
    for (int i=0;i<mesh.Nx;i++) {
      int id = i+j*mesh.Nx;
      h_U[id] = 0.0; //initial guess
      h_RHS[id] = 0.0; //forcing
    }
  }

  // add boundary contributions
  uint64_t offsetX = grid.mycol*mesh.Nx+1;
  uint64_t offsetY = grid.myrow*mesh.Ny+1;
  dfloat totalX  = grid.Ncol*mesh.Nx+2;
  dfloat totalY  = grid.Nrow*mesh.Ny+2;

  if (grid.Neighbor[SIDE_DOWN]==-1) {
    for (int i=0;i<mesh.Nx;i++) {
      dfloat bc = sin(M_PI*(i+offsetX)/(totalX));
      h_RHS[i] += bc/(mesh.dx*mesh.dx);
    }
  }
  if (grid.Neighbor[SIDE_UP]==-1) {
    for (int i=0;i<mesh.Nx;i++) {
      dfloat bc = sin(M_PI*(i+offsetX)/(totalX));
      h_RHS[i+(mesh.Ny-1)*mesh.Nx] += bc/(mesh.dx*mesh.dx);
    }
  }

  if (grid.Neighbor[SIDE_LEFT]==-1) {
    for (int j=0;j<mesh.Ny;j++) {
      dfloat bc = sin(M_PI*(j+offsetY)/(totalY));
      h_RHS[j*mesh.Nx] += bc/(mesh.dy*mesh.dy);
    }
  }
  if (grid.Neighbor[SIDE_RIGHT]==-1) {
    for (int j=0;j<mesh.Ny;j++) {
      dfloat bc = sin(M_PI*(j+offsetY)/(totalY));
      h_RHS[(mesh.Nx-1)+j*mesh.Nx] += bc/(mesh.dy*mesh.dy);
    }
  }

  for (int j=0;j<mesh.Ny;j++) {
    for (int i=0;i<mesh.Nx;i++) {
      int id = i+j*mesh.Nx;
      h_AU[id] = 0.0;
      h_RES[id] = h_RHS[id];
    }
  }
}

// Display a number in a pretty format
char * FormatNumber(double value, const char * suffix, char * printBuf) {
  std::string magnitude = " kMGT";
  size_t orderIdx = 0;

  value = fabs(value);
  while ((value > 1000.0) && (orderIdx < magnitude.length() - 1))
  {
    ++orderIdx;
    value /= 1000.0;
  }

  sprintf(printBuf, "%.2lf %c%s", value, magnitude[orderIdx], suffix);

  return printBuf;
}

// Print a performance counter in a specific format
void PrintPerfCounter(const char * counterDesc, const char * counterUnit,
                      double counter, double elapsedTime, int size) {
  char printBuf[256];
  double avgCounter = counter / elapsedTime;
  double rankAvgCounter = avgCounter / size;

  printf("%s: %s (total), ", counterDesc, FormatNumber(avgCounter, counterUnit, printBuf));
  printf("%s (per process)\n", FormatNumber(rankAvgCounter, counterUnit, printBuf));
}

void Jacobi_t::PrintResults() {
  double lattUpdates = 0.0, flops = 0.0, bandWidth = 0.0;

  MPI_Barrier(grid.comm);

  double localCommTime = totalCommTime;
  double localComputeTime = totalLocalComputeTime;

  MPI_Allreduce(&localCommTime, &totalCommTime, 1, MPI_DOUBLE, MPI_SUM, grid.comm);
  MPI_Allreduce(&localComputeTime, &totalLocalComputeTime, 1, MPI_DOUBLE, MPI_SUM, grid.comm);

  // Show the performance counters
  if (!grid.rank) {
    printf("Total Jacobi run time: %.4lf sec.\n", elasped);
    // printf("Average per-process communication time: %.4lf sec.\n", avgTransferTime);

    // Compute the performance counters over all MPI processes
    lattUpdates = 1.0 * (grid.Ncol * mesh.Nx) * (grid.Nrow * mesh.Ny) * iterations;
    flops = 17.0 * (lattUpdates);              // Operations per Jacobi kernel run
    bandWidth = 12.0 * (lattUpdates) * sizeof(dfloat);     // Transfers per Jacobi kernel run

    double extraCommTime = totalCommTime - totalLocalComputeTime;
    extraCommTime = (extraCommTime<0) ? 0.0 : extraCommTime;

    double fractionCommHidden = 100.0 - 100.0*extraCommTime/totalCommTime;

    PrintPerfCounter("Measured lattice updates", "LU/s", lattUpdates, elasped, grid.size);
    PrintPerfCounter("Measured FLOPS", "FLOPS", flops, elasped, grid.size);
    PrintPerfCounter("Measured device bandwidth", "B/s", bandWidth, elasped, grid.size);

    if (grid.size>1)
      printf("Percentage of MPI traffic hidden by computation: %3.1f\n", fractionCommHidden);
  }
}
