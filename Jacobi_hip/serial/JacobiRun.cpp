//**************************************************************************
//* Copyright (c) 2019, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "Jacobi.hpp"

void Jacobi_t::Run() {

  OnePrintf((!grid.rank), "Starting Jacobi run.\n");

  iterations = 0;

  //compute initial residual (assumes zero initial guess)
  dfloat residual = Norm(grid, mesh, h_RES);

  OnePrintf((!grid.rank),
      "Iteration:   %d - Residual: %.6f\n", iterations, residual);

  MPI_Barrier(grid.comm);
  timerStart = MPI_Wtime();

  totalCommTime = 0.0;
  totalLocalComputeTime = 0.0;

  while ((iterations < JACOBI_MAX_LOOPS) && (residual > JACOBI_TOLERANCE)) {

    //record when comms start
    double commStart = MPI_Wtime();

    //Extract data off GPU exchange Halo with MPI
    HaloExchangeStart(grid, mesh, h_U);

    //record when local compute starts
    double localComputeStart = MPI_Wtime();

    //compute local part of Laplacian
    LocalLaplacian(grid, mesh, h_U, h_AU);

    //record when local compute ends
    double localComputeElapsed = MPI_Wtime() - localComputeStart;

    //finish Halo exchange with MPI
    HaloExchangeFinish(grid, mesh, h_U);

    //comms are complete at this point, record the elapsed time
    double commElapsed = MPI_Wtime() - commStart;

    //use halo data to complete Laplacian computation
    HaloLaplacian(grid, mesh, h_U, h_AU);

    //Jacobi iterative method
    // U = U + D^{-1}*(RHS - AU)
    JacobiIteration(grid, mesh, h_RHS, h_AU, h_RES, h_U);

    //residual = ||U||
    residual = Norm(grid, mesh, h_RES);

    //query the completed events to find the time the local laplacian took
    totalCommTime += commElapsed;
    totalLocalComputeTime += localComputeElapsed;

    ++iterations;
    OnePrintf((!grid.rank) && ((iterations) % 100 == 0),
      "Iteration: %d - Residual: %.6f\n", iterations, residual);
  }

  MPI_Barrier(grid.comm);
  timerStop = MPI_Wtime();
  elasped = timerStop - timerStart;

  OnePrintf((!grid.rank), "Stopped after %d iterations with residue %.6f\n", iterations, residual);

  PrintResults();
}

