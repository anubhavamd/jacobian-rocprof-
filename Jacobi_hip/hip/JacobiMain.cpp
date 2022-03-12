//**************************************************************************
//* Copyright (c) 2019, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "Jacobi.hpp"

/**
 * @file JacobiMain.cpp
 * @brief This contains the application entry point
 */

/**
 * @brief The application entry point
 *
 * @param[in] argc  The number of command-line arguments
 * @param[in] argv  The list of command-line arguments
 */
int main(int argc, char ** argv)
{
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  grid_t grid;
  mesh_t mesh;

  // Extract topology and domain dimensions from the command-line arguments
  ParseCommandLineArguments(argc, argv,
                            comm,
                            grid,
                            mesh);

  Jacobi_t Jacobi(grid, mesh);

  Jacobi.Run();

  // Finalize the MPI process
  MPI_Finalize();
  return STATUS_OK;
}
