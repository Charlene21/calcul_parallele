/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   newmain.cpp
 * Author: paviotch
 *
 * Created on December 7, 2016, 4:46 PM
 */

#include <cstdlib>
#include <iostream>
#include <mpi.h>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {


    MPI_Init(&argc, &argv);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int number;
    if (world_rank == 0) {
        number = -1;
        cout << "send " << endl;
        MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        cout << "send " << endl;
    } else if (world_rank == 1) {
        cout << "rcv " << endl;
        MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        cout << "Process 1 received number " << number << "from process 0" << endl;
    }

    MPI_Finalize();
    return 0;
}

