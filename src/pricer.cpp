/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   pricer.cpp
 * Author: paviotch
 *
 * Created on September 26, 2016, 11:05 AM
 */

#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#include "Simulation.hpp"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int size, rank;
    //affecte à rank le numéro du processus
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double time1;

    if (rank == 0) {
        time1 = MPI_Wtime();
    }

    //affecte à size le nombre de processus qui exécutent le programme
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Initialisation de la graine
    int seed = time(NULL);

    //Création des générateurs aléatoires indépendants
    PnlRng * rng = pnl_rng_dcmt_create_id(rank, seed);
    pnl_rng_sseed(rng, seed);

    
    //Lecture du fichier donné en entrée
    if (argv[1] == NULL) {
        cout << "Attention : Vous devez fournir un fichier en entrée" << endl;
    } else {

        char *infile = argv[1];
        Param *P = new Parser(infile);
        Simulation *sim;

        if (rank == 0){
            sim = new Simulation(P, rng, size);
        }

        if (rank != 0) {
            sim = new Simulation();
        }

        double prix = 0;
        double ic = 0;

        if (argv[2] == NULL) {

            if (size > 1) {


                double variance = 0;
                double memorized_sum;
                double memorized_sum_square;

                sim->monte_carlo->price_parallelisation(variance, size, rank, prix, ic, false, 0, 0, memorized_sum, memorized_sum_square);
                if (rank == 0) {
                    cout << "Nombre de tirages : " << sim->monte_carlo->nbSamples_ << endl;
                    cout << "prix en 0 lance par " << rank << ": " << prix << endl;
                    cout << "largeur de l'intervalle de confiance en 0 pour le prix lance par : " << ic << endl;

                }
            } else {

                sim->monte_carlo->price(prix, ic);
                cout << "Nombre de tirages : " << sim->monte_carlo->nbSamples_ << endl;
                cout << "prix en 0 : " << prix << endl;
                cout << "largeur de l'intervalle de confiance en 0 pour le prix : " << ic << endl;
            }

        } else {

            if (size <= 1) {
                cout << "Pour trouver le nombre de tirages, il est nécessaire d'utiliser au moins 2 processus" << endl;
            } else {

                char * prec = argv[2];
                double precision = atof(prec);
                uintmax_t nb_tirages = 1000;
              
                double variance, sigma_n, coeff;
                double memorized_sum = 0; 
                double memorized_sum_square = 0;
                int nb_tirages_previous = 0;

                sim->monte_carlo->price_parallelisation(variance, size, rank, prix, ic, true, nb_tirages_previous, nb_tirages, memorized_sum, memorized_sum_square);

                nb_tirages_previous = nb_tirages; 

                sigma_n = sqrt(variance);

                nb_tirages = ceil(pow((2 * 1.96 * sigma_n) / precision, 2));

                sim->monte_carlo->price_parallelisation(variance, size, rank, prix, ic, true, nb_tirages_previous, nb_tirages, memorized_sum, memorized_sum_square);
                nb_tirages_previous = nb_tirages;
 
                sigma_n = sqrt(variance);

                coeff = 2 * 1.96 * sigma_n / sqrt(nb_tirages);
               
                while (coeff > precision) {
        
                    nb_tirages = ceil(pow((2 * 1.96 * sigma_n) / precision, 2));
                    sim->monte_carlo->price_parallelisation(variance, size, rank, prix, ic, true, nb_tirages_previous, nb_tirages, memorized_sum, memorized_sum_square);
                    nb_tirages_previous = nb_tirages;
                    sigma_n = sqrt(variance);
                    coeff = 2 * 1.96 * sigma_n / sqrt(nb_tirages);


                }
                
                
                if (rank == 0) {
                    printf("nb tirages : %d\n", nb_tirages);
                    cout << "prix en 0 : " << prix << endl;
                    cout << "largeur de l'intervalle de confiance en 0 pour le prix : " << ic << endl;
                }

            }

        }


    }
    double time2;
    if (rank == 0) {
        time2 = MPI_Wtime();
        cout << "Temps d'exécution du programme : " << time2 - time1 << endl;
    }
    MPI_Finalize();
    return 0;
}

