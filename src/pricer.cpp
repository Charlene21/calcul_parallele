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

struct data {
    double correlation;
    double fdstep;
    int hedgingDatesNumber;
    double interestRate;
    double maturity;
    int optionSize;
    string optionType;
    PnlVect *payoffCoefficients;
    double sampleNumber;
    PnlVect *spot;
    double strike;
    int TimeStepNumber;
    PnlVect trend;
    PnlVect volatility;
};

/*
 * 
 */
int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int size, rank;
    //affecte à rank le numéro du processus
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double time1;

    cout << "rank : " << rank << endl;
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

        //MPI_Barrier(MPI_COMM_WORLD);

        if (rank != 0) {
            sim = new Simulation(rank);
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
                    //sim->monte_carlo->price_master(prix, ic, tmp_sum, tmp_sum_square); 
                    cout << "prix en 0 lance par " << rank << ": " << prix << endl;
                    cout << "largeur de l'intervalle de confiance en 0 pour le prix lance par : " << ic << endl;

                }
            } else {

                sim->monte_carlo->price(prix, ic);
                cout << "prix en 0 : " << prix << endl;
                cout << "largeur de l'intervalle de confiance en 0 pour le prix : " << ic << endl;
            }

        } else {

            if (size <= 1) {
                cout << "Pour trouver le nombre de tirages, il est nécessaire d'utiliser au moins 2 processus" << endl;
            } else {

                char * prec = argv[2];
                double precision = atof(prec);
                cout << "precision : " << precision << endl;
                int nb_tirages = 1000;
                double variance, sigma_n, coeff;
                double memorized_sum = 0; 
                double memorized_sum_square = 0;
                int nb_tirages_previous = 0;

                cout << ">>NB_TIRAGES 1 : " << nb_tirages << endl;
                sim->monte_carlo->price_parallelisation(variance, size, rank, prix, ic, true, nb_tirages_previous, nb_tirages, memorized_sum, memorized_sum_square);

                nb_tirages_previous = nb_tirages;
                cout << "nb_tirages_previous (1) : " << nb_tirages_previous << endl;

                sigma_n = sqrt(variance);

                nb_tirages = ceil(pow((2 * 1.96 * sigma_n) / precision, 2));
                cout << "nb_tirages_previous (1) : " << nb_tirages_previous << endl;
                cout << "nb_tirages (1) : " << nb_tirages<< endl;

                sim->monte_carlo->price_parallelisation(variance, size, rank, prix, ic, true, nb_tirages_previous, nb_tirages, memorized_sum, memorized_sum_square);
                nb_tirages_previous = nb_tirages;
                cout << "nb_tirages_previous (2) : " << nb_tirages_previous << endl;
                sigma_n = sqrt(variance);
                //cout << ">>>>>sigma_n : " << sigma_n << endl;

                coeff = 2 * 1.96 * sigma_n / sqrt(nb_tirages);
                //cout << "coeff : " << coeff << endl;

                //cout << "BOUCLE >>>>>>>>>>>>>>>>>>>> " << endl;
                while (coeff > precision) {
                  //  cout << ">>NB_TIRAGES : " << nb_tirages << endl;
                    //sigma_n = sim->monte_carlo->estimation_variance(nb_tirages); 
                    nb_tirages = ceil(pow((2 * 1.96 * sigma_n) / precision, 2));
                    //cout << ">>NB_TIRAGES : " << nb_tirages << endl;
                    sim->monte_carlo->price_parallelisation(variance, size, rank, prix, ic, true, nb_tirages_previous, nb_tirages, memorized_sum, memorized_sum_square);
                    nb_tirages_previous = nb_tirages;
                    sigma_n = sqrt(variance);
                    coeff = 2 * 1.96 * sigma_n / sqrt(nb_tirages);


                }


               // cout << "END BOUCLE >>>>>>>>>>>>>>>>>>>> " << endl;
                printf("nb tirages : %d\n", nb_tirages);

              //  sim->monte_carlo->price_parallelisation(variance, size, rank, prix, ic, true, nb_tirages);
                if (rank == 0) {
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

