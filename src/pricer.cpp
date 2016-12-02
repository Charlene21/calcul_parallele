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
    if (rank == 0){
    time1 = MPI_Wtime();
    }

    //affecte à size le nombre de processus qui exécutent le programme
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int seed = time(NULL);
   
    PnlRng * rng = pnl_rng_dcmt_create_id(rank, seed);
    pnl_rng_sseed(rng, seed); 

    const char* extension = "-c";
    
    
    if (strcmp(argv[1], extension) == 0) {
        char *infile = argv[2];
        Param *P = new Parser(infile);
        Simulation *sim = new Simulation(P,rng);
        PnlVect * val_pf = pnl_vect_create(sim->nbTimeStepH + 1);
        PnlVect * price = pnl_vect_create(sim->nbTimeStepH + 1);
        double err = 0;

        sim->simu_couverture(val_pf, err, price);

        cout << "erreur P&L : " << err << endl;


    } else {
        char *infile = argv[1];          
        Param *P = new Parser(infile);
        Simulation *sim;
        sim = new Simulation(P,rng);
        double prix = 0;
        double ic = 0;


        if(argv[2] == NULL){
            
            if (size > 1){

//                double sum = 0;
//                double sum_square = 0;
                double variance = 0;
//                double tmp_sum = 0;
//                double tmp_sum_square = 0;

//                if(rank != 0)
//                sim->monte_carlo->price_slave(sum, sum_square, size, rank);

                //MPI_Barrier(MPI_COMM_WORLD);

//                MPI_Reduce(&sum, &tmp_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//                MPI_Reduce(&sum_square, &tmp_sum_square, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                //printf("END_REDUCE appele par n°%d : valeur de tmp_sum : %d\n", rank, tmp_sum);
                sim->monte_carlo->price_parallelisation(variance, size, rank, prix, ic, false, 0);
                
                cout << "VARIANCE : " << variance << endl;

                if (rank==0){
                //sim->monte_carlo->price_master(prix, ic, tmp_sum, tmp_sum_square); 
                cout << "prix en 0 lance par " << rank << ": " << prix << endl;
                cout << "largeur de l'intervalle de confiance en 0 pour le prix lance par : " << rank << " : " << ic << endl;

                }
            } else {

                sim->monte_carlo->price(prix, ic);
                cout << "prix en 0 : " << prix << endl;
                cout << "largeur de l'intervalle de confiance en 0 pour le prix : " << ic << endl;
            }
            
        }else{
//            if(size<=1){
//                cout << "Attention : il n"
//            }
            char * prec = argv[2];
            double precision = atof(prec);
            cout << "precision : " << precision << endl;
            int nb_tirages = 1000;
            double variance, sigma_n, coeff;
            
            cout << ">>NB_TIRAGES 1 : " << nb_tirages << endl;
            sim->monte_carlo->price_parallelisation(variance, size, rank, prix, ic, true, nb_tirages);
           
            cout << ">>>>>>>>>>>>>>>>>>>>>sigma_n : " << variance << endl;
            
            sigma_n = sqrt(variance);
            
            nb_tirages = floor(pow((2*1.96*sigma_n)/precision,2));
            cout << ">>NB_TIRAGES 2: " << nb_tirages << endl;
   
            sim->monte_carlo->price_parallelisation(variance, size, rank, prix, ic, true, nb_tirages);
          
            sigma_n = sqrt(variance);
            cout << ">>>>>sigma_n : " << sigma_n << endl;
            
            coeff = 2*1.96*sigma_n/sqrt(nb_tirages);
            cout << "coeff : " << coeff << endl;
            
            cout << "BOUCLE >>>>>>>>>>>>>>>>>>>> " << endl;
            while(coeff > precision)
            {    
                cout << ">>NB_TIRAGES : " << nb_tirages << endl;
                //sigma_n = sim->monte_carlo->estimation_variance(nb_tirages); 
                nb_tirages = floor(pow((2*1.96*sigma_n)/precision,2));
                printf("nb tirages boucle : %d\n",nb_tirages);
                cout << ">>NB_TIRAGES : " << nb_tirages << endl;
                sim->monte_carlo->price_parallelisation(variance, size, rank, prix, ic, true, nb_tirages);
                sigma_n = sqrt(variance);
                coeff = 2*1.96*sigma_n/sqrt(nb_tirages);


            } 
            
            
             cout << "END BOUCLE >>>>>>>>>>>>>>>>>>>> " << endl;
             printf("nb tirages : %d\n",nb_tirages);
             
             sim->monte_carlo->price_parallelisation(variance, size, rank, prix, ic, true, nb_tirages);
             if (rank == 0){
                sim->monte_carlo->price(prix, ic);
                cout << "prix en 0 : " << prix << endl;
                cout << "largeur de l'intervalle de confiance en 0 pour le prix : " << ic << endl;
             }

            }
            
        }

    double time2;
    if (rank==0){
    time2 = MPI_Wtime();
    cout << "Temps d'exécution du programme : " << time2-time1 << endl;
    }


    MPI_Finalize();
    
    
    return 0;
}

