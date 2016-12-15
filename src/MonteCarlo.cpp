/* 
 * File:   MonteCarlo.cpp
 * Author: cpaviot
 * 
 * Created on 17 septembre 2016, 13:49
 */

#include "OptionBasket.hpp"
#include "OptionAsiatique.hpp"
#include "OptionPerformance.hpp"
#include "MonteCarlo.hpp"
#include "parser.hpp"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdexcept> 

using namespace std;

MonteCarlo::MonteCarlo() {
    mod_ = new BlackScholesModel();

    int bufsize, pos = 0;
    char *buf;
    MPI_Status status;
    char * type;
    int sizeType;
    double maturity;
    double strike;
    int nbTimeSteps;
    PnlVect* lambda = pnl_vect_create(mod_->size_);

    int tag = 0;

    MPI_Probe(0, tag, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_PACKED, &bufsize);
    buf = (char *) malloc(bufsize);

    MPI_Recv(buf, bufsize, MPI_PACKED, 0, tag, MPI_COMM_WORLD, &status);

    MPI_Unpack(buf, bufsize, &pos, lambda->array, lambda->size, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buf, bufsize, &pos, &maturity, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buf, bufsize, &pos, &nbTimeSteps, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buf, bufsize, &pos, &strike, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_Unpack(buf, bufsize, &pos, &sizeType, 1, MPI_INT, MPI_COMM_WORLD);

    type = (char*) malloc(sizeType);
    MPI_Unpack(buf, bufsize, &pos, type, sizeType, MPI_CHAR, MPI_COMM_WORLD);
    MPI_Unpack(buf, bufsize, &pos, &fdStep_, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buf, bufsize, &pos, &nbSamples_, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    if (((string) type).compare("asian") == 0) {
        opt_ = new OptionAsiatique(maturity, nbTimeSteps, mod_->size_, strike, lambda);
    } else if (((string) type).compare("basket") == 0) {
        opt_ = new OptionBasket(maturity, nbTimeSteps, mod_->size_, strike, lambda);
    } else if (((string) type).compare("performance") == 0) {
        opt_ = new OptionPerformance(maturity, nbTimeSteps, mod_->size_, lambda);
    }

    rng_ = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng_, time(NULL));
    shiftPlus_ = pnl_mat_new();
    shiftMoins_ = pnl_mat_new();
    path_ = pnl_mat_create_from_zero(this->opt_->nbTimeSteps_ + 1, this->mod_->size_);
}

MonteCarlo::MonteCarlo(Param *P) {

    mod_ = new BlackScholesModel(P);
    P->extract("fd step", fdStep_);

    //Option
    double maturity = 0;
    int nbTimeSteps = 0;
    double strike = 0;
    string type = "";
    P->extract("maturity", maturity);
    P->extract("TimeStep Number", nbTimeSteps);
    P->extract("strike", strike);
    P->extract("option type", type);
    PnlVect* lambda = pnl_vect_create(mod_->size_);
    P->extract("payoff coefficients", lambda, mod_->size_);

    P->extract("Sample Number", nbSamples_);

    if (type.compare("asian") == 0) {
        opt_ = new OptionAsiatique(maturity, nbTimeSteps, mod_->size_, strike, lambda);
    } else if (type.compare("basket") == 0) {
        opt_ = new OptionBasket(maturity, nbTimeSteps, mod_->size_, strike, lambda);
    } else if (type.compare("performance") == 0) {
        opt_ = new OptionPerformance(maturity, nbTimeSteps, mod_->size_, lambda);
    }

    rng_ = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng_, time(NULL));

    shiftPlus_ = pnl_mat_new();
    shiftMoins_ = pnl_mat_new();
    path_ = pnl_mat_create_from_zero(this->opt_->nbTimeSteps_ + 1, this->mod_->size_);
}

MonteCarlo::MonteCarlo(Param *P, PnlRng *rng, int size) {


    mod_ = new BlackScholesModel(P, size);
    P->extract("fd step", fdStep_);

    //Option
    double maturity = 0;
    int nbTimeSteps = 0;
    double strike = 0;
    string type = "";
    P->extract("maturity", maturity);
    P->extract("TimeStep Number", nbTimeSteps);
    P->extract("strike", strike);
    P->extract("option type", type);
    PnlVect* lambda = pnl_vect_create(mod_->size_);
    P->extract("payoff coefficients", lambda, mod_->size_);

    P->extract("Sample Number", nbSamples_);

    if (type.compare("asian") == 0) {
        opt_ = new OptionAsiatique(maturity, nbTimeSteps, mod_->size_, strike, lambda);
    } else if (type.compare("basket") == 0) {
        opt_ = new OptionBasket(maturity, nbTimeSteps, mod_->size_, strike, lambda);
    } else if (type.compare("performance") == 0) {
        opt_ = new OptionPerformance(maturity, nbTimeSteps, mod_->size_, lambda);
    }

    char * buf;
    int bufsize = 0;
    int count, pos = 0;

    //maturity, strike, fdstep et le nombre de doubles dans lambda
    MPI_Pack_size(4 + lambda->size, MPI_DOUBLE, MPI_COMM_WORLD, &count);
    bufsize += count;
    //nbTimeStep, nbSamples_ la taille de lambda est déjà transmise dans la classe BlackScholes
    MPI_Pack_size(2, MPI_INT, MPI_COMM_WORLD, &count);
    bufsize += count;
    MPI_Pack_size(type.size() + 1, MPI_CHAR, MPI_COMM_WORLD, &count);
    bufsize += count;

    buf = (char *) malloc(bufsize);

    //on pack les valeurs de lambda
    MPI_Pack(lambda->array, lambda->size, MPI_DOUBLE, buf, bufsize, &pos, MPI_COMM_WORLD);
    //on pack maturity
    MPI_Pack(&maturity, 1, MPI_DOUBLE, buf, bufsize, &pos, MPI_COMM_WORLD);
    //on pack nbTimeSteps
    MPI_Pack(&nbTimeSteps, 1, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
    //on pack strike
    MPI_Pack(&strike, 1, MPI_DOUBLE, buf, bufsize, &pos, MPI_COMM_WORLD);
    //on pack la taille de type
    int s = (int) type.size() + 1;

    MPI_Pack(&(s), 1, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
    //on pack ce que contient type
    MPI_Pack(type.c_str(), s, MPI_CHAR, buf, bufsize, &pos, MPI_COMM_WORLD);
    MPI_Pack(&fdStep_, 1, MPI_DOUBLE, buf, bufsize, &pos, MPI_COMM_WORLD);
    MPI_Pack(&nbSamples_, 1, MPI_DOUBLE, buf, bufsize, &pos, MPI_COMM_WORLD);


    for (int i = 1; i < size; i++) {
        MPI_Send(buf, bufsize, MPI_PACKED, i, 0, MPI_COMM_WORLD);
    }

    rng_ = rng;
    shiftPlus_ = pnl_mat_new();
    shiftMoins_ = pnl_mat_new();
    path_ = pnl_mat_create_from_zero(this->opt_->nbTimeSteps_ + 1, this->mod_->size_);

}

void MonteCarlo::price(double &prix, double &ic) {
    double sum = 0;
    double tmp = sum;
    double sum_square = 0;
    double variance = 0;
    double payoff = 0;

    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);


    for (int i = 0; i < nbSamples_; i++) {
        pnl_mat_set_all(path, 0);
        mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
        payoff = opt_->payoff(path);
        sum += payoff;
        sum_square += pow(payoff, 2);
    }

    pnl_mat_free(&path);

    tmp = sum;
    variance = getVariance(tmp, sum_square, 0);
    prix = getPrice(sum, 0);
    //std::cout << "VARIANCE : " << variance << std::endl;
    ic = getIntervalleConfiance(variance);

}

void MonteCarlo::price_parallelisation(double &variance, int size, int rank,
        double &prix, double &ic, bool cond, int nb_tirages_previous,
        int nb_tirages, double &memorized_sum, double &memorized_sum_square) {

    double sum = 0;
    double sum_square = 0;
    double tmp_sum = 0;
    double tmp_sum_square = 0;


    if (rank != 0) {
        if (!cond)
            this->price_slave(sum, sum_square, size, rank);
        else {

            if (nb_tirages - nb_tirages_previous > 0) {
                this->monte_carlo_slave(sum, sum_square, size, nb_tirages - nb_tirages_previous);
            } else {
                this->monte_carlo_slave(sum, sum_square, size, nb_tirages);
            }
        }
    }


    //On agrège les résultats pour ensuite les utiliser dans le processus 0
    MPI_Reduce(&sum, &tmp_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_square, &tmp_sum_square, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    if (rank == 0) {
        if (nb_tirages - nb_tirages_previous > 0) {
            tmp_sum += memorized_sum;
            tmp_sum_square += memorized_sum_square;
        }

        memorized_sum = tmp_sum;
        memorized_sum_square = tmp_sum_square;

        //si la précision n'est pas donnée
        if (!cond)
            price_master(prix, ic, tmp_sum, tmp_sum_square, variance);
        //si elle est donnée
        else {
            monte_carlo_master(prix, ic, tmp_sum, tmp_sum_square, variance, nb_tirages);
        }
    }


    if (cond) {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&variance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

}

void MonteCarlo::price_slave(double &sum, double &sum_square, int size, int rank) {

    double variance = 0;
    double payoff = 0;

    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);

    for (int i = 0; i < nbSamples_ / (size - 1); i++) {

        pnl_mat_set_all(path, 0);
        mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);

        payoff = opt_->payoff(path);
        sum += payoff;
        sum_square += pow(payoff, 2);
    }

    pnl_mat_free(&path);

}

void MonteCarlo::monte_carlo_slave(double &sum, double &sum_square, int size, int nb_tirages) {
    double payoff = 0;

    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);

    for (int i = 0; i < nb_tirages / (size - 1); i++) {
        pnl_mat_set_all(path, 0);
        mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
        payoff = opt_->payoff(path);
        sum += payoff;
        sum_square += pow(payoff, 2);
    }

    pnl_mat_free(&path);

}

void MonteCarlo::price_master(double &prix, double &ic, double sum, double sum_square, double &variance) {

    variance = getVariance(sum, sum_square, 0);
    prix = getPrice(sum, 0);
    //    std::cout << "VARIANCE : " << variance << std::endl;
    ic = getIntervalleConfiance(variance);

}

void MonteCarlo::monte_carlo_master(double &prix, double &ic, double sum, double sum_square, double &variance, int nb_tirages) {

    variance = getVarianceNbTirages(sum, sum_square, 0, nb_tirages);
    prix = getPriceNbTirages(sum, 0, nb_tirages);
    //    std::cout << "VARIANCE : " << variance << std::endl;
    ic = getIntervalleConfianceNbTirages(variance, nb_tirages);

}

void MonteCarlo::price(const PnlMat *past, double t, double &prix, double &ic) {
    double sum = 0;
    double tmp = sum;
    double sum_square = 0;
    double variance = 0;
    double payoff = 0;

    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);
    for (int i = 0; i < nbSamples_; i++) {
        mod_->asset(path, t, opt_->T_, opt_->nbTimeSteps_, rng_, past);
        payoff = opt_->payoff(path);
        sum += payoff;
        sum_square += pow(payoff, 2);
    }

    pnl_mat_free(&path);

    tmp = sum;
    variance = getVariance(tmp, sum_square, t);
    prix = getPrice(sum, t);
    // std::cout << "VARIANCE : " << variance << std::endl;
    ic = getIntervalleConfiance(variance);
}

void MonteCarlo::delta(const PnlMat *past, double t, PnlVect *delta) {
    double h = fdStep_;
    double payoff = 0;
    double prec = 0;

    PnlVect *sum_square = pnl_vect_create_from_zero(delta->size);
    PnlVect *ic = pnl_vect_create_from_zero(delta->size);

    pnl_vect_set_all(delta, 0);
    // Moyenne des payoffs
    for (int j = 0; j < nbSamples_; j++) {
        // Simulation du path
        if (t == 0) {
            this->mod_->asset(path_, this->opt_->T_, this->opt_->nbTimeSteps_, rng_);
        } else {
            this->mod_->asset(path_, t, this->opt_->T_, this->opt_->nbTimeSteps_, rng_, past);
        }

        // Shift_path
        for (int i = 0; i < this->mod_->size_; i++) {
            // Création des trajectoires shiftées
            this->mod_->shiftAsset(shiftPlus_, path_, i, h, t, this->opt_->T_ / this->opt_->nbTimeSteps_);
            this->mod_->shiftAsset(shiftMoins_, path_, i, -h, t, this->opt_->T_ / this->opt_->nbTimeSteps_);
            payoff = this->opt_->payoff(shiftPlus_) - this->opt_->payoff(shiftMoins_);
            prec = pnl_vect_get(delta, i);
            pnl_vect_set(delta, i, prec + payoff);

            //pour l'intervalle de confiance
            pnl_vect_set(sum_square, i, pnl_vect_get(sum_square, i) + pow(payoff, 2));
        }
    }
    //Pour l'intervalle de confiance en 0
    if (t == 0) {
        for (int i = 0; i < this->mod_->size_; i++) {
            pnl_vect_set(ic, i, sqrt((pnl_vect_get(sum_square, i) / nbSamples_) - (pow(pnl_vect_get(delta, i), 2) / pow(nbSamples_, 2))));
        }
    }

    double coeff = exp(-mod_->r_ * (opt_->T_ - t)) / (2 * nbSamples_ * h);

    pnl_vect_mult_scalar(ic, 2 * coeff * nbSamples_ * 1.96 / sqrt(nbSamples_));
    pnl_vect_mult_scalar(delta, coeff);
    PnlVect *copy = pnl_vect_create_from_zero(past->n);
    pnl_mat_get_row(copy, past, past->m - 1);

    pnl_vect_div_vect_term(delta, copy);

    pnl_vect_div_vect_term(ic, copy);
    if (t == 0) {
        cout << "largeur des intervalles de confiance pour DELTA en t=0 : " << endl;
        pnl_vect_print(ic);
    }

    pnl_vect_free(&copy);
    pnl_vect_free(&sum_square);
    pnl_vect_free(&ic);

}

MonteCarlo::~MonteCarlo() {
    pnl_rng_free(&rng_);
    pnl_mat_free(&shiftMoins_);
    pnl_mat_free(&shiftPlus_);
    pnl_mat_free(&path_);
}

/*  ----------- fonctions auxiliaires de factorisation du code ----------  */

double MonteCarlo::getVariance(double sum, double sum_square, double t) {
    sum /= nbSamples_;
    sum = pow(sum, 2);

    sum_square /= nbSamples_;
    return exp(-2 * mod_->r_ * (opt_->T_ - t))*(sum_square - sum);
}

double MonteCarlo::getVarianceNbTirages(double sum, double sum_square, double t, int nb_tirages) {
    sum /= nb_tirages;
    sum = pow(sum, 2);

    sum_square /= nb_tirages;
    return exp(-2 * mod_->r_ * (opt_->T_ - t))*(sum_square - sum);
}

double MonteCarlo::getIntervalleConfiance(double variance) {
    return 2 * 1.96 * sqrt(variance / nbSamples_);
}

double MonteCarlo::getPrice(double sum, double t) {
    sum *= exp(-mod_->r_ * (opt_->T_ - t)) / nbSamples_;
    return sum;
}

double MonteCarlo::getPriceNbTirages(double sum, double t, int nbTirages) {
    sum *= exp(-mod_->r_ * (opt_->T_ - t)) / nbTirages;
    return sum;
}

double MonteCarlo::getIntervalleConfianceNbTirages(double variance, int nbTirages) {
    return 2 * 1.96 * sqrt(variance / nbTirages);
}

/*  ----------- fonctions déterministes pour les tests ----------  */

void MonteCarlo::price(double &prix, double &ic, PnlVect *G) {

    double sum = 0;
    double tmp = sum;
    double sum_square = 0;
    double variance = 0;
    double payoff = 0;
    prix = 0;
    ic = 0;

    PnlMat *path = pnl_mat_create_from_ptr(1, mod_->size_, mod_->spot_->array);

    for (int i = 0; i < nbSamples_; i++) {

        mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, G);
        payoff = opt_->payoff(path);
        sum += payoff;
        sum_square += pow(payoff, 2);
    }
    pnl_mat_free(&path);

    tmp = sum;
    tmp /= nbSamples_;
    tmp = pow(tmp, 2);
    sum_square /= nbSamples_;

    variance = exp(-2 * mod_->r_ * opt_->T_)*(sum_square - tmp);
    //std::cout << "VARIANCE : " << variance << std::endl;
    sum *= exp(-mod_->r_ * opt_->T_) / nbSamples_;
    prix = sum;
    ic = 2 * 1.96 * sqrt(variance / nbSamples_);
}

void MonteCarlo::price(const PnlMat *past, double t, double &prix, double &ic, PnlVect *G) {
    double sum = 0;
    double tmp = sum;
    double sum_square = 0;
    double variance = 0;
    double payoff;
    prix = 0;
    ic = 0;


    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);

    for (int i = 0; i < nbSamples_; i++) {
        mod_->asset(path, t, opt_->T_, opt_->nbTimeSteps_, G, past);
        payoff = opt_->payoff(path);
        sum += payoff;
        sum_square += pow(payoff, 2);

    }

    pnl_mat_free(&path);
    tmp = sum;
    variance = getVariance(tmp, sum_square, t);
    prix = getPrice(sum, t);
    //    std::cout << "VARIANCE : " << variance << std::endl;
    ic = getIntervalleConfiance(variance);
}

void MonteCarlo::delta(const PnlMat *past, double t, PnlVect *delta, PnlVect *vect) {

    double h = fdStep_;
    double payoff = 0;
    double prec = 0;
    PnlMat *path = pnl_mat_create_from_zero(this->opt_->nbTimeSteps_ + 1, this->mod_->size_);
    PnlMat *shiftPlus = pnl_mat_new();
    PnlMat *shiftMoins = pnl_mat_new();


    pnl_vect_set_all(delta, 0);
    // Moyenne des payoffs
    for (int j = 0; j < nbSamples_; j++) {
        // Simulation du path
        if (t == 0) {
            this->mod_->asset(path, this->opt_->T_, this->opt_->nbTimeSteps_, vect);
        } else {
            this->mod_->asset(path, t, this->opt_->T_, this->opt_->nbTimeSteps_, vect, past);
        }

        // Shift_path
        for (int i = 0; i < this->mod_->size_; i++) {
            // Création des trajectoires shiftées

            this->mod_->shiftAsset(shiftPlus, path, i, h, t, this->opt_->T_ / this->opt_->nbTimeSteps_);
            this->mod_->shiftAsset(shiftMoins, path, i, -h, t, this->opt_->T_ / this->opt_->nbTimeSteps_);
            payoff = this->opt_->payoff(shiftPlus) - this->opt_->payoff(shiftMoins);
            prec = pnl_vect_get(delta, i);
            pnl_vect_set(delta, i, prec + payoff);
        }
    }

    double coeff = exp(-mod_->r_ * (opt_->T_ - t)) / (2 * nbSamples_ * h);
    pnl_vect_mult_scalar(delta, coeff);
    PnlVect *copy = pnl_vect_create_from_zero(past->n);
    pnl_mat_get_row(copy, past, past->m - 1);
    pnl_vect_div_vect_term(delta, copy);

    pnl_vect_free(&copy);
    pnl_mat_free(&shiftPlus);
    pnl_mat_free(&shiftMoins);
    pnl_mat_free(&path);
}

