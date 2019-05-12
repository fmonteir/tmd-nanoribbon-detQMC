//
//  simulation_eq_time.cpp
//
//
//  Created by Francisco Brito on 06/05/2019.
//
//  This program simulates the Hubbard model for a tmd nanoribbon
//  using auxiliary field (or determinant) Quantum Monte Carlo: in particular,
//  the BSS algorithm. This is version is optimized to compute observables that
//  require measuring equal-time Green's functions
//  The used notation is based on the lecture notes "Numerical Methods for
//  Quantum Monte Carlo Simulations of the Hubbard Model by Zhaojun Bai,
//  Wenbin Chen, Richard Scalettar, and Ichitaro Yamazaki (2009)
//

//  Number of threads
#ifndef NTH
#define NTH 4
#endif

//  Total number of "sites" (actual spatial sites plus orbitals)
#ifndef NX
#define NX 10
#endif

//  Width of the ribbon
#ifndef NY
#define NY 5
#endif

//  Inverse Trotter error
#ifndef DT_INV
#define DT_INV 8
#endif

//  Inverse temperature
#ifndef BETA
#define BETA 8
#endif

//  How often to calculate Green's functions afresh
//  (measured in number of imaginary-time slices)
#ifndef GREEN_AFRESH_FREQ
#define GREEN_AFRESH_FREQ 4
#endif

//  Output information about the progress of the run
#ifndef VERBOSE
#define VERBOSE 1
#endif

//   Number of orbitals
#ifndef NORB
#define NORB 3
#endif

//  Includes

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"
#include "model.hpp"
#include "green.hpp"
#include "aux.hpp"

int main(int argc, char **argv)
{
    if ( argc != 7) //  tmd, U, mu,
                    //  # sweeps, # warm-up sweeps, # auto-correlation sweeps
    {
        std::cout << "Not enough arguments given to simulation."
        << std::endl << std::endl;

        return -1;
    }

    double tmd = atoi(argv[1]);  //1,...,6 (MoS2, WS2, MoSe2, WSe2, MoTe2, WTe2)
    double U = atof(argv[2]);  //on-site interaction
    double mu = atof(argv[3]);  //chemical potential
    int totalMCSweeps = atof(argv[4]);  //number of sweeps
    int W = atof(argv[5]);  //number of warm-up sweeps
    int A = atof(argv[6]);  //number of auto-correlation sweeps

    double dt = 1. / DT_INV;  //  Trotter error. The error scales as dt^2
    const int NSITES = (int)(NORB * NX * NY);  //  # sites (real + orbital)
    const int L = (int)(BETA * DT_INV);  //  # slices
    //  Lbda = # intervals in which the product of B's is divided to stabilize.
    const int Lbda = (int)(L / GREEN_AFRESH_FREQ);
    //  HS transformation parameter (to order dtau^2)
    double nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12;

    double elDens = 0; double elDoubleOc = 0; double kineticEnergy = 0;
    Eigen::MatrixXd spin_corr =
      Eigen::Matrix<double, NORB * NY, NSITES>::Zero();
    double * av_weights = new double[W * L];
    double av_sign = 0;

    #pragma omp parallel num_threads(NTH)
    {
        //  RANDOM NUMBER GENERATION AND MONTE CARLO-RELATED VARIABLES.
        std::mt19937 gen;  //  mt19937 algorithm to generate random numbers
        std::uniform_real_distribution<> dis(0.0, 1.0);
        //  Uncorrelated random number generation
        std::random_device r;
        std::array<int, 624> seed_data;
        std::generate(seed_data.begin(), seed_data.end(), std::ref(r));
        std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
        gen.seed(seq);
        double decisionMaker; int totalMCSteps = totalMCSweeps * NSITES * L;


        // -- INITIALIZATION ---


        //  TIGHT-BINDING MODEL OF A TMD NANORIBBON
        Geometry< NSITES > K;
        K.setParams(tmd);
        K.HoppingMatrix();
        K.computeExponential(dt);

        //  INITIALIZE THE HS MATRIX WITH +1 AND -1 RANDOMLY.
        Configuration< L , NSITES > * h = new Configuration< L , NSITES >;
        h->genHsMatrix();

        //  GENERATE THE B-MATRICES.
        OneParticlePropagators< NSITES, L > * Bup =
          new OneParticlePropagators< NSITES, L >;
        OneParticlePropagators< NSITES, L > * Bdown=
          new OneParticlePropagators< NSITES, L >;
        Bup->fillMatrices( true, nu, h->matrix(), K.BpreFactor(dt, mu) );
        Bdown->fillMatrices( false, nu, h->matrix(), K.BpreFactor(dt, mu) );

        //  GENERATE THE SPIN-UP AND SPIN-DOWN GREEN FUNCTIONS.
        Green< NSITES, L, Lbda> * Gup = new Green< NSITES, L, Lbda>;
        Green< NSITES, L, Lbda> * Gdown = new Green< NSITES, L, Lbda>;
        Gup->storeTDQ( Bup->list() ); Gdown->storeTDQ( Bdown->list() );
        Gup->computeGreenFromTDQ(); Gdown->computeGreenFromTDQ();

        //  INITIALIZE RANK-ONE UPDATE-RELATED QUANTITIES AND ACCEPTANCE RATIO.
        double alphaUp; double alphaDown; double dUp; double dDown; double accRatio;

        //  INITIALIZE ARRAYS TO STORE MEASUREMENTS.
        // double * weights = new double[W * L];
        double LOGweight = 0.;

        double electronDensities = 0;
        double doubleOcs = 0;
        double energies = 0;
        Eigen::MatrixXd magCorrZZs =
          Eigen::Matrix<double, NORB * NY, NSITES>::Zero();

        double sign = 1; double meanSign = 0;

        double electronDensity; double doubleOc; double energy;
        Eigen::MatrixXd magCorrZZ = Eigen::Matrix<double, NORB * NY,NSITES>::Zero();

        double nEl = 0; double nUp_nDw = 0; double Hkin = 0;
        Eigen::MatrixXd SiSjZ =
          Eigen::Matrix<double, NORB * NY, NSITES>::Zero();

        //  INITIALIZE (l, i) <- (0, 0). INITIATIALIZE SPATIAL SWEEP COUNTER.
        //  FOR EACH IMAGINARY TIME SLICE l, LOOP OVER ALL SPATIAL LATTICE,
        //  THEN CHANGE SLICE, AND SO ON UNTIL l=L. REPEAT.
        int l = 0; int i = 0; int latticeSweepUntilAfresh = 0; int sweep = 0;


        // --- MC LOOP ---

        if (VERBOSE == 1)
        {
          std::cout << "\nMC loop started. Progress:\n";
        }

        for (int step = 0; step < totalMCSteps; step++)
        {
            if (VERBOSE == 1)
            {
                 // DISPLAY PROGRESS OF THE RUN.
                if ( (step + 1)  % (totalMCSteps/8) == 0 )
                {
                    std::cout << (step + 1) * 1. / totalMCSteps * 100 << " %"
                      << std::endl << std::endl;
                    std::cout << "Average Sign: " << meanSign
                      << std::endl << std::endl;
                    std::cout << "Log Weight: " << LOGweight
                      << std::endl << std::endl;
                }
            }
            //  COMPUTE THE ACCEPTANCE RATIO.
            alphaUp = ( exp( -2 * h->get(l, i) * nu ) - 1 );
            alphaDown = ( exp( 2 * h->get(l, i) * nu ) - 1 );
            dUp = ( 1 + alphaUp  * ( 1 - Gup->get(i, i) ) );
            dDown = ( 1 + alphaDown  * ( 1 - Gdown->get(i, i) ) );
            //  SAMPLING: METROPOLIS
            accRatio = fabs( dUp * dDown );

            //  DECIDE WHETHER OR NOT TO ACCEPT THE STEP.
            decisionMaker = dis(gen);

            if (decisionMaker <= accRatio )
            {
                //  KEEP TRACK OF WEIGHT
                LOGweight += log( fabs( dUp ) ) + log( fabs ( dDown ) );
                sign *= std::copysign(1, dUp * dDown );
                //  FLIP A SPIN
                h->flip(l, i);
                //  UPDATE Bs
                Bup->update(l, i, alphaUp); Bdown->update(l, i, alphaDown);
                //  RANK-ONE UPDATE -> O(N^2)
                Gup->update(alphaUp, dUp, i); Gdown->update(alphaDown, dDown, i);
            }


            // --- COMPUTE WRAPPED GREEN'S FUNCTIONS. ---


            if (i < NSITES - 1)
            {   //  CONTINUE LOOPING THROUGH THE SPATIAL LATTICE.
                i += 1;
            }
            else
            {
                //  EITHER WRAP OR COMPUTE GREEN'S FUNCTIONS FROM SCRATCH.
                latticeSweepUntilAfresh += 1;
                //  --- MEASUREMENTS ---
                if ( sweep < W )
                {
                    //  STORE WEIGHT OF ACCEPTED CONFIGURATIONS
                    // weights[sweep * L + l] = LOGweight;
                    av_weights[sweep * L + l] += LOGweight / 4;
                }

                //  STORE ELECTRON DENSITY, DOUBLE OCCUPANCY
                //  AND SPIN-SPIN CORRELATIONS.
                electronDensity = 0.; doubleOc = 0.; energy = 0.;
                magCorrZZ = Eigen::Matrix<double, NORB * NY, NSITES>::Zero();
                for (int a = 0; a < NORB; a++)
                for (int x1 = 0; x1 < NX; x1++)
                for (int y1 = 0; y1 < NY; y1++)
                for (int b = 0; b < NORB; b++)
                for (int x2 = 0; x2 < NX; x2++)
                for (int y2 = 0; y2 < NY; y2++)
                {
                    if (NORB*( NX * y1 + x1 ) + a == NORB*( NX * y2 + x2 ) + b)
                    {
                        electronDensity -= (Gup->get(NORB*( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        + Gdown->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b) );
                        doubleOc += - Gup->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        - Gdown->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        + Gup->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        * Gdown->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b);
                        energy -= 2 * ( Gup->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        + Gdown->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b) )
                        * (K.get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b) + mu);

                        magCorrZZ( NORB * y1 + a ,
                        NORB * ( NX * y2 + abs(x2 - x1) ) + b ) +=
                        ( Gup->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        + Gdown->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        - 2 * Gup->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        * Gdown->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b) ) / 2 / NX;

                        magCorrZZ( NORB * y1 + a ,
                        NORB * ( NX * y2 + ( NX - abs(x2 - x1) ) % NX ) + b ) +=
                        ( Gup->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        + Gdown->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        - 2 * Gup->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        * Gdown->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b) ) / 2 / NX;
                    }
                    else
                    {
                        energy -= ( Gup->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        + Gdown->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        + Gup->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y1 + x1 ) + a)
                        + Gdown->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y1 + x1 ) + a) )
                        * K.get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b);

                        magCorrZZ( NORB * y1 + a ,
                        NORB * ( NX * y2 + abs(x2 - x1) ) + b ) +=
                        ( - ( 1 - Gup->get(NORB * ( NX * y1 + x1 ) + a,
                        NORB * ( NX * y1 + x1 ) + a) )
                        * ( 1 - Gdown->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y2 + x2 ) + b) )
                        - ( 1 - Gdown->get(NORB * ( NX * y1 + x1 ) + a,
                        NORB * ( NX * y1 + x1 ) + a) )
                        * ( 1 - Gup->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y2 + x2 ) + b) )
                        + ( 1 - Gup->get(NORB * ( NX * y1 + x1 ) + a,
                        NORB * ( NX * y1 + x1 ) + a) )
                        * ( 1 - Gup->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y2 + x2 ) + b) )
                        + ( 1 - Gdown->get(NORB * ( NX * y1 + x1 ) + a,
                        NORB * ( NX * y1 + x1 ) + a) )
                        * ( 1 - Gdown->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y2 + x2 ) + b) )
                        - Gup->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y1 + x1 ) + a)
                        * Gup->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        - Gdown->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y1 + x1 ) + a)
                        * Gdown->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b) ) / 2 / NX;

                        magCorrZZ( NORB * y1 + a ,
                        NORB * ( NX * y2 + ( NX - abs(x2 - x1) ) % NX ) + b ) +=
                        ( - ( 1 - Gup->get(NORB * ( NX * y1 + x1 ) + a,
                        NORB * ( NX * y1 + x1 ) + a) )
                        * ( 1 - Gdown->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y2 + x2 ) + b) )
                        - ( 1 - Gdown->get(NORB * ( NX * y1 + x1 ) + a,
                        NORB * ( NX * y1 + x1 ) + a) )
                        * ( 1 - Gup->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y2 + x2 ) + b) )
                        + ( 1 - Gup->get(NORB * ( NX * y1 + x1 ) + a,
                        NORB * ( NX * y1 + x1 ) + a) )
                        * ( 1 - Gup->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y2 + x2 ) + b) )
                        + ( 1 - Gdown->get(NORB * ( NX * y1 + x1 ) + a,
                        NORB * ( NX * y1 + x1 ) + a) )
                        * ( 1 - Gdown->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y2 + x2 ) + b) )
                        - Gup->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y1 + x1 ) + a)
                        * Gup->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b)
                        - Gdown->get(NORB * ( NX * y2 + x2 ) + b,
                        NORB * ( NX * y1 + x1 ) + a)
                        * Gdown->get(NORB * ( NX * y1 + x1 ) + a ,
                        NORB * ( NX * y2 + x2 ) + b) ) / 2 / NX;
                    }
                }
                electronDensity /= NSITES; electronDensity += 2;
                doubleOc /= NSITES; doubleOc += 1;
                energy /= NSITES;

                electronDensities +=
                  ( electronDensity * sign - electronDensities ) / ( l + 1 ) ;
                doubleOcs +=
                  ( doubleOc * sign - doubleOcs ) / ( l + 1 ) ;
                magCorrZZs +=
                  (magCorrZZ * sign - magCorrZZs ) / ( l + 1 );
                energies +=
                  ( energy * sign - energies ) / ( l + 1 ) ;

                //  DEAL WITH THE GREEN'S FUNCTIONS.


                //  DECIDE WHETHER TO COMPUTE GREEN'S FUNCTIONS AFRESH OR TO WRAP.
                if (latticeSweepUntilAfresh == GREEN_AFRESH_FREQ)
                {   //  COMPUTE SPIN-UP AND SPIN-DOWN GREEN'S FUNCTIONS AFRESH.
                    if (l != ( L - 1 ) )
                    {
                        Gup->storeQDT(Bup->list(), l, GREEN_AFRESH_FREQ);
                        Gdown->storeQDT(Bdown->list(), l, GREEN_AFRESH_FREQ);
                        //  This is the standard way described in
                        //  "Stable simulations of models of interacting electrons"
                        Gup->computeStableGreen(l, GREEN_AFRESH_FREQ);
                        Gdown->computeStableGreen(l, GREEN_AFRESH_FREQ);
                    }
                    else
                    {
                        Gup->storeTDQ(Bup->list()); Gdown->storeTDQ(Bdown->list());
                        Gup->computeGreenFromTDQ(); Gdown->computeGreenFromTDQ();
                    }
                    latticeSweepUntilAfresh = 0;
                }
                else
                {   //  WRAPPING.
                    Gup->wrap( Bup->matrix(l) ); Gdown->wrap( Bdown->matrix(l) );
                }
                if (l < L - 1)
                {
                    l += 1; i = 0;
                }
                else
                {
                    if ( (sweep >= W) )
                    {
                        if ( sweep % A == 0 )
                        {
                          meanSign += ( sign - meanSign )
                           / ( ( sweep - W ) / A + 1 );
                          nEl += ( electronDensities - nEl )
                           / ( (sweep - W)/A + 1 ) ;
                          nUp_nDw += ( doubleOcs - nUp_nDw )
                           / ( (sweep - W)/A + 1 ) ;
                          SiSjZ += ( magCorrZZs - SiSjZ )
                           / ( (sweep - W)/A + 1 ) ;
                          Hkin += ( energies - Hkin )
                           / ( (sweep - W)/A + 1 ) ;

                        }
                        electronDensities = 0.; doubleOcs = 0.;
                        magCorrZZs=Eigen::Matrix<double, NORB * NY, NSITES>::Zero();
                        energies = 0.;

                    }
                    //  MOVE SWEEP COUNTER
                    sweep += 1;
                    l = 0; i = 0;

                }
            }
        }   //  END OF MC LOOP.

        //  Normalize to mean sign
        nEl /= meanSign; nUp_nDw /= meanSign; SiSjZ /= meanSign;
        Hkin /= meanSign;

        elDens += nEl;
        elDoubleOc += nUp_nDw;
        kineticEnergy += Hkin;
        spin_corr += SiSjZ;
        av_sign += abs(meanSign);

        delete Gup; delete Gdown; delete h; delete Bup; delete Bdown;
    }

    elDens /= NTH;
    elDoubleOc /= NTH;
    kineticEnergy /= NTH;
    spin_corr /= NTH;
    av_sign /= NTH;

    write(L, totalMCSweeps, W, A, av_sign, av_weights,
      elDens, U, elDoubleOc, kineticEnergy, spin_corr);

    delete[] av_weights;

    return 0;
}
