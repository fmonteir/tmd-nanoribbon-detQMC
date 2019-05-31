//
//  model.hpp
//
//
//  Created by Francisco Brito on 06/05/2019.
//

#ifndef model_hpp
#define model_hpp

struct node {
    //  This structure defines the model. A point or node on the M-atom
    //  triangular lattice has 3 d-orbitals from which the electrons hop.
    //  There are hoppings to the 3 orbitals of each of the 6 nearest neighbors.
    //  There are also on-site terms. See Phys. Rev. B 88, 085433 (2013)

    //  Number of orbitals in the tight-binding model
    unsigned int Norb;
    //  Number of hoppings from a point on the lattice
    unsigned int *NHoppings;
    //  Longitudinal distance between sites
    std::vector<int> *x_idxs;
    //  Transverse distance between sites
    std::vector<int> *y_idxs;
    //  d-orbital index (0, 1,..., Norb) of the orbital to which we hop
    std::vector<unsigned int> *end_orbs;
    //  Hopping
    std::vector<double> *hoppings;

    node(int number_of_orbitals) {
      Norb = number_of_orbitals;
      NHoppings = new unsigned[Norb];
      for (unsigned o = 0; o < Norb; o++)
      {
          NHoppings[o] = 0;
      }
      x_idxs = new std::vector<int>[Norb];
      y_idxs = new std::vector<int>[Norb];
      end_orbs = new std::vector<unsigned int>[Norb];
      hoppings = new std::vector<double>[Norb];
    }

    void add_hopping(unsigned start_orb, unsigned end_orb,
      int x_idx, int y_idx, double hopping)
    {
        NHoppings[start_orb]++; //  increment neighbor
        end_orbs[start_orb].push_back(end_orb); //  write d-orbital
        x_idxs[start_orb].push_back(x_idx); //  write longitudinal distance
        y_idxs[start_orb].push_back(y_idx); //  write transverse distance
        hoppings[start_orb].push_back(hopping); //  write hopping
    }
};

template<int N>
class Geometry
{
    Eigen::MatrixXd Bup;
    Eigen::MatrixXd Bdw;
    Eigen::MatrixXd Hoppings;
    double e0; double e1; double t0;
    double abs_t0; double t1; double t2;
    double t11; double t12; double t22;
public:
    void setParams(int tmd);
    void HoppingMatrix();
    void HoppingMatrixMF(double U_MF, Eigen::VectorXd nUp, Eigen::VectorXd nDw);
    void computeExponential(double dt);
    Eigen::MatrixXd BpreFactorUp(double dt, double mu);
    Eigen::MatrixXd BpreFactorDw(double dt, double mu);
    double get_t0();
    double get(int x, int y);
    Eigen::MatrixXd matrix();
    Geometry() : Bup(N, N), Bdw(N, N) {
    };
};

template<int N>
void Geometry<N>::computeExponential(double dt)
{
    Bup = (dt * Bup).exp();
    Bdw = (dt * Bdw).exp();
}

template<int N>
void Geometry<N>::setParams(int tmd)
{
    if (tmd == 1) //  MoS2
    {abs_t0 = 0.184; e0 = 1.046 / abs_t0; e1 = 2.104 / abs_t0;
     t1 = 0.401 / abs_t0; t2 = 0.507 / abs_t0; t11 = 0.218 / abs_t0;
     t12 = 0.338 / abs_t0; t22 = 0.057 / abs_t0;}
    if (tmd == 2) //  WS2
    {abs_t0 = 0.206; e0 = 1.130 / abs_t0; e1 = 2.275 / abs_t0;
     t1 = 0.567 / abs_t0; t2 = 0.536 / abs_t0; t11 = 0.286 / abs_t0;
     t12 = 0.384 / abs_t0; t22 = -0.061 / abs_t0;}
    if (tmd == 3) //  MoSe2
    {abs_t0 = 0.188;   e0 = 0.919 / abs_t0; e1 = 2.065 / abs_t0;
     t1 = 0.317 / abs_t0; t2 = 0.456 / abs_t0; t11 = 0.211 / abs_t0;
     t12 = 0.290 / abs_t0; t22 = 0.130 / abs_t0;}
    if (tmd == 4) //  WSe2
    {abs_t0 = 0.207; e0 = 0.943 / abs_t0; e1 = 2.179 / abs_t0;
     t1 = 0.457 / abs_t0; t2 = 0.486 / abs_t0; t11 = 0.263 / abs_t0;
     t12 = 0.329 / abs_t0;t22 = 0.034 / abs_t0;}
    if (tmd == 5) //  MoTe2
    {abs_t0 = 0.169; e0 = 0.605 / abs_t0; e1 = 1.972 / abs_t0;
     t1 = 0.228 / abs_t0; t2 = 0.390 / abs_t0; t11 = 0.207 / abs_t0;
     t12 = 0.239 / abs_t0; t22 = 0.252 / abs_t0;}
    if (tmd == 6) //  WTe2
    {abs_t0 = 0.175; e0 = 0.606 / abs_t0; e1 = 2.102 / abs_t0;
     t1 = 0.342 / abs_t0; t2 = 0.410 / abs_t0; t11 = 0.233 / abs_t0;
     t12 = 0.270 / abs_t0; t22 = 0.190 / abs_t0;}
}

template<int N>
void Geometry<N>::HoppingMatrix()
{
    t0 = -1;

    node n(NORB);
    // Add hoppings
    n.add_hopping(0, 0, 0, 0, e0);
    n.add_hopping(1, 1, 0, 0, e1);
    n.add_hopping(2, 2, 0, 0, e1);
    // R1
    n.add_hopping(0, 0, 1, 0, t0);
    n.add_hopping(1, 1, 1, 0, t11);
    n.add_hopping(2, 2, 1, 0, t22);
    n.add_hopping(0, 1, 1, 0, t1);
    n.add_hopping(0, 2, 1, 0, t2);
    n.add_hopping(1, 2, 1, 0, t12);
    n.add_hopping(1, 0, 1, 0, -t1);
    n.add_hopping(2, 0, 1, 0, t2);
    n.add_hopping(2, 1, 1, 0, -t12);
    // // R4
    // n.add_hopping(0, 0, -1, 0, t0);
    // n.add_hopping(1, 1, -1, 0, t11);
    // n.add_hopping(2, 2, -1, 0, t22);
    // n.add_hopping(0, 1, -1, 0, -t1);
    // n.add_hopping(0, 2, -1, 0, t2);
    // n.add_hopping(1, 2, -1, 0, -t12);
    // n.add_hopping(1, 0, -1, 0, t1);
    // n.add_hopping(2, 0, -1, 0, t2);
    // n.add_hopping(2, 1, -1, 0, t12);
    // R2
    n.add_hopping(0, 0, 1, -1, t0);
    n.add_hopping(1, 1, 1, -1, (t11 + 3 * t22)/4);
    n.add_hopping(2, 2, 1, -1, ( 3 * t11 + t22 ) / 4);
    n.add_hopping(0, 1, 1, -1, t1/2 - sqrt(3)*t2/2);
    n.add_hopping(0, 2, 1, -1, -sqrt(3)*t1/2 - t2/2);
    n.add_hopping(1, 2, 1, -1, sqrt(3)*(t22-t11)/4-t12);
    n.add_hopping(1, 0, 1, -1, -t1/2-sqrt(3)*t2/2);
    n.add_hopping(2, 0, 1, -1, sqrt(3)*t1/2-t2/2);
    n.add_hopping(2, 1, 1, -1, sqrt(3)*(t22-t11)/4+t12);
    // // R5
    // n.add_hopping(0, 0, -1, 1, t0);
    // n.add_hopping(1, 1, -1, 1, (t11 + 3 * t22 ) / 4);
    // n.add_hopping(2, 2, -1, 1, (3 * t11 + t22 ) / 4);
    // n.add_hopping(0, 1, -1, 1, -t1/2-sqrt(3)*t2/2);
    // n.add_hopping(0, 2, -1, 1, sqrt(3)*t1/2-t2/2);
    // n.add_hopping(1, 2, -1, 1, sqrt(3)*(t22-t11)/4+t12);
    // n.add_hopping(1, 0, -1, 1, t1/2-sqrt(3)*t2/2);
    // n.add_hopping(2, 0, -1, 1, -sqrt(3)*t1/2-t2/2);
    // n.add_hopping(2, 1, -1, 1, sqrt(3)*(t22-t11)/4-t12);
    // R3
    n.add_hopping(0, 0, 0, -1, t0);
    n.add_hopping(1, 1, 0, -1, ( t11 + 3 * t22 ) / 4);
    n.add_hopping(2, 2, 0, -1, ( 3 * t11 + t22 ) / 4);
    n.add_hopping(0, 1, 0, -1, -t1/2+sqrt(3)*t2/2);
    n.add_hopping(0, 2, 0, -1, -sqrt(3)*t1/2-t2/2);
    n.add_hopping(1, 2, 0, -1, -sqrt(3)*(t22-t11)/4+t12);
    n.add_hopping(1, 0, 0, -1, t1/2+sqrt(3)*t2/2);
    n.add_hopping(2, 0, 0, -1, sqrt(3)*t1/2-t2/2);
    n.add_hopping(2, 1, 0, -1, -sqrt(3)*(t22-t11)/4-t12);
    // // R6
    // n.add_hopping(0, 0, 0, 1, t0);
    // n.add_hopping(1, 1, 0, 1, ( t11 + 3 * t22 ) / 4);
    // n.add_hopping(2, 2, 0, 1, ( 3 * t11 + t22 ) / 4);
    // n.add_hopping(0, 1, 0, 1, t1/2+sqrt(3)*t2/2);
    // n.add_hopping(0, 2, 0, 1, sqrt(3)*t1/2-t2/2);
    // n.add_hopping(1, 2, 0, 1, -sqrt(3)*(t22-t11)/4-t12);
    // n.add_hopping(1, 0, 0, 1, -t1/2+sqrt(3)*t2/2);
    // n.add_hopping(2, 0, 0, 1, -sqrt(3)*t1/2-t2/2);
    // n.add_hopping(2, 1, 0, 1, -sqrt(3)*(t22-t11)/4+t12);

    Bup = Eigen::MatrixXd::Zero(NX * NY * NORB, NX * NY * NORB);
    Bdw = Eigen::MatrixXd::Zero(NX * NY * NORB, NX * NY * NORB);

    for (int x = 0; x < NX; x++ )
        for (int y = 0; y < NY; y++ )
          for (int orb = 0; orb < NORB; orb++ )
          {
              unsigned start_idx = orb + NORB * ( NX * y + x );
              for (unsigned neighbor_idx = 0;
                  neighbor_idx < n.NHoppings[orb]; neighbor_idx++)
              {
                  unsigned end_x =
                    ( x + n.x_idxs[orb].at(neighbor_idx) + NX ) % NX;
                  int end_y =
                    y + n.y_idxs[orb].at(neighbor_idx);
                  unsigned end_idx =
                    n.end_orbs[orb].at(neighbor_idx) + NORB
                    * ( NX * end_y + end_x );
                  if ( end_y >= 0 && end_y < NY )
                  {
                      Bup(start_idx, end_idx)
                         = -n.hoppings[orb].at(neighbor_idx);
                      Bup(end_idx, start_idx)
                         = -n.hoppings[orb].at(neighbor_idx);
                      Bdw(start_idx, end_idx)
                         = -n.hoppings[orb].at(neighbor_idx);
                      Bdw(end_idx, start_idx)
                         = -n.hoppings[orb].at(neighbor_idx);
                  }
              }
          }
    Hoppings = -Bup;
}

template<int N>
void Geometry<N>::HoppingMatrixMF(double U_MF, Eigen::VectorXd nUp, Eigen::VectorXd nDw)
{
    t0 = -1;

    node n(NORB);
    // Add hoppings
    n.add_hopping(0, 0, 0, 0, e0);
    n.add_hopping(1, 1, 0, 0, e1);
    n.add_hopping(2, 2, 0, 0, e1);
    // R1
    n.add_hopping(0, 0, 1, 0, t0);
    n.add_hopping(1, 1, 1, 0, t11);
    n.add_hopping(2, 2, 1, 0, t22);
    n.add_hopping(0, 1, 1, 0, t1);
    n.add_hopping(0, 2, 1, 0, t2);
    n.add_hopping(1, 2, 1, 0, t12);
    n.add_hopping(1, 0, 1, 0, -t1);
    n.add_hopping(2, 0, 1, 0, t2);
    n.add_hopping(2, 1, 1, 0, -t12);
    // // R4
    // n.add_hopping(0, 0, -1, 0, t0);
    // n.add_hopping(1, 1, -1, 0, t11);
    // n.add_hopping(2, 2, -1, 0, t22);
    // n.add_hopping(0, 1, -1, 0, -t1);
    // n.add_hopping(0, 2, -1, 0, t2);
    // n.add_hopping(1, 2, -1, 0, -t12);
    // n.add_hopping(1, 0, -1, 0, t1);
    // n.add_hopping(2, 0, -1, 0, t2);
    // n.add_hopping(2, 1, -1, 0, t12);
    // R2
    n.add_hopping(0, 0, 1, -1, t0);
    n.add_hopping(1, 1, 1, -1, (t11 + 3 * t22)/4);
    n.add_hopping(2, 2, 1, -1, ( 3 * t11 + t22 ) / 4);
    n.add_hopping(0, 1, 1, -1, t1/2 - sqrt(3)*t2/2);
    n.add_hopping(0, 2, 1, -1, -sqrt(3)*t1/2 - t2/2);
    n.add_hopping(1, 2, 1, -1, sqrt(3)*(t22-t11)/4-t12);
    n.add_hopping(1, 0, 1, -1, -t1/2-sqrt(3)*t2/2);
    n.add_hopping(2, 0, 1, -1, sqrt(3)*t1/2-t2/2);
    n.add_hopping(2, 1, 1, -1, sqrt(3)*(t22-t11)/4+t12);
    // // R5
    // n.add_hopping(0, 0, -1, 1, t0);
    // n.add_hopping(1, 1, -1, 1, (t11 + 3 * t22 ) / 4);
    // n.add_hopping(2, 2, -1, 1, (3 * t11 + t22 ) / 4);
    // n.add_hopping(0, 1, -1, 1, -t1/2-sqrt(3)*t2/2);
    // n.add_hopping(0, 2, -1, 1, sqrt(3)*t1/2-t2/2);
    // n.add_hopping(1, 2, -1, 1, sqrt(3)*(t22-t11)/4+t12);
    // n.add_hopping(1, 0, -1, 1, t1/2-sqrt(3)*t2/2);
    // n.add_hopping(2, 0, -1, 1, -sqrt(3)*t1/2-t2/2);
    // n.add_hopping(2, 1, -1, 1, sqrt(3)*(t22-t11)/4-t12);
    // R3
    n.add_hopping(0, 0, 0, -1, t0);
    n.add_hopping(1, 1, 0, -1, ( t11 + 3 * t22 ) / 4);
    n.add_hopping(2, 2, 0, -1, ( 3 * t11 + t22 ) / 4);
    n.add_hopping(0, 1, 0, -1, -t1/2+sqrt(3)*t2/2);
    n.add_hopping(0, 2, 0, -1, -sqrt(3)*t1/2-t2/2);
    n.add_hopping(1, 2, 0, -1, -sqrt(3)*(t22-t11)/4+t12);
    n.add_hopping(1, 0, 0, -1, t1/2+sqrt(3)*t2/2);
    n.add_hopping(2, 0, 0, -1, sqrt(3)*t1/2-t2/2);
    n.add_hopping(2, 1, 0, -1, -sqrt(3)*(t22-t11)/4-t12);
    // // R6
    // n.add_hopping(0, 0, 0, 1, t0);
    // n.add_hopping(1, 1, 0, 1, ( t11 + 3 * t22 ) / 4);
    // n.add_hopping(2, 2, 0, 1, ( 3 * t11 + t22 ) / 4);
    // n.add_hopping(0, 1, 0, 1, t1/2+sqrt(3)*t2/2);
    // n.add_hopping(0, 2, 0, 1, sqrt(3)*t1/2-t2/2);
    // n.add_hopping(1, 2, 0, 1, -sqrt(3)*(t22-t11)/4-t12);
    // n.add_hopping(1, 0, 0, 1, -t1/2+sqrt(3)*t2/2);
    // n.add_hopping(2, 0, 0, 1, -sqrt(3)*t1/2-t2/2);
    // n.add_hopping(2, 1, 0, 1, -sqrt(3)*(t22-t11)/4+t12);

    Bup = Eigen::MatrixXd::Zero(NX * NY * NORB, NX * NY * NORB);
    Bdw = Eigen::MatrixXd::Zero(NX * NY * NORB, NX * NY * NORB);

    for (int x = 0; x < NX; x++ )
        for (int y = 0; y < NY; y++ )
          for (int orb = 0; orb < NORB; orb++ )
          {
              unsigned start_idx = orb + NORB * ( NX * y + x );
              for (unsigned neighbor_idx = 0;
                  neighbor_idx < n.NHoppings[orb]; neighbor_idx++)
              {
                  unsigned end_x =
                    ( x + n.x_idxs[orb].at(neighbor_idx) + NX ) % NX;
                  int end_y =
                    y + n.y_idxs[orb].at(neighbor_idx);
                  unsigned end_idx =
                    n.end_orbs[orb].at(neighbor_idx) + NORB
                    * ( NX * end_y + end_x );
                  if ( end_y >= 0 && end_y < NY )
                  {
                      Bup(start_idx, end_idx)
                         = -n.hoppings[orb].at(neighbor_idx);
                      Bdw(start_idx, end_idx)
                         = -n.hoppings[orb].at(neighbor_idx);
                      Bup(end_idx, start_idx)
                         = -n.hoppings[orb].at(neighbor_idx);
                      Bdw(end_idx, start_idx)
                         = -n.hoppings[orb].at(neighbor_idx);
                  }
              }
          }
    Bup -= U_MF * Eigen::MatrixXd(nDw.asDiagonal());
    Bdw -= U_MF * Eigen::MatrixXd(nUp.asDiagonal());
    Hoppings = -Bup;
}

template<int N>
Eigen::MatrixXd Geometry<N>::BpreFactorUp(double dt, double mu)
{
    // add magnetic field via spin dependence
    return exp(dt * mu) * Bup;
}

template<int N>
Eigen::MatrixXd Geometry<N>::BpreFactorDw(double dt, double mu)
{
    // add magnetic field via spin dependence
    return exp(dt * mu) * Bdw;
}

template<int N>
double Geometry<N>::get_t0()
{
    return abs_t0;
}

template<int N>
double Geometry<N>::get(int x, int y)
{
    return Hoppings(x, y);
}

template<int N>
Eigen::MatrixXd Geometry<N>::matrix()
{
    return Hoppings;
}

template<int L, int N>
class Configuration
{
    Eigen::MatrixXd HSfield;
public:
    void genHsMatrix();
    Eigen::MatrixXd matrix();
    double get(int x, int y);
    void flip(int l, int i);
    Configuration() : HSfield(L, N) {
    };
};

template<int L, int N>
void Configuration<L, N>::genHsMatrix()
{
    //  Generate the HS field matrix
    int l;
    int i;
    HSfield = Eigen::Matrix<double, L, N>::Random(L,N);
    for (l = 0; l < L; l++)
    {
        for (i = 0; i < N; i++)
        {
            if ( HSfield(l, i) < 0 )
            {
                HSfield(l, i) = -1;
            }
            else
            {
                HSfield(l, i) = 1;
            }
        }
    }
}

template<int L, int N>
Eigen::MatrixXd Configuration<L, N>::matrix()
{
    return HSfield;
}

template<int L, int N>
double Configuration<L, N>::get(int x, int y)
{
    return HSfield(x, y);
}

template<int L, int N>
void Configuration<L, N>::flip(int l, int i)
{
    HSfield(l, i) *= -1;
}

template<int N, int L>
class OneParticlePropagators
{
    Eigen::MatrixXd B[L];
public:
    void fillMatrices(bool spin, double nu, Eigen::MatrixXd h,
      Eigen::MatrixXd BpreFactor);
    Eigen::MatrixXd matrix(int l);
    Eigen::MatrixXd * list();
    void update(int l, int i, double alpha);
};

template<int N, int L>
void OneParticlePropagators<N, L>::fillMatrices(bool spin, double nu,
   Eigen::MatrixXd h, Eigen::MatrixXd BpreFactor)
{
    int i;
    int l;
    for (l = 0; l < L; l++)
    {
        B[l] = Eigen::Matrix<double, N, N>::Zero();

        if (spin == true)
        {
            for (i = 0; i < N; i++)
            {
                B[l](i, i) = exp ( nu * h(l, i) ) ;
            }
        }
        else
        {
            for (i = 0; i < N; i++)
            {
                B[l](i, i) = exp ( - 1. * nu * h(l, i) );
            }
        }
        B[l] = BpreFactor * B[l] ;
    }
}

template<int N, int L>
void OneParticlePropagators<N, L>::update(int l, int i, double alpha)
{
    B[l].col(i) *= ( alpha + 1 );
}

template<int N, int L>
Eigen::MatrixXd OneParticlePropagators<N, L>::matrix(int l)
{
    return B[l];
}

template<int N, int L>
Eigen::MatrixXd * OneParticlePropagators<N, L>::list()
{
    return B;
}

#endif /* model_hpp */
