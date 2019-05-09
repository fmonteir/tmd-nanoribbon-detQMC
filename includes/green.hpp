//
//  green.hpp
//
//
//  Created by Francisco Brito on 06/05/2019.
//

#ifndef green_hpp
#define green_hpp

#include "QDT.hpp"
#include "TDQ.hpp"

template<int N, int L, int Lbda>
class Green
{
    Eigen::Matrix<double, -1, -1> G;
    Eigen::Matrix<double, -1, -1> u;
    Eigen::Matrix<double, -1, -1> w;
    Eigen::Matrix<double, -1, -1> Q;
    Eigen::Matrix<double, -1, -1> D;
    Eigen::Matrix<double, -1, -1> T;
    //  ALLOCATE MEMORY TO STORE THE PARTIAL PRODUCTS INVOLVED IN
    //  SPEEDING UP THE LOW TEMPERATURE STABILIZATION.
    Eigen::MatrixXd Qs[Lbda];
    Eigen::MatrixXd Ds[Lbda];
    Eigen::MatrixXd Ts[Lbda];
public:
    //  COMPUTING ELEMENTS OF THE GREEN'S FUNCTION
    void update(double alpha, double d, int i);
    void wrap(Eigen::MatrixXd B);
    void computeGreenFromTDQ();
    void storeTDQ(Eigen::MatrixXd* Bs);
    void computeStableGreen(int l, int greenAfreshFreq);
    void storeQDT(Eigen::MatrixXd* Bs, int l, int greenAfreshFreq);
    double get(int x, int y);
    Green() : G(N, N), u(N, 1), w(1, N), Q(N, N), D(N, N), T(N, N) {
    };
};

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::update(double alpha, double d, int i)
{
    u = ( Eigen::Matrix<double, N, N>::Identity() - G ).col(i);
    w = G.row(i);
    for (int x = 0; x < N; x++)
    {
        for (int y = 0; y < N; y++)
        {
            G(x, y) -= alpha / d * u(x) * w(y);
        }
    }
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::wrap(Eigen::MatrixXd B)
{
    G = B * G * B.inverse();
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::computeGreenFromTDQ()
{
    Q = Qs[Lbda - 1];
    D = Ds[Lbda - 1];
    T = Ts[Lbda - 1];
    //  COMPUTE GREEN'S FUNCTION.
    //  See eq. (4.3) of "Stable simulations of models of interacting electrons"
    //  by Loh Jr and Gubernatis
    TDQ< N > middleMatrix;
    Q = middleMatrix.QR_and_getQ( T.inverse() * Q.inverse() + D ) * Q; //  U' U
    D = middleMatrix.getD();
    //  Compute inverse of the diagonal matrix in the decomposition
    for (int i = 0; i < N; i++)
    {
        D(i, i) = 1 / D(i, i);
    }
    T = T * middleMatrix.getT();    //  V V'
    //  Final form of eq. (4.3)
    G = Q.inverse() * D * T.inverse();
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::storeTDQ(Eigen::MatrixXd* Bs)
{
    //  This is very similar to the function above. The only difference is that
    //  the Green's function is not computed, and this is to be used at the
    //  beginning of each imaginary time sweep, i.e. M = I + B_{L-1} ... B_{0}
    Q = Eigen::MatrixXd::Identity(N, N);
    D = Eigen::MatrixXd::Identity(N, N);
    T = Eigen::MatrixXd::Identity(N, N);

    int slice;
    int sliceCounter = 0;
    int lbda = 0;

    //  INITIALIZE PARTIAL PRODUCTS.
    Eigen::MatrixXd partialProdBs
    = Eigen::Matrix<double, N, N>::Identity();  //  an array of partial products

    //  COMPUTE UDVs.
    for (slice = L - 1; slice >= 0; slice--)
    {
        partialProdBs *= Bs[slice];
        sliceCounter += 1;
        if (sliceCounter == L / Lbda)
        {
            //  Householder QR applied to (previous partial product) * U * D
            //  where U, D are from the current partial product's decomposition
            TDQ< N > vduLeft;
            Q = vduLeft.QR_and_getQ( D * Q * partialProdBs );
            Qs[lbda] = Q;
            D = vduLeft.getD();
            Ds[lbda] = D;
            T = T * vduLeft.getT();
            Ts[lbda] = T;
            sliceCounter = 0;
            lbda += 1;
            partialProdBs = Eigen::Matrix<double, N, N>::Identity();
        }
    }
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::storeQDT(Eigen::MatrixXd* Bs, int l, int greenAfreshFreq)
{
    //  NOTE THAT THE ARGUMENT l will never be zero, since at l = 0,
    //  we compute the VDUs from the right and the Green's function
    int lbda = (l + 1) / greenAfreshFreq - 1;
    int slice;
    //  an array of partial products
    Eigen::MatrixXd partialProdBs = Eigen::Matrix<double, N, N>::Identity();

    //  INITIALIZE PARTIAL PRODUCTS.
    if ( lbda == 0 )
    {
        Q = Eigen::Matrix<double, N, N>::Identity();
        D = Eigen::Matrix<double, N, N>::Identity();
        T = Eigen::Matrix<double, N, N>::Identity();
    }
    else
    {
        Q = Qs[Lbda - lbda] ;
        D = Ds[Lbda - lbda] ;
        T = Ts[Lbda - lbda] ;
    }

    //  COMPUTE UDVs.
    for (slice = ( l - (greenAfreshFreq - 1) ) ; slice <= l ; slice++)
    {
        partialProdBs = Bs[slice] * partialProdBs;
    }

    //  Householder QR applied to (previous partial product) * U * D
    //  where U, and D are from the current partial product's decomposition
    QDT< N > udvRight;
    Q = udvRight.QR_and_getQ( partialProdBs * Q * D );
    D = udvRight.getD();
    T = udvRight.getT() * T;
    //  STORE THEM. NOTE THAT WE SHOULD NEVER REACH lbda = Lbda.
    //  At that point, we must store the VDUs to restart.
    //  Hence, this is safe! Otherwise, it wouldn't be.
    Qs[Lbda - lbda - 1] = Q;
    Ds[Lbda - lbda - 1] = D;
    Ts[Lbda - lbda - 1] = T;
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::computeStableGreen(int l, int greenAfreshFreq)
{
    //  NOTE THAT THE ARGUMENT l will never be zero, since at l = 0,
    //  we compute the VDUs from the right and the Green's function
    int lbda = (l + 1) / greenAfreshFreq - 1;
    QDT< N > udvIllConditionedSum;
    //  U_R.inv * U_L.inv + D_R * V_R * V_L * D_L
    Q = udvIllConditionedSum.QR_and_getQ(
      Qs[Lbda - lbda - 1].inverse() * Qs[Lbda - lbda - 2].inverse()
      + Ds[Lbda - lbda - 1] * Ts[Lbda - lbda - 1]
      * Ts[Lbda - lbda - 2] * Ds[Lbda - lbda - 2] );
    D = udvIllConditionedSum.getD();
    T = udvIllConditionedSum.getT();
    G = Qs[Lbda - lbda - 2].inverse() * T.inverse()
      * D.inverse() * Q.inverse() * Qs[Lbda - lbda - 1].inverse();
    //  U_L.inv V.inv D.inv U.inv U_R.inv
}

template<int N, int L, int Lbda>
double Green<N, L, Lbda>::get(int x, int y)
{
    return G(x, y);
}

#endif /* green_hpp */
