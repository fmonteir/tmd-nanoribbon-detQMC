//
//  aux.hpp
//
//
//  Created by Francisco Brito on 06/05/2019.
//

#ifndef aux_hpp
#define aux_hpp

void write(int L, int sweep, int W, int A, double meanSign, double * weights,
  double nEl, double U, double nUp_nDw, double Hkin, Eigen::MatrixXd SiSjZ)
{
    //  Normalize to mean sign
    nEl /= meanSign; nUp_nDw /= meanSign; SiSjZ /= meanSign;
    Hkin /= meanSign;

    int precision = 10;
    Eigen::IOFormat CleanFmt(precision, 0, ", ", "\n", "", "");

    if (VERBOSE == 1)
    {
        std::cout << "Writing results" << std::endl << std::endl;
        std::cout << "<s>: " << meanSign << std::endl << std::endl;
        std::cout << "ds / <s>: " << sqrt( 1 - pow(meanSign, 2) )
         / sqrt( ( (sweep - W) / A - 1 ) ) / meanSign
          << std::endl << std::endl;
        std::cout << "nEl: " << nEl << std::endl << std::endl;
        std::cout << "Hkin: " << Hkin << std::endl << std::endl;
        std::cout << "U nUp_nDw: " << U * nUp_nDw << std::endl << std::endl;
    }

    //  STORE MEASUREMENTS
    std::ofstream file_log_weights("temp-data/Log-weights.csv");
    std::ofstream file_meas("temp-data/MeasurementsScalars.csv");
    std::ofstream file_spin_corr("temp-data/EqTimeSzCorrelations.csv");
    if ( file_log_weights.is_open() and file_meas.is_open()
     and file_spin_corr.is_open() )
    {
        file_log_weights << std::left << std::setw(50)
        << "Configuration log weight" << '\n';
        for (int s = 0; s < W; s++)
        {
            for (int slice = 0; slice < L; slice++)
            {
                file_log_weights << std::left << std::setw(50)
                << weights[s * L + slice] << '\n';
            }
        }
        file_meas << std::left << std::setw(50) << "Electron density <n>,";
        file_meas << std::left << std::setw(50) << std::setprecision(10)
        << nEl << '\n';
        file_meas << std::left << std::setw(50) << "Double occupancy <n+ n->,";
        file_meas << std::left << std::setw(50) << std::setprecision(10)
        << nUp_nDw << '\n';
        file_meas << std::left << std::setw(50) << "Hkin,";
        file_meas << std::left << std::setw(50) << std::setprecision(10)
        << Hkin << '\n';
        file_meas << std::left << std::setw(50) << "Average sign,";
        file_meas << std::left << std::setw(50) << std::setprecision(10)
        << meanSign << '\n';
        file_spin_corr << std::left << std::setw(50) << "<SiSj>" << '\n';
        file_spin_corr << std::setprecision(10)
        << SiSjZ.format(CleanFmt) << '\n';
    }
    file_log_weights.close();
    file_meas.close();
    file_spin_corr.close();
}

#endif /* aux_hpp */
