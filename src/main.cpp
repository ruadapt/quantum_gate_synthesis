/** \file main.cpp
 */
#include "types.h"
#include "gridSynth.h"
#include "multiQubitSynthesis.h"
#include "cliffordT.h"
#include "matrix.h"
#include <boost/lexical_cast.hpp>

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage: ./main <angle in radians> <precision in bits> [effort = 25]" << std::endl;
        return 1;
    }
    const char *theta_str = argv[1];
    const char *prec_str = argv[2];
    Real theta(theta_str);
    Real prec(prec_str);
    const char *effort_str = (argc >= 4) ? argv[3] : "25";
    int effort = boost::lexical_cast<int>(effort_str);
    std::cout << "Using " << REAL_DIGITS << " digits of precision for Reals" << std::endl;
    U2<DOmega> uU;
    Maybe<double> error;
    List<std::tuple<DOmega, Integer, gridsynth::DStatus>> candidate_info;
    std::tie(uU, error, candidate_info) = gridsynth::gridsynth_stats<Real>(prec, theta, effort);
    List<Gate> gates = clifford_t::synthesis_u2(uU);
    std::cout << "Gates: " << clifford_t::to_string(gates) << std::endl;
    std::cout << "Error: " << (error.has_value() ? std::to_string(error.value()) : "Nothing") << std::endl;
}
