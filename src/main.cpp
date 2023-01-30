#include "types.h"
#include "gridSynth.h"
#include "multiQubitSynthesis.h"
#include "cliffordT.h"
#include "matrix.h"

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage: ./main <angle in radians> <precision in bits>" << std::endl;
        return 1;
    }
    char *prec_str = argv[2];
    Real prec(prec_str);
    char *theta_str = argv[1];
    Real theta(theta_str);
    int effort = 25; // TODO make this configurable?
    U2<DOmega> uU;
    Maybe<double> error;
    List<std::tuple<DOmega, Integer, gridsynth::DStatus>> candidate_info;
    std::tie(uU, error, candidate_info) = gridsynth::gridsynth_internal<Real>(prec, theta, effort);
    List<Gate> gates = clifford_t::synthesis_u2(uU);
    std::cout << "Gates: " << clifford_t::to_string(gates) << std::endl;
    std::cout << "Error: " << error.value_or(0) << std::endl;
}
