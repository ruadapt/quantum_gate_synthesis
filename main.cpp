#include "types.h"
#include "gridSynth.h"
#include "multiQubitSynthesis.h"
#include "cliffordT.h"

int main()
{
    Real prec = 5;
    Real theta = 1.7;
    int effort = 25;
    U2<DOmega> uU;
    Maybe<double> error;
    List<std::tuple<DOmega, Integer, gridsynth::DStatus>> candidate_info;
    std::tie(uU, error, candidate_info) = gridsynth::gridsynth_internal<Real>(prec, theta, effort);
    List<Gate> gates = clifford_t::synthesis_u2(uU);
    std::cout << "gates:" << std::endl;
    std::cout << clifford_t::to_string(gates) << std::endl;
}
