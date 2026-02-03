#include "simulation.hpp"
#include <emscripten.h>
#include <vector>

int main() { return 0; }

extern "C" EMSCRIPTEN_KEEPALIVE double* simulation(double* params_pointer) {
    static std::vector<double> results;
    
    results.clear();

    // params pointer: [m1, m2, l1, l2, t_start, t_end, steps, g, state[0], state[1], state[2], state[3], fps]
    // simulation args: (state = {}, g, m1, m2, l1, l2, t_start, t_end, steps)
    full_simulation(results, params_pointer[12], {params_pointer[8], params_pointer[9], params_pointer[10], params_pointer[11]}, params_pointer[7], params_pointer[0], params_pointer[1], params_pointer[2], params_pointer[3], params_pointer[4], params_pointer[5], params_pointer[6]);
    
    return results.data(); // t, th1, w1, th2, w2, l1, l2, l3, l4, h
}

