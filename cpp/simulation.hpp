#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <boost/numeric/odeint.hpp>


double g = 9.81;
double m1 = 2;
double m2 = 1;
double l1 = 2;
double l2 = 1;
double t_start = 0.0;
double t_end = 40.0;
double steps = 1001;
double dt = (t_end - t_start) / (steps - 1);

// Manual integration code data types

struct State {
    double th1; // Ang pos
    double w1; // Ang vel
    double th2; // Ang pos
    double w2; // Ang vel
};

struct LStates { // L for lyapunov
    struct State l_s1 {1.0, 0.0, 0.0, 0.0};
    struct State l_s2 {0.0, 1.0, 0.0, 0.0};
    struct State l_s3 {0.0, 0.0, 1.0, 0.0};
    struct State l_s4 {0.0, 0.0, 0.0, 1.0};
    double log_sum1 = 0;
    double log_sum2 = 0;
    double log_sum3 = 0;
    double log_sum4 = 0;
};

// Boost integration data types
using boost_state = std::vector<double>;
const int N_LYAP = 4;

inline double compute_energy(const State& s) {

    double T1 = 0.5 * (m1 + m2) * l1 * l1 * s.w1 * s.w1; // KE of m1
    double T2 = 0.5 * m2 * l2 * l2 * s.w2 * s.w2; // KE of m2
    double T3 = m2 * l1 * l2 * s.w1 * s.w2 * std::cos(s.th1 - s.th2); // KE of m2
    double T = T1 + T2 + T3; // Total KE
    
    double V1 = -(m1 + m2) * g * l1 * std::cos(s.th1); // PE of m1
    double V2 = -m2 * g * l2 * std::cos(s.th2); // PE of m2
    double V = V1 + V2; // Total PE
    
    return T + V; // Total Energy
}

// Helpers
inline void push_buffer(std::vector<double> &vec, const double t, const double th1, const double w1,const double th2, const double w2,const double l1,const double l2,const double l3,const double l4) {
    vec.push_back(t);
    vec.push_back(th1);
    vec.push_back(w1);
    vec.push_back(th2);
    vec.push_back(w2);
    vec.push_back(l1);
    vec.push_back(l2);
    vec.push_back(l3);
    vec.push_back(l4);
    vec.push_back(compute_energy({th1, w1, th2, w2}));
    // t, th1, w1, th2, w2, l1, l2, l3, l4, h
}

// Calculate derivatives (Equations of motion)
inline State compute_derivatives(const State& s) {
    double delta = s.th1 - s.th2;
    double den = 2 * m1 + m2 - m2 * std::cos(2* s.th1 - 2 * s.th2);

    double num1 = -g * (2 * m1 + m2) * std::sin(s.th1) - m2 * g * std::sin(s.th1 - 2 * s.th2) - 2 * std::sin(delta) * m2 * (s.w2 * s.w2 * l2 + s.w1 * s.w1 * l1 * std::cos(delta));
    double alpha1 = num1 / (l1 * den);

    double num2 = 2 * std::sin(delta) * (s.w1 * s.w1 * l1 * (m1 + m2) + g * (m1 + m2) * std::cos(s.th1) + s.w2 * s.w2 * l2 * m2 * std::cos(delta));
    double alpha2 = num2 / (l2 * den);

    return State{s.w1, alpha1, s.w2, alpha2};
}


// Helper function to add states
inline State add_states(const State& a, const State& b, double scale = 1.0) {
    return { 
        a.th1 + b.th1 * scale,
        a.w1 + b.w1 * scale,
        a.th2 + b.th2 * scale,
        a.w2 + b.w2 * scale
    };
}

// Helper to subtract states
inline State subtract_states(const State& a, const State& b, double scale = 1.0) {
    return { 
        a.th1 - b.th1 * scale,
        a.w1 - b.w1 * scale,
        a.th2 - b.th2 * scale,
        a.w2 - b.w2 * scale
    };
}

inline LStates add_lstates(const LStates& a, const LStates& b, double scale) {
    return { 
        add_states(a.l_s1, b.l_s1, scale),
        add_states(a.l_s2, b.l_s2, scale),
        add_states(a.l_s3, b.l_s3, scale),
        add_states(a.l_s4, b.l_s4, scale),
        a.log_sum1, a.log_sum2, a.log_sum3, a.log_sum4
    };
}

inline State multiply_state(const State& a, double scale) {
    return { 
        a.th1 * scale,
        a.w1 * scale,
        a.th2 * scale,
        a.w2 * scale
    };
}


// Converts from state struct to boost vector format, derives then converts back 
inline void boost_compute_derivatives(const boost_state &x, boost_state &dxdt) {
    State current = {x[0], x[1], x[2], x[3]};
    State current_derived = compute_derivatives(current);

    dxdt[0] = current_derived.th1;
    dxdt[1] = current_derived.w1;
    dxdt[2] = current_derived.th2;
    dxdt[3] = current_derived.w2;


}

inline State compute_sensitivity(const State& real, const State& fake, const double epsilon = 1e-6) {
    return multiply_state(subtract_states(compute_derivatives(add_states(real, fake, epsilon)), compute_derivatives(real)), 1/epsilon);
}

// Calculates Lyapunov System for Boost
inline void boost_lyap_sys(const boost_state &x, boost_state &dxdt, const double t) {
    State real = {x[0], x[1], x[2], x[3]};
    State d_real = compute_derivatives(real);

    dxdt[0] = d_real.th1;
    dxdt[1] = d_real.w1;
    dxdt[2] = d_real.th2;
    dxdt[3] = d_real.w2;


    for (int k = 0; k < N_LYAP; ++k) {
        int offset = 4 + (k * 4);
        
        State fake_vec = {x[offset], x[offset + 1], x[offset + 2], x[offset + 3]};

        State d_fake = compute_sensitivity(real, fake_vec);

        dxdt[offset] = d_fake.th1;
        dxdt[offset+1] = d_fake.w1;
        dxdt[offset+2] = d_fake.th2;
        dxdt[offset+3] = d_fake.w2;
    }

    
}
// Lyapunov Exponent
inline LStates fake_derivatives(const State& real, const LStates& fake) {
    return {
        compute_sensitivity(real, fake.l_s1),
        compute_sensitivity(real, fake.l_s2),
        compute_sensitivity(real, fake.l_s3),
        compute_sensitivity(real, fake.l_s4)
    };
}

inline double state_mag(const State& s) {
    return std::sqrt(s.th1*s.th1 + s.w1*s.w1 + s.th2*s.th2 + s.w2*s.w2);
}

inline double dot_product(const State& a, const State& b) {
    return (a.th1*b.th1 + a.w1*b.w1 + a.th2*b.th2 + a.w2*b.w2);
}

inline void gram_shmidt_shuffle(LStates& fake) {
    double mag1 = state_mag(fake.l_s1);
    fake.log_sum1 += std::log(mag1);

    fake.l_s1 = multiply_state(fake.l_s1, 1.0/mag1);

    double overlap12 = dot_product(fake.l_s2, fake.l_s1);
    fake.l_s2 = subtract_states(fake.l_s2, multiply_state(fake.l_s1, overlap12));
    
    double mag2 = state_mag(fake.l_s2);
    fake.log_sum2 += std::log(mag2);
    fake.l_s2 = multiply_state(fake.l_s2, 1.0/mag2);

    double overlap13 = dot_product(fake.l_s3, fake.l_s1);
    double overlap23 = dot_product(fake.l_s3, fake.l_s2);
    
    fake.l_s3 = subtract_states(fake.l_s3, multiply_state(fake.l_s1, overlap13));
    fake.l_s3 = subtract_states(fake.l_s3, multiply_state(fake.l_s2, overlap23));

    double mag3 = state_mag(fake.l_s3);
    fake.log_sum3 += std::log(mag3);
    fake.l_s3 = multiply_state(fake.l_s3, 1.0/mag3);

    double overlap14 = dot_product(fake.l_s4, fake.l_s1);
    double overlap24 = dot_product(fake.l_s4, fake.l_s2);
    double overlap34 = dot_product(fake.l_s4, fake.l_s3);

    fake.l_s4 = subtract_states(fake.l_s4, multiply_state(fake.l_s1, overlap14));
    fake.l_s4 = subtract_states(fake.l_s4, multiply_state(fake.l_s2, overlap24));
    fake.l_s4 = subtract_states(fake.l_s4, multiply_state(fake.l_s3, overlap34));
    
    double mag4 = state_mag(fake.l_s4);
    fake.log_sum4 += std::log(mag4);
    fake.l_s4 = multiply_state(fake.l_s4, 1.0/mag4);
}

inline std::vector<double> boost_gram_schmidt(boost_state &x) {
    LStates ghosts;

    ghosts.l_s1 = {x[4],  x[5],  x[6],  x[7]};
    ghosts.l_s2 = {x[8],  x[9],  x[10], x[11]};
    ghosts.l_s3 = {x[12], x[13], x[14], x[15]};
    ghosts.l_s4 = {x[16], x[17], x[18], x[19]};

    ghosts.log_sum1 = 0; ghosts.log_sum2 = 0; ghosts.log_sum3 = 0; ghosts.log_sum4 = 0;

    gram_shmidt_shuffle(ghosts);
    
    x[4]  = ghosts.l_s1.th1; x[5]  = ghosts.l_s1.w1; x[6]  = ghosts.l_s1.th2; x[7]  = ghosts.l_s1.w2;
    x[8]  = ghosts.l_s2.th1; x[9]  = ghosts.l_s2.w1; x[10] = ghosts.l_s2.th2; x[11] = ghosts.l_s2.w2;
    x[12] = ghosts.l_s3.th1; x[13] = ghosts.l_s3.w1; x[14] = ghosts.l_s3.th2; x[15] = ghosts.l_s3.w2;
    x[16] = ghosts.l_s4.th1; x[17] = ghosts.l_s4.w1; x[18] = ghosts.l_s4.th2; x[19] = ghosts.l_s4.w2;


    return {ghosts.log_sum1, ghosts.log_sum2, ghosts.log_sum3, ghosts.log_sum4};
}

inline LStates rk4_lyapunov_step(const LStates& fake, const State& current, double step_size) {

    State rk1 = compute_derivatives(current);
    State rmid1 = add_states(current, rk1, step_size / 2.0);
    State rk2 = compute_derivatives(rmid1);
    State rmid2 = add_states(current, rk2, step_size / 2.0);
    State rk3 = compute_derivatives(rmid2);
    State rmid3 = add_states(current, rk3, step_size);

    LStates k1 = fake_derivatives(current, fake);
    LStates k2 = fake_derivatives(rmid1, add_lstates(fake, k1, step_size / 2.0));
    LStates k3 = fake_derivatives(rmid2, add_lstates(fake, k2, step_size / 2.0));
    LStates k4 = fake_derivatives(rmid3, add_lstates(fake, k3, step_size));
    
    LStates next = fake;
    next.l_s1.th1 = fake.l_s1.th1 + (step_size / 6.0) * (k1.l_s1.th1 + 2 * k2.l_s1.th1 + 2 * k3.l_s1.th1 + k4.l_s1.th1);
    next.l_s1.w1 = fake.l_s1.w1 + (step_size / 6.0) * (k1.l_s1.w1 + 2 * k2.l_s1.w1 + 2 * k3.l_s1.w1 + k4.l_s1.w1);
    next.l_s1.th2 = fake.l_s1.th2 + (step_size / 6.0) * (k1.l_s1.th2 + 2 * k2.l_s1.th2 + 2 * k3.l_s1.th2 + k4.l_s1.th2);
    next.l_s1.w2 = fake.l_s1.w2 + (step_size / 6.0) * (k1.l_s1.w2 + 2 * k2.l_s1.w2 + 2 * k3.l_s1.w2 + k4.l_s1.w2);

    next.l_s2.th1 = fake.l_s2.th1 + (step_size / 6.0) * (k1.l_s2.th1 + 2 * k2.l_s2.th1 + 2 * k3.l_s2.th1 + k4.l_s2.th1);
    next.l_s2.w1 = fake.l_s2.w1 + (step_size / 6.0) * (k1.l_s2.w1 + 2 * k2.l_s2.w1 + 2 * k3.l_s2.w1 + k4.l_s2.w1);
    next.l_s2.th2 = fake.l_s2.th2 + (step_size / 6.0) * (k1.l_s2.th2 + 2 * k2.l_s2.th2 + 2 * k3.l_s2.th2 + k4.l_s2.th2);
    next.l_s2.w2 = fake.l_s2.w2 + (step_size / 6.0) * (k1.l_s2.w2 + 2 * k2.l_s2.w2 + 2 * k3.l_s2.w2 + k4.l_s2.w2);

    next.l_s3.th1 = fake.l_s3.th1 + (step_size / 6.0) * (k1.l_s3.th1 + 2 * k2.l_s3.th1 + 2 * k3.l_s3.th1 + k4.l_s3.th1);
    next.l_s3.w1 = fake.l_s3.w1 + (step_size / 6.0) * (k1.l_s3.w1 + 2 * k2.l_s3.w1 + 2 * k3.l_s3.w1 + k4.l_s3.w1);
    next.l_s3.th2 = fake.l_s3.th2 + (step_size / 6.0) * (k1.l_s3.th2 + 2 * k2.l_s3.th2 + 2 * k3.l_s3.th2 + k4.l_s3.th2);
    next.l_s3.w2 = fake.l_s3.w2 + (step_size / 6.0) * (k1.l_s3.w2 + 2 * k2.l_s3.w2 + 2 * k3.l_s3.w2 + k4.l_s3.w2);

    next.l_s4.th1 = fake.l_s4.th1 + (step_size / 6.0) * (k1.l_s4.th1 + 2 * k2.l_s4.th1 + 2 * k3.l_s4.th1 + k4.l_s4.th1);
    next.l_s4.w1 = fake.l_s4.w1 + (step_size / 6.0) * (k1.l_s4.w1 + 2 * k2.l_s4.w1 + 2 * k3.l_s4.w1 + k4.l_s4.w1);
    next.l_s4.th2 = fake.l_s4.th2 + (step_size / 6.0) * (k1.l_s4.th2 + 2 * k2.l_s4.th2 + 2 * k3.l_s4.th2 + k4.l_s4.th2);
    next.l_s4.w2 = fake.l_s4.w2 + (step_size / 6.0) * (k1.l_s4.w2 + 2 * k2.l_s4.w2 + 2 * k3.l_s4.w2 + k4.l_s4.w2);

    return next;
}

// Runge-Kutta 4 Integrator
inline State rk4_step(const State& current, double step_size) {
    State k1 = compute_derivatives(current);
    State k2 = compute_derivatives(add_states(current, k1, step_size / 2.0));
    State k3 = compute_derivatives(add_states(current, k2, step_size / 2.0));
    State k4 = compute_derivatives(add_states(current, k3, step_size));
    
    State next;
    next.th1 = current.th1 + (step_size / 6.0) * (k1.th1 + 2 * k2.th1 + 2 * k3.th1 + k4.th1);
    next.w1 = current.w1 + (step_size / 6.0) * (k1.w1 + 2 * k2.w1 + 2 * k3.w1 + k4.w1);
    next.th2 = current.th2 + (step_size / 6.0) * (k1.th2 + 2 * k2.th2 + 2 * k3.th2 + k4.th2);
    next.w2 = current.w2 + (step_size / 6.0) * (k1.w2 + 2 * k2.w2 + 2 * k3.w2 + k4.w2);
    
    return next;
}

inline void run_naive_rk4_simulation(std::vector<double>& output_buffer, bool verbose = true, State state = {1.0, -3.0, -1.0, 5.0}, const double saveEvery = 0.1) {
    LStates ghosts;

    double t = t_start;
    double lastSaveTime = t_start;
    
    if (verbose) std::cout << "Simulating fixed step size Runge-Kutta 4..." << std::endl;
    push_buffer(output_buffer, t, state.th1, state.w1, state.th2, state.w2, ghosts.log_sum1, ghosts.log_sum2, ghosts.log_sum3, ghosts.log_sum4);
    while (t < t_end) {
        // Write to buffer
        if (t >= (lastSaveTime + saveEvery)) {
            push_buffer(output_buffer, t, state.th1, state.w1, state.th2, state.w2, ghosts.log_sum1, ghosts.log_sum2, ghosts.log_sum3, ghosts.log_sum4);
            lastSaveTime += std::floor((t - lastSaveTime) / saveEvery) * saveEvery;
        }
        ghosts = rk4_lyapunov_step(ghosts, state, dt);
        state = rk4_step(state, dt);
        gram_shmidt_shuffle(ghosts);
        t += dt;
    }


    if (verbose) std::cout << "Simulation complete." << std::endl;

};


inline void run_adaptive_rk4_simulation(std::vector<double>& output_buffer, bool verbose = true, State state = {1.0, -3.0, -1.0, 5.0}, const double saveEvery = 0.1) {
    // Init cond
    LStates ghosts;

    if (verbose) std::cout << "Simulating adaptive step size Runge-Kutta 4..." << std::endl;
    
    double t = t_start;
    double step_size = dt;
    double error_tolerance = 1e-6;
    double error;
    double lastSaveTime = t_start;

    State full_step;
    State partial_step;
    int n = 2; // For later experimentation
    push_buffer(output_buffer, t, state.th1, state.w1, state.th2, state.w2, ghosts.log_sum1, ghosts.log_sum2, ghosts.log_sum3, ghosts.log_sum4);
    while (t < t_end) {
        // Write to file
        
        full_step = rk4_step(state, step_size);

        partial_step = rk4_step(state, step_size/n);
        for (int i = 0; i < n - 1; ++i) partial_step = rk4_step(partial_step, step_size/n);
        
        error = std::max({
            std::abs(full_step.th1 - partial_step.th1),
            std::abs(full_step.w1 - partial_step.w1),
            std::abs(full_step.th2 - partial_step.th2),
            std::abs(full_step.w2 - partial_step.w2)
        });
        
        if (error > error_tolerance) {
            step_size *= 0.5;
        } else {
            if (t >= (lastSaveTime + saveEvery)) {
                push_buffer(output_buffer, t, state.th1, state.w1, state.th2, state.w2, ghosts.log_sum1, ghosts.log_sum2, ghosts.log_sum3, ghosts.log_sum4);
                lastSaveTime += std::floor((t - lastSaveTime) / saveEvery) * saveEvery;
            }
            t += step_size;
            ghosts = rk4_lyapunov_step(ghosts, state, step_size);
            state = partial_step;
            gram_shmidt_shuffle(ghosts);
            if (error < error_tolerance/10 && step_size < dt) step_size *= 2;
        }
    }

    if (verbose) std::cout << "Adaptive simulation complete. Data saved to " << std::endl;
}

// Simulation parameters


inline void run_boost_rkd5_simulation(std::vector<double>& output_buffer, bool verbose = true, State state = {1.0, -3.0, -1.0, 5.0}, const double saveEvery = 0.1) {
    using namespace boost::numeric::odeint;
    
    boost_state x(20, 0.0);

    x[0] = state.th1; x[1] = state.w1; x[2] = state.th2; x[3] = state.w2;

    x[4] = 1.0; x[9] = 1.0; x[14] = 1.0; x[19] = 1.0;

    std::vector<double> lyap_sums = {0.0,0.0,0.0,0.0};

    if (verbose) std::cout << "Simulating with Boost Runge-Kutta-Dopri5" << std::endl;

    double absolute_error = 1e-6;
    double relative_error = 1e-6;
    double lastSaveTime = t_start;

    auto stepper = make_controlled(absolute_error, relative_error, runge_kutta_dopri5<boost_state>());
    
    double t = t_start;
    double normalization_interval = dt; 

    push_buffer(output_buffer, t, x[0], x[1], x[2], x[3], lyap_sums[0], lyap_sums[1], lyap_sums[2], lyap_sums[3]);
    while (t < t_end) {
        if (t >= (lastSaveTime + saveEvery)) {
            push_buffer(output_buffer, t, x[0], x[1], x[2], x[3], lyap_sums[0], lyap_sums[1], lyap_sums[2], lyap_sums[3]);
            lastSaveTime += std::floor((t - lastSaveTime) / saveEvery) * saveEvery;
        }
        size_t steps_taken = integrate_adaptive(
            stepper,
            boost_lyap_sys,
            x,
            t,
            t + normalization_interval,
            normalization_interval / 10.0
        );

        t += normalization_interval;

        std::vector<double> step_logs = boost_gram_schmidt(x);

        for(int i=0; i<4; i++) lyap_sums[i] += step_logs[i];
    }

    if (verbose) std::cout << "Boost simulation complete. Saved to: " << std::endl;
}

inline void full_simulation(std::vector<double>& output_buffer, const double fps, State _state = {1.0, -3.0, -1.0, 5.0}, const double _g = 9.81, const double _m1 = 2, const double _m2 = 1, const double _l1 = 2, const double _l2 = 1, const double _t_start = 0.0, const double _t_end = 40.0, const double _steps = 1001) {
    g = _g;
    m1 = _m1;
    m2 = _m2;
    l1 = _l1;
    l2 = _l2;
    t_start = _t_start;
    t_end = _t_end;
    steps = _steps;
    dt = (_t_end - _t_start) / (_steps - 1);

    const double saveEvery = (1.0/fps);

    output_buffer.push_back(0); output_buffer.push_back(0); output_buffer.push_back(0); // Reserve first slots

    run_naive_rk4_simulation(output_buffer, false, _state, saveEvery);
    output_buffer[1] = output_buffer.size();
    run_adaptive_rk4_simulation(output_buffer, false, _state, saveEvery);
    output_buffer[2] = output_buffer.size();
    run_boost_rkd5_simulation(output_buffer, false, _state, saveEvery);

    output_buffer[0] = output_buffer.size(); // Get vec size and set first slot to that
}


#endif
