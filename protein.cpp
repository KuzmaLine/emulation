#include "emulation_state.hpp"

namespace {
    constexpr int STEPS_COUNT = 100000;
}

int main(void) {
    using namespace QComputations;
    using OpType = Operator<EmulationState>;  

    auto H_op = OpType(exc_relax_matter);
    OpType A_dec(term_dec);

    std::vector<double> alpha_vec = {0.0001, 0.001, 0.05}; //0.1, 0.2, 0.4, 0.8};

    for (auto alpha: alpha_vec) {
        double dt = find_best_dt({M_PI/0.4, M_PI/(0.2)}) / 10;

        std::cout << "DT = " << dt << std::endl;

        EmulationState st(6, "|0;0>");
        st.set_alpha(0, 1, alpha);
        st.set_alpha(1, 2, alpha);
        st.set_g(2, 3, 0.2);
        st.set_g(3, 4, 0.1);
        st.set_g(4, 5, 0.1);
        //st.set_alpha(5, 6, alpha);
        st.set_alpha(5, 0, alpha);
        //st.set_f({6});

        //st.set_energy({0.3, 1, 0.6, 0.1, 0.3, 0.6, 0.6}); 
        st.set_energy({0.3, 1, 0.6, 0.1, 0.3, 0.6}); 

        auto basis = State_Graph<EmulationState>(st, H_op, {A_dec}).get_basis();
        show_basis(basis);

        H_by_Operator<EmulationState> H(st, H_op, {std::make_pair(1, A_dec)});

        H.show();

        //operator_to_matrix(A_dec, basis).show();

        auto basis_size = basis.size();

        Matrix<double> probs_simple(C_STYLE, basis_size, STEPS_COUNT + 1);

        State<EmulationState> state_simple(basis);
        state_simple[0] = COMPLEX(1, 0);

        for (size_t i = 0; i < basis_size; i++) {
            probs_simple[i][0] = std::abs(state_simple[i] * std::conj(state_simple[i]));
        }
        
        for (size_t step = 1; step <= STEPS_COUNT; step++) {
            //state = schrodinger_step(state, H, dt, basis);
            state_simple = schrodinger_step(state_simple, H, dt, basis);

            //state = A_dec.run(state);
            //state.normalize();

            simple_term(state_simple);

            double sum_probs = 0;
            for (size_t i = 0; i < basis_size; i++) {
                //probs[i][step] = std::abs(state[i] * std::conj(state[i]));
                probs_simple[i][step] = std::abs(state_simple[i] * std::conj(state_simple[i]));
                sum_probs += probs_simple[i][step];
            }

            //std::cout << "SUM_PROBS = " << sum_probs << std::endl;
        }

        auto time_vec = linspace(0, STEPS_COUNT, STEPS_COUNT + 1);
        //auto probs_qme = quantum_master_equation(st, H_qme, time_vec);

        //make_probs_files(H_qme, probs_qme, time_vec, H_qme.get_basis(), "Protein/Protein. КОУ. Разные энергии. Вероятности. gamma = 0.3");
        //make_probs_files(H, probs, time_vec, H.get_basis(), "Protein/Protein. Через нормализацию. Разные энергии. Вероятности. alpha = " + std::to_string(alpha));
        make_probs_files(H, probs_simple, time_vec, H.get_basis(), "Protein/Protein. Вероятности. alpha = " + std::to_string(alpha));
    }

    return 0;
}