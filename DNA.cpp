#include "emulation_state.hpp"

int main(void) {
    using namespace QComputations;
    using OpType = Operator<EmulationState>;  

    auto H_op = OpType(exc_relax_matter);
    OpType A_dec(term_dec);

    std::vector<double> alpha_vec = {0, 0.05, 0.1, 0.2, 0.4, 0.8, 0.99999999};

    for (auto alpha: alpha_vec) {
        double dt = find_best_dt({M_PI/0.4, M_PI/(0.2)}) / 10;

        std::cout << "DT = " << dt << std::endl;

        EmulationState st(12, "|0;0>");
        st.set_g(0, 2, 0.1);
        st.set_alpha(2, 3, alpha);
        st.set_g(2, 4, 0.1);
        st.set_g(4, 6, 0.1);
        st.set_alpha(6, 8, alpha);
        st.set_alpha(8, 10, alpha);
        st.set_g(0, 1, 0.1);
        st.set_g(1, 5, 0.2);
        st.set_g(5, 7, 0.2);
        st.set_alpha(7, 9, alpha);
        st.set_alpha(9, 11, alpha);
        st.set_f({3, 10, 11});

        st.set_energy({1, 0.9, 0.9, 0.95, 0.8, 0.8, 0.6, 0.6, 0.2, 0.2, 0.2, 0.2}); 

        auto basis = State_Graph<EmulationState>(st, H_op, {A_dec}).get_basis();
        show_basis(basis);

        H_by_Operator<EmulationState> H(st, H_op, {std::make_pair(1, A_dec)});
        H_by_Operator<EmulationState> H_qme(st, H_op, {std::make_pair(0.3, OpType(true_term))});

        H.show();

        operator_to_matrix(A_dec, basis).show();

        auto basis_size = basis.size();

        Matrix<double> probs(C_STYLE, basis_size, STEPS_COUNT + 1);
        Matrix<double> probs_simple(C_STYLE, basis_size, STEPS_COUNT + 1);
        //Matrix<double> probs_H(C_STYLE, basis_size, STEPS_COUNT + 1);

        State<EmulationState> state(basis);
        State<EmulationState> state_qme(basis);
        State<EmulationState> state_simple(basis);
        state[0] = COMPLEX(1, 0);
        state_qme[0] = COMPLEX(1, 0);
        state_simple[0] = COMPLEX(1, 0);

        for (size_t i = 0; i < basis_size; i++) {
            probs[i][0] = std::abs(state[i] * std::conj(state[i]));
            probs_simple[i][0] = std::abs(state[i] * std::conj(state[i]));
        }
        
        for (size_t step = 1; step <= STEPS_COUNT; step++) {
            state = schrodinger_step(state, H, dt, basis);
            state_simple = schrodinger_step(state_simple, H, dt, basis);

            state = A_dec.run(state);
            state.normalize();

            simple_term(state_simple);

            double sum_probs = 0;
            for (size_t i = 0; i < basis_size; i++) {
                probs[i][step] = std::abs(state[i] * std::conj(state[i]));
                probs_simple[i][step] = std::abs(state_simple[i] * std::conj(state_simple[i]));
                sum_probs += probs_simple[i][step];
            }

            std::cout << "SUM_PROBS = " << sum_probs << std::endl;
        }

        auto time_vec = linspace(0, STEPS_COUNT, STEPS_COUNT + 1);
        auto probs_qme = quantum_master_equation(st, H_qme, time_vec);

        make_probs_files(H_qme, probs_qme, time_vec, H_qme.get_basis(), "DNA/DNA. КОУ. Разные энергии. Вероятности. gamma = 0.3");
        make_probs_files(H, probs, time_vec, H.get_basis(), "DNA/DNA. Через нормализацию. Разные энергии. Вероятности. alpha = " + std::to_string(alpha));
        make_probs_files(H, probs_simple, time_vec, H.get_basis(), "DNA/DNA. Упрощённый метод. Разные энергии. Вероятности. alpha = " + std::to_string(alpha));
    }

    return 0;
}