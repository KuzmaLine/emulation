#include "QComputations_SINGLE.hpp"

namespace QComputations {

namespace {
    constexpr int QUDITS_COUNT = 2;

    int STEPS_COUNT = 100;

    // Во сколько раз dt должен быть меньше минимального полупериода
    constexpr int DT_DECREASE = 10;
    // Точность при нахождении общего полупериода
    constexpr double ACCURACY = 1e-15;
    const std::string LOWER_NUMBERS[] = {"₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"};

    // Добавить строчные маленькие цифры
    void add_lower(std::string& res, int num) {
        std::string tmp_res;
        do {
            tmp_res = LOWER_NUMBERS[num % 10] + tmp_res;
            num /= 10;
        } while(num != 0);

        res += tmp_res;
    }

    // Найти общий полупериод делённый на DT_DEACREASE
    double find_best_dt(const std::vector<double>& periods) {
        auto period = periods[0];
        for (size_t i = 1; i < periods.size(); i++) {
            int n = 1, k = 1;
            while(std::abs((period * n) / (periods[i] * k) - 1) >= ACCURACY) {
                if (period * n < periods[i] * k) {
                    while (period * n < periods[i] * k) {
                        n++;
                    }
                } else {
                    while (period * n > periods[i] * k) {
                        k++;
                    }
                }
            }

            period = period * n;
        }

        auto min_period = (*std::min(periods.begin(), periods.end()));
        while (period > min_period / DT_DECREASE) {
            period /= 2;
        }
        return period;
    }
}

// Реализация нашего состояния с дополненным понятием терминального состояния
class JC_State : public Basis_State {
    public:
        JC_State(size_t g_count): Basis_State(QUDITS_COUNT), g_map_(g_count + 1) {
            this->set_max_val(1000, 0);
        }
        JC_State(size_t g_count, const std::string& state_str): Basis_State(state_str), g_map_(g_count + 1) {
            this->set_max_val(1000, 0);
        }

        // Сделать шаг по i-тому g 
        JC_State JC_Step(int i) const {
            JC_State res(*this);
            int target = g_map_[this->get_qudit(0)][i].first;
            res.set_qudit(target, 0);
            res.ph_ = (ph_ + 1) % 2;
            res.set_qudit(ph_, 1);
            res.group_ = target;

            return res;
        }

        // Перевести текущее состояние в терминальное (если возможно)
        JC_State Terminal_Step() const {
            JC_State res(*this);

            int i = -1;
            while(f_[++i] != this->get_qudit(0)) {}

            res.set_qudit(g_map_.size() + i, 0);
            res.ph_ = 0;

            return res;
        }

        // Указать состояния, переходящие в терминальные
        void set_f(const std::vector<int>& f) { 
            f_ = f;
        }

        // Есть ли переход в терминальное состояние
        bool is_f() const {
            bool res = false;
            int i = -1;
            for (int i = 0; i < f_.size(); i++) {
                if (f_[i] == this->get_qudit(0)) return true;
            }

            return false;
        }

        // Является ли состояние терминальным
        bool is_terminal() const {
            return (this->get_qudit(0) >= g_map_.size());
        }

        //void set_energy(double energy) { state_energy_ = energy; }
        void set_energy(const std::vector<double>& energy) {state_energy_ = energy;}
        double get_energy() const { return state_energy_[this->get_qudit(0)]; }
        //double get_energy() const { return state_energy_; }

        // Получить фотон
        int get_photon() const { return ph_; }
        // Получить значение matter
        int get_matter() const { return this->get_qudit(1); }

        // Установить g
        void set_g(int i, int j, COMPLEX g) {
            for (auto& p: g_map_[i]) {
                if (p.first == j) {
                    p.second = g;

                    for (auto& p_second: g_map_[j]) {
                        if (p_second.first == i) {
                            p_second.second = g;
                        }
                    }

                    return;
                }
            }

            g_map_[i].emplace_back(std::make_pair(j, g));
            g_map_[j].emplace_back(std::make_pair(i, g));
        }

        // Число возможных переходов
        size_t get_neighbours_count() const { return (this->get_qudit(0) < g_map_.size() ? g_map_[this->get_qudit(0)].size() : 0); }
        // Получить g_i
        COMPLEX get_g(int i) const { return g_map_[this->get_qudit(0)][i].second; }

        // Строковое представление состояния
        std::string to_string() const override {
            std::string res = "|";

            if (this->is_terminal()) {
                res += "f";
            } else {
                res += std::to_string(ph_);
            }

            //res += std::to_string(this->get_qudit(0));

            add_lower(res, group_);

            res += ";" + std::to_string(this->get_matter()) + ">";

            return res;
        }

        // Перегрузка сортировки
        bool operator<(const Basis_State& other) const override {
            return this->to_string() < other.to_string();
        }
    private:
        std::vector<std::vector<std::pair<int, COMPLEX>>> g_map_;
        std::vector<int> f_;
        int ph_ = 0;
        int group_ = 0;
        std::vector<double> state_energy_;
        //double state_energy_;
        double sum_energy = QConfig::instance().h() * QConfig::instance().w();
};

// Реализация всевозможных переходов
State<JC_State> exc_relax_matter(const JC_State& st) {
    State<JC_State> res;
    for (size_t i = 0; i < st.get_neighbours_count(); i++) {
        res += State<JC_State>(st.JC_Step(i)) * st.get_g(i);
    }

    res += State<JC_State>(st) * st.get_energy();

    return res;
}

// Переход в терминальное состояние
State<JC_State> term_dec(const JC_State& st) {
    if(st.is_f()) {
        auto f = st.Terminal_Step();
         State<JC_State> res(f);
         res.insert(st, 0);

         return res;
    }

    return State<JC_State>(st);
}

}

int main(void) {
    using namespace QComputations;
    using OpType = Operator<JC_State>;  
    auto sum_energy = QConfig::instance().h() * QConfig::instance().w();
    auto vacuum_energy = QConfig::instance().h() * QConfig::instance().w() / 4;

    auto H_op = OpType(exc_relax_matter);
    OpType A_dec(term_dec);

    // Перебор g_3.
    std::vector<double> g_vec = {0.05, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1};

    for (auto g: g_vec) {
        STEPS_COUNT = 500;
        JC_State st(3, "|0;1>");
        st.set_g(0, 1, 0.1);
        st.set_g(0, 2, 0.4);
        st.set_g(0, 3, g);
        st.set_f({1, 2, 3});

        
        st.set_energy({sum_energy + vacuum_energy,           // |0_0, 1>
                       sum_energy + vacuum_energy,       // |1_1, 1>
                       sum_energy + vacuum_energy,   // |1_2, 1>
                       sum_energy + vacuum_energy,   // |1_3, 1>
                       sum_energy + vacuum_energy,       // |F_1, 1>
                       sum_energy + vacuum_energy,   // |F_2, 1>
                       sum_energy + vacuum_energy}); // |F_3, 1>
        
        //st.set_energy(sum_energy + vacuum_energy);

        double dt = find_best_dt({M_PI/0.2, M_PI/0.8, M_PI/(g*2)});


        // Находим базис
        auto basis = State_Graph<JC_State>(st, H_op, {A_dec}).get_basis();

        // Генерация самого гамильтониана
        H_by_Operator<JC_State> H(st, H_op, {std::make_pair(1, A_dec)});

        show_basis(H.get_basis());

        H.show();

        auto basis_size = basis.size();

        Matrix<double> probs(C_STYLE, basis_size, STEPS_COUNT + 1);

        State<JC_State> state(st, basis);

        for (size_t i = 0; i < basis_size; i++) {
            probs[i][0] = std::abs(state[i] * std::conj(state[i]));
        }
        
        for (size_t step = 1; step <= STEPS_COUNT; step++) {
            // Шаг моделирования уравнения Шрёдингера
            state = schrodinger_step(state, H, dt, basis);

            // Прогон оператора декогеренции. (Переход в терминальные состояния)
            state = A_dec.run(state);

            // Нормировка
            state.normalize();


            for (size_t i = 0; i < basis_size; i++) {
                probs[i][step] = std::abs(state[i] * std::conj(state[i]));
            }
        }


        // Вектор времени
        auto time_vec = linspace(0, STEPS_COUNT, STEPS_COUNT + 1);

        make_probs_files(H, probs, time_vec, H.get_basis(), "res/3 альтернативы. Вероятности. g" + LOWER_NUMBERS[1] + "=0.1, g" + LOWER_NUMBERS[2] + "=0.4, g" + LOWER_NUMBERS[3] + "=" + to_string_double_with_precision(g, 2, 4));
    }


    g_vec = {0.05, 0.1, 0.2, 0.4, 0.8, 1};

    for (auto g: g_vec) {
        JC_State st(2, "|0;1>");
        st.set_g(0, 1, 0.2);
        st.set_g(0, 2, g);
        st.set_f({1, 2});

        st.set_energy({sum_energy + vacuum_energy,
                sum_energy / 3 + vacuum_energy,
                sum_energy / 2 + vacuum_energy / 2,
                sum_energy * 2 / 3 + vacuum_energy,
                sum_energy / 2 + vacuum_energy / 2});

        //st.set_energy(sum_energy + vacuum_energy);


        double dt = find_best_dt({M_PI/0.4, M_PI/(g*2)});

        auto basis = State_Graph<JC_State>(st, H_op, {A_dec}).get_basis();

        H_by_Operator<JC_State> H(st, H_op, {std::make_pair(1, A_dec)});

        show_basis(H.get_basis());

        H.show();

        auto basis_size = basis.size();

        Matrix<double> phases(C_STYLE, basis_size, STEPS_COUNT + 1);

        Matrix<double> probs(C_STYLE, basis_size, STEPS_COUNT + 1);

        State<JC_State> state(st, basis);

        for (size_t i = 0; i < basis_size; i++) {
            probs[i][0] = std::abs(state[i] * std::conj(state[i]));
            phases[i][0] = std::arg(state[i]);
        }
        
        for (size_t step = 1; step <= STEPS_COUNT; step++) {
            state = schrodinger_step(state, H, dt, basis);
            state = A_dec.run(state);
            state.normalize();


            for (size_t i = 0; i < basis_size; i++) {
                phases[i][step] = std::arg(state[i]);
                probs[i][step] = std::abs(state[i] * std::conj(state[i]));
            }
        }

        auto time_vec = linspace(0, STEPS_COUNT, STEPS_COUNT + 1);

        make_probs_files(H, probs, time_vec, H.get_basis(), "res/2 альтернативы. Вероятности. g" + LOWER_NUMBERS[1] + "=0.2, g" + LOWER_NUMBERS[2] + "=" + to_string_double_with_precision(g, 2, 4));

        make_probs_files(H, phases, time_vec, H.get_basis(), "res/2 альтернативы. Фазы. g" + LOWER_NUMBERS[1] + "=0.2, g" + LOWER_NUMBERS[2] + "=" + to_string_double_with_precision(g, 2, 4));
    }

    {
        JC_State st(2, "|0;1>");
        st.set_g(0, 1, 0.2);
        st.set_g(0, 2, 0.5);
        st.set_f({1, 2});


        
        st.set_energy({sum_energy + vacuum_energy,
                sum_energy / 2 + vacuum_energy,
                sum_energy / 2 + vacuum_energy,
                sum_energy / 2  + vacuum_energy,
                sum_energy / 2 + vacuum_energy});
        
        //st.set_energy(sum_energy + vacuum_energy);
        double dt = find_best_dt({M_PI/0.4, M_PI});

        auto basis = State_Graph<JC_State>(st, H_op, {A_dec}).get_basis();

        H_by_Operator<JC_State> H(st, H_op, {std::make_pair(1, A_dec)});

        show_basis(H.get_basis());

        H.show();

        auto basis_size = basis.size();

        Matrix<double> probs(C_STYLE, basis_size, STEPS_COUNT + 1);

        State<JC_State> state(st, basis);

        for (size_t i = 0; i < basis_size; i++) {
            probs[i][0] = std::abs(state[i] * std::conj(state[i]));
        }
        
        for (size_t step = 1; step <= STEPS_COUNT; step++) {
            state = schrodinger_step(state, H, dt, basis);
            state = A_dec.run(state);
            state.normalize();


            for (size_t i = 0; i < basis_size; i++) {
                probs[i][step] = std::abs(state[i] * std::conj(state[i]));
            }
        }

        auto time_vec = linspace(0, STEPS_COUNT, STEPS_COUNT + 1);

        make_probs_files(H, probs, time_vec, H.get_basis(), "res/2 альтернативы. Вероятности. Равные соотношения энергий");
    }

    //!!!!!!!!!!!!
    return 0;

    {
        JC_State st(1, "|0;1>");

        double g = 0.1;
        double dt = M_PI/0.2/DT_DECREASE;

        st.set_g(0, 1, g);
        st.set_f({1});
        //STEPS_COUNT = 1000;

        //st.set_energy(sum_energy + vacuum_energy);

        auto basis = State_Graph<JC_State>(st, H_op, {A_dec}).get_basis();

        H_by_Operator<JC_State> H(st, H_op, {std::make_pair(1, A_dec)});

        show_basis(H.get_basis());

        H.show();

        auto basis_size = basis.size();

        Matrix<double> phases(C_STYLE, basis_size, STEPS_COUNT + 1);

        Matrix<double> probs(C_STYLE, basis_size, STEPS_COUNT + 1);

        State<JC_State> state(st, basis);

        for (size_t i = 0; i < basis_size; i++) {
            probs[i][0] = std::abs(state[i] * std::conj(state[i]));
            phases[i][0] = std::arg(state[i]);
        }
        
        for (size_t step = 1; step <= STEPS_COUNT; step++) {
            state = schrodinger_step(state, H, dt, basis);
            state = A_dec.run(state);
            state.normalize();

            for (size_t i = 0; i < basis_size; i++) {
                phases[i][step] = std::arg(state[i]);
                probs[i][step] = std::abs(state[i] * std::conj(state[i]));
            }
        }

        auto time_vec = linspace(0, STEPS_COUNT, STEPS_COUNT + 1);

        make_probs_files(H, probs, time_vec, H.get_basis(), "res/1 альтернатива. Вероятности. g=" + to_string_double_with_precision(g, 2, 4));

        make_probs_files(H, phases, time_vec, H.get_basis(), "res/1 альтернатива. Фазы. g=" + to_string_double_with_precision(g, 2, 4));
    }

    return 0;
}