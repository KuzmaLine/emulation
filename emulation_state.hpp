#pragma once
#include "QComputations_SINGLE.hpp"
#include <algorithm>

namespace QComputations {

namespace {
    auto sum_energy = QConfig::instance().h() * QConfig::instance().w();
    auto vacuum_energy = QConfig::instance().h() * QConfig::instance().w() / 4;
    //double sum_energy = 0;
    //double vacuum_energy = 0;
    constexpr int QUDITS_COUNT = 2;

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
class EmulationState : public Basis_State {
    public:
        EmulationState(size_t states_count): Basis_State(QUDITS_COUNT), g_map_(states_count), alpha_map_(states_count) {
            this->set_max_val(1000, 0);
        }
        EmulationState(size_t states_count, const std::string& state_str): Basis_State(state_str), g_map_(states_count), alpha_map_(states_count) {
            this->set_max_val(1000, 0);
        }

        // Сделать шаг по i-тому g 
        EmulationState emulation_step(int i) const {
            EmulationState res(*this);
            int target = g_map_[this->get_qudit(0)][i].first;
            res.set_qudit(target, 0);
            //res.ph_ = (ph_ + 1) % 2;
            res.ph_ = 1;
            res.set_qudit(0, 1);
            res.group_ = target;

            return res;
        }

        // Перевести текущее состояние в терминальное (если возможно)
        EmulationState alpha_step(int i) const {
            EmulationState res(*this);

            int target = alpha_map_[this->get_qudit(0)][i].first;
            res.set_qudit(target, 0);
            res.ph_ = (ph_ + 1) % 2;
            res.ph_ = 1;
            res.set_qudit(0, 1);
            
            if (!res.is_terminal()) {
                res.group_ = target;
            }

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
            return std::find(f_.begin(), f_.end(), this->get_qudit(0)) != f_.end();
            //return (this->get_qudit(0) >= g_map_.size());
        }

        void set_energy(double energy) { state_energy_ = std::vector<double>(g_map_.size(), energy); }
        void set_energy(const std::vector<double>& energy) {state_energy_ = energy;}
        double get_energy() const { return state_energy_[this->get_qudit(0)]; }
        //double get_energy() const { return state_energy_; }

        // Получить фотон
        int get_photon() const { return ph_; }
        // Получить значение matter
        int get_matter() const { return (this->get_qudit(0) == 0 ? 1 : this->get_qudit(1)); }

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

        void set_alpha(int i, int j, double alpha) {
            for (auto& p: g_map_[i]) {
                if (p.first == j) {
                    p.second = alpha;

                    return;
                }
            }

            alpha_map_[i].emplace_back(std::make_pair(j, alpha));
        }

        // Число возможных переходов
        size_t get_neighbours_count() const { return g_map_[this->get_qudit(0)].size(); }
        size_t get_alpha_count() const { return alpha_map_[this->get_qudit(0)].size(); }
        // Получить g_i
        COMPLEX get_g(int i) const { return g_map_[this->get_qudit(0)][i].second; }
        double get_alpha(int i) const { return alpha_map_[this->get_qudit(0)][i].second; }
        int get_alpha_target(int i) const { return alpha_map_[this->get_qudit(0)][i].first; }

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
            //return this->get_qudit(0) < other.get_qudit(0);
        }
    private:
        std::vector<std::vector<std::pair<int, COMPLEX>>> g_map_;
        std::vector<std::vector<std::pair<int, double>>> alpha_map_;
        std::vector<int> f_;
        int ph_ = 0;
        int group_ = 0;
        std::vector<double> state_energy_;
};

// Реализация всевозможных переходов
State<EmulationState> exc_relax_matter(const EmulationState& st) {
    State<EmulationState> res;
    for (size_t i = 0; i < st.get_neighbours_count(); i++) {
        res += State<EmulationState>(st.emulation_step(i)) * st.get_g(i);
    }

    res += State<EmulationState>(st) * st.get_energy();

    return res;
}

// Переход в терминальное состояние
State<EmulationState> term_dec(const EmulationState& st) {
    State<EmulationState> res(st, 1);
    for (size_t i = 0; i < st.get_alpha_count(); i++) {
        res += State<EmulationState>(st.alpha_step(i)) * st.get_alpha(i);
        res[st] -= st.get_alpha(i);
    }

    return res;
}

State<EmulationState> true_term(const EmulationState& st) {
    State<EmulationState> res;
    for (size_t i = 0; i < st.get_alpha_count(); i++) {
        res += State<EmulationState>(st.alpha_step(i));
    }

    return res;
}


void simple_term_alpha(State<EmulationState>& state, int i, int j, double alpha) {
    int index_source = 0;
    int index_target = 0;
    for (auto st: state.state_components()) {
        if (st->get_qudit(0) == i) {
            break;
        }

        index_source++;
    }

    for (auto st: state.state_components()) {
        if (st->get_qudit(0) == j) {
            break;
        }

        index_target++;
    }

    COMPLEX source_ampl = state[index_source];
    COMPLEX target_ampl = state[index_target];


    if (!is_zero(std::abs(source_ampl))) {
        /*
        COMPLEX new_source_ampl = std::sqrt(1 - alpha) * source_ampl;
        COMPLEX new_target_ampl = std::sqrt(target_ampl * target_ampl + alpha * source_ampl * source_ampl);

        state[index_source] = new_source_ampl;
        state[index_target] = new_target_ampl;
        */

        COMPLEX new_source_ampl = (1 - alpha) * source_ampl;
        COMPLEX new_target_ampl = target_ampl + alpha * source_ampl;

        double coef = std::sqrt((std::abs(source_ampl) * std::abs(source_ampl) + std::abs(target_ampl) * std::abs(target_ampl)) / (std::abs(new_source_ampl) * std::abs(new_source_ampl) + std::abs(new_target_ampl) * std::abs(new_target_ampl)));

        state[index_source] = coef * new_source_ampl;
        state[index_target] = coef * new_target_ampl;
    }
}

void simple_term(State<EmulationState>& state) {
    int index = 0;
    for (auto& st: state.state_components()) {
        for (int i = 0; i < st->get_alpha_count(); i++) {
            simple_term_alpha(state, st->get_qudit(0), st->get_alpha_target(i), st->get_alpha(i));
        }
        index++;
    }
}

}