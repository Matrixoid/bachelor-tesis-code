#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <random>
#include <time.h>

enum command_type {
    Encode,
    Decode,
    Simulate
};

struct command {

    command() = default;
    command(command_type ct, std::vector<double> vector) : com(ct), vec(vector) {}

    void clear() {
        vec.clear();
    }

    friend std::ostream& operator<<(std::ostream& os, const command& cmd) {
        if (cmd.com == command_type::Encode)
            os << "Encode ";
        else if (cmd.com == command_type::Decode)
            os << "Decode ";
        else if (cmd.com == command_type::Simulate)
            os << "Simulate ";
        for (auto a : cmd.vec) {
            os << a << " ";
        }
        std::cout << std::endl;
        return os;
    }

    friend std::istream& operator>>(std::istream& is, command& cmd) {
        std::string cm;
        std::stringstream ss;
        std::getline(is, cm);
        ss.str(cm);
        std::string ct;
        ss >> ct;
        if (ct == "Encode")
            cmd.com = command_type::Encode;
        else if (ct == "Decode")
            cmd.com = command_type::Decode;
        else if (ct == "Simulate")
            cmd.com = command_type::Simulate;

        while (!ss.eof()) {
            double a;
            ss >> a;
            cmd.vec.push_back(a);
        }
        return is;
    }

    command_type com;
    std::vector<double> vec;
};

struct PolarCodes {
private:

    void print_permutations_to_file() {
        std::ofstream file;
        file.open("permutations.txt");
        for (int i = 0; i < 64; i++) {
            for (int j = 0; j < n; j++) {
                file << permutations[i][j] << " ";
            }
            file << permutations_prob[i] << "\n";
        }
        file.close();
    }

    int* bit_permutation(int* permutation) {
        int* res = new int[n];
        for (int i = 0; i < n; i++) {
            res[i] = i;
        }
        return res;
    }

    double permutation_prob_compute(int* permutation) {
        double prob = 1;
        for (int i = 0; i < n; i++) {
            if (!F[i]) {
                prob *= (1 - erasure_probabilities[permutation[i]]);
            }
        }
        return 1 - prob;
    }

    void swap(int* arr, int a, int b)
    {
        int t = arr[a];
        arr[a] = arr[b];
        arr[b] = t;
    }
    void generate_permutations(int* cur_perm, int k, int max_permutations)
    {
        if (k == m)
        {
            int* bit_perm = bit_permutation(cur_perm);
            double prob = permutation_prob_compute(bit_perm);
            delete[] bit_perm;
            if (permutations_index >= max_permutations) {
                int min_index = 0;
                int min_prob = permutations_prob[0];
                for (int i = 1; i < max_permutations; i++) {
                    if (permutations_prob[i] < min_prob) {
                        min_index = i;
                        min_prob = permutations_prob[i];
                    }
                }
                if (prob > min_prob) {
                    if (permutations[min_index]) {
                        delete[] permutations[min_index];
                    }
                    permutations[min_index] = cur_perm;
                    permutations_prob[min_index] = prob;
                }
            }
            else {
                permutations[permutations_index] = cur_perm;
                permutations_prob[permutations_index] = prob;
            }
        }
        else
        {
            for (int j = k; j < m; j++)
            {
                swap(cur_perm, k, j);
                generate_permutations(cur_perm, k + 1, max_permutations);
                swap(cur_perm, k, j);
            }
        }
    }

    void generate_permutations(int max_permutations) {
        int* arr = new int[m];
        permutations_prob = new double[max_permutations];
        permutations = new int* [max_permutations];
        permutations_index = 0;
        for (int i = 0; i < m; i++) {
            arr[i] = i;
        }

        generate_permutations(arr, 0, max_permutations);
        delete[] arr;
    }

    double* get_permutation(double* L, int* permutation) {
        double* res = new double[n];
        for (int i = 0; i < n; i++) {
            res[i] = L[permutation[i]];
        }

        return res;
    }

    double* erasure_probabilities_compute(int n, double p) {
        double* res = new double[(m + 1) * n];
        int l = n;

        int layer = 0;
        for (int i = 0; i < l; i++) {
            res[layer * n + i] = p;
        }
        layer++;

        l /= 2;
        while (l > 0) {
            for (int i = 0; i < n;) {
                for (int j = i; j < i + l; j++) {
                    res[layer * n + j] = 1 - (1 - res[(layer - 1) * n + j]) * (1 - res[(layer - 1) * n + j]);
                }
                i += l;
                for (int j = i; j < i + l; j++) {
                    res[layer * n + j] = res[(layer - 1) * n + j] * res[(layer - 1) * n + j];
                }
                i += l;
            }
            l /= 2;
            layer++;
        }

        return res;
    }

    bool* find_k_mins(int k, const double* probabilities) {
        std::vector<std::pair<int, double>> tmp;
        bool* res = new bool[n];

        for (int i = 0; i < n; i++) {
            tmp.push_back(std::make_pair(i, probabilities[i]));
            res[i] = true;
        }

        std::sort(tmp.begin(), tmp.end(), [](std::pair<int, double> a, std::pair<int, double> b) { return a.second < b.second; });

        for (int i = 0; i < k; i++) {
            res[tmp[i].first] = false;
        }

        return res;
    }

    int* multiply(int* vec, int matrix[2][2]) {
        int* res = new int[2] { (vec[0] + vec[1]) % 2, vec[1] };
        return res;
    }

    int* encoding(int* vec, int size, int m) {
        int matrix[2][2] = { {1, 0}, {1, 1} };

        int* part1 = new int[size / 2];
        int* part2 = new int[size / 2];
        if (m != 1) {
            int* tmp1 = new int[size / 2];
            int* tmp2 = new int[size / 2];

            for (int i = 0; i < size / 2; i++) {
                tmp1[i] = vec[i];
            }
            for (int i = 0; i < size / 2; i++) {
                tmp1[i] = (tmp1[i] + vec[i + size / 2]) % 2;
                tmp2[i] = vec[i + size / 2];
            }
            int* e1 = encoding(tmp1, size / 2, m - 1);
            int* e2 = encoding(tmp2, size / 2, m - 1);
            for (int i = 0; i < size / 2; i++) {
                part1[i] = e1[i];
            }
            for (int i = 0; i < size / 2; i++) {
                part2[i] = e2[i];
            }
            delete[] e1;
            delete[] e2;
            delete[] tmp1;
            delete[] tmp2;
        }
        else {
            int* tmp1 = new int[size / 2];
            int* tmp2 = new int[size / 2];

            for (int i = 0; i < size / 2; i++) {
                tmp1[i] = vec[i];
            }
            for (int i = 0; i < size / 2; i++) {
                tmp1[i] = tmp1[i] = (tmp1[i] + vec[i + size / 2]) % 2;
                tmp2[i] = vec[i + size / 2];
            }
            part1 = multiply(tmp1, matrix);
            part2 = multiply(tmp2, matrix);
            delete[] tmp1;
            delete[] tmp2;
        }
        int* res = new int[size];
        for (int i = 0; i < size / 2; i++) {
            res[i] = part1[i];
        }
        for (int i = 0; i < size / 2; i++) {
            res[size / 2 + i] = part2[i];
        }
        delete[] part1;
        delete[] part2;
        return res;
    }

    void initialize_data_structures() {
        inactive_path_indices = new int[list_size];
        inactive_path_indices_size = list_size;
        active_path = new bool[list_size];
        array_pointer_P = new double** [m + 1];
        array_pointer_C = new bool** [m + 1];
        path_index_to_array_index = new int* [m + 1];
        inactive_array_indices = new int* [m + 1];
        inactive_array_indices_sizes = new int[m + 1];
        array_references_count = new int* [m + 1];

        for (int lambda = 0; lambda < m + 1; lambda++) {
            inactive_array_indices[lambda] = new int[list_size];
            inactive_array_indices_sizes[lambda] = list_size;
            array_pointer_P[lambda] = new double* [list_size];
            array_pointer_C[lambda] = new bool* [list_size];
            for (int s = 0; s < list_size; s++) {
                array_pointer_P[lambda][s] = new double[2 * (1 << (m - lambda))];
                array_pointer_C[lambda][s] = new bool[2 * (1 << (m - lambda))];
                for (int i = 0; i < 2 * (1 << (m - lambda)); i++) {
                    array_pointer_P[lambda][s][i] = 0;
                    array_pointer_C[lambda][s][i] = false;
                }
                inactive_array_indices[lambda][s] = s;
            }
            path_index_to_array_index[lambda] = new int[list_size];
            array_references_count[lambda] = new int[list_size];
        }
        for (int l = 0; l < list_size; l++) {
            inactive_path_indices[l] = l;
            active_path[l] = false;
        }
    }


    int assign_initial_path() {
        int l = inactive_path_indices[inactive_path_indices_size - 1];
        inactive_path_indices_size--;
        active_path[l] = true;
        for (int lambda = 0; lambda < m + 1; lambda++) {
            int s = inactive_array_indices[lambda][inactive_array_indices_sizes[lambda] - 1];
            inactive_array_indices_sizes[lambda]--;
            path_index_to_array_index[lambda][l] = s;
            array_references_count[lambda][s] = 1;
        }

        return l;
    }

    int clone_path(int l) {
        int l1 = inactive_path_indices[inactive_path_indices_size - 1];
        inactive_path_indices_size--;
        active_path[l1] = true;
        for (int lambda = 0; lambda < m + 1; lambda++) {
            int s = path_index_to_array_index[lambda][l];
            path_index_to_array_index[lambda][l1] = s;
            array_references_count[lambda][s]++;
        }

        return l1;
    }

    void kill_path(int l) {
        active_path[l] = false;
        inactive_path_indices[inactive_path_indices_size] = l;
        inactive_path_indices_size++;
        for (int lambda = 0; lambda < m + 1; lambda++) {
            int s = path_index_to_array_index[lambda][l];
            array_references_count[lambda][s]--;
            if (array_references_count[lambda][s] == 0) {
                inactive_array_indices[lambda][inactive_array_indices_sizes[lambda]] = s;
                inactive_array_indices_sizes[lambda]++;
            }
        }
    }

    double* get_array_pointer_P(int lambda, int l) {
        int s = path_index_to_array_index[lambda][l];
        int s1;
        if (array_references_count[lambda][s] == 1) {
            s1 = s;
        }
        else {
            s1 = inactive_array_indices[lambda][inactive_array_indices_sizes[lambda] - 1];
            inactive_array_indices_sizes[lambda]--;
            for (int i = 0; i < 2 * (1 << (m - lambda)); i++) {
                array_pointer_P[lambda][s1][i] = array_pointer_P[lambda][s][i];
                array_pointer_C[lambda][s1][i] = array_pointer_C[lambda][s][i];
            }
            array_references_count[lambda][s]--;
            array_references_count[lambda][s1] = 1;
            path_index_to_array_index[lambda][l] = s1;
        }

        return array_pointer_P[lambda][s1];
    }

    bool* get_array_pointer_C(int lambda, int l) {
        int s = path_index_to_array_index[lambda][l];
        int s1;
        if (array_references_count[lambda][s] == 1) {
            s1 = s;
        }
        else {
            s1 = inactive_array_indices[lambda][inactive_array_indices_sizes[lambda] - 1];
            inactive_array_indices_sizes[lambda]--;
            for (int i = 0; i < 2 * (1 << (m - lambda)); i++) {
                array_pointer_P[lambda][s1][i] = array_pointer_P[lambda][s][i];
                array_pointer_C[lambda][s1][i] = array_pointer_C[lambda][s][i];
            }
            array_references_count[lambda][s]--;
            array_references_count[lambda][s1] = 1;
            path_index_to_array_index[lambda][l] = s1;
        }

        return array_pointer_C[lambda][s1];
    }

    int compute_index(int p, int b, int lambda) {
        return p + (1 << lambda) * b;
    }

    void recursively_calc_P(int lambda, int phase) {
        if (lambda == 0) {
            return;
        }
        int psi = phase / 2;
        if (phase % 2 == 0) {
            recursively_calc_P(lambda - 1, psi);
        }
        double sig = 0;
        for (int l = 0; l < list_size; l++) {
            if (!active_path[l]) {
                continue;
            }
            double* Pl = get_array_pointer_P(lambda, l);
            double* Pl1 = get_array_pointer_P(lambda - 1, l);
            bool* Cl = get_array_pointer_C(lambda, l);

            for (int b = 0; b < (1 << (m - lambda)); b++) {
                if (phase % 2 == 0) {
                    Pl[2 * b] = (Pl1[2 * (2 * b)] * Pl1[2 * (2 * b + 1)] + Pl1[2 * (2 * b) + 1] * Pl1[2 * (2 * b + 1) + 1]) / 2;
                    sig = std::max(sig, Pl[2 * b]);
                    Pl[2 * b + 1] = (Pl1[2 * (2 * b) + 1] * Pl1[2 * (2 * b + 1)] + Pl1[2 * (2 * b)] * Pl1[2 * (2 * b + 1) + 1]) / 2;
                    sig = std::max(sig, Pl[2 * b + 1]);
                }
                else {
                    int u = (Cl[2 * b]) ? 1 : 0;
                    for (int uu = 0; uu < 2; uu++) {
                        Pl[2 * b + uu] = 0.5 * Pl1[2 * (2 * b) + (u xor uu)] * Pl1[2 * (2 * b + 1) + uu];
                        sig = std::max(sig, Pl[2 * b + uu]);
                    }
                }
            }
        }
        for (int l = 0; l < list_size; l++) {
            if (!active_path[l]) {
                continue;
            }
            double* Pl = get_array_pointer_P(lambda, l);
            for (int b = 0; b < (1 << (m - lambda)); b++) {
                Pl[2 * b] /= sig;
                Pl[2 * b + 1] /= sig;
            }
        }
    }

    void recursively_update_C(int lambda, int phase) {
        int psi = phase / 2;

        for (int l = 0; l < list_size; l++) {
            if (!active_path[l]) {
                continue;
            }
            bool* Cl = get_array_pointer_C(lambda, l);
            bool* Cl1 = get_array_pointer_C(lambda - 1, l);
            for (int b = 0; b < (1 << (m - lambda)); b++) {
                Cl1[2 * 2 * b + psi % 2] = Cl[2 * b] xor Cl[2 * b + 1];
                Cl1[2 * (2 * b + 1) + psi % 2] = Cl[2 * b + 1];
            }
        }
        if (psi % 2 == 1) {
            recursively_update_C(lambda - 1, psi);
        }
    }

    void continue_paths_frozen_bits(int phase) {
        for (int l = 0; l < list_size; l++) {
            if (!active_path[l]) {
                continue;
            }
            bool* Cm = get_array_pointer_C(m, l);
            Cm[phase % 2] = (frozen_bits_values[frozen_bits_values_index] == 0) ? false : true;
        }
        frozen_bits_values_index++;
    }

    void continue_paths_unfrozen_bit(int phase) {

        std::vector<std::pair<double, int>> tmp;
        int i = 0;
        for (int l = 0; l < list_size; l++) {
            if (active_path[l]) {
                double* Pm = get_array_pointer_P(m, l);
                tmp.emplace_back(Pm[0], l + 1);
                tmp.emplace_back(Pm[1], -(l + 1));
                i++;
            }
        }
        int p = std::min(2 * i, list_size);
        std::sort(tmp.rbegin(), tmp.rend());

        bool** continue_forks = new bool* [list_size];
        for (int i = 0; i < list_size; i++) {
            continue_forks[i] = new bool[2] { false, false };
        }
        for (int j = 0; j < p; j++) {
            if (tmp[j].second > 0) {
                continue_forks[(tmp[j].second - 1)][0] = true;
            }
            else {
                continue_forks[(-tmp[j].second - 1)][1] = true;
            }
        }

        for (int l = 0; l < list_size; l++) {
            if (!active_path[l]) {
                continue;
            }
            if (!continue_forks[l][0] && !continue_forks[l][1]) {
                kill_path(l);
            }
        }

        for (int l = 0; l < list_size; l++) {
            if (!continue_forks[l][0] && !continue_forks[l][1]) {
                continue;
            }
            bool* Cm = get_array_pointer_C(m, l);
            if (continue_forks[l][0] && continue_forks[l][1]) {
                Cm[phase % 2] = false;

                int lp = clone_path(l);
                bool* Cmp = get_array_pointer_C(m, lp);
                Cmp[phase % 2] = true;
            }
            else {
                if (continue_forks[l][0]) {
                    Cm[phase % 2] = false;
                }
                else {
                    Cm[phase % 2] = true;
                }
            }
        }
        for (int i = 0; i < list_size; i++) {
            delete[] continue_forks[i];
        }
        delete[] continue_forks;
    }

    int get_reversed(int ind) {
        int newInd = 0;
        for (int i = 0; i < m; i++) {
            newInd = newInd * 2 + ind % 2;
            ind /= 2;
        }
        return newInd;
    }

    double compute_metric1(double x, double y) {
        return (x * y) / std::abs(x * y) * std::min(std::abs(x), std::abs(y));
    }

    double compute_metric2(double x, double y, int u) {
        return (1 - 2 * u) * x + y;
    }

    double compute_metric(int* dec) {
        double res = 0;
        double* y = new double[n];
        for (int i = 0; i < n; i++) {
            y[i] = 1 - 2 * dec[i];
        }
        for (int l = m; l >= 0; l--) {
            for (int i = 0; i < std::pow(2, l - 1); i++) {
                y[i] = compute_metric1(y[i], y[i + (int)std::pow(2, l - 1)]);
            }
        }
        res = std::min(0.0, y[0]);
        for (int i = 0; i < n; i++) {
            y[i] = 1 - 2 * dec[i];
        }
        for (int l = m; l >= 0; l--) {
            for (int i = 0; i < std::pow(2, l - 1); i++) {
                y[i] = compute_metric2(y[i], y[i + (int)std::pow(2, l - 1)], dec[i]);
            }
        }
        res += std::min(0.0, y[0]);
        delete[] y;
        return res;
    }

    int* reverse_permutation(int* dec, int* permutation) {
        int* res = new int[n];
        for (int i = 0; i < n; i++) {
            res[i] = dec[permutation[i]];
        }
        return res;
    }

    void delete_all_arrays() {
        delete[] inactive_path_indices;
        delete[] active_path;
        for (int i = 0; i < m + 1; i++) {
            delete[] path_index_to_array_index[i];
            delete[] array_references_count[i];
            delete[] inactive_array_indices[i];
            for (int j = 0; j < list_size; j++) {
                delete[] array_pointer_P[i][j];
                delete[] array_pointer_C[i][j];
            }
            delete[] array_pointer_P[i];
            delete[] array_pointer_C[i];
        }
        delete[] path_index_to_array_index;
        delete[] array_references_count;
        delete[] inactive_array_indices_sizes;
        delete[] inactive_array_indices;
        delete[] array_pointer_P;
        delete[] array_pointer_C;
    }

    double* erasure_probabilities;
    int** frozen_bits_restrictions;
    bool* F;
    int* frozen_bits_values;
    int frozen_bits_values_index;
    int* inactive_path_indices;
    int inactive_path_indices_size;
    bool* active_path;
    double*** array_pointer_P;
    bool*** array_pointer_C;
    int** path_index_to_array_index;
    int** inactive_array_indices;
    int* inactive_array_indices_sizes;
    int** array_references_count;
    int** permutations;
    double* permutations_prob;
    int permutations_index;
public:
    PolarCodes() = default;

    PolarCodes(int n, int k, double prob, int list_size) : n(n),
        k(k),
        prob(prob),
        list_size(list_size),
        m((int)log2(n)),
        frozen_bits_values_index(0)
    {
        erasure_probabilities = erasure_probabilities_compute(n, prob);
        F = find_k_mins(k, erasure_probabilities + m * n);
        frozen_bits_restrictions = new int* [n];
        for (int i = 0; i < n; i++) {
            if (F[i]) {
                frozen_bits_restrictions[i] = new int[2];
                frozen_bits_restrictions[i][0] = 1;
                frozen_bits_restrictions[i][1] = i;
            }
            else {
                frozen_bits_restrictions[i] = new int[1];
                frozen_bits_restrictions[i][0] = 0;
            }
        }

        int a = 0;
    }

    PolarCodes(int m, int k, int d, int n, int s, int p, int list_size, int** frozen_bits_restriction, double prob) : m(m),
        k(k),
        prob(prob),
        d(d),
        n(n),
        s(s),
        p(p),
        list_size(list_size),
        frozen_bits_restrictions(frozen_bits_restriction),
        frozen_bits_values_index(0)
    {
        erasure_probabilities = erasure_probabilities_compute(n, prob);
        F = new bool[n];
        for (int i = 0; i < n; i++) {
            if (frozen_bits_restriction[i][0] != 0) {
                F[i] = true;
            }
            else {
                F[i] = false;
            }
        }
        generate_permutations(64);
    }

    int* get_frozen_bits() {
        int* frozen_bits = new int[n - k];
        int fb = 0;
        for (int i = 0; i < n; i++) {
            if (F[i]) {
                frozen_bits[fb] = i;
                fb++;
            }
        }
        return frozen_bits;
    }

    int* get_unfrozen_bits() {
        int* unfrozen_bits = new int[k];
        int fb = 0;
        for (int i = 0; i < n; i++) {
            if (!F[i]) {
                unfrozen_bits[fb] = i;
                fb++;
            }
        }
        return unfrozen_bits;
    }

    int* encoding(int* word) {
        frozen_bits_values = new int[n - k];
        int* to_encode = new int[n];
        int j = 0;
        int z = 0;
        for (int i = 0; i < n; i++) {
            if (j < n - k && F[j + z]) {
                to_encode[j + z] = 0;
                for (int ii = 2; ii < frozen_bits_restrictions[j][0]; ii++) {
                    to_encode[j + z] += to_encode[frozen_bits_restrictions[j][ii]];
                }
                to_encode[j + z] %= 2;
                frozen_bits_values[j] = to_encode[j + z];
                j++;
            }
            else {
                to_encode[j + z] = word[z];
                z++;
            }
        }
        int* res = new int[n];
        int* e = encoding(to_encode, n, m - 1);
        for (int i = 0; i < n; i++) {
            res[i] = e[i];
        }
        delete[] e;
        delete[] to_encode;

        return res;
    }

    int* decoding(double* L) {
        initialize_data_structures();
        int l = assign_initial_path();
        double* P0 = get_array_pointer_P(0, l);

        for (int b = 0; b < n; b++) {
            int idX = get_reversed(b);
            P0[2 * b] = exp(L[idX]) / (exp(L[idX]) + 1);
            P0[2 * b + 1] = 1 - P0[2 * b];
        }
        for (int phase = 0; phase < n; phase++) {
            recursively_calc_P(m, phase);
            if (F[phase]) {
                continue_paths_frozen_bits(phase);
            }
            else {
                continue_paths_unfrozen_bit(phase);
            }
            if (phase % 2 == 1) {
                recursively_update_C(m, phase);
            }
        }
        frozen_bits_values_index = 0;

        int l1 = 0;
        double pp = 0;
        for (int l = 0; l < list_size; l++) {
            if (!active_path[l]) {
                continue;
            }

            bool* Cm = get_array_pointer_C(m, l);
            double* Pm = get_array_pointer_P(m, l);

            if (pp < Pm[Cm[1]]) {
                l1 = l;
                pp = Pm[Cm[1]];
            }
        }
        bool* C0 = get_array_pointer_C(0, l1);

        int* res = new int[n];
        for (int b = 0; b < n; b++) {
            res[b] = (C0[2 * get_reversed(b)]) ? 1 : 0;
        }
        delete_all_arrays();
        return res;
    }

    int* decoding(double* L, int permutation_count) {
        int metric = -10000000;
        int* res = new int[n];
        for (int i = 0; i < permutation_count; i++) {
            double* L_tmp = get_permutation(L, bit_permutation(permutations[i]));
            int* dec = decoding(L_tmp);
            delete[] L_tmp;
            double cur_metric = compute_metric(dec);
            if (cur_metric > metric) {
                metric = cur_metric;
                delete[] res;
                res = reverse_permutation(dec, bit_permutation(permutations[i]));
                delete[] dec;
            }
        }
        return res;
    }

    ~PolarCodes() {
        delete[] F;
        delete[] erasure_probabilities;
        delete[] frozen_bits_values;
    }

    int n;
    int k;
    double prob;
    int m;
    int d;
    int s;
    int p;
    int list_size;
};

void print_res(int* arr, int n) {
    for (int j = 0; j < n; j++) {
        std::cout << arr[j] << " ";
    }
}


void encode_command_execute(command com, PolarCodes& pc, std::ofstream& output_file) {
    int* to_encode = new int[pc.k];
    for (int i = 0; i < pc.k; i++) {
        to_encode[i] = com.vec[i];
    }
    int* enc = pc.encoding(to_encode);
    for (int i = 0; i < pc.n; i++) {
        output_file << enc[i] << " ";
    }
    output_file << "\n";
    delete[] enc;
    delete[] to_encode;
}

void decode_command_execute(command com, int permutation_count, PolarCodes& pc, std::ofstream& output_file) {
    double* L = new double[pc.n];
    for (int i = 0; i < pc.n; i++) {
        L[i] = com.vec[i];
    }
    int* dec = pc.decoding(L, permutation_count);
    for (int i = 0; i < pc.n; i++) {
        output_file << dec[i] << " ";
    }
    output_file << "\n";
    delete[] L;
    delete[] dec;
}

void simulate_command_execute(double snrb, int num_of_simulations, int max_errors, int permutation_count, PolarCodes& pc, std::ostream& output_file) {
    std::random_device rd{};
    std::mt19937 gen{ rd() };
    auto gen_b = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());

    double snr = pow(10, -snrb / 10);
    snr = snr * pc.n / (pc.k + 50);
    double sigma = sqrt(0.5 * snr);
    std::normal_distribution<> nd{ 0, sigma };

    int errs = 0;
    int iters = 0;
    for (int i = 0; i < num_of_simulations; i++) {
        int* b = new int[pc.k];
        for (int j = 0; j < pc.k; j++)
            b[j] = gen_b();

        int* enc = pc.encoding(b);
        delete[] b;

        double* L = new double[pc.n];
        for (int i = 0; i < pc.n; i++) {
            L[i] = 1 - 2 * enc[i] + nd(gen);
        }
        int* dec = pc.decoding(L, permutation_count);
        delete[] L;

        for (int i = 0; i < pc.n; i++) {
            if (enc[i] != dec[i]) {
                errs++;
                break;
            }
        }
        delete[] enc;
        delete[] dec;
        iters = i + 1;
        if (errs >= max_errors)
            break;
    }
    output_file << (double)errs / iters << "\n";
}

int main()
{
    bool dynamic_froze = true;
    std::random_device rd{};
    std::mt19937 gen{ rd() };
    auto gen_b = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());

    if (dynamic_froze) {
        std::string input_file_name = "1024_512_24_5.spec";

        std::ifstream input_file;
        input_file.open(input_file_name);
        int m = 0, k = 0, d = 0, n = 0, s = 0, p = 0;
        if (input_file.is_open())
            input_file >> m >> k >> d >> n >> s >> p;

        int** frozen_bits_restrictions = new int* [n];
        for (int i = 0; i < n; i++) {
            frozen_bits_restrictions[i] = new int[1] {0};
        }
        for (int i = 0; i < n - k; i++) {
            int restrict_count;
            input_file >> restrict_count;
            int* restrict = new int[restrict_count];
            for (int j = 0; j < restrict_count; j++) {
                input_file >> restrict[j];
            }
            delete[] frozen_bits_restrictions[restrict[restrict_count - 1]];
            frozen_bits_restrictions[restrict[restrict_count - 1]] = new int[restrict_count + 1];
            frozen_bits_restrictions[restrict[restrict_count - 1]][0] = restrict_count;
            frozen_bits_restrictions[restrict[restrict_count - 1]][1] = restrict[restrict_count - 1];
            for (int j = 0; j < restrict_count - 1; j++) {
                frozen_bits_restrictions[restrict[restrict_count - 1]][j + 2] = restrict[j];
            }
            delete[] restrict;
        }

        for (int lst = 16; lst <= 64; lst += 16) {
            for (int prm_cnt = 16; prm_cnt <= 64; prm_cnt += 16) {
                if (lst == 48 || prm_cnt == 48)
                    continue;
                std::cout << "List size = " << lst << "\n";
                PolarCodes pc = PolarCodes(m, k, d, n, s, p, lst, frozen_bits_restrictions, 0.5);

                for (double snrb = 0.5; snrb <= 2.5; snrb += 0.25) {
                    clock_t start, end;
                    start = clock();
                    std::cout << "snrb = " << snrb << ": ";
                    int num_of_simulations = 100000;
                    int max_errors = 100;

                    simulate_command_execute(snrb, num_of_simulations, max_errors, prm_cnt, pc, std::cout);
                    end = clock();
                    std::cout << (double)(end - start) / ((double)CLOCKS_PER_SEC) << "\n";
                }
                std::cout << "\n\n";
            }
        }
        for (int i = 0; i < n; i++) {
            delete[] frozen_bits_restrictions[i];
        }
        delete[] frozen_bits_restrictions;
    }
    else {
        std::string input_file_name = "input.txt";
        std::string output_file_name = "output.txt";
        std::ifstream input_file;
        input_file.open(input_file_name);
        int n, k;
        double p;
        int list_size;
        if (input_file.is_open())
            input_file >> n >> k >> p >> list_size;

        std::string cd;
        std::getline(input_file, cd);
        std::vector<command> vc;
        command cmd;
        while (input_file >> cmd) {
            vc.push_back(cmd);
            cmd.clear();
        }
        input_file.close();

        PolarCodes pc = PolarCodes(n, k, p, list_size);

        std::ofstream output_file;
        output_file.open(output_file_name);
        int* frozen_bits = pc.get_frozen_bits();
        for (int i = 0; i < k; i++) {
            output_file << frozen_bits[i] << "\n";
        }
        output_file << "\n";
        delete[] frozen_bits;

        for (int com_num = 0; com_num < vc.size(); com_num++) {
            if (vc[com_num].com == command_type::Encode) {
                encode_command_execute(vc[com_num], pc, output_file);
            }
            else if (vc[com_num].com == command_type::Decode) {
                decode_command_execute(vc[com_num], 16, pc, output_file);
            }
            else if (vc[com_num].com == command_type::Simulate) {
                simulate_command_execute(vc[com_num].vec[0], vc[com_num].vec[1], vc[com_num].vec[2], 16, pc, output_file);
            }
        }
    }
}