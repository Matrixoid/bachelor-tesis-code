#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <random>
#include <time.h> // Потом удалить!

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
                    array_pointer_P[lambda][s][i] = false;
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
                if (phase % 2) {
                    int up = Cl[2 * b];
                    for (int upp = 0; upp < 2; upp++) {
                        Pl[2 * b + upp] = (Pl1[2 * 2 * b + up ^ upp] * Pl1[2 * (2 * b + 1) + upp]) / 2;
                    }
                }
                else {
                    for (int up = 0; up < 2; up++) {
                        Pl[2 * b + up] = 0;
                        for (int upp = 0; upp < 2; upp++) {
                            Pl[2 * b + up] += (Pl1[2 * 2 * b + up ^ upp] * Pl1[2 * (2 * b + 1) + upp]) / 2;
                        }
                    }
                }
                sig = std::max(sig, Pl[2 * b]);
                sig = std::max(sig, Pl[2 * b + 1]);
            }
        }
        //sig = 1;
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
                Cl1[2 * 2 * b + psi % 2] = Cl[2 * b] ^ Cl[2 * b + 1];
                Cl1[2 * (2 * b + 1) + psi % 2] = Cl[2 * b + 1];
            }
        }
        if (psi % 2 == 1) {
            recursively_update_C(lambda - 1, psi);
        }
    }

    void continue_paths_unfrozen_bit(int phase) {
        std::vector<std::pair<double, int>> tmp;
        double* prob_forks = new double[2 * list_size];
        for (int i = 0; i < 2 * list_size; i++) {
            prob_forks[i] = 0;
        }
        int i = 0;
        for (int l = 0; l < list_size; l++) {
            if (!active_path[l]) {
                prob_forks[2 * l] = -1;
                prob_forks[2 * l + 1] = -1;
            }
            else {
                double* Pm = get_array_pointer_P(m, l);
                prob_forks[2 * l] = Pm[0];
                prob_forks[2 * l + 1] = Pm[1];
                i++;
                tmp.emplace_back(prob_forks[2 * l], l + 1);
                tmp.emplace_back(prob_forks[2 * l + 1], -(l + 1));

            }
        }
        delete[] prob_forks;
        int p = std::min(2 * i, list_size);
        std::sort(tmp.rbegin(), tmp.rend());

        bool* continue_forks = new bool[2 * list_size];
        for (int i = 0; i < 2 * list_size; i++) {
            continue_forks[i] = false;
        }
        for (int j = 0; j < p; j++) {
            if (tmp[j].second > 0) {
                continue_forks[2 * (tmp[j].second - 1)] = true;
            }
            else {
                continue_forks[2 * (-tmp[j].second - 1) + 1] = true;
            }
        }

        for (int l = 0; l < list_size; l++) {
            if (!active_path[l]) {
                continue;
            }
            if (!continue_forks[2 * l] && !continue_forks[2 * l + 1]) {
                kill_path(l);
            }
        }

        for (int l = 0; l < list_size; l++) {
            if (!continue_forks[2 * l] && !continue_forks[2 * l + 1]) {
                continue;
            }
            bool* Cm = get_array_pointer_C(m, l);
            if (continue_forks[2 * l] && continue_forks[2 * l + 1]) {
                Cm[phase % 2] = false;

                int lp = clone_path(l);
                bool* Cmp = get_array_pointer_C(m, lp);
                Cmp[phase % 2] = true;
            }
            else {
                if (continue_forks[2 * l]) {
                    Cm[phase % 2] = false;
                }
                else {
                    Cm[phase % 2] = true;
                }
            }
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
        delete[] frozen_bits_values;
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

    PolarCodes(int m, int k, int d, int n, int s, int p, int list_size, int** frozen_bits_restriction) : m(m),
        k(k),
        d(d),
        n(n),
        s(s),
        p(p),
        list_size(list_size),
        frozen_bits_restrictions(frozen_bits_restriction),
        frozen_bits_values_index(0)
    {
        F = new bool[n];
        for (int i = 0; i < n; i++) {
            if (frozen_bits_restriction[i][0] != 0) {
                F[i] = true;
            }
            else {
                F[i] = false;
            }
        }
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
                for (int l = 0; l < list_size; l++) {
                    if (!active_path[l]) {
                        continue;
                    }
                    bool* Cm = get_array_pointer_C(m, l);
                    Cm[phase % 2] = (frozen_bits_values[frozen_bits_values_index] == 0) ? false : true;
                }
                frozen_bits_values_index++;
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

    ~PolarCodes() {
        delete[] F;
        //delete[] erasure_probabilities;
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

/*
0 1 2 4
0 1 1 0 1 0 0 1
0 0 0 0 1 1 1 1
0.151286
0.00890789
*/

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

void decode_command_execute(command com, PolarCodes& pc, std::ofstream& output_file) {
    double* L = new double[pc.n];
    for (int i = 0; i < pc.n; i++) {
        L[i] = com.vec[i];
    }
    int* dec = pc.decoding(L);
    for (int i = 0; i < pc.n; i++) {
        output_file << dec[i] << " ";
    }
    output_file << "\n";
    delete[] L;
    delete[] dec;
}

void simulate_command_execute(double snrb, int num_of_simulations, int max_errors, PolarCodes& pc, std::ostream& output_file) {
    std::random_device rd{};
    std::mt19937 gen{ rd() };
    auto gen_b = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());

    double snr = pow(10, -snrb / 10);
    snr = snr * pc.n / pc.k;
    double sigma = sqrt(1.0 / 2.0 * snr);
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
        int* dec = pc.decoding(L);
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
        int m, k, d, n, s, p;
        if (input_file.is_open())
            input_file >> m >> k >> d >> n >> s >> p;

        int r = (int)(k * 0.1);

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

        for (int lst = 32; lst <= 64; lst += 16) {
            if (lst == 48)
                continue;
            std::cout << "List size = " << lst << "\n";
            PolarCodes pc = PolarCodes(m, k, d, n, s, p, lst, frozen_bits_restrictions);

            int* unfrozen_bits = pc.get_unfrozen_bits();

            for (double snrb = 0.5; snrb <= 2.5; snrb += 0.25) {
                clock_t start, end;
                start = clock();
                std::cout << "snrb = " << snrb << ": ";
                int num_of_simulations = 100000;
                int max_errors = 100;

                simulate_command_execute(snrb, num_of_simulations, max_errors, pc, std::cout);
                end = clock();
                std::cout << (double)(end - start) / ((double)CLOCKS_PER_SEC) << "\n";
            }
            std::cout << "\n\n";
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
            output_file << frozen_bits[i] << " ";
        }
        output_file << "\n";
        delete[] frozen_bits;

        for (int com_num = 0; com_num < vc.size(); com_num++) {
            if (vc[com_num].com == command_type::Encode) {
                encode_command_execute(vc[com_num], pc, output_file);
            }
            else if (vc[com_num].com == command_type::Decode) {
                decode_command_execute(vc[com_num], pc, output_file);
            }
            else if (vc[com_num].com == command_type::Simulate) {
                simulate_command_execute(vc[com_num].vec[0], vc[com_num].vec[1], vc[com_num].vec[2], pc, output_file);
            }
        }
    }
}