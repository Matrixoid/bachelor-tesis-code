#include <iostream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <random>
#include <functional>

int n, p, d;
std::vector<std::vector<int>> field;
std::unordered_map<int, int> from_number_to_alpha_table;
std::unordered_map<int, int> from_alpha_to_number_table;
std::vector<std::vector<int>> multiply_table;
std::vector<std::set<int>> cyclotomic_classes;
std::vector<int> generating_poly;

std::vector<int> roots;

std::vector<int> mul_poly(std::vector<int>& a, std::vector<int>& b) {
    std::vector<int> r(a.size() + b.size() - 1);
    for (int i = b.size() - 1; i >= 0; i--) {
        for (int j = a.size() - 1; j >= 0; j--) {
            r[i + j] = r[i + j] ^ (b[i] & a[j]);
        }
    }
    return r;
}

void xor_rows(std::vector<int>& a, const std::vector<int>& b) {
    for (int i = 0; a.size() == b.size() && i < a.size(); i++) {
        a[i] ^= b[i];
    }
}

void xor_rows(std::vector<int> a, std::vector<int>& b, int n) {
    for (int i = 0; i < n; i++)
        b[i] ^= a[i];
}

int compute_power(int num) {
    int res = 0;
    while (num != 1) {
        res++;
        num /= 2;
    }
    return res;
}

void from_int_to_binary(int p, std::vector<int>& poly) {
    while (p != 0) {
        poly.push_back(p % 2);
        p /= 2;
    }
    std::reverse(poly.begin(), poly.end());
}

int from_binary_to_int(const std::vector<int>& poly) {
    int p = 1;
    int res = 0;
    for (int i = poly.size() - 1; i >= 0; i--) {
        res += poly[i] * p;
        p *= 2;
    }
    return res;
}

void next_vec(std::vector<int>& cur_vec) {
    for (int i = 0; i < cur_vec.size() - 1; i++) {
        cur_vec[i] = cur_vec[i + 1];
    }
    cur_vec[cur_vec.size() - 1] = 0;
}

void build_field(int m) {
    std::vector<int> poly;
    from_int_to_binary(p, poly);
    poly.erase(poly.begin());

    std::unordered_set<int> ss;

    std::vector<int> cur_vec;
    for (int i = 0; i < m; i++) {
        cur_vec.push_back(0);
    }
    ss.insert(0);

    field.push_back(cur_vec);
    cur_vec[m - 1] = 1;
    field.push_back(cur_vec);
    ss.insert(1);

    do {
        bool f = cur_vec[0] == 1;
        next_vec(cur_vec);
        if (f)
            xor_rows(cur_vec, poly);
        if (ss.find(from_binary_to_int(cur_vec)) == ss.end()) {
            field.push_back(cur_vec);
            ss.insert(from_binary_to_int(cur_vec));
        }
        else {
            break;
        }
    } while (true);

}

void build_cyclotomic_classes() {
    std::set<int> cyclotomic_class;
    cyclotomic_class.insert(0);
    cyclotomic_classes.push_back(cyclotomic_class);
    cyclotomic_class.clear();

    std::set<int> c;
    for (int i = 1; i < n - 1; i++) {
        c.insert(i);
    }

    while (!c.empty()) {
        int new_c = *c.begin();
        if (new_c >= d)
            break;

        while (cyclotomic_class.find(new_c) == cyclotomic_class.end()) {
            c.erase(new_c);
            cyclotomic_class.insert(new_c);
            new_c *= 2;
            new_c %= n;
        }

        cyclotomic_classes.push_back(cyclotomic_class);
        cyclotomic_class.clear();
    }
}

void build_from_number_to_alpha_table() {
    from_number_to_alpha_table[0] = -1;
    for (int i = 1; i < field.size(); i++) {
        from_number_to_alpha_table[from_binary_to_int(field[i])] = i - 1;
    }
}

void build_from_alpha_to_number_table() {
    for (int i = 1; i < field.size(); i++) {
        from_alpha_to_number_table[i - 1] = from_binary_to_int(field[i]);
    }
}

void build_multiply_table() {
    multiply_table = std::vector<std::vector<int>>(n + 1, std::vector<int>(n + 1));

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            multiply_table[i][j] = from_binary_to_int(field[(from_number_to_alpha_table[i] + from_number_to_alpha_table[j]) % n + 1]);
        }
    }
}

void build_generating_poly() {
    generating_poly.push_back(1);
    for (int i = 1; i < cyclotomic_classes.size(); i++) {
        std::set<int> current_cyclotomic_class = cyclotomic_classes[i];
        std::vector<int> polynom;
        polynom.push_back(0);

        for (int num : current_cyclotomic_class) {
            std::vector<int> tmp_poly;
            tmp_poly.push_back(-1);
            for (int j = 0; j < polynom.size(); j++)
                tmp_poly.push_back(polynom[j]);
            for (int l = 0; l < polynom.size(); l++) {
                if (polynom[l] >= 0) {
                    polynom[l] += num;
                    polynom[l] %= n;
                }
                if (polynom[l] >= 0 && tmp_poly[l] >= 0) {
                    polynom[l] = from_number_to_alpha_table[from_binary_to_int(field[polynom[l] + 1]) ^ from_binary_to_int(field[tmp_poly[l] + 1])];
                }
                else if (tmp_poly[l] >= 0)
                    polynom[l] = tmp_poly[l];
            }
            polynom.push_back(tmp_poly[tmp_poly.size() - 1]);
        }
        for (int k = 0; k < polynom.size(); k++) {
            if (polynom[k] == 0)
                polynom[k] = 1;
            else
                polynom[k] = 0;
        }

        generating_poly = mul_poly(generating_poly, polynom);
    }
}

std::vector<int> encoding(const std::vector<int>& code_word) {
    std::vector<int> res(generating_poly.size() - 1);
    res.insert(res.end(), code_word.begin(), code_word.end());
    std::vector<int> tmp = res;
    for (int i = res.size() - 1; i >= generating_poly.size() - 1; i--) {
        if (tmp[i] == 1) {
            for (int j = 0; j < generating_poly.size(); j++) {
                tmp[i - j] ^= generating_poly[generating_poly.size() - 1 - j];
            }
        }
    }

    for (int i = 0; i < generating_poly.size() - 1; i++)
        res[i] = tmp[i];

    return res;
}

void decoding(std::vector<int>& dec) {
    std::vector<int> s;
    int count_zero_s = 0;
    for (int j = 1; j <= d - 1; j++) {
        int t = 0;
        for (int i = 0; i < dec.size(); i++) {
            if (dec[i] == 0)
                t ^= from_alpha_to_number_table[(i * j) % n];
        }
        if (t == 0) count_zero_s++;
        s.push_back(t);
    }

    if (count_zero_s == d - 1) return;

    int L = 0;
    int m = 0;
    std::vector<int> lambda = { 1 };
    std::vector<int> b = { 1 };
    for (int rr = 1; rr <= d - 1; rr++) {
        int delta = 0;
        for (int j = 0; j <= std::min(L, (int)lambda.size()); j++) {
            delta ^= multiply_table[lambda[j]][s[rr - j - 1]];
        }

        if (delta != 0) {
            std::vector<int> temp = lambda;
            temp.resize(std::max(temp.size(), rr - m + b.size()), 0);
            for (int i = 0; i < b.size(); i++) temp[rr - m + i] ^= multiply_table[delta][b[i]];

            if (2 * L <= rr - 1) {
                b.clear();
                int dd = from_alpha_to_number_table[(n - from_number_to_alpha_table[delta]) % n];
                for (int i : lambda) {
                    b.push_back(multiply_table[dd][i]);
                }
                L = rr - L;
                m = rr;
            }

            lambda = temp;
        }
    }

    if (L != lambda.size() - 1) return;

    std::vector<int> err;
    for (int i = 0; i < n; i++) {
        int res = lambda[0];
        for (int j = 1; j < lambda.size(); j++) {
            res ^= multiply_table[from_alpha_to_number_table[(j * i) % n]][lambda[j]];
        }
        if (res == 0) {
            err.push_back((n - i) % n);
            if (err.size() == L) break;
        }
    }

    for (int e : err) {
        if (dec[e] == 1)
            dec[e] = 0;
        else
            dec[e] = 1;
    }
}

std::vector<std::vector<int>> multiply(std::vector<std::vector<int>> arikan_kernel, std::vector<std::vector<int>> kernel) {
    std::vector<std::vector<int>> res = std::vector<std::vector<int>>(kernel.size() * arikan_kernel.size(), std::vector<int>(kernel[0].size() * arikan_kernel[0].size()));

    for (int i = 0; i < arikan_kernel.size(); i++) {
        for (int j = 0; j < arikan_kernel[0].size(); j++) {
            for (int ii = 0; ii < kernel.size(); ii++) {
                for (int jj = 0; jj < kernel[0].size(); jj++) {
                    res[i * kernel.size() + ii][j * kernel[0].size() + jj] = kernel[ii][jj] * arikan_kernel[i][j];
                }
            }
        }
    }

    return res;
}

std::vector<std::vector<int>> Kronecker_multiply(std::vector<std::vector<int>> kernel, int m) {
    std::vector<std::vector<int>> res = kernel;
    std::vector<std::vector<int>> arikan_kernel = { {1, 0}, {1, 1} };

    while (m > 1) {
        res = multiply(arikan_kernel, res);
        m--;
    }

    return res;
}

void subtitude_to_g() {
    int pp = field[0].size();
    int cnt = 0;
    for (int i = 1; i < field.size(); i++) {
        std::vector<int> res = std::vector<int>(pp);
        for (int j = 0; j < generating_poly.size(); j++) {
            if (generating_poly[j] != 0) {
                int degree = ((i - 1) * j) % (field.size() - 1);
                for (int l = 0; l < pp; l++) {
                    res[l] = (res[l] + field[degree + 1][l]) % 2;
                }
            }
        }

        bool flag = true;
        for (int j = 0; j < res.size(); j++) {
            if (res[j] != 0) {
                flag = false;
                break;
            }
        }

        if (flag) {
            cnt++;
            std::string s = "";
            for (int j = 0; j < field[i].size(); j++) {
                if (field[i][j] == 0) continue;
                s += "x^" + std::to_string(field[i].size() - j - 1) + " + ";
            }
            if (s.size() > 2) {
                s.pop_back();
                s.pop_back();
            }
            std::cout << "i = " << i - 1 << ": " << s << std::endl;
            roots.push_back(i - 1);
        }
    }
    std::cout << cnt << "\n\n";
}

void generate_bit_vector(std::vector<int>& bit_vector) {
    if (bit_vector[bit_vector.size() - 1] == 0) {
        bit_vector[bit_vector.size() - 1] = 1;
    }
    else {
        int l = bit_vector.size() - 1;
        while (bit_vector[l] == 1) {
            bit_vector[l] = 0;
            l--;
        }
        bit_vector[l] = 1;
    }
}

bool all_zero(std::vector<int> t) {
    bool flag = true;
    for (int i = 0; i < t.size(); i++) {
        if (t[i] == 1) {
            flag = false;
            break;
        }
    }
    return flag;
}

void minimal_span_form(std::vector<std::vector<int>>& G) {
    //Гаус в одну сторону
    int j = 0;
    int i = 0;
    std::vector<int> first_one(G.size());
    for (int k = 0; k < G.size(); k++) {
        for (int l = 0; l < G[0].size(); l++) {
            if (G[k][l] == 1) {
                first_one[k] = l;
                break;
            }
            if (l == G[0].size() - 1 && G[k][l] == 0) {
                first_one[k] = -1;
            }
        }
    }

    while (i <= G.size() && j <= G[0].size()) {
        if (i < G.size() && first_one[i] == -1) {
            i++;
            continue;
        }
        int tmp = -1;
        for (int k = i; k < G.size(); k++) {
            if (first_one[k] == j) {
                tmp = k;
                if (tmp != i) {
                    std::swap(G[tmp], G[i]);
                    std::swap(first_one[tmp], first_one[i]);
                    tmp = i;
                }
                break;
            }
        }
        if (tmp == -1) {
            j++;
            continue;
        }

        for (int k = i + 1; k < G.size(); k++) {
            if (first_one[k] == j) {
                xor_rows(G[tmp], G[k], G[0].size());
                for (int l = 0; l < G[0].size(); l++) {
                    if (G[k][l] == 1) {
                        first_one[k] = l;
                        break;
                    }
                    if (l == G[0].size() - 1 && G[k][l] == 0) {
                        first_one[k] = -1;
                    }
                }
            }
        }
        j++;
        if (tmp != -1)
            i++;
    }

    //Гаус в другую сторону
    i = G.size() - 1;
    j = G[0].size() - 1;
    std::vector<int> last_one(i + 1);
    for (int k = 0; k < i + 1; k++) {
        for (int l = j; l >= 0; l--) {
            if (G[k][l] == 1) {
                last_one[k] = l;
                break;
            }
        }
    }

    while (j > 0) {

        int tmp = -1;
        for (int k = i; k >= 0; k--) {
            if (last_one[k] == j) {
                tmp = k;
                break;
            }
        }
        if (tmp == -1) {
            j--;
            continue;
        }

        for (int k = i; k >= 0; k--) {
            if (k != tmp && last_one[k] == last_one[tmp]) {
                xor_rows(G[tmp], G[k], G[0].size());
                for (int p = G[0].size() - 1; p >= 0; p--) {
                    if (G[k][p] == 1) {
                        last_one[k] = p;
                        break;
                    }
                }
            }
        }

        j--;
    }

    int a = 15;
}

std::vector<std::vector<int>> get_check_matrix_16() {
    std::vector<std::vector<int>> tmp = std::vector<std::vector<int>>(14, std::vector<int>(2));
    tmp[0] = { 1, 3 };
    std::vector<int> vec_start = std::vector<int>(4);
    vec_start[vec_start.size() - 2] = 1;
    for (int i = 1; i < 14; i++) {
        generate_bit_vector(vec_start);
        tmp[i][0] = from_number_to_alpha_table[from_binary_to_int(vec_start)];
        for (int j = 1; j < 2; j++) {
            tmp[i][j] = (tmp[i][0] * tmp[0][j]) % 15;
        }
    }
    std::vector<std::vector<int>> res = std::vector<std::vector<int>>(16, std::vector<int>(16 - 7));
    res[0][0] = 1;
    for (int i = 0; i < res.size(); i++) {
        res[i][0] = 1;
    }

    int pos = 1;
    for (int j = 0; j < 2; j++) {
        std::vector<int> v = field[1];
        std::reverse(v.begin(), v.end());
        for (int vi = 0; vi < v.size(); vi++) {
            res[1][pos] = v[vi];
            pos++;
        }
    }
    for (int i = 2; i < res.size(); i++) {
        pos = 1;
        for (int j = 0; j < 2; j++) {
            std::vector<int> v = field[tmp[i - 2][j] + 1];
            std::reverse(v.begin(), v.end());
            for (int vi = 0; vi < v.size(); vi++) {
                res[i][pos] = v[vi];
                pos++;
            }
        }
    }
    return res;
}

std::vector<std::vector<int>> get_check_matrix() {
    std::vector<std::vector<int>> tmp = std::vector<std::vector<int>>(1022, std::vector<int>(13));
    tmp[0] = { 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25 };
    std::vector<int> vec_start = std::vector<int>(10);
    vec_start[vec_start.size() - 2] = 1;
    for (int i = 1; i < 1022; i++) {
        generate_bit_vector(vec_start);
        tmp[i][0] = from_number_to_alpha_table[from_binary_to_int(vec_start)];
        for (int j = 1; j < 13; j++) {
            tmp[i][j] = (tmp[i][0] / tmp[0][0] * tmp[0][j]) % 1023;
        }
    }
    std::vector<std::vector<int>> res = std::vector<std::vector<int>>(1024, std::vector<int>(1024 - 893));
    res[0][0] = 1;
    for (int i = 0; i < res.size(); i++) {
        res[i][0] = 1;
    }

    int pos = 1;
    for (int i = 0; i < 13; i++) {
        std::vector<int> v = field[1];
        std::reverse(v.begin(), v.end());
        for (int vi = 0; vi < v.size(); vi++) {
            res[1][pos] = v[vi];
            pos++;
        }
    }
    for (int i = 2; i < res.size(); i++) {
        pos = 1;
        for (int j = 0; j < 13; j++) {
            std::vector<int> v = field[tmp[i - 2][j] + 1];
            std::reverse(v.begin(), v.end());
            for (int vi = 0; vi < v.size(); vi++) {
                res[i][pos] = v[vi];
                pos++;
            }
        }
    }
    return res;
}

std::vector<std::vector<int>> matrix_multiply(std::vector<std::vector<int>> kernel, std::vector<std::vector<int>> H) {
    std::vector<std::vector<int>> res = std::vector<std::vector<int>>(kernel.size(), std::vector<int>(H[0].size()));
    for (int i = 0; i < kernel.size(); i++) {
        for (int j = 0; j < H[0].size(); j++) {
            for (int k = 0; k < kernel[0].size(); k++) {
                res[i][j] += kernel[i][k] * H[k][j];
                res[i][j] %= 2;
            }
        }
    }

    return res;
}

std::vector<std::vector<int>> matrix_transpose(std::vector<std::vector<int>> a) {
    std::vector<std::vector<int>> res = std::vector<std::vector<int>>(a[0].size(), std::vector<int>(a.size()));
    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j < a[i].size(); j++) {
            res[j][i] = a[i][j];
        }
    }
    return res;
}

void print_field() {
    for (int i = 0; i < field.size(); i++) {
        std::string s = "";
        for (int j = 0; j < field[i].size(); j++) {
            if (field[i][j] == 0) continue;
            s += "x^" + std::to_string(field[i].size() - j - 1) + " + ";
        }
        if (s.size() > 2) {
            s.pop_back();
            s.pop_back();
        }
        std::cout << s << std::endl;
    }
}

void print_matrix_into_file(std::ostream& output_file, std::vector<std::vector<int>> matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            output_file << matrix[i][j] << " ";
        }
        output_file << "\n";
    }
    output_file << "\n";
}

void code_to_file(std::string output_file_name, std::vector<std::vector<int>> matrix, int m, int k, int d, int n) {
    std::ofstream output_file;
    output_file.open(output_file_name);

    output_file << m << " " << k << " " << d << " " << n << " " << 0 << " " << 0 << "\n\n";

    for (int i = 0; i < matrix.size(); i++) {
        int sum = 0;
        std::string s = "";
        for (int j = 0; j < matrix[i].size(); j++) {
            if (matrix[i][j] == 1) {
                sum++;
                s += std::to_string(j) + " ";
            }
        }
        if (sum != 0) {
            output_file << std::to_string(sum) << " " << s << std::endl;
        }
    }
}

void try_min(std::vector<std::vector<int>>& matrix, int index) {
    for (int i = 0; i < index; i++) {
        int first_sum = 0;
        for (int j = 0; j < matrix[i].size(); j++) {
            if (matrix[index][j] == 1) first_sum++;
        }
        xor_rows(matrix[index], matrix[i]);
        int second_sum = 0;
        for (int j = 0; j < matrix[i].size(); j++) {
            if (matrix[index][j] == 1) second_sum++;
        }
        if (first_sum < second_sum) xor_rows(matrix[index], matrix[i]);
    }
}

std::vector<int> get_bit_reversal_sequence(std::vector<int> bit_sequence) {
    if (bit_sequence.size() == 1) return bit_sequence;

    std::vector<int> res = std::vector<int>(bit_sequence.size());
    std::vector<int> tmp1;
    for (int i = 0; i < bit_sequence.size() / 2; i++) {
        tmp1.push_back(bit_sequence[i]);
    }
    std::vector<int> tmp2;
    for (int i = bit_sequence.size() / 2; i < bit_sequence.size(); i++) {
        tmp2.push_back(bit_sequence[i]);
    }

    tmp1 = get_bit_reversal_sequence(tmp1);
    tmp2 = get_bit_reversal_sequence(tmp2);
    for (int i = 0; i < bit_sequence.size() / 2; i++) {
        res[2 * i] = tmp1[i];
        res[2 * i + 1] = tmp2[i];
    }
    return res;
}

std::vector<int> get_bit_reversal_sequence(int n) {
    std::vector<int> res;
    for (int i = 0; i < n; i++) {
        res.push_back(i);
    }
    res = get_bit_reversal_sequence(res);
    return res;
}

std::vector<std::vector<int>> get_reversal_matrix(std::vector<std::vector<int>> matrix, std::vector<int> bit_reversal_sequence) {
    std::vector<std::vector<int>> res = std::vector<std::vector<int>>(matrix.size(), std::vector<int>(matrix[0].size()));
    for (int j = 0; j < bit_reversal_sequence.size(); j++) {
        for (int i = 0; i < matrix.size(); i++) {
            res[i][j] = matrix[i][bit_reversal_sequence[j]];
        }
    }

    return res;
}

std::vector<std::vector<int>> get_check_matrix(int n) {
    if (n == 16) {
        return get_check_matrix_16();
    }
    return get_check_matrix();
}

std::vector<std::vector<int>> sort_by_end(std::vector<std::vector<int>> matrix) {
    int i = matrix.size() - 1;
    int j = matrix[0].size() - 1;
    std::vector<std::pair<int, int>> last_one(i + 1);
    for (int k = 0; k < i + 1; k++) {
        for (int l = j; l >= 0; l--) {
            if (matrix[k][l] == 1) {
                last_one[k] = std::make_pair(k, l);
                break;
            }
        }
    }
    std::sort(last_one.begin(), last_one.end(), [](std::pair<int, int> a, std::pair<int, int> b) { return a.second < b.second || (a.second == b.second && a.first < b.first);});

    std::vector<std::vector<int>> res = std::vector<std::vector<int>>(matrix.size());
    for (int i = 0; i < last_one.size(); i++) {
        res[i] = matrix[last_one[i].first];
    }
    return res;
}

std::vector<std::vector<int>> get_constraint_matrix(std::vector<std::vector<int>> matrix, int frozen_bits_count) {
    std::vector<std::vector<int>> res = std::vector<std::vector<int>>(frozen_bits_count, std::vector<int>(matrix[0].size()));
    std::ifstream frozen_bits_file;
    frozen_bits_file.open("frozen_bits.txt");
    std::vector<int> frozen_bits;
    for (int i = 0; i < frozen_bits_count; i++) {
        int frozen_bit;
        frozen_bits_file >> frozen_bit;
        frozen_bits.push_back(frozen_bit);
    }

    int i = matrix.size() - 1;
    int j = matrix[0].size() - 1;
    std::vector<int> last_one(i + 1);
    for (int k = 0; k < i + 1; k++) {
        for (int l = j; l >= 0; l--) {
            if (matrix[k][l] == 1) {
                last_one[k] = l;
                break;
            }
        }
    }

    std::vector<int> f_bits;
    int pos = 0;
    while (f_bits.size() < frozen_bits_count - matrix.size()) {
        bool flag = false;
        for (int s = 0; s < last_one.size(); s++) {
            if (last_one[s] == frozen_bits[pos]) flag = true;
        }

        if (!flag) {
            f_bits.push_back(frozen_bits[pos]);
        }
        pos++;
    }

    for (int s = 0; s < f_bits.size(); s++) {
        std::vector<int> tmp = std::vector<int>(1024);
        tmp[f_bits[s]] = 1;
        matrix.push_back(tmp);
    }

    return sort_by_end(matrix);
}

int main()
{
    std::string input_file_name = "input.txt";

    std::ifstream input_file;
    input_file.open(input_file_name);
    if (input_file.is_open())
        input_file >> n >> p >> d;

    int m = compute_power(n + 1);
    build_field(m);
    build_from_number_to_alpha_table();
    build_multiply_table();
    build_cyclotomic_classes();
    build_generating_poly();
    build_from_alpha_to_number_table();
    subtitude_to_g();
    std::vector<std::vector<int>> kernel = Kronecker_multiply({ {1, 0}, {1, 1} }, m);
    std::vector<int> brs = get_bit_reversal_sequence(n + 1);
    std::vector<std::vector<int>> kernel4 = get_reversal_matrix(kernel, brs);
    std::vector<std::vector<int>> matr = get_check_matrix(n + 1);
    std::vector<std::vector<int>> aaa = matrix_multiply(kernel4, matr);

    aaa = matrix_transpose(aaa);
    minimal_span_form(aaa);

    if (n + 1 == 1024) {
        aaa = get_constraint_matrix(aaa, 512);
    }

    for (int i = 1; i < aaa.size(); i++) {
        try_min(aaa, i);
    }

    code_to_file("1024_512_28.spec", aaa, m, 512, 28, 1024);
}