// Born-subtracted determinant of the uniform radial lattice. All matrices
// act in the real rotating basis. Neighboring lattice blocks are -I.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

struct Matrix {
    std::array<double, 9> a{};
    double &operator()(int i, int j) { return a[3 * i + j]; }
    double operator()(int i, int j) const { return a[3 * i + j]; }
};

Matrix operator+(Matrix a, const Matrix &b) {
    for (int i = 0; i < 9; ++i) a.a[i] += b.a[i];
    return a;
}

Matrix operator-(Matrix a, const Matrix &b) {
    for (int i = 0; i < 9; ++i) a.a[i] -= b.a[i];
    return a;
}

Matrix operator*(const Matrix &a, const Matrix &b) {
    Matrix c;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            c(i, j) = a(i, 0) * b(0, j) + a(i, 1) * b(1, j) + a(i, 2) * b(2, j);
    return c;
}

double trace(const Matrix &a) { return a(0, 0) + a(1, 1) + a(2, 2); }

double trace_square(const Matrix &a) {
    return a(0, 0) * a(0, 0) + a(1, 1) * a(1, 1) + a(2, 2) * a(2, 2) +
           2 * (a(0, 1) * a(1, 0) + a(0, 2) * a(2, 0) + a(1, 2) * a(2, 1));
}

double determinant(const Matrix &m) {
    return m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) +
           m(0, 1) * (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) +
           m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
}

Matrix inverse(const Matrix &m) {
    Matrix a;
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3, k = (i + 2) % 3;
        for (int n = 0; n < 3; ++n) {
            const int p = (n + 1) % 3, q = (n + 2) % 3;
            a(n, i) = m(j, p) * m(k, q) - m(j, q) * m(k, p);
        }
    }
    const double invdet = 1 / (m(0, 0) * a(0, 0) + m(0, 1) * a(1, 0) + m(0, 2) * a(2, 0));
    for (double &value : a.a) value *= invdet;
    return a;
}

struct Site {
    double inverse_index_squared, reporter;
    Matrix scaled_potential;
};

struct Probe {
    int ell;
    double zeta;
};

struct Point {
    double coupled, reporter, minimum_pivot_ratio;
};

Point remainder(const std::vector<Site> &sites, double omega, double h, Probe probe) {
    const double zeta = probe.zeta, h2 = h * h;
    const double d_x = (zeta * zeta + 1 - omega * omega) * h2;
    const double d_phi = (zeta * zeta + 2) * h2;
    const double d_d = (zeta * zeta + 1.01) * h2;
    const double gyro = 2 * zeta * omega * h2;
    const double ell_factor = probe.ell * (probe.ell + 1);
    Matrix vacuum_inverse, full_inverse, delta, born1, born2;
    double scalar_vacuum = 0, scalar_full = 0, scalar_delta = 0, scalar_b1 = 0, scalar_b2 = 0;
    Point result{0, 0, 1};
    for (const Site &site : sites) {
        const double diag = 2 + ell_factor * site.inverse_index_squared;
        Matrix vacuum;
        vacuum(0, 0) = vacuum(1, 1) = diag + d_x - vacuum_inverse(0, 0);
        vacuum(2, 2) = diag + d_phi - vacuum_inverse(2, 2);
        vacuum(0, 1) = gyro - vacuum_inverse(0, 1);
        vacuum(1, 0) = -vacuum(0, 1);

        delta = site.scaled_potential + vacuum_inverse * delta * full_inverse;
        const Matrix first_product = vacuum_inverse * born1;
        born2 = (vacuum_inverse * born2 - first_product * first_product) * vacuum_inverse;
        born1 = site.scaled_potential + first_product * vacuum_inverse;
        full_inverse = inverse(vacuum + delta);
        const double inverse_x = 1 / (vacuum(0, 0) * vacuum(0, 0) + vacuum(0, 1) * vacuum(0, 1));
        vacuum_inverse(0, 0) = vacuum_inverse(1, 1) = vacuum(0, 0) * inverse_x;
        vacuum_inverse(0, 1) = -vacuum(0, 1) * inverse_x;
        vacuum_inverse(1, 0) = -vacuum_inverse(0, 1);
        vacuum_inverse(2, 2) = 1 / vacuum(2, 2);
        const Matrix a = vacuum_inverse * delta, b = vacuum_inverse * born1;
        const double tr = trace(a);
        const double ratio_change = tr + 0.5 * (tr * tr - trace_square(a)) + determinant(a);
        result.coupled += std::log1p(ratio_change) - trace(b) -
                          trace(vacuum_inverse * born2) + 0.5 * trace_square(b);

        const double scalar_pivot = diag + d_d - scalar_vacuum;
        scalar_delta = site.reporter + scalar_vacuum * scalar_delta * scalar_full;
        scalar_b2 = scalar_vacuum * scalar_vacuum * (scalar_b2 - scalar_vacuum * scalar_b1 * scalar_b1);
        scalar_b1 = site.reporter + scalar_vacuum * scalar_vacuum * scalar_b1;
        scalar_full = 1 / (scalar_pivot + scalar_delta);
        scalar_vacuum = 1 / scalar_pivot;
        const double a_d = scalar_vacuum * scalar_delta, b_d = scalar_vacuum * scalar_b1;
        result.reporter += std::log1p(a_d) - b_d - scalar_vacuum * scalar_b2 + 0.5 * b_d * b_d;
        result.minimum_pivot_ratio = std::min({result.minimum_pivot_ratio, 1 + ratio_change, 1 + a_d});
    }
    return result;
}

int main(int argc, char **argv) {
    if (argc != 2) return 2;
    std::ifstream input(argv[1]);
    int n, count;
    double h, omega;
    if (!(input >> n >> count >> h >> omega) || n <= 0 || count <= 0 || !(h > 0)) return 2;
    std::vector<Probe> probes(count);
    for (auto &probe : probes) {
        if (!(input >> probe.ell >> probe.zeta) || probe.ell < 2 || probe.zeta <= 0) return 2;
    }
    std::vector<Site> sites(n);
    const double h2 = h * h;
    for (int i = 0; i < n; ++i) {
        double f, s;
        if (!(input >> f >> s)) return 2;
        Site &site = sites[i];
        site.inverse_index_squared = 1.0 / ((i + 1.0) * (i + 1.0));
        const double f2 = f * f, s2m1 = (s - 1) * (s + 1);
        site.scaled_potential(0, 0) = (s2m1 + .3 * f2) * h2;
        site.scaled_potential(1, 1) = (s2m1 + .1 * f2) * h2;
        site.scaled_potential(2, 2) = (3 * s2m1 + f2) * h2;
        site.scaled_potential(0, 2) = site.scaled_potential(2, 0) = 2 * f * s * h2;
        site.reporter = .01 * s2m1 * h2;
    }
    std::vector<Point> output(count);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < count; ++i) output[i] = remainder(sites, omega, h, probes[i]);

    std::cout << std::setprecision(17);
    for (const Point &point : output) {
        std::cout << point.coupled << ' ' << point.reporter << ' ' << point.minimum_pivot_ratio << '\n';
        if (!std::isfinite(point.coupled) || !std::isfinite(point.reporter) || !(point.minimum_pivot_ratio > 0))
            return 1;
    }
    return 0;
}
