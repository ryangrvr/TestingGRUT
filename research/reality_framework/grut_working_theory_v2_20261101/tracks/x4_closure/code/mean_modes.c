/* Homogeneous Majorana helicity modes on a supplied H55 mean waveform.
 * The chiral-rotated Hamiltonian has axis b=(M,0,h*q/a+theta_dot/2).
 * This code neither evolves the scalar source nor computes a stress tensor.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

enum { NSPECIES = 4, NODE_FIELDS = 3, OBS_FIELDS = 6 };

typedef struct {
    double x, y, z;
} Bloch;

static void read_exact(FILE *stream, void *data, size_t size, size_t count) {
    if (fread(data, size, count, stream) != count) {
        fputs("Incomplete mean-mode input\n", stderr);
        exit(1);
    }
}

static inline Bloch rotate(Bloch v, double bx, double bz, double dt) {
    const double norm = sqrt(bx * bx + bz * bz);
    const double ux = bx / norm, uz = bz / norm;
    const double angle = 2 * norm * dt;
    const double cosine = cos(angle), sine = sin(angle);
    const double along = (1 - cosine) * (ux * v.x + uz * v.z);
    return (Bloch){cosine * v.x - sine * uz * v.y + along * ux,
                   cosine * v.y + sine * (uz * v.x - ux * v.z),
                   cosine * v.z + sine * ux * v.y + along * uz};
}

static inline Bloch vacuum(const double *obs, double hq, double mass, int parent, int order) {
    double k = hq * obs[0];
    double bz = k + obs[2 + parent];
    double energy = sqrt(bz * bz + mass * mass);
    double x = -mass / energy, z = -bz / energy, y = 0;
    if (order == 1) {
        double bz_dot = -obs[1] * k + obs[4 + parent];
        y = -mass * bz_dot / (2 * energy * energy * energy);
    }
    double norm = sqrt(1 + y * y);
    return (Bloch){x / norm, y / norm, z / norm};
}

static inline double occupation(Bloch state, Bloch reference) {
    double x = state.x - reference.x, y = state.y - reference.y, z = state.z - reference.z;
    return (x * x + y * y + z * z) / 4;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fputs("Usage: mean_modes INPUT OUTPUT\n", stderr);
        return 1;
    }
    FILE *input = fopen(argv[1], "rb");
    if (!input) return 1;
    int32_t shape[3];
    double physical[2 + NSPECIES];
    read_exact(input, shape, sizeof *shape, 3);
    read_exact(input, physical, sizeof *physical, 2 + NSPECIES);
    const int nstep = shape[0], nq = shape[1], nobs = shape[2];
    const double dt = physical[0], qmax = physical[1], dq = qmax / nq;
    if (nstep < 1 || nq < 1 || nobs < 2 || dt <= 0 || qmax <= 0) return 1;
    const int nmodes = 2 * NSPECIES * nq;
    int32_t *obs_steps = malloc(nobs * sizeof *obs_steps);
    double *obs = malloc((size_t)nobs * OBS_FIELDS * sizeof *obs);
    double *nodes = malloc((size_t)nstep * 2 * NODE_FIELDS * sizeof *nodes);
    int ncol = 1 + 2 * (nobs - 1);
    double *results = malloc((size_t)nmodes * ncol * sizeof *results);
    if (!obs_steps || !obs || !nodes || !results) return 1;
    read_exact(input, obs_steps, sizeof *obs_steps, nobs);
    read_exact(input, obs, sizeof *obs, (size_t)nobs * OBS_FIELDS);
    read_exact(input, nodes, sizeof *nodes, (size_t)nstep * 2 * NODE_FIELDS);
    fclose(input);
    if (obs_steps[0] != 0 || obs_steps[nobs-1] != nstep) return 1;

    /* Fourth-order two-exponential commutator-free Magnus method. */
    const double early = (3 + 2 * sqrt(3.0)) / 12;
    const double late = (3 - 2 * sqrt(3.0)) / 12;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < nmodes; i++) {
        int species = i / (2 * nq), parent = species != 0;
        double mass = physical[2 + species], q = dq * (i % nq + 0.5);
        double hq = ((i / nq) % 2 ? 1 : -1) * q;
        double *out = results + (size_t)i * ncol;
        Bloch v = vacuum(obs, hq, mass, parent, 1);
        int snapshot = 1;
        out[0] = 0;
        for (int step = 0; step < nstep; step++) {
            const double *first = nodes + (size_t)step * 2 * NODE_FIELDS;
            const double *second = first + NODE_FIELDS;
            double z1 = hq * first[0] + first[1 + parent];
            double z2 = hq * second[0] + second[1 + parent];
            v = rotate(v, mass / 2, early * z1 + late * z2, dt);
            v = rotate(v, mass / 2, late * z1 + early * z2, dt);
            if (step + 1 != obs_steps[snapshot]) continue;
            double norm = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
            out[0] = fmax(out[0], fabs(norm - 1));
            v = (Bloch){v.x / norm, v.y / norm, v.z / norm};
            const double *at = obs + (size_t)snapshot * OBS_FIELDS;
            for (int order = 0; order < 2; order++)
                out[1 + 2 * (snapshot - 1) + order] = occupation(v, vacuum(at, hq, mass, parent, order));
            snapshot++;
        }
    }
    FILE *output = fopen(argv[2], "w");
    if (!output) return 1;
    for (int i = 0; i < nmodes; i++) {
        fprintf(output, "%d %d %.17g", i / (2 * nq) + 1, (i / nq) % 2 ? 1 : -1, dq * (i % nq + 0.5));
        for (int j = 0; j < ncol; j++) fprintf(output, " %.17g", results[(size_t)i * ncol + j]);
        fputc('\n', output);
    }
    fclose(output);
    free(obs_steps);
    free(obs);
    free(nodes);
    free(results);
    return 0;
}
