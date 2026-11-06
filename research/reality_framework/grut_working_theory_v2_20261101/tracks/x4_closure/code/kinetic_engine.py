"""Extend the inherited dephased engine in scratch with the full energy ledger."""
import subprocess
from pathlib import Path

from trace import TRACKS, replace_once

ENGINE = TRACKS / 'e2e_kinetic/code/ek_modes.c'

SPECTATOR = r'''
/* Exact evolution for constant H over this expansion step, with
 * ds/dt=-(3H+Gamma)s and dDR/dt=-4H DR+Gamma s. */
static void spectator_step(double dt, double aold, double anew) {
    double decay = gamma_s * dt;
    double x = log(anew / aold) - decay;
    double exprel = fabs(x) < 1e-7 ? 1 + x / 2 + x * x / 6 : expm1(x) / x;
    double r = aold / anew, r3 = r * r * r;
    rho_dr = (rho_dr + rho_s * decay * exprel) * r3 * r;
    rho_s *= r3 * exp(-decay);
}

static double pressure(void) {
    double result = (rho_h + rho_th + rho_dr) / 3.0;
    for (int j = 0; j < NSP; j++) result += S[j].pressure;
    return result;
}
'''


def source(dt_factor):
    generated = ENGINE.read_text()
    patches = [
        ('    double E, Nn;', '    double E, Nn, pressure;'),
        ('static FILE *frec;', '''
static FILE *frec, *fledger;
static double rho_s, rho_dr, gamma_s;
static double exchange_error, budget_error, occ_min = 1.0, occ_max = 0.0;'''),
        ('    fclose(fp);\n    double rho = rhoP', '''
    rho_s = rd(fp); rho_dr = rd(fp); gamma_s = rd(fp);
    if (coherent || rho_s < 0 || rho_dr < 0 || gamma_s < 0) die("X4 dephased ledger input");
    fclose(fp);
    double rho = rho_s + rho_dr + rhoP'''),
        ('    double E = 0, N = 0, ia3', '    double E = 0, N = 0, P3 = 0, ia3'),
        ('            E += s->w[i] * om * nn; N += s->w[i] * nn;',
         '''
            E += s->w[i] * om * nn; N += s->w[i] * nn;
            P3 += s->w[i] * k * k / om * nn;
            if (nn < occ_min) occ_min = nn;
            if (nn > occ_max) occ_max = nn;'''),
        ('    s->E = E * ia3; s->Nn = N * ia3;',
         '    s->E = E * ia3; s->Nn = N * ia3; s->pressure = P3 * ia3 / 3.0;'),
        ('    double r = rhoP[0] + rhoP[1] + rho_h + rho_th;',
         '    double r = rhoP[0] + rhoP[1] + rho_h + rho_th + rho_s + rho_dr;'),
        ('/* ------------------------------------------------------------------ record */', SPECTATOR),
        ('    fflush(frec);', r'''
    fprintf(fledger, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
            t, a, rho_s, rho_dr, rho_h + rho_th, rho, exchange_error,
            budget_error, occ_min, occ_max);
    fflush(fledger);
    fflush(frec);'''),
        ('    double dt = 1e-4, trec = t;', '''
    double dt = 1e-4, trec = t;
    const double initial_energy = rho * a * a * a;
    double pressure_work = 0.0;'''),
        ('        double lim = 0.02 / H;', f'        double lim = {dt_factor:.17g} * 0.02 / H;'),
        ('0.01 * rhoP[P] / Pc[P]', f'{dt_factor:.17g} * 0.01 * rhoP[P] / Pc[P]'),
        ('0.2 * tk)', f'{dt_factor:.17g} * 0.2 * tk)'),
        ('''        /* expansion */
        double aold = a, anew = noexp ? a : a * exp(H * dt);''', '''
        /* All transfers occur before expansion and at the same volume. */
        sm_update(dt, a, a);
        for (int j = 0; j < NSP; j++) species_energy(&S[j]);
        double exchanged = rho_tot();
        exchange_error = fmax(exchange_error, fabs(exchanged / rho - 1.0));
        double pressure_before = pressure();

        double aold = a, anew = noexp ? a : a * exp(H * dt);'''),
        ('''        sm_update(dt, aold, anew);
        for (int j = 0; j < NSP; j++) species_energy(&S[j]);
        rho = rho_tot(); H = noexp ? H : sqrt(rho / (3 * mP * mP));''', '''
        sm_update(0, aold, anew);
        spectator_step(dt, aold, anew);
        for (int j = 0; j < NSP; j++) species_energy(&S[j]);
        rho = rho_tot(); H = noexp ? H : sqrt(rho / (3 * mP * mP));
        pressure_work += 0.5 * (pressure_before + pressure()) * (a*a*a - aold*aold*aold);
        budget_error = fmax(budget_error, fabs((rho*a*a*a + pressure_work) / initial_energy - 1.0));'''),
        ('''        if ((rhoP[0] + rhoP[1]) < 1e-9 * rho && EN < 1e-9 * rho && t > 10) break;''',
         '''        if ((rhoP[0] + rhoP[1]) < 1e-9 * rho && EN < 1e-9 * rho &&
            rho_s < 1e-9 * rho && t > 10) break;'''),
        ('    if (argc < 3) die("usage");', r'''
    if (argc == 8 && !strcmp(argv[1], "spectator-test")) {
        rho_s = atof(argv[2]); rho_dr = atof(argv[3]); gamma_s = atof(argv[4]);
        spectator_step(atof(argv[5]), atof(argv[6]), atof(argv[7]));
        printf("%.17g %.17g\n", rho_s, rho_dr);
        return 0;
    }
    if (argc < 3) die("usage");'''),
        ('''    frec = fopen(fn, "w");
    if (coherent) run_coherent(); else run_dephased();
    fclose(frec);''', r'''
    frec = fopen(fn, "w");
    snprintf(fn, sizeof fn, "%s_ledger.txt", argv[3]);
    fledger = fopen(fn, "w");
    if (!frec || !fledger) die("open records");
    fprintf(fledger, "# t a rho_s rho_dr rho_sm rho exchange_error budget_error f_min f_max\n");
    run_dephased();
    fclose(fledger);
    fclose(frec);'''),
    ]
    for old, new in patches:
        generated = replace_once(generated, old, new)
    return generated


def compile_engine(work: Path, dt_factor):
    path = work / f'engine_dt{dt_factor:g}.c'
    path.write_text(source(dt_factor))
    binary = path.with_suffix('')
    subprocess.run(['gcc', '-O3', '-march=native', '-o', str(binary), str(path), '-lm'], check=True)
    return binary
