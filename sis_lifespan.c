/* ============================================================================
 *  sis_lifespan.c
 *
 *  Dinamica SIS sobre redes del configuration model con grado power-law,
 *  simulada con el algoritmo de Gillespie y analizada con el "lifespan method"
 *  (Mata, Boguna, Castellano, Pastor-Satorras, Phys. Rev. E 91, 052117 (2015)).
 *
 *  - delta (tasa de recuperacion) = 1  ->  el control es lambda = lambda/delta.
 *  - Parametro de orden:       Pend(lambda,N)  = P(realizacion endemica).
 *  - Susceptibilidad:          <tau(lambda,N)> = vida media de las NO endemicas.
 *  - Endemica := la cobertura (nodos distintos infectados) alcanza cth*N.
 *
 *  Estructuras de datos clave (todas con insercion/borrado O(1)):
 *    inf_list[]  : nodos infectados (los primeros n_inf).
 *    act[]       : aristas activas, identificadas por su id de arista DIRIGIDA
 *                  (infectado -> susceptible). act_pos[] guarda la posicion de
 *                  cada arista dirigida dentro de act[] para poder borrar en O(1)
 *                  (esto es el "additional data structure" del enunciado).
 *    rev[]       : para cada arista dirigida, el id de su arista inversa.
 *
 *  Compilar:  gcc -O3 -march=native -o sis_lifespan sis_lifespan.c -lm
 * ==========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

/* ----------------------------- RNG: xoshiro256** --------------------------- */
static uint64_t RNG[4];

static inline uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static void rng_seed(uint64_t seed) {
    uint64_t s = seed ? seed : 0x123456789ABCDEFULL;
    for (int i = 0; i < 4; i++) RNG[i] = splitmix64(&s);
}
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rng_next(void) {
    const uint64_t result = rotl(RNG[1] * 5, 7) * 9;
    const uint64_t t = RNG[1] << 17;
    RNG[2] ^= RNG[0]; RNG[3] ^= RNG[1]; RNG[1] ^= RNG[2]; RNG[0] ^= RNG[3];
    RNG[2] ^= t; RNG[3] = rotl(RNG[3], 45);
    return result;
}
/* doble en [0,1) */
static inline double rng_unit(void) {
    return (double)(rng_next() >> 11) * (1.0 / 9007199254740992.0);
}
/* Multiplicacion 64x64 -> mitad alta, portable (sin __uint128_t).
 * Funciona en MSVC, MinGW y GCC/Linux. */
static inline uint64_t mul_hi64(uint64_t x, uint64_t y) {
    uint64_t xlo = (uint32_t)x, xhi = x >> 32;
    uint64_t ylo = (uint32_t)y, yhi = y >> 32;
    uint64_t ll = xlo * ylo;
    uint64_t lh = xlo * yhi;
    uint64_t hl = xhi * ylo;
    uint64_t hh = xhi * yhi;
    uint64_t mid = (ll >> 32) + (uint32_t)lh + (uint32_t)hl;
    return hh + (lh >> 32) + (hl >> 32) + (mid >> 32);
}

/* entero uniforme sin sesgo en [0,n)  (metodo de Lemire, portable) */
static inline uint64_t rng_below(uint64_t n) {
    uint64_t x  = rng_next();
    uint64_t lo = x * n;          /* mitad baja (desbordamiento intencionado) */
    uint64_t hi = mul_hi64(x, n);
    if (lo < n) {
        uint64_t t = (-n) % n;
        while (lo < t) { x = rng_next(); lo = x * n; hi = mul_hi64(x, n); }
    }
    return hi;
}

/* ----------------------------- Estado global ------------------------------- */
static int      N;            /* numero de nodos                               */
static long     n_edges;      /* aristas no dirigidas                          */
static long     n_dir;        /* aristas dirigidas = 2*n_edges                 */

static long    *adj_start;    /* CSR: offsets, tam N+1                         */
static int     *adj;          /* CSR: vecinos, tam n_dir                       */
static long    *rev;          /* id de la arista dirigida inversa, tam n_dir   */

static uint8_t *state;        /* 0=S, 1=I, tam N                               */
static int     *inf_list;     /* nodos infectados, tam N                       */
static long     n_inf;

static long    *act;          /* aristas activas (id dirigido), tam n_edges    */
static long    *act_pos;      /* posicion en act[] de cada arista dirigida     */
static long     n_act;

static uint8_t *covered;      /* nodo ya infectado alguna vez en la realizacion*/
static int     *touched;      /* lista de nodos cubiertos (para resetear O(cov))*/
static long     n_touch;
static long     cov;          /* cobertura actual                              */

static int     *seed_nodes;   /* nodos candidatos a semilla (grado == kmin)    */
static long     n_seed;

/* --------------------- helpers de listas dinamicas O(1) -------------------- */
static inline void add_active(long de) {
    act[n_act] = de; act_pos[de] = n_act; n_act++;
}
static inline void remove_active(long de) {
    long p = act_pos[de];
    long last = act[--n_act];
    act[p] = last; act_pos[last] = p; act_pos[de] = -1;
}
static inline void add_inf(int v) {
    inf_list[n_inf++] = v; state[v] = 1;
}

/* ---- actualizacion de aristas activas cuando un nodo cambia de estado ----- */
/* v acaba de INFECTARSE (S->I). state[v] ya = 1. */
static inline void update_SI(int v) {
    long a0 = adj_start[v], a1 = adj_start[v + 1];
    for (long de = a0; de < a1; de++) {
        int x = adj[de];
        if (state[x] == 0) add_active(de);        /* v->x pasa a activa        */
        else               remove_active(rev[de]);/* x->v estaba activa: fuera */
    }
}
/* v acaba de RECUPERARSE (I->S). state[v] ya = 0. */
static inline void update_IS(int v) {
    long a0 = adj_start[v], a1 = adj_start[v + 1];
    for (long de = a0; de < a1; de++) {
        int x = adj[de];
        if (state[x] == 0) remove_active(de);      /* v->x estaba activa: fuera */
        else               add_active(rev[de]);    /* x->v pasa a activa        */
    }
}

static int cmp_u64(const void *p, const void *q) {
    uint64_t x = *(const uint64_t *)p, y = *(const uint64_t *)q;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* =========================== Configuration model ========================== */
/* Muestrea grados ~ k^-gamma en [kmin,kmax], empareja stubs, elimina
 * autoenlaces y multienlaces (grafo simple) y construye la CSR + rev[].       */
static void build_network(int Nn, double gamma, int kmin, int kmax,
                          uint64_t net_seed,
                          double *mean_k, double *mean_k2) {
    rng_seed(net_seed);
    N = Nn;

    /* --- CDF discreta de P(k) ~ k^-gamma --- */
    int K = kmax - kmin + 1;
    double *cdf = (double *)malloc((size_t)K * sizeof(double));
    double acc = 0.0;
    for (int i = 0; i < K; i++) { acc += pow((double)(kmin + i), -gamma); cdf[i] = acc; }
    for (int i = 0; i < K; i++) cdf[i] /= acc;

    int *deg = (int *)malloc((size_t)N * sizeof(int));
    long long stub_sum = 0;
    for (int v = 0; v < N; v++) {
        double u = rng_unit();
        /* busqueda binaria del primer cdf>=u */
        int lo = 0, hi = K - 1;
        while (lo < hi) { int mid = (lo + hi) >> 1; if (cdf[mid] < u) lo = mid + 1; else hi = mid; }
        deg[v] = kmin + lo;
        stub_sum += deg[v];
    }
    free(cdf);
    if (stub_sum & 1LL) { int v = (int)rng_below((uint64_t)N); deg[v]++; stub_sum++; }

    /* --- lista de stubs y barajado de Fisher-Yates --- */
    long long nstub = stub_sum;
    int *stubs = (int *)malloc((size_t)nstub * sizeof(int));
    long long idx = 0;
    for (int v = 0; v < N; v++) for (int j = 0; j < deg[v]; j++) stubs[idx++] = v;
    for (long long i = nstub - 1; i > 0; i--) {
        long long j = (long long)rng_below((uint64_t)(i + 1));
        int tmp = stubs[i]; stubs[i] = stubs[j]; stubs[j] = tmp;
    }
    free(deg);

    /* --- aristas como clave a*N+b (a<b), sin autoenlaces --- */
    long long npair = nstub / 2;
    uint64_t *keys = (uint64_t *)malloc((size_t)npair * sizeof(uint64_t));
    long long cnt = 0;
    for (long long m = 0; m < npair; m++) {
        int a = stubs[2 * m], b = stubs[2 * m + 1];
        if (a == b) continue;                 /* fuera autoenlaces            */
        if (a > b) { int t = a; a = b; b = t; }
        keys[cnt++] = (uint64_t)a * (uint64_t)N + (uint64_t)b;
    }
    free(stubs);

    /* ordenar + unique  -> elimina multienlaces */
    qsort(keys, (size_t)cnt, sizeof(uint64_t), cmp_u64);
    long long E = 0;
    for (long long i = 0; i < cnt; i++) if (i == 0 || keys[i] != keys[i - 1]) keys[E++] = keys[i];
    n_edges = E;
    n_dir   = 2 * E;

    /* --- construir CSR + rev[] --- */
    adj_start = (long *)calloc((size_t)N + 1, sizeof(long));
    for (long long i = 0; i < E; i++) {
        int a = (int)(keys[i] / (uint64_t)N);
        int b = (int)(keys[i] % (uint64_t)N);
        adj_start[a + 1]++; adj_start[b + 1]++;
    }
    for (int v = 0; v < N; v++) adj_start[v + 1] += adj_start[v];

    adj = (int  *)malloc((size_t)n_dir * sizeof(int));
    rev = (long *)malloc((size_t)n_dir * sizeof(long));
    long *cur = (long *)malloc(((size_t)N) * sizeof(long));
    for (int v = 0; v < N; v++) cur[v] = adj_start[v];
    for (long long i = 0; i < E; i++) {
        int a = (int)(keys[i] / (uint64_t)N);
        int b = (int)(keys[i] % (uint64_t)N);
        long pa = cur[a]++; long pb = cur[b]++;
        adj[pa] = b; adj[pb] = a;
        rev[pa] = pb; rev[pb] = pa;
    }
    free(cur);
    free(keys);

    /* --- momentos del grado y lista de semillas (grado == kmin) --- */
    double s1 = 0.0, s2 = 0.0;
    int min_deg = N;
    for (int v = 0; v < N; v++) {
        long d = adj_start[v + 1] - adj_start[v];
        s1 += (double)d; s2 += (double)d * (double)d;
        if ((int)d < min_deg) min_deg = (int)d;
    }
    *mean_k  = s1 / N;
    *mean_k2 = s2 / N;

    int target = kmin;                 /* semilla de grado minimo nominal kmin */
    long c0 = 0;
    for (int v = 0; v < N; v++) if ((adj_start[v + 1] - adj_start[v]) == target) c0++;
    if (c0 == 0) { target = min_deg; for (int v = 0; v < N; v++) if ((adj_start[v+1]-adj_start[v])==target) c0++; }
    seed_nodes = (int *)malloc((size_t)c0 * sizeof(int));
    n_seed = 0;
    for (int v = 0; v < N; v++) if ((adj_start[v + 1] - adj_start[v]) == target) seed_nodes[n_seed++] = v;

    fprintf(stderr,
        "# red construida: N=%d  E=%ld  <k>=%.4f  <k^2>=%.4f  kmin_real=%d  "
        "semillas(grado=%d)=%ld  lambda_c^MF=<k>/<k^2>=%.5f\n",
        N, n_edges, *mean_k, *mean_k2, min_deg, target, n_seed, (*mean_k) / (*mean_k2));
}

/* ============================== una realizacion =========================== */
/* Devuelve 1 si endemica (cobertura>=Cth), 0 si finita. Si finita, *tau.     */
static int run_realization(double lambda, double delta, long Cth, double *tau_out) {
    int v0 = seed_nodes[rng_below((uint64_t)n_seed)];
    double t = 0.0;

    add_inf(v0);
    covered[v0] = 1; touched[n_touch++] = v0; cov = 1;
    update_SI(v0);

    int endemic = 0;
    while (n_inf > 0) {
        double rate_inf = lambda * (double)n_act;
        double rate_rec = delta  * (double)n_inf;
        double ltot = rate_inf + rate_rec;

        /* avance temporal: xi en (0,1] para evitar log(0) */
        double xi = 1.0 - rng_unit();
        t += -log(xi) / ltot;

        if (rng_unit() * ltot < rate_inf) {
            /* --- INFECCION: arista activa al azar -> infectar el extremo S --- */
            long a = (long)rng_below((uint64_t)n_act);
            long de = act[a];
            int w = adj[de];                 /* extremo susceptible             */
            add_inf(w);
            if (!covered[w]) { covered[w] = 1; touched[n_touch++] = w; cov++; }
            update_SI(w);
            if (cov >= Cth) { endemic = 1; break; }
        } else {
            /* --- RECUPERACION: nodo infectado al azar --- */
            long r = (long)rng_below((uint64_t)n_inf);
            int v = inf_list[r];
            inf_list[r] = inf_list[--n_inf];
            state[v] = 0;
            update_IS(v);
        }
    }
    if (!endemic) *tau_out = t;

    /* --- reset O(cobertura): solo lo tocado --- */
    for (long i = 0; i < n_act; i++) act_pos[act[i]] = -1;
    n_act = 0;
    for (long i = 0; i < n_touch; i++) { int v = touched[i]; state[v] = 0; covered[v] = 0; }
    n_touch = 0; n_inf = 0; cov = 0;
    return endemic;
}

/* ================================== main ================================== */
static void usage(const char *p) {
    fprintf(stderr,
      "uso: %s --N <int> --gamma <f> [--kmin 4] [--kmax -1] [--delta 1.0]\n"
      "        [--cth 0.5] --realizations <int>\n"
      "        (--lambda-list \"l1,l2,...\"  |  --lambda-min <f> --lambda-max <f> --lambda-steps <int>)\n"
      "        [--net-seed <u64>] [--dyn-seed <u64>] [--out file.csv] [--save-degrees file]\n"
      "  kmax=-1 -> usa N-1 (sin corte). Para UCM pasar kmax = floor(sqrt(N)).\n", p);
}

int main(int argc, char **argv) {
    int    Nn = 0, kmin = 4, kmax = -1, M = 0, lsteps = 0;
    double gamma = 0.0, delta = 1.0, cth = 0.5, lmin = 0, lmax = 0;
    uint64_t net_seed = 12345ULL, dyn_seed = 67890ULL;
    char  *out = NULL, *llist = NULL, *degfile = NULL;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--N"))            Nn = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--gamma"))        gamma = atof(argv[++i]);
        else if (!strcmp(argv[i], "--kmin"))         kmin = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--kmax"))         kmax = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--delta"))        delta = atof(argv[++i]);
        else if (!strcmp(argv[i], "--cth"))          cth = atof(argv[++i]);
        else if (!strcmp(argv[i], "--realizations")) M = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--lambda-list"))  llist = argv[++i];
        else if (!strcmp(argv[i], "--lambda-min"))   lmin = atof(argv[++i]);
        else if (!strcmp(argv[i], "--lambda-max"))   lmax = atof(argv[++i]);
        else if (!strcmp(argv[i], "--lambda-steps")) lsteps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--net-seed"))     net_seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--dyn-seed"))     dyn_seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--out"))          out = argv[++i];
        else if (!strcmp(argv[i], "--save-degrees")) degfile = argv[++i];
        else { fprintf(stderr, "arg desconocido: %s\n", argv[i]); usage(argv[0]); return 1; }
    }
    if (Nn <= 0 || gamma <= 0 || M <= 0 || (!llist && lsteps <= 0)) { usage(argv[0]); return 1; }
    if (kmax < 0) kmax = Nn - 1;
    if (kmax > Nn - 1) kmax = Nn - 1;

    /* lista de lambdas */
    double *lam = NULL; int nlam = 0;
    if (llist) {
        int cap = 8; lam = malloc(cap * sizeof(double));
        char *s = strdup(llist), *tok = strtok(s, ",");
        while (tok) { if (nlam == cap) { cap *= 2; lam = realloc(lam, cap * sizeof(double)); }
                      lam[nlam++] = atof(tok); tok = strtok(NULL, ","); }
        free(s);
    } else {
        nlam = lsteps; lam = malloc(nlam * sizeof(double));
        for (int i = 0; i < nlam; i++)
            lam[i] = (nlam == 1) ? lmin : lmin + (lmax - lmin) * i / (nlam - 1);
    }

    /* construir red */
    double mean_k, mean_k2;
    build_network(Nn, gamma, kmin, kmax, net_seed, &mean_k, &mean_k2);

    if (degfile) {
        FILE *fd = fopen(degfile, "w");
        for (int v = 0; v < N; v++) fprintf(fd, "%ld\n", adj_start[v + 1] - adj_start[v]);
        fclose(fd);
    }

    /* reservar estructuras de la dinamica (una vez) */
    state    = (uint8_t *)calloc((size_t)N, 1);
    inf_list = (int     *)malloc((size_t)N * sizeof(int));
    covered  = (uint8_t *)calloc((size_t)N, 1);
    touched  = (int     *)malloc((size_t)N * sizeof(int));
    act      = (long    *)malloc((size_t)n_edges * sizeof(long));
    act_pos  = (long    *)malloc((size_t)n_dir   * sizeof(long));
    for (long i = 0; i < n_dir; i++) act_pos[i] = -1;
    n_inf = n_act = n_touch = cov = 0;

    long Cth = (long)(cth * (double)N);
    if (Cth < 1) Cth = 1;

    rng_seed(dyn_seed);

    FILE *f = out ? fopen(out, "w") : stdout;
    fprintf(f, "# SIS lifespan  N=%d gamma=%.3f kmin=%d kmax=%d delta=%.3f cth=%.3f "
               "realizations=%d net_seed=%llu dyn_seed=%llu\n",
               N, gamma, kmin, kmax, delta, cth, M,
               (unsigned long long)net_seed, (unsigned long long)dyn_seed);
    fprintf(f, "# E=%ld mean_k=%.6f mean_k2=%.6f lambda_c_MF=%.6f Cth=%ld\n",
               n_edges, mean_k, mean_k2, mean_k / mean_k2, Cth);
    fprintf(f, "lambda,p_end,p_end_err,tau_mean,tau_err,tau2_mean,tau2_err,n_finite,n_endemic,n_real,mean_k,mean_k2\n");

    for (int il = 0; il < nlam; il++) {
        double lambda = lam[il];
        long n_end = 0, n_fin = 0;
        double S1 = 0, S2 = 0, S3 = 0, S4 = 0;   /* sumas de potencias de tau */
        for (int r = 0; r < M; r++) {
            double tau = 0.0;
            int e = run_realization(lambda, delta, Cth, &tau);
            if (e) n_end++;
            else { n_fin++; double t2 = tau * tau; S1 += tau; S2 += t2; S3 += t2 * tau; S4 += t2 * t2; }
        }
        double p = (double)n_end / (double)M;
        double p_err = sqrt(p * (1.0 - p) / (double)M);
        double nf = (double)n_fin;
        double tau_mean  = (n_fin > 0) ? S1 / nf : 0.0;
        double tau2_mean = (n_fin > 0) ? S2 / nf : 0.0;
        double tau_err   = (n_fin > 1) ? sqrt((S2 / nf - tau_mean  * tau_mean ) / nf) : 0.0;
        double tau2_err  = (n_fin > 1) ? sqrt((S4 / nf - tau2_mean * tau2_mean) / nf) : 0.0;
        fprintf(f, "%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%ld,%ld,%d,%.6f,%.6f\n",
                lambda, p, p_err, tau_mean, tau_err, tau2_mean, tau2_err,
                n_fin, n_end, M, mean_k, mean_k2);
        fflush(f);
        fprintf(stderr, "  lambda=%.5f  Pend=%.4f  <tau>=%.3f  <tau2>=%.3f  (fin=%ld end=%ld)\n",
                lambda, p, tau_mean, tau2_mean, n_fin, n_end);
    }
    if (out) fclose(f);
    return 0;
}
