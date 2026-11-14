/*
 * _fast_enum.c: Fast C99 kernel for enumerating candidate minimal admissible ch-diagrams on a
 * 4-regular planar shadow using:
 *   1. Early marked-marked bigon (Omega_7 + Omega_6-) pruning on (m_mask, yz_mask),
 *   2. Linking-parity pruning on (S_+, S_-) before classical sign assignment,
 *   3. Lazy 2-stage Fast Kauffman Bracket Butterfly Transform (O(c * 2^c)) modulo P = 10^9 + 7
 *      at two independent evaluation points A = 3 and A = 5 (evaluating S_- and A = 5 ONLY when S_+
 *      at A = 3 has surviving sign assignments!),
 *   4. Exact pairwise linking-number == 0 verification on surviving sign masks,
 *   5. Optional exact 1-step Yoshikawa reduction pruning (Omega_3 / Omega_4 / Omega_5 / Omega_7
 *      followed by Omega_1- / Omega_6- / Omega_2-).
 *
 * Every candidate returned by this kernel is subsequently certified in Python by explicit
 * Reidemeister move unlinking sequences (zero false negatives AND zero false positives).
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MOD_P 1000000007LL
#define MAX_N 16
#define MAX_STATES (1 << 12)

typedef struct {
    uint16_t m_mask;
    uint16_t yz_mask;
    uint16_t sign_mask;
} CandidateDiagram;

static int64_t mod_pow(int64_t base, int exp) {
    int64_t res = 1;
    base %= MOD_P;
    if (base < 0) base += MOD_P;
    while (exp > 0) {
        if (exp & 1) res = (res * base) % MOD_P;
        base = (base * base) % MOD_P;
        exp >>= 1;
    }
    return res;
}

static int64_t mod_inv(int64_t x) {
    return mod_pow(x, (int)(MOD_P - 2));
}

/* Vertex types: 0 = X0, 1 = X1, 2 = Y, 3 = Z */
static inline int c_is_over(int vt, int slot) {
    return (vt == 0) ? ((slot & 1) == 1) : ((slot & 1) == 0);
}

static int has_immediate_reduction(int n, const int *opp, const int *vt) {
    /* 1. Omega_1- (classical kink) or Omega_6- (marked kink) */
    for (int v = 0; v < n; v++) {
        for (int s = 0; s < 4; s++) {
            if (opp[4 * v + s] == 4 * v + ((s + 1) & 3)) {
                return 1;
            }
        }
    }
    /* 2. Omega_2- (classical bigon) */
    for (int u = 0; u < n; u++) {
        if (vt[u] >= 2) continue;
        for (int su = 0; su < 4; su++) {
            int od1 = opp[4 * u + su];
            int w = od1 >> 2;
            if (w <= u || vt[w] >= 2) continue;
            int sw = od1 & 3;
            int su2 = (su + 1) & 3;
            int sw2 = (sw - 1) & 3;
            if (opp[4 * u + su2] != ((w << 2) | sw2)) continue;
            if (c_is_over(vt[u], su) == c_is_over(vt[w], sw)) {
                return 1;
            }
        }
    }
    return 0;
}

static void rewire_ports_c(
    int n,
    const int *opp,
    const int *vt,
    int num_ws,
    const int *ws,
    const int *new_ws_vt,
    int num_int_edges,
    const int int_edges[][2],
    int num_ports,
    const int *ports,
    const int port_assign[][2], /* {dart_local, port_idx} */
    int *out_opp,
    int *out_vt
) {
    memcpy(out_opp, opp, sizeof(int) * 4 * n);
    memcpy(out_vt, vt, sizeof(int) * n);

    int ext_targets[8];
    int port_to_new_dart[4 * MAX_N];
    for (int i = 0; i < 4 * n; i++) port_to_new_dart[i] = -1;

    for (int p = 0; p < num_ports; p++) {
        ext_targets[p] = opp[ports[p]];
    }
    for (int a = 0; a < num_ports; a++) {
        int dart_local = port_assign[a][0];
        int p_idx = port_assign[a][1];
        port_to_new_dart[ports[p_idx]] = dart_local;
    }

    for (int e = 0; e < num_int_edges; e++) {
        int da = int_edges[e][0];
        int db = int_edges[e][1];
        out_opp[da] = db;
        out_opp[db] = da;
    }

    for (int a = 0; a < num_ports; a++) {
        int dart_local = port_assign[a][0];
        int p_idx = port_assign[a][1];
        int target = ext_targets[p_idx];
        if (port_to_new_dart[target] >= 0) {
            out_opp[dart_local] = port_to_new_dart[target];
        } else {
            out_opp[dart_local] = target;
            out_opp[target] = dart_local;
        }
    }

    for (int i = 0; i < num_ws; i++) {
        out_vt[ws[i]] = new_ws_vt[i];
    }
}

static int has_0_or_2_edge_cut_c(int n, const int *opp);

static int has_1step_yoshikawa_reduction(int n, const int *opp, const int *vt, int check_2cut) {
    if (n < 2) return 0;
    int tmp_opp[4 * MAX_N];
    int tmp_vt[MAX_N];

    /* 1. Omega_7 neighbors along marked-marked edges */
    for (int u = 0; u < n; u++) {
        int ut = vt[u];
        if (ut < 2) continue;
        for (int su = 0; su < 4; su++) {
            int d1 = (u << 2) | su;
            int d2 = opp[d1];
            if (d2 < d1) continue;
            int w = d2 >> 2;
            int sw = d2 & 3;
            if (w == u || vt[w] < 2) continue;
            int wt = vt[w];
            int same_type = (ut == wt) ? 1 : 0;
            int same_parity = ((su & 1) == (sw & 1)) ? 1 : 0;
            if (same_type == same_parity) continue;

            int ports[6] = {
                (u << 2) | ((su + 1) & 3),
                (u << 2) | ((su + 2) & 3),
                (u << 2) | ((su + 3) & 3),
                (w << 2) | ((sw + 1) & 3),
                (w << 2) | ((sw + 2) & 3),
                (w << 2) | ((sw + 3) & 3)
            };
            int lplus_pairs_p2_p3 = ((ut == 2 && (su & 1) == 0) || (ut == 3 && (su & 1) == 1)) ? 1 : 0;
            int ws[2] = {u, w};
            int int_edges[1][2] = {{4 * u + 0, 4 * w + 0}};

            for (int shift = 1; shift <= 2; shift++) {
                int want_yz = lplus_pairs_p2_p3 ? (shift == 2) : (shift == 1);
                int new_wsvt[2] = {want_yz ? 2 : 3, want_yz ? 3 : 2};
                int port_assign[6][2] = {
                    {4 * u + 1, (shift + 0) % 6},
                    {4 * u + 2, (shift + 1) % 6},
                    {4 * u + 3, (shift + 2) % 6},
                    {4 * w + 1, (shift + 3) % 6},
                    {4 * w + 2, (shift + 4) % 6},
                    {4 * w + 3, (shift + 5) % 6}
                };
                rewire_ports_c(n, opp, vt, 2, ws, new_wsvt, 1, int_edges, 6, ports, port_assign, tmp_opp, tmp_vt);
                if (has_immediate_reduction(n, tmp_opp, tmp_vt)) return 1;
                if (check_2cut && n >= 4 && has_0_or_2_edge_cut_c(n, tmp_opp)) return 1;
            }
        }
    }

    /* 2. Omega_5 neighbors along classical-marked 2-gons */
    for (int u = 0; u < n; u++) {
        int ut = vt[u];
        if (ut >= 2) continue;
        for (int su = 0; su < 4; su++) {
            int su2 = (su + 1) & 3;
            int od1 = opp[4 * u + su];
            int v = od1 >> 2;
            int sv = od1 & 3;
            if (vt[v] < 2) continue;
            int sv2 = (sv - 1) & 3;
            if (opp[4 * u + su2] != ((v << 2) | sv2)) continue;

            int ports[4] = {
                (u << 2) | ((su + 2) & 3),
                (u << 2) | ((su + 3) & 3),
                (v << 2) | ((sv + 1) & 3),
                (v << 2) | ((sv + 2) & 3)
            };
            int vt_v = vt[v];
            int lplus_is_thru = ((vt_v == 3 && (sv & 1) == 0) || (vt_v == 2 && (sv & 1) == 1)) ? 1 : 0;
            int ac_over_bd = c_is_over(ut, su);
            int ws[2] = {v, u};
            int new_wsvt[2] = {lplus_is_thru ? 2 : 3, ac_over_bd ? 1 : 0};
            int int_edges[2][2] = {{4 * v + 0, 4 * u + 3}, {4 * v + 1, 4 * u + 2}};
            int port_assign[4][2] = {
                {4 * v + 2, 0},
                {4 * v + 3, 1},
                {4 * u + 0, 2},
                {4 * u + 1, 3}
            };
            rewire_ports_c(n, opp, vt, 2, ws, new_wsvt, 2, int_edges, 4, ports, port_assign, tmp_opp, tmp_vt);
            if (has_immediate_reduction(n, tmp_opp, tmp_vt)) return 1;
        }
    }

    /* 3. Omega_3 and Omega_4 / Omega_4' triangle neighbors */
    if (n >= 3) {
        for (int start_d = 0; start_d < 4 * n; start_d++) {
            int d0 = start_d;
            int od0 = opp[d0];
            int v1 = od0 >> 2, t1 = od0 & 3;
            int s1 = (t1 - 1) & 3;
            int d1 = (v1 << 2) | s1;

            int od1 = opp[d1];
            int v2 = od1 >> 2, t2 = od1 & 3;
            int s2 = (t2 - 1) & 3;
            int d2 = (v2 << 2) | s2;

            int od2 = opp[d2];
            int v0 = od2 >> 2, t0 = od2 & 3;
            int s0 = (t0 - 1) & 3;
            if (((v0 << 2) | s0) != d0 || v0 == v1 || v1 == v2 || v2 == v0) continue;
            /* Visit each oriented triangle once: require d0 < d1 && d0 < d2 */
            if (d0 > d1 || d0 > d2) continue;

            int vs[3] = {v0, v1, v2};
            int ss[3] = {s0, s1, s2};
            int m_cnt = 0, m_idx = -1;
            for (int i = 0; i < 3; i++) {
                if (vt[vs[i]] >= 2) {
                    m_cnt++;
                    m_idx = i;
                }
            }

            if (m_cnt == 0) {
                /* Omega_3 */
                int a1 = c_is_over(vt[v0], s0), a2 = c_is_over(vt[v1], t1);
                int b1 = c_is_over(vt[v1], s1), b2 = c_is_over(vt[v2], t2);
                int c1 = c_is_over(vt[v2], s2), c2 = c_is_over(vt[v0], t0);
                int ha = (a1 && a2) ? 2 : ((!a1 && !a2) ? 0 : 1);
                int hb = (b1 && b2) ? 2 : ((!b1 && !b2) ? 0 : 1);
                int hc = (c1 && c2) ? 2 : ((!c1 && !c2) ? 0 : 1);
                if (((1 << ha) | (1 << hb) | (1 << hc)) != 7) continue;

                int ports[6] = {
                    (v0 << 2) | ((s0 + 2) & 3),
                    (v0 << 2) | ((s0 + 3) & 3),
                    (v1 << 2) | ((s1 + 2) & 3),
                    (v1 << 2) | ((s1 + 3) & 3),
                    (v2 << 2) | ((s2 + 2) & 3),
                    (v2 << 2) | ((s2 + 3) & 3)
                };
                int ws[3] = {v0, v1, v2};
                int new_wsvt[3] = {
                    (hb > hc) ? 0 : 1,
                    (hc > ha) ? 0 : 1,
                    (ha > hb) ? 0 : 1
                };
                int int_edges[3][2] = {
                    {4 * ws[0] + 0, 4 * ws[1] + 1},
                    {4 * ws[1] + 0, 4 * ws[2] + 1},
                    {4 * ws[2] + 0, 4 * ws[0] + 1}
                };
                int port_assign[6][2] = {
                    {4 * ws[0] + 2, 1},
                    {4 * ws[0] + 3, 2},
                    {4 * ws[1] + 2, 3},
                    {4 * ws[1] + 3, 4},
                    {4 * ws[2] + 2, 5},
                    {4 * ws[2] + 3, 0}
                };
                rewire_ports_c(n, opp, vt, 3, ws, new_wsvt, 3, int_edges, 6, ports, port_assign, tmp_opp, tmp_vt);
                if (has_immediate_reduction(n, tmp_opp, tmp_vt)) return 1;
                if (check_2cut && n >= 4 && has_0_or_2_edge_cut_c(n, tmp_opp)) return 1;
            } else if (m_cnt == 1) {
                /* Omega_4 / Omega_4' */
                int mv = vs[m_idx], cv1 = vs[(m_idx + 1) % 3], cv2 = vs[(m_idx + 2) % 3];
                int ms = ss[m_idx], cs1 = ss[(m_idx + 1) % 3], cs2 = ss[(m_idx + 2) % 3];
                int ct2 = (cs2 + 1) & 3;
                int b_over1 = c_is_over(vt[cv1], cs1);
                int b_over2 = c_is_over(vt[cv2], ct2);
                if (b_over1 != b_over2) continue;

                int ports[6] = {
                    (mv << 2) | ((ms + 2) & 3),
                    (mv << 2) | ((ms + 3) & 3),
                    (cv1 << 2) | ((cs1 + 2) & 3),
                    (cv1 << 2) | ((cs1 + 3) & 3),
                    (cv2 << 2) | ((cs2 + 2) & 3),
                    (cv2 << 2) | ((cs2 + 3) & 3)
                };
                int vt_m = vt[mv];
                int lplus_pairs_ms_ms3 = ((vt_m == 2 && (ms & 1) == 0) || (vt_m == 3 && (ms & 1) == 1)) ? 1 : 0;
                int ws[3] = {cv2, mv, cv1};
                int new_wsvt[3] = {
                    b_over1 ? 0 : 1,
                    lplus_pairs_ms_ms3 ? 2 : 3,
                    b_over1 ? 1 : 0
                };
                int int_edges[3][2] = {
                    {4 * ws[0] + 0, 4 * ws[1] + 1},
                    {4 * ws[1] + 0, 4 * ws[2] + 1},
                    {4 * ws[2] + 0, 4 * ws[0] + 1}
                };
                int port_assign[6][2] = {
                    {4 * ws[0] + 2, 1},
                    {4 * ws[0] + 3, 2},
                    {4 * ws[1] + 2, 3},
                    {4 * ws[1] + 3, 4},
                    {4 * ws[2] + 2, 5},
                    {4 * ws[2] + 3, 0}
                };
                rewire_ports_c(n, opp, vt, 3, ws, new_wsvt, 3, int_edges, 6, ports, port_assign, tmp_opp, tmp_vt);
                if (has_immediate_reduction(n, tmp_opp, tmp_vt)) return 1;
                if (check_2cut && n >= 4 && has_0_or_2_edge_cut_c(n, tmp_opp)) return 1;
            }
        }
    }
    return 0;
}

/* Reduced classical shadow built from MarkedDiagram shadow after splicing markers */
typedef struct {
    int c;
    int free_loops;
    int opp[4 * MAX_N];
    int num_circles;
    int base_sign[MAX_N];      /* +1 or -1 when crossing i is X0 */
    int comp_a[MAX_N];         /* circle id of strand (0,2) */
    int comp_b[MAX_N];         /* circle id of strand (1,3) */
    int odd_linking_fail;      /* 1 if any pair of distinct circles crosses an odd # of times */
    int num_lk_pairs;
    int lk_pair_mask[MAX_N * MAX_N]; /* bitmask of crossings between pair p */
    int lk_pair_sign[MAX_N * MAX_N][MAX_N];
} ReducedShadow;

static void build_reduced_shadow(
    int n,
    const int *orig_opp,
    uint16_t m_mask,
    uint16_t yz_mask,
    int sign_res, /* +1 for L_+, -1 for L_- */
    ReducedShadow *out
) {
    int old_to_c[MAX_N];
    int splice[4 * MAX_N];
    int is_marked_dart[4 * MAX_N];
    memset(is_marked_dart, 0, sizeof(int) * 4 * n);

    int c = 0;
    int m_idx = 0;
    for (int v = 0; v < n; v++) {
        if ((m_mask >> v) & 1) {
            old_to_c[v] = -1;
            int is_z = (yz_mask >> m_idx) & 1;
            m_idx++;
            int pair_03 = ((!is_z && sign_res == 1) || (is_z && sign_res == -1)) ? 1 : 0;
            for (int s = 0; s < 4; s++) {
                is_marked_dart[4 * v + s] = 1;
                if (pair_03) {
                    splice[4 * v + s] = 4 * v + (3 - s);
                } else {
                    splice[4 * v + s] = 4 * v + (s ^ 1);
                }
            }
        } else {
            old_to_c[v] = c++;
        }
    }

    out->c = c;
    int visited_md[4 * MAX_N];
    memset(visited_md, 0, sizeof(int) * 4 * n);

    for (int v = 0; v < n; v++) {
        int cv = old_to_c[v];
        if (cv < 0) continue;
        for (int s = 0; s < 4; s++) {
            int cur = orig_opp[4 * v + s];
            while (is_marked_dart[cur]) {
                visited_md[cur] = 1;
                int nxt = splice[cur];
                visited_md[nxt] = 1;
                cur = orig_opp[nxt];
            }
            int w_old = cur >> 2;
            int ws = cur & 3;
            out->opp[4 * cv + s] = (old_to_c[w_old] << 2) | ws;
        }
    }

    int extra_loops = 0;
    for (int d = 0; d < 4 * n; d++) {
        if (!is_marked_dart[d] || visited_md[d]) continue;
        extra_loops++;
        int cur = d;
        while (is_marked_dart[cur] && !visited_md[cur]) {
            visited_md[cur] = 1;
            int nxt = splice[cur];
            visited_md[nxt] = 1;
            cur = orig_opp[nxt];
        }
    }
    out->free_loops = extra_loops;

    /* Trace transverse circles of reduced classical shadow */
    int seen[4 * MAX_N];
    int dart_dir[4 * MAX_N];
    int dart_comp[4 * MAX_N];
    memset(seen, 0, sizeof(int) * 4 * c);

    int num_cyc = 0;
    for (int sd = 0; sd < 4 * c; sd++) {
        if (seen[sd]) continue;
        int cur = sd;
        while (!seen[cur]) {
            seen[cur] = 1;
            int ex = cur ^ 2;
            seen[ex] = 1;
            dart_dir[cur] = 1;
            dart_dir[ex] = -1;
            dart_comp[cur] = num_cyc;
            dart_comp[ex] = num_cyc;
            cur = out->opp[ex];
        }
        num_cyc++;
    }
    out->num_circles = num_cyc + extra_loops;

    int pair_cnt[MAX_N][MAX_N];
    memset(pair_cnt, 0, sizeof(pair_cnt));

    for (int i = 0; i < c; i++) {
        int s_under_in = (dart_dir[4 * i + 0] == 1) ? 0 : 2;
        int s_over_in = (dart_dir[4 * i + 1] == 1) ? 1 : 3;
        int diff = (s_over_in - s_under_in) & 3;
        out->base_sign[i] = (diff == 3) ? 1 : -1;
        int ca = dart_comp[4 * i + 0];
        int cb = dart_comp[4 * i + 1];
        out->comp_a[i] = ca;
        out->comp_b[i] = cb;
        if (ca != cb) {
            int p1 = (ca < cb) ? ca : cb;
            int p2 = (ca < cb) ? cb : ca;
            pair_cnt[p1][p2]++;
        }
    }

    out->odd_linking_fail = 0;
    out->num_lk_pairs = 0;
    for (int p1 = 0; p1 < num_cyc; p1++) {
        for (int p2 = p1 + 1; p2 < num_cyc; p2++) {
            if (pair_cnt[p1][p2] & 1) {
                out->odd_linking_fail = 1;
                return;
            }
            if (pair_cnt[p1][p2] > 0) {
                int pid = out->num_lk_pairs++;
                int mask = 0;
                for (int i = 0; i < c; i++) {
                    int ca = out->comp_a[i];
                    int cb = out->comp_b[i];
                    if ((ca == p1 && cb == p2) || (ca == p2 && cb == p1)) {
                        mask |= (1 << i);
                        out->lk_pair_sign[pid][i] = out->base_sign[i];
                    } else {
                        out->lk_pair_sign[pid][i] = 0;
                    }
                }
                out->lk_pair_mask[pid] = mask;
            }
        }
    }
}

static void compute_loop_counts(const ReducedShadow *rs, uint8_t *loops_out) {
    int c = rs->c;
    int total_states = 1 << c;
    const int *opp = rs->opp;
    uint8_t seen[4 * MAX_N];

    for (int st = 0; st < total_states; st++) {
        memset(seen, 0, 4 * c);
        int loops = rs->free_loops;
        for (int sd = 0; sd < 4 * c; sd++) {
            if (seen[sd]) continue;
            loops++;
            int cur = sd;
            while (!seen[cur]) {
                seen[cur] = 1;
                int v = cur >> 2;
                int s = cur & 3;
                int partner = ((st >> v) & 1) ? ((v << 2) | (3 - s)) : ((v << 2) | (s ^ 1));
                seen[partner] = 1;
                cur = opp[partner];
            }
        }
        loops_out[st] = (uint8_t)loops;
    }
}

static void kauffman_butterfly(
    int c,
    const uint8_t *loops_arr,
    int64_t A_val,
    int64_t A_inv,
    const int64_t *d_pows,
    int64_t *buf
) {
    int total = 1 << c;
    for (int s = 0; s < total; s++) {
        buf[s] = d_pows[loops_arr[s] - 1];
    }
    for (int i = 0; i < c; i++) {
        int step = 1 << i;
        for (int blk = 0; blk < total; blk += (step << 1)) {
            for (int j = 0; j < step; j++) {
                int i0 = blk + j;
                int i1 = i0 + step;
                int64_t v0 = buf[i0];
                int64_t v1 = buf[i1];
                buf[i0] = (A_val * v0 + A_inv * v1) % MOD_P;
                buf[i1] = (A_inv * v0 + A_val * v1) % MOD_P;
            }
        }
    }
}

int enumerate_shadow_candidates(
    int n,
    const int *orig_opp,
    int m_min,
    int m_max,
    int prune_level, /* 0: none, 1: classical R-II bigons, 2: + 1-step Yoshikawa reductions */
    CandidateDiagram *out_buf,
    int max_out
) {
    const int64_t A_list[2] = {3LL, 5LL};
    int64_t A_inv_list[2];
    int64_t d_pows[2][2 * MAX_N + 2];
    int64_t neg_A3_pos[2][4 * MAX_N + 5];

    for (int t = 0; t < 2; t++) {
        int64_t A = A_list[t];
        int64_t Ainv = mod_inv(A);
        A_inv_list[t] = Ainv;
        int64_t A2 = (A * A) % MOD_P;
        int64_t Ainv2 = (Ainv * Ainv) % MOD_P;
        int64_t d_val = (MOD_P - ((A2 + Ainv2) % MOD_P)) % MOD_P;
        d_pows[t][0] = 1;
        for (int k = 1; k <= 2 * MAX_N + 1; k++) {
            d_pows[t][k] = (d_pows[t][k - 1] * d_val) % MOD_P;
        }
        int64_t neg_A3 = (MOD_P - ((A2 * A) % MOD_P)) % MOD_P;
        int64_t neg_A3_inv = mod_inv(neg_A3);
        int center = 2 * MAX_N;
        neg_A3_pos[t][center] = 1;
        for (int w = 1; w <= 2 * MAX_N; w++) {
            neg_A3_pos[t][center + w] = (neg_A3_pos[t][center + w - 1] * neg_A3) % MOD_P;
            neg_A3_pos[t][center - w] = (neg_A3_pos[t][center - w + 1] * neg_A3_inv) % MOD_P;
        }
    }

    /* Pre-detect planar 2-gons in `orig_opp` */
    int num_bigons = 0;
    int bigon_u[MAX_N * 2], bigon_w[MAX_N * 2], bigon_same_parity[MAX_N * 2];
    if (prune_level >= 1) {
        for (int u = 0; u < n; u++) {
            for (int su = 0; su < 4; su++) {
                int su2 = (su + 1) & 3;
                int od1 = orig_opp[4 * u + su];
                int w = od1 >> 2;
                if (w <= u) continue;
                int sw = od1 & 3;
                int sw2 = (sw - 1) & 3;
                if (orig_opp[4 * u + su2] == ((w << 2) | sw2)) {
                    bigon_u[num_bigons] = u;
                    bigon_w[num_bigons] = w;
                    bigon_same_parity[num_bigons] = ((su & 1) == (sw & 1)) ? 1 : 0;
                    num_bigons++;
                }
            }
        }
    }

    int out_count = 0;
    uint8_t loops_p[MAX_STATES];
    uint8_t loops_m[MAX_STATES];
    int64_t kp1[MAX_STATES], km1[MAX_STATES];
    int64_t kp2[MAX_STATES], km2[MAX_STATES];
    int surv_smask[MAX_STATES];
    int surv_wp[MAX_STATES];
    int surv_wm[MAX_STATES];

    int total_masks = 1 << n;
    for (int m_mask = 0; m_mask < total_masks; m_mask++) {
        int m = __builtin_popcount((unsigned int)m_mask);
        if (m < m_min || m > m_max) continue;
        int c = n - m;

        int orig_to_c[MAX_N];
        int orig_to_m[MAX_N];
        int cc = 0, mm = 0;
        for (int v = 0; v < n; v++) {
            if ((m_mask >> v) & 1) {
                orig_to_c[v] = -1;
                orig_to_m[v] = mm++;
            } else {
                orig_to_c[v] = cc++;
                orig_to_m[v] = -1;
            }
        }

        /* Separate classical-classical bigons and marked-marked bigons */
        int act_bg = 0;
        int act_ci[MAX_N * 2], act_cj[MAX_N * 2], act_forbid_xor[MAX_N * 2];
        int mm_bg = 0;
        int mm_mi[MAX_N * 2], mm_mj[MAX_N * 2], mm_forbid_xor[MAX_N * 2];

        if (prune_level >= 1) {
            for (int b = 0; b < num_bigons; b++) {
                int u = bigon_u[b], w = bigon_w[b];
                int ci = orig_to_c[u], cj = orig_to_c[w];
                if (ci >= 0 && cj >= 0) {
                    act_ci[act_bg] = ci;
                    act_cj[act_bg] = cj;
                    act_forbid_xor[act_bg] = bigon_same_parity[b] ? 0 : 1;
                    act_bg++;
                } else if (prune_level >= 2 && n >= 3 && orig_to_m[u] >= 0 && orig_to_m[w] >= 0) {
                    /* Marked-marked bigon on n >= 3: Omega_7 applies iff same_type != same_parity,
                     * turning the second edge of the bigon into a 1-gon (Omega_6-)!
                     * Forbid (yz_bit_u ^ yz_bit_w) == (same_parity ? 1 : 0)! */
                    mm_mi[mm_bg] = orig_to_m[u];
                    mm_mj[mm_bg] = orig_to_m[w];
                    mm_forbid_xor[mm_bg] = bigon_same_parity[b] ? 1 : 0;
                    mm_bg++;
                }
            }
        }

        int yz_limit = (m > 0) ? (1 << (m - 1)) : 1;
        int sign_limit = (c > 0) ? (1 << (c - 1)) : 1;

        for (int yz_half = 0; yz_half < yz_limit; yz_half++) {
            uint16_t yz_mask = (uint16_t)(yz_half << 1);

            /* Early marked-marked bigon filter BEFORE building ReducedShadow! */
            int mm_hit = 0;
            for (int b = 0; b < mm_bg; b++) {
                int bi = (yz_mask >> mm_mi[b]) & 1;
                int bj = (yz_mask >> mm_mj[b]) & 1;
                if ((bi ^ bj) == mm_forbid_xor[b]) {
                    mm_hit = 1;
                    break;
                }
            }
            if (mm_hit) continue;

            ReducedShadow rsp, rsm;
            build_reduced_shadow(n, orig_opp, (uint16_t)m_mask, yz_mask, 1, &rsp);
            if (rsp.odd_linking_fail) continue;
            build_reduced_shadow(n, orig_opp, (uint16_t)m_mask, yz_mask, -1, &rsm);
            if (rsm.odd_linking_fail) continue;

            if (c == 0) {
                int vt[MAX_N];
                for (int v = 0; v < n; v++) {
                    vt[v] = ((yz_mask >> orig_to_m[v]) & 1) ? 3 : 2;
                }
                if (prune_level >= 2 && has_1step_yoshikawa_reduction(n, orig_opp, vt, prune_level >= 3 ? 1 : 0)) {
                    continue;
                }
                if (out_count < max_out) {
                    out_buf[out_count].m_mask = (uint16_t)m_mask;
                    out_buf[out_count].yz_mask = yz_mask;
                    out_buf[out_count].sign_mask = 0;
                }
                out_count++;
                continue;
            }

            /* Stage 1: compute L_+ loops + A1 = 3 butterfly ONLY! */
            compute_loop_counts(&rsp, loops_p);
            kauffman_butterfly(c, loops_p, A_list[0], A_inv_list[0], d_pows[0], kp1);

            int center = 2 * MAX_N;
            int64_t dp1 = d_pows[0][rsp.num_circles - 1];
            int surv_cnt = 0;

            for (int shalf = 0; shalf < sign_limit; shalf++) {
                int smask = shalf << 1;

                int r2_hit = 0;
                for (int b = 0; b < act_bg; b++) {
                    int bi = (smask >> act_ci[b]) & 1;
                    int bj = (smask >> act_cj[b]) & 1;
                    if ((bi ^ bj) == act_forbid_xor[b]) {
                        r2_hit = 1;
                        break;
                    }
                }
                if (r2_hit) continue;

                int wp = 0, wm = 0;
                for (int i = 0; i < c; i++) {
                    int mult = ((smask >> i) & 1) ? -1 : 1;
                    wp += mult * rsp.base_sign[i];
                    wm += mult * rsm.base_sign[i];
                }

                int64_t exp_kp1 = (neg_A3_pos[0][center + wp] * dp1) % MOD_P;
                if (kp1[smask] != exp_kp1) continue;

                /* Pairwise linking numbers == 0 in L_+ and L_- */
                int lk_fail = 0;
                for (int pid = 0; pid < rsp.num_lk_pairs; pid++) {
                    int sum = 0;
                    int pmask = rsp.lk_pair_mask[pid];
                    for (int i = 0; i < c; i++) {
                        if ((pmask >> i) & 1) {
                            sum += (((smask >> i) & 1) ? -1 : 1) * rsp.lk_pair_sign[pid][i];
                        }
                    }
                    if (sum != 0) { lk_fail = 1; break; }
                }
                if (lk_fail) continue;

                for (int pid = 0; pid < rsm.num_lk_pairs; pid++) {
                    int sum = 0;
                    int pmask = rsm.lk_pair_mask[pid];
                    for (int i = 0; i < c; i++) {
                        if ((pmask >> i) & 1) {
                            sum += (((smask >> i) & 1) ? -1 : 1) * rsm.lk_pair_sign[pid][i];
                        }
                    }
                    if (sum != 0) { lk_fail = 1; break; }
                }
                if (lk_fail) continue;

                surv_smask[surv_cnt] = smask;
                surv_wp[surv_cnt] = wp;
                surv_wm[surv_cnt] = wm;
                surv_cnt++;
            }

            if (surv_cnt == 0) continue;

            /* Stage 2: ONLY when L_+ has survivors, evaluate L_- at A1 = 3 and both at A2 = 5! */
            compute_loop_counts(&rsm, loops_m);
            kauffman_butterfly(c, loops_m, A_list[0], A_inv_list[0], d_pows[0], km1);

            int64_t dm1 = d_pows[0][rsm.num_circles - 1];
            int64_t dp2 = d_pows[1][rsp.num_circles - 1];
            int64_t dm2 = d_pows[1][rsm.num_circles - 1];
            int computed_a2 = 0;

            for (int k = 0; k < surv_cnt; k++) {
                int smask = surv_smask[k];
                int wp = surv_wp[k];
                int wm = surv_wm[k];

                int64_t exp_km1 = (neg_A3_pos[0][center + wm] * dm1) % MOD_P;
                if (km1[smask] != exp_km1) continue;

                if (!computed_a2) {
                    kauffman_butterfly(c, loops_p, A_list[1], A_inv_list[1], d_pows[1], kp2);
                    kauffman_butterfly(c, loops_m, A_list[1], A_inv_list[1], d_pows[1], km2);
                    computed_a2 = 1;
                }

                int64_t exp_kp2 = (neg_A3_pos[1][center + wp] * dp2) % MOD_P;
                if (kp2[smask] != exp_kp2) continue;
                int64_t exp_km2 = (neg_A3_pos[1][center + wm] * dm2) % MOD_P;
                if (km2[smask] != exp_km2) continue;

                if (prune_level >= 2) {
                    int vt[MAX_N];
                    for (int v = 0; v < n; v++) {
                        if (orig_to_m[v] >= 0) {
                            vt[v] = ((yz_mask >> orig_to_m[v]) & 1) ? 3 : 2;
                        } else {
                            vt[v] = ((smask >> orig_to_c[v]) & 1) ? 1 : 0;
                        }
                    }
                    if (has_1step_yoshikawa_reduction(n, orig_opp, vt, prune_level >= 3 ? 1 : 0)) {
                        continue;
                    }
                }

                if (out_count < max_out) {
                    out_buf[out_count].m_mask = (uint16_t)m_mask;
                    out_buf[out_count].yz_mask = yz_mask;
                    out_buf[out_count].sign_mask = (uint16_t)smask;
                }
                out_count++;
            }
        }
    }
    return out_count;
}

/* =========================================================================
 * Lossless Canonical Key & Global C Yoshikawa BFS Deduplicator
 * ========================================================================= */

typedef struct {
    uint32_t words[MAX_N];
} CanonKey;

static inline int canon_key_less(int n, const CanonKey *a, const CanonKey *b) {
    for (int i = 0; i < n; i++) {
        if (a->words[i] != b->words[i]) {
            return a->words[i] < b->words[i];
        }
    }
    return 0;
}

static inline int canon_key_eq(int n, const CanonKey *a, const CanonKey *b) {
    for (int i = 0; i < n; i++) {
        if (a->words[i] != b->words[i]) return 0;
    }
    return 1;
}

static inline uint64_t canon_key_hash(int n, const CanonKey *k) {
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < n; i++) {
        h ^= (uint64_t)k->words[i];
        h *= 1099511628211ULL;
    }
    h ^= (h >> 33);
    h *= 0xff51afd7ed558ccdULL;
    h ^= (h >> 33);
    return h;
}

/* Compute exact canonical spherical code across all 8 O(4)/reflection variants for connected (n, opp, vt) */
static void canonical_sym_key_c(int n, const int *opp, const int *vt, CanonKey *best) {
    int opp_refl[4 * MAX_N];
    static const int slot_refl[4] = {0, 3, 2, 1};
    int has_classical = 0;
    for (int v = 0; v < n; v++) {
        if (vt[v] < 2) has_classical = 1;
        for (int s = 0; s < 4; s++) {
            int od = opp[4 * v + s];
            opp_refl[4 * v + slot_refl[s]] = ((od >> 2) << 2) | slot_refl[od & 3];
        }
    }

    const int *opps[2] = {opp, opp_refl};
    int has_best = 0;
    memset(best, 0, sizeof(CanonKey));

    int ov2n[MAX_N];
    int slot_shift[MAX_N];
    int q[MAX_N];
    uint32_t cur_words[MAX_N];

    for (int r_idx = 0; r_idx < 2; r_idx++) {
        const int *cur_opp = opps[r_idx];
        for (int v0 = 0; v0 < n; v0++) {
            /* Any classical start gives vt_rel(q[0]) == 0, strictly beating any marked start (>= 2) */
            if (has_classical && vt[v0] >= 2) continue;
            for (int s0 = 0; s0 < 4; s0++) {
                for (int i = 0; i < n; i++) ov2n[i] = -1;
                ov2n[v0] = 0;
                slot_shift[v0] = s0;
                q[0] = v0;
                int qlen = 1;
                int cmp_eq = 1;
                int xor_x = -1;
                int xor_y = -1;

                for (int qi = 0; qi < n; qi++) {
                    int v = q[qi];
                    int sh = slot_shift[v];
                    int vt_orig = vt[v];
                    int vt_rel;
                    if (vt_orig < 2) {
                        int raw_b = vt_orig ^ (sh & 1);
                        if (xor_x < 0) {
                            xor_x = raw_b;
                            vt_rel = 0;
                        } else {
                            vt_rel = raw_b ^ xor_x;
                        }
                    } else {
                        int raw_b = (vt_orig - 2) ^ (sh & 1);
                        if (xor_y < 0) {
                            xor_y = raw_b;
                            vt_rel = 2;
                        } else {
                            vt_rel = 2 + (raw_b ^ xor_y);
                        }
                    }

                    int d_enc[4];
                    for (int k = 0; k < 4; k++) {
                        int cd = (v << 2) | ((sh + k) & 3);
                        int od = cur_opp[cd];
                        int ov = od >> 2;
                        int os = od & 3;
                        if (ov2n[ov] < 0) {
                            ov2n[ov] = qlen;
                            slot_shift[ov] = os;
                            q[qlen++] = ov;
                        }
                        int rel_os = (os - slot_shift[ov]) & 3;
                        d_enc[k] = (ov2n[ov] << 2) | rel_os;
                    }
                    uint32_t w = ((uint32_t)vt_rel << 24)
                               | ((uint32_t)d_enc[0] << 18)
                               | ((uint32_t)d_enc[1] << 12)
                               | ((uint32_t)d_enc[2] << 6)
                               | (uint32_t)d_enc[3];
                    if (has_best && cmp_eq) {
                        if (w > best->words[qi]) goto next_start;
                        if (w < best->words[qi]) cmp_eq = 0;
                    }
                    cur_words[qi] = w;
                }
                if (!has_best || !cmp_eq) {
                    for (int i = 0; i < n; i++) best->words[i] = cur_words[i];
                    has_best = 1;
                }
            next_start:;
            }
        }
    }
}

static int has_0_or_2_edge_cut_c(int n, const int *opp) {
    if (n <= 1) return 0;
    uint8_t mult[MAX_N][MAX_N];
    uint32_t adj[MAX_N];
    for (int i = 0; i < n; i++) {
        adj[i] = 0u;
        for (int j = 0; j < n; j++) mult[i][j] = 0;
    }

    int eu[2 * MAX_N], ew[2 * MAX_N], ne = 0;
    for (int d = 0; d < 4 * n; d++) {
        int od = opp[d];
        if (d < od) {
            int u = d >> 2;
            int w = od >> 2;
            if (u != w) {
                mult[u][w]++;
                mult[w][u]++;
                adj[u] |= (1u << w);
                adj[w] |= (1u << u);
                eu[ne] = u;
                ew[ne] = w;
                ne++;
            }
        }
    }

    uint32_t vis = 1u, todo = 1u;
    while (todo) {
        int u = __builtin_ctz(todo);
        todo &= todo - 1u;
        uint32_t nw = adj[u] & ~vis;
        vis |= nw;
        todo |= nw;
    }
    if (__builtin_popcount(vis) < n) return 1;
    if (n < 4) return 0;

    for (int i = 0; i < ne; i++) {
        int u1 = eu[i], w1 = ew[i];
        uint8_t rem1 = --mult[u1][w1];
        mult[w1][u1] = rem1;
        if (rem1 >= 2) {
            mult[u1][w1] = ++rem1;
            mult[w1][u1] = rem1;
            continue;
        }
        if (rem1 == 0) {
            adj[u1] &= ~(1u << w1);
            adj[w1] &= ~(1u << u1);
        }

        uint32_t w1_bit = 1u << w1;
        for (int j = i + 1; j < ne; j++) {
            int u2 = eu[j], w2 = ew[j];
            if (rem1 == 1 && !((u2 == u1 && w2 == w1) || (u2 == w1 && w2 == u1))) {
                continue;
            }
            uint8_t rem2 = --mult[u2][w2];
            mult[w2][u2] = rem2;
            if (rem2 == 0) {
                adj[u2] &= ~(1u << w2);
                adj[w2] &= ~(1u << u2);
            } else {
                mult[u2][w2] = ++rem2;
                mult[w2][u2] = rem2;
                continue;
            }

            vis = 1u << u1;
            todo = vis;
            int reached_w1 = 0;
            while (todo) {
                int u = __builtin_ctz(todo);
                todo &= todo - 1u;
                uint32_t nw = adj[u] & ~vis;
                if (nw & w1_bit) {
                    reached_w1 = 1;
                    break;
                }
                vis |= nw;
                todo |= nw;
            }

            adj[u2] |= (1u << w2);
            adj[w2] |= (1u << u2);
            mult[u2][w2] = 1;
            mult[w2][u2] = 1;

            if (!reached_w1) {
                int cnt = __builtin_popcount(vis);
                if (cnt >= 2 && cnt <= n - 2) {
                    return 1;
                }
            }
            if (rem1 == 1) break;
        }

        if (rem1 == 0) {
            adj[u1] |= (1u << w1);
            adj[w1] |= (1u << u1);
        }
        mult[u1][w1] = rem1 + 1;
        mult[w1][u1] = rem1 + 1;
    }
    return 0;
}

static int list_all_level_neighbors_c(
    int n,
    const int *opp,
    const int *vt,
    int out_opps[][4 * MAX_N],
    int out_vts[][MAX_N],
    int max_nbrs
) {
    if (n < 2) return 0;
    int cnt = 0;

    /* 1. Omega_7 */
    for (int u = 0; u < n; u++) {
        int ut = vt[u];
        if (ut < 2) continue;
        for (int su = 0; su < 4; su++) {
            int d1 = (u << 2) | su;
            int d2 = opp[d1];
            if (d2 < d1) continue;
            int w = d2 >> 2;
            int sw = d2 & 3;
            if (w == u || vt[w] < 2) continue;
            int wt = vt[w];
            if (((ut == wt) ? 1 : 0) == (((su & 1) == (sw & 1)) ? 1 : 0)) continue;

            int ports[6] = {
                (u << 2) | ((su + 1) & 3),
                (u << 2) | ((su + 2) & 3),
                (u << 2) | ((su + 3) & 3),
                (w << 2) | ((sw + 1) & 3),
                (w << 2) | ((sw + 2) & 3),
                (w << 2) | ((sw + 3) & 3)
            };
            int lplus_pairs_p2_p3 = ((ut == 2 && (su & 1) == 0) || (ut == 3 && (su & 1) == 1)) ? 1 : 0;
            int ws[2] = {u, w};
            int int_edges[1][2] = {{4 * u + 0, 4 * w + 0}};

            for (int shift = 1; shift <= 2; shift++) {
                if (cnt >= max_nbrs) return cnt;
                int want_yz = lplus_pairs_p2_p3 ? (shift == 2) : (shift == 1);
                int new_wsvt[2] = {want_yz ? 2 : 3, want_yz ? 3 : 2};
                int port_assign[6][2] = {
                    {4 * u + 1, (shift + 0) % 6},
                    {4 * u + 2, (shift + 1) % 6},
                    {4 * u + 3, (shift + 2) % 6},
                    {4 * w + 1, (shift + 3) % 6},
                    {4 * w + 2, (shift + 4) % 6},
                    {4 * w + 3, (shift + 5) % 6}
                };
                rewire_ports_c(n, opp, vt, 2, ws, new_wsvt, 1, int_edges, 6, ports, port_assign, out_opps[cnt], out_vts[cnt]);
                cnt++;
            }
        }
    }

    /* 2. Omega_5 */
    for (int u = 0; u < n; u++) {
        int ut = vt[u];
        if (ut >= 2) continue;
        for (int su = 0; su < 4; su++) {
            int su2 = (su + 1) & 3;
            int od1 = opp[4 * u + su];
            int v = od1 >> 2;
            int sv = od1 & 3;
            if (vt[v] < 2) continue;
            int sv2 = (sv - 1) & 3;
            if (opp[4 * u + su2] != ((v << 2) | sv2)) continue;
            if (cnt >= max_nbrs) return cnt;

            int ports[4] = {
                (u << 2) | ((su + 2) & 3),
                (u << 2) | ((su + 3) & 3),
                (v << 2) | ((sv + 1) & 3),
                (v << 2) | ((sv + 2) & 3)
            };
            int vt_v = vt[v];
            int lplus_is_thru = ((vt_v == 3 && (sv & 1) == 0) || (vt_v == 2 && (sv & 1) == 1)) ? 1 : 0;
            int ac_over_bd = c_is_over(ut, su);
            int ws[2] = {v, u};
            int new_wsvt[2] = {lplus_is_thru ? 2 : 3, ac_over_bd ? 1 : 0};
            int int_edges[2][2] = {{4 * v + 0, 4 * u + 3}, {4 * v + 1, 4 * u + 2}};
            int port_assign[4][2] = {
                {4 * v + 2, 0},
                {4 * v + 3, 1},
                {4 * u + 0, 2},
                {4 * u + 1, 3}
            };
            rewire_ports_c(n, opp, vt, 2, ws, new_wsvt, 2, int_edges, 4, ports, port_assign, out_opps[cnt], out_vts[cnt]);
            cnt++;
        }
    }

    /* 3. Omega_3 & Omega_4 */
    if (n >= 3) {
        for (int start_d = 0; start_d < 4 * n; start_d++) {
            int d0 = start_d;
            int od0 = opp[d0];
            int v1 = od0 >> 2, t1 = od0 & 3;
            int s1 = (t1 - 1) & 3;
            int d1 = (v1 << 2) | s1;

            int od1 = opp[d1];
            int v2 = od1 >> 2, t2 = od1 & 3;
            int s2 = (t2 - 1) & 3;
            int d2 = (v2 << 2) | s2;

            int od2 = opp[d2];
            int v0 = od2 >> 2, t0 = od2 & 3;
            int s0 = (t0 - 1) & 3;
            if (((v0 << 2) | s0) != d0 || v0 == v1 || v1 == v2 || v2 == v0) continue;
            if (d0 > d1 || d0 > d2) continue;

            int vs[3] = {v0, v1, v2};
            int ss[3] = {s0, s1, s2};
            int m_cnt = 0, m_idx = -1;
            for (int i = 0; i < 3; i++) {
                if (vt[vs[i]] >= 2) {
                    m_cnt++;
                    m_idx = i;
                }
            }

            if (m_cnt == 0) {
                int a1 = c_is_over(vt[v0], s0), a2 = c_is_over(vt[v1], t1);
                int b1 = c_is_over(vt[v1], s1), b2 = c_is_over(vt[v2], t2);
                int c1 = c_is_over(vt[v2], s2), c2 = c_is_over(vt[v0], t0);
                int ha = (a1 && a2) ? 2 : ((!a1 && !a2) ? 0 : 1);
                int hb = (b1 && b2) ? 2 : ((!b1 && !b2) ? 0 : 1);
                int hc = (c1 && c2) ? 2 : ((!c1 && !c2) ? 0 : 1);
                if (((1 << ha) | (1 << hb) | (1 << hc)) != 7) continue;
                if (cnt >= max_nbrs) return cnt;

                int ports[6] = {
                    (v0 << 2) | ((s0 + 2) & 3),
                    (v0 << 2) | ((s0 + 3) & 3),
                    (v1 << 2) | ((s1 + 2) & 3),
                    (v1 << 2) | ((s1 + 3) & 3),
                    (v2 << 2) | ((s2 + 2) & 3),
                    (v2 << 2) | ((s2 + 3) & 3)
                };
                int ws[3] = {v0, v1, v2};
                int new_wsvt[3] = {(hb > hc) ? 0 : 1, (hc > ha) ? 0 : 1, (ha > hb) ? 0 : 1};
                int int_edges[3][2] = {
                    {4 * ws[0] + 0, 4 * ws[1] + 1},
                    {4 * ws[1] + 0, 4 * ws[2] + 1},
                    {4 * ws[2] + 0, 4 * ws[0] + 1}
                };
                int port_assign[6][2] = {
                    {4 * ws[0] + 2, 1},
                    {4 * ws[0] + 3, 2},
                    {4 * ws[1] + 2, 3},
                    {4 * ws[1] + 3, 4},
                    {4 * ws[2] + 2, 5},
                    {4 * ws[2] + 3, 0}
                };
                rewire_ports_c(n, opp, vt, 3, ws, new_wsvt, 3, int_edges, 6, ports, port_assign, out_opps[cnt], out_vts[cnt]);
                cnt++;
            } else if (m_cnt == 1) {
                int mv = vs[m_idx], cv1 = vs[(m_idx + 1) % 3], cv2 = vs[(m_idx + 2) % 3];
                int ms = ss[m_idx], cs1 = ss[(m_idx + 1) % 3], cs2 = ss[(m_idx + 2) % 3];
                int ct2 = (cs2 + 1) & 3;
                int b_over1 = c_is_over(vt[cv1], cs1);
                int b_over2 = c_is_over(vt[cv2], ct2);
                if (b_over1 != b_over2) continue;
                if (cnt >= max_nbrs) return cnt;

                int ports[6] = {
                    (mv << 2) | ((ms + 2) & 3),
                    (mv << 2) | ((ms + 3) & 3),
                    (cv1 << 2) | ((cs1 + 2) & 3),
                    (cv1 << 2) | ((cs1 + 3) & 3),
                    (cv2 << 2) | ((cs2 + 2) & 3),
                    (cv2 << 2) | ((cs2 + 3) & 3)
                };
                int vt_m = vt[mv];
                int lplus_pairs_ms_ms3 = ((vt_m == 2 && (ms & 1) == 0) || (vt_m == 3 && (ms & 1) == 1)) ? 1 : 0;
                int ws[3] = {cv2, mv, cv1};
                int new_wsvt[3] = {b_over1 ? 0 : 1, lplus_pairs_ms_ms3 ? 2 : 3, b_over1 ? 1 : 0};
                int int_edges[3][2] = {
                    {4 * ws[0] + 0, 4 * ws[1] + 1},
                    {4 * ws[1] + 0, 4 * ws[2] + 1},
                    {4 * ws[2] + 0, 4 * ws[0] + 1}
                };
                int port_assign[6][2] = {
                    {4 * ws[0] + 2, 1},
                    {4 * ws[0] + 3, 2},
                    {4 * ws[1] + 2, 3},
                    {4 * ws[1] + 3, 4},
                    {4 * ws[2] + 2, 5},
                    {4 * ws[2] + 3, 0}
                };
                rewire_ports_c(n, opp, vt, 3, ws, new_wsvt, 3, int_edges, 6, ports, port_assign, out_opps[cnt], out_vts[cnt]);
                cnt++;
            }
        }
    }

    /* 4. Omega_8 (band swim) */
    if (n >= 6) {
        for (int start_d = 0; start_d < 4 * n; start_d++) {
            int vs[4], ss[4];
            int cur = start_d;
            int ok = 1;
            for (int k = 0; k < 4; k++) {
                int v = cur >> 2, s = cur & 3;
                if (vt[v] >= 2) { ok = 0; break; }
                for (int j = 0; j < k; j++) {
                    if (vs[j] == v) { ok = 0; break; }
                }
                if (!ok) break;
                vs[k] = v;
                ss[k] = s;
                int od = opp[cur];
                cur = ((od >> 2) << 2) | (((od & 3) - 1) & 3);
            }
            if (!ok || cur != start_d) continue;

            int c_sw = vs[0], c_se = vs[1], c_ne = vs[2], c_nw = vs[3];
            int s_sw = ss[0], s_se = ss[1], s_ne = ss[2], s_nw = ss[3];

            int d_sw_s = (c_sw << 2) | ((s_sw + 3) & 3);
            int d_se_s = (c_se << 2) | ((s_se + 2) & 3);
            int od_sw_s = opp[d_sw_s];
            int od_se_s = opp[d_se_s];
            int w = od_sw_s >> 2;
            if ((od_se_s >> 2) != w || vt[w] < 2) continue;
            int sw_w = od_se_s & 3;
            if (((sw_w + 1) & 3) != (od_sw_s & 3)) continue;

            int hb0 = c_is_over(vt[c_sw], s_sw);
            int hb1 = c_is_over(vt[c_se], (s_se + 1) & 3);
            int ht0 = c_is_over(vt[c_ne], s_ne);
            int ht1 = c_is_over(vt[c_nw], (s_nw + 1) & 3);
            if (hb0 != hb1 || ht0 != ht1 || hb0 == ht0) continue;

            int w_lplus_thru = ((vt[w] == 2 && (sw_w & 1) == 0) || (vt[w] == 3 && (sw_w & 1) == 1)) ? 1 : 0;
            int mk = w_lplus_thru ? 2 : 3;
            int cb = hb0 ? 0 : 1;
            int ct = hb0 ? 1 : 0;

            /* Case A: u on West -> move u to East (same level) or w to North (opposite level) */
            int d_sw_w = (c_sw << 2) | ((s_sw + 2) & 3);
            int d_nw_w = (c_nw << 2) | ((s_nw + 3) & 3);
            int od_sw_w = opp[d_sw_w];
            int od_nw_w = opp[d_nw_w];
            int u = od_sw_w >> 2;
            if (u != w && (od_nw_w >> 2) == u && vt[u] >= 2) {
                int su_u = od_sw_w & 3;
                if (((su_u + 1) & 3) == (od_nw_w & 3)) {
                    if (cnt >= max_nbrs) return cnt;
                    int u_lplus_thru = ((vt[u] == 2 && (su_u & 1) == 0) || (vt[u] == 3 && (su_u & 1) == 1)) ? 1 : 0;
                    int ports[8] = {
                        (w << 2) | ((sw_w + 2) & 3),
                        (w << 2) | ((sw_w + 3) & 3),
                        (c_se << 2) | ((s_se + 3) & 3),
                        (c_ne << 2) | ((s_ne + 2) & 3),
                        (c_ne << 2) | ((s_ne + 3) & 3),
                        (c_nw << 2) | ((s_nw + 2) & 3),
                        (u << 2) | ((su_u + 2) & 3),
                        (u << 2) | ((su_u + 3) & 3)
                    };
                    int ws[6] = {u, w, c_sw, c_se, c_ne, c_nw};
                    if (w_lplus_thru == u_lplus_thru) {
                        int new_wsvt[6] = {mk, mk, cb, cb, ct, ct};
                        int int_edges[8][2] = {
                            {4 * c_sw + 1, 4 * c_se + 3},
                            {4 * c_se + 2, 4 * c_ne + 0},
                            {4 * c_ne + 3, 4 * c_nw + 1},
                            {4 * c_nw + 0, 4 * c_sw + 2},
                            {4 * w + 3, 4 * c_sw + 0},
                            {4 * w + 2, 4 * c_se + 0},
                            {4 * u + 3, 4 * c_se + 1},
                            {4 * u + 2, 4 * c_ne + 1}
                        };
                        int port_assign[8][2] = {
                            {4 * w + 0, 0},
                            {4 * w + 1, 1},
                            {4 * u + 0, 2},
                            {4 * u + 1, 3},
                            {4 * c_ne + 2, 4},
                            {4 * c_nw + 2, 5},
                            {4 * c_nw + 3, 6},
                            {4 * c_sw + 3, 7}
                        };
                        rewire_ports_c(n, opp, vt, 6, ws, new_wsvt, 8, int_edges, 8, ports, port_assign, out_opps[cnt], out_vts[cnt]);
                        cnt++;
                    } else {
                        int mk_u = u_lplus_thru ? 2 : 3;
                        int new_wsvt[6] = {mk_u, mk, cb, cb, ct, ct};
                        int int_edges[8][2] = {
                            {4 * c_sw + 1, 4 * c_se + 3},
                            {4 * c_se + 2, 4 * c_ne + 0},
                            {4 * c_ne + 3, 4 * c_nw + 1},
                            {4 * c_nw + 0, 4 * c_sw + 2},
                            {4 * w + 1, 4 * c_ne + 2},
                            {4 * w + 0, 4 * c_nw + 2},
                            {4 * u + 0, 4 * c_sw + 3},
                            {4 * u + 1, 4 * c_nw + 3}
                        };
                        int port_assign[8][2] = {
                            {4 * c_sw + 0, 0},
                            {4 * c_se + 0, 1},
                            {4 * c_se + 1, 2},
                            {4 * c_ne + 1, 3},
                            {4 * w + 2, 4},
                            {4 * w + 3, 5},
                            {4 * u + 2, 6},
                            {4 * u + 3, 7}
                        };
                        rewire_ports_c(n, opp, vt, 6, ws, new_wsvt, 8, int_edges, 8, ports, port_assign, out_opps[cnt], out_vts[cnt]);
                        cnt++;
                    }
                }
            }

            /* Case B: u2 on East -> move u2 to West (same level) or w to North (opposite level) */
            int d_se_e = (c_se << 2) | ((s_se + 3) & 3);
            int d_ne_e = (c_ne << 2) | ((s_ne + 2) & 3);
            int od_se_e = opp[d_se_e];
            int od_ne_e = opp[d_ne_e];
            int u2 = od_ne_e >> 2;
            if (u2 != w && (od_se_e >> 2) == u2 && vt[u2] >= 2) {
                int su_u2 = od_ne_e & 3;
                if (((su_u2 + 1) & 3) == (od_se_e & 3)) {
                    if (cnt >= max_nbrs) return cnt;
                    int u2_lplus_thru = ((vt[u2] == 2 && (su_u2 & 1) == 0) || (vt[u2] == 3 && (su_u2 & 1) == 1)) ? 1 : 0;
                    int ports[8] = {
                        (w << 2) | ((sw_w + 2) & 3),
                        (w << 2) | ((sw_w + 3) & 3),
                        (u2 << 2) | ((su_u2 + 2) & 3),
                        (u2 << 2) | ((su_u2 + 3) & 3),
                        (c_ne << 2) | ((s_ne + 3) & 3),
                        (c_nw << 2) | ((s_nw + 2) & 3),
                        (c_nw << 2) | ((s_nw + 3) & 3),
                        (c_sw << 2) | ((s_sw + 2) & 3)
                    };
                    int ws[6] = {u2, w, c_sw, c_se, c_ne, c_nw};
                    if (w_lplus_thru == u2_lplus_thru) {
                        int new_wsvt[6] = {mk, mk, cb, cb, ct, ct};
                        int int_edges[8][2] = {
                            {4 * c_sw + 1, 4 * c_se + 3},
                            {4 * c_se + 2, 4 * c_ne + 0},
                            {4 * c_ne + 3, 4 * c_nw + 1},
                            {4 * c_nw + 0, 4 * c_sw + 2},
                            {4 * w + 3, 4 * c_sw + 0},
                            {4 * w + 2, 4 * c_se + 0},
                            {4 * u2 + 0, 4 * c_sw + 3},
                            {4 * u2 + 1, 4 * c_nw + 3}
                        };
                        int port_assign[8][2] = {
                            {4 * w + 0, 0},
                            {4 * w + 1, 1},
                            {4 * c_se + 1, 2},
                            {4 * c_ne + 1, 3},
                            {4 * c_ne + 2, 4},
                            {4 * c_nw + 2, 5},
                            {4 * u2 + 2, 6},
                            {4 * u2 + 3, 7}
                        };
                        rewire_ports_c(n, opp, vt, 6, ws, new_wsvt, 8, int_edges, 8, ports, port_assign, out_opps[cnt], out_vts[cnt]);
                        cnt++;
                    } else {
                        int mk_u2 = u2_lplus_thru ? 2 : 3;
                        int new_wsvt[6] = {mk_u2, mk, cb, cb, ct, ct};
                        int int_edges[8][2] = {
                            {4 * c_sw + 1, 4 * c_se + 3},
                            {4 * c_se + 2, 4 * c_ne + 0},
                            {4 * c_ne + 3, 4 * c_nw + 1},
                            {4 * c_nw + 0, 4 * c_sw + 2},
                            {4 * w + 1, 4 * c_ne + 2},
                            {4 * w + 0, 4 * c_nw + 2},
                            {4 * u2 + 3, 4 * c_se + 1},
                            {4 * u2 + 2, 4 * c_ne + 1}
                        };
                        int port_assign[8][2] = {
                            {4 * c_sw + 0, 0},
                            {4 * c_se + 0, 1},
                            {4 * u2 + 0, 2},
                            {4 * u2 + 1, 3},
                            {4 * w + 2, 4},
                            {4 * w + 3, 5},
                            {4 * c_nw + 3, 6},
                            {4 * c_sw + 3, 7}
                        };
                        rewire_ports_c(n, opp, vt, 6, ws, new_wsvt, 8, int_edges, 8, ports, port_assign, out_opps[cnt], out_vts[cnt]);
                        cnt++;
                    }
                }
            }
        }
    }

    /* 5. Marked 2-gon tube involution (certified +3-vertex Omega_6+ / Omega_2+ / Omega_7 composite) */
    for (int u = 0; u < n; u++) {
        if (vt[u] < 2) continue;
        for (int su = 0; su < 4; su++) {
            int od1 = opp[4 * u + su];
            int w = od1 >> 2;
            if (w <= u || vt[w] < 2) continue;
            int sw = od1 & 3;
            int su2 = (su + 1) & 3;
            int sw2 = (sw - 1) & 3;
            if (opp[4 * u + su2] != ((w << 2) | sw2)) continue;
            if (((vt[u] == vt[w]) ? 1 : 0) != (((su & 1) == (sw & 1)) ? 1 : 0)) continue;
            if (cnt >= max_nbrs) return cnt;
            memcpy(out_opps[cnt], opp, sizeof(int) * 4 * n);
            memcpy(out_vts[cnt], vt, sizeof(int) * n);
            out_vts[cnt][u] = (vt[u] == 2) ? 3 : 2;
            out_vts[cnt][w] = (vt[w] == 2) ? 3 : 2;
            cnt++;
        }
    }
    return cnt;
}

typedef struct {
    uint8_t used;
    CanonKey key;
} HashSlot;

typedef struct {
    size_t cap;
    size_t count;
    HashSlot *slots;
} CanonHashSet;

static void hashset_init(CanonHashSet *hs, size_t initial_cap_pow2) {
    hs->cap = initial_cap_pow2;
    hs->count = 0;
    hs->slots = (HashSlot *)calloc(hs->cap, sizeof(HashSlot));
}

static void hashset_free(CanonHashSet *hs) {
    free(hs->slots);
    hs->slots = NULL;
    hs->cap = 0;
    hs->count = 0;
}

/* Returns 0 if newly inserted (with `init_status`), else returns existing `used` status (1 or 2) */
static int hashset_insert_or_get(CanonHashSet *hs, int n, const CanonKey *k, uint8_t init_status) {
    if ((hs->count + 1) * 2 >= hs->cap) {
        size_t old_cap = hs->cap;
        HashSlot *old_slots = hs->slots;
        hs->cap = old_cap * 2;
        hs->slots = (HashSlot *)calloc(hs->cap, sizeof(HashSlot));
        size_t mask = hs->cap - 1;
        for (size_t i = 0; i < old_cap; i++) {
            if (old_slots[i].used) {
                size_t idx = (size_t)canon_key_hash(n, &old_slots[i].key) & mask;
                while (hs->slots[idx].used) {
                    idx = (idx + 1) & mask;
                }
                hs->slots[idx] = old_slots[i];
            }
        }
        free(old_slots);
    }

    size_t mask = hs->cap - 1;
    size_t idx = (size_t)canon_key_hash(n, k) & mask;
    while (hs->slots[idx].used) {
        if (canon_key_eq(n, &hs->slots[idx].key, k)) {
            if (init_status == 2) hs->slots[idx].used = 2;
            return (int)hs->slots[idx].used;
        }
        idx = (idx + 1) & mask;
    }
    hs->slots[idx].used = init_status;
    hs->slots[idx].key = *k;
    hs->count++;
    return 0;
}

static void hashset_mark_discarded(CanonHashSet *hs, int n, const CanonKey *k) {
    size_t mask = hs->cap - 1;
    size_t idx = (size_t)canon_key_hash(n, k) & mask;
    while (hs->slots[idx].used) {
        if (canon_key_eq(n, &hs->slots[idx].key, k)) {
            hs->slots[idx].used = 2;
            return;
        }
        idx = (idx + 1) & mask;
    }
}

static int hashset_lookup(const CanonHashSet *hs, int n, const CanonKey *k) {
    size_t mask = hs->cap - 1;
    size_t idx = (size_t)canon_key_hash(n, k) & mask;
    while (hs->slots[idx].used) {
        if (canon_key_eq(n, &hs->slots[idx].key, k)) {
            return (int)hs->slots[idx].used;
        }
        idx = (idx + 1) & mask;
    }
    return 0;
}

static int explore_single_component_c(
    int n,
    const int *start_opp,
    const int *start_vt,
    const CanonKey *k0,
    int max_bfs_states,
    int start_is_prime_shadow,
    CanonHashSet *visited,
    int (*q_opp)[4 * MAX_N],
    int (*q_vt)[MAX_N],
    CanonKey *q_keys,
    int *best_opp_out,
    int *best_vt_out
) {
    int nbr_opps[128][4 * MAX_N];
    int nbr_vts[128][MAX_N];

    int discarded = 0;
    if (!start_is_prime_shadow) {
        if (has_immediate_reduction(n, start_opp, start_vt)) {
            discarded = 1;
        } else if (n >= 4 && has_0_or_2_edge_cut_c(n, start_opp)) {
            discarded = 1;
        }
    }

    hashset_insert_or_get(visited, n, k0, discarded ? 2 : 1);
    if (discarded) {
        return 1;
    }

    CanonKey best_k = *k0;
    memcpy(best_opp_out, start_opp, sizeof(int) * 4 * n);
    memcpy(best_vt_out, start_vt, sizeof(int) * n);

    memcpy(q_opp[0], start_opp, sizeof(int) * 4 * n);
    memcpy(q_vt[0], start_vt, sizeof(int) * n);
    q_keys[0] = *k0;
    int qh = 0, qt = 1;

    while (qh < qt && !discarded) {
        const int *cur_opp = q_opp[qh];
        const int *cur_vt = q_vt[qh];
        qh++;

        int nnbr = list_all_level_neighbors_c(n, cur_opp, cur_vt, nbr_opps, nbr_vts, 128);
        for (int ni = 0; ni < nnbr; ni++) {
            if (has_immediate_reduction(n, nbr_opps[ni], nbr_vts[ni])) {
                discarded = 1;
                break;
            }
        }
        if (discarded) break;

        for (int ni = 0; ni < nnbr; ni++) {
            if (n >= 4) {
                if (has_0_or_2_edge_cut_c(n, nbr_opps[ni])) {
                    discarded = 1;
                    break;
                }
            } else {
                uint32_t cvis = 1u;
                int cq[MAX_N], cqh = 0, cqt = 0;
                cq[cqt++] = 0;
                while (cqh < cqt) {
                    int u = cq[cqh++];
                    for (int s = 0; s < 4; s++) {
                        int w = nbr_opps[ni][4 * u + s] >> 2;
                        if (!((cvis >> w) & 1u)) {
                            cvis |= (1u << w);
                            cq[cqt++] = w;
                        }
                    }
                }
                if (cqt < n) {
                    discarded = 1;
                    break;
                }
            }
        }
        if (discarded) break;

        for (int ni = 0; ni < nnbr; ni++) {
            CanonKey nk;
            canonical_sym_key_c(n, nbr_opps[ni], nbr_vts[ni], &nk);
            int prev_st = hashset_lookup(visited, n, &nk);
            if (prev_st != 0) {
                if (prev_st == 2) {
                    discarded = 1;
                    break;
                }
                continue;
            }

            if (canon_key_less(n, &nk, &best_k)) {
                best_k = nk;
                memcpy(best_opp_out, nbr_opps[ni], sizeof(int) * 4 * n);
                memcpy(best_vt_out, nbr_vts[ni], sizeof(int) * n);
            }
            if (qt < max_bfs_states) {
                hashset_insert_or_get(visited, n, &nk, 1);
                memcpy(q_opp[qt], nbr_opps[ni], sizeof(int) * 4 * n);
                memcpy(q_vt[qt], nbr_vts[ni], sizeof(int) * n);
                q_keys[qt] = nk;
                qt++;
            }
        }
    }

    if (discarded) {
        for (int i = 0; i < qt; i++) {
            hashset_mark_discarded(visited, n, &q_keys[i]);
        }
        return 1;
    }
    return 0;
}

int enumerate_census_level_c(
    int n,
    int num_shadows,
    const int *all_shadows_opp,
    int m_min,
    int m_max,
    int max_bfs_states,
    int *out_opp,
    int *out_vt,
    int max_kept,
    int64_t *out_stats /* [0]=raw_cands, [1]=components, [2]=discarded, [3]=kept */
) {
    CanonHashSet visited;
    hashset_init(&visited, 1ULL << 20);

    size_t cand_cap = 65536;
    CandidateDiagram *cand_buf = (CandidateDiagram *)malloc(cand_cap * sizeof(CandidateDiagram));

    int (*q_opp)[4 * MAX_N] = malloc((size_t)max_bfs_states * sizeof(*q_opp));
    int (*q_vt)[MAX_N] = malloc((size_t)max_bfs_states * sizeof(*q_vt));
    CanonKey *q_keys = (CanonKey *)malloc((size_t)max_bfs_states * sizeof(CanonKey));

    int64_t total_cands = 0;
    int64_t total_comps = 0;
    int64_t total_discarded = 0;
    int kept_cnt = 0;

    for (int sh_idx = 0; sh_idx < num_shadows; sh_idx++) {
        const int *sh_opp = all_shadows_opp + sh_idx * (4 * n);
        int cnt = enumerate_shadow_candidates(n, sh_opp, m_min, m_max, 3, cand_buf, (int)cand_cap);
        if (cnt > (int)cand_cap) {
            cand_cap = (size_t)cnt + 1024;
            free(cand_buf);
            cand_buf = (CandidateDiagram *)malloc(cand_cap * sizeof(CandidateDiagram));
            cnt = enumerate_shadow_candidates(n, sh_opp, m_min, m_max, 3, cand_buf, (int)cand_cap);
        }
        total_cands += cnt;

        for (int ci = 0; ci < cnt; ci++) {
            uint16_t m_mask = cand_buf[ci].m_mask;
            uint16_t yz_mask = cand_buf[ci].yz_mask;
            uint16_t smask = cand_buf[ci].sign_mask;

            int d0_vt[MAX_N];
            int mi = 0, ki = 0;
            for (int v = 0; v < n; v++) {
                if ((m_mask >> v) & 1) {
                    d0_vt[v] = ((yz_mask >> mi) & 1) ? 3 : 2;
                    mi++;
                } else {
                    d0_vt[v] = ((smask >> ki) & 1) ? 1 : 0;
                    ki++;
                }
            }

            CanonKey k0;
            canonical_sym_key_c(n, sh_opp, d0_vt, &k0);
            if (hashset_lookup(&visited, n, &k0) != 0) {
                continue;
            }

            total_comps++;
            int best_opp[4 * MAX_N];
            int best_vt[MAX_N];
            int disc = explore_single_component_c(
                n, sh_opp, d0_vt, &k0, max_bfs_states, 1, &visited, q_opp, q_vt, q_keys, best_opp, best_vt
            );
            if (disc) {
                total_discarded++;
            } else {
                if (kept_cnt < max_kept) {
                    memcpy(out_opp + kept_cnt * (4 * n), best_opp, sizeof(int) * 4 * n);
                    memcpy(out_vt + kept_cnt * n, best_vt, sizeof(int) * n);
                }
                kept_cnt++;
            }
        }
    }

    out_stats[0] = total_cands;
    out_stats[1] = total_comps;
    out_stats[2] = total_discarded;
    out_stats[3] = kept_cnt;

    free(q_keys);
    free(q_opp);
    free(q_vt);
    free(cand_buf);
    hashset_free(&visited);
    return kept_cnt;
}

int merge_and_filter_diagrams_c(
    int n,
    int num_in,
    const int *in_opp,
    const int *in_vt,
    int max_bfs_states,
    int *out_opp,
    int *out_vt,
    int max_kept
) {
    CanonHashSet visited;
    hashset_init(&visited, 1ULL << 18);

    int (*q_opp)[4 * MAX_N] = malloc((size_t)max_bfs_states * sizeof(*q_opp));
    int (*q_vt)[MAX_N] = malloc((size_t)max_bfs_states * sizeof(*q_vt));
    CanonKey *q_keys = (CanonKey *)malloc((size_t)max_bfs_states * sizeof(CanonKey));

    int kept_cnt = 0;
    for (int i = 0; i < num_in; i++) {
        const int *d_opp = in_opp + i * (4 * n);
        const int *d_vt = in_vt + i * n;
        CanonKey k0;
        canonical_sym_key_c(n, d_opp, d_vt, &k0);
        if (hashset_lookup(&visited, n, &k0) != 0) {
            continue;
        }
        int best_opp[4 * MAX_N];
        int best_vt[MAX_N];
        int disc = explore_single_component_c(
            n, d_opp, d_vt, &k0, max_bfs_states, 0, &visited, q_opp, q_vt, q_keys, best_opp, best_vt
        );
        if (!disc) {
            if (kept_cnt < max_kept) {
                memcpy(out_opp + kept_cnt * (4 * n), best_opp, sizeof(int) * 4 * n);
                memcpy(out_vt + kept_cnt * n, best_vt, sizeof(int) * n);
            }
            kept_cnt++;
        }
    }

    free(q_keys);
    free(q_opp);
    free(q_vt);
    hashset_free(&visited);
    return kept_cnt;
}

/* =========================================================================
 * Greedy Reduction & +1/+2 Stabilization Mountain-Pass Tunneling in C
 * ========================================================================= */

static int remove_spliced_vertices_c(
    int n,
    const int *opp,
    const int *vt,
    uint32_t rem_mask,
    const int *internal_splice, /* -1 if dart not spliced */
    int *out_opp,
    int *out_vt,
    int *extra_loops
) {
    int o2n[MAX_N];
    int m = 0;
    for (int v = 0; v < n; v++) {
        if ((rem_mask >> v) & 1u) {
            o2n[v] = -1;
        } else {
            o2n[v] = m;
            out_vt[m] = vt[v];
            m++;
        }
    }

    uint8_t vis[4 * MAX_N];
    memset(vis, 0, 4 * n);
    int new_opp[4 * MAX_N];

    for (int v = 0; v < n; v++) {
        int vn = o2n[v];
        if (vn < 0) continue;
        for (int s = 0; s < 4; s++) {
            int cur = opp[4 * v + s];
            while (internal_splice[cur] >= 0) {
                vis[cur] = 1;
                int nxt = internal_splice[cur];
                vis[nxt] = 1;
                cur = opp[nxt];
            }
            new_opp[4 * vn + s] = (o2n[cur >> 2] << 2) | (cur & 3);
        }
    }

    int loops = 0;
    for (int d = 0; d < 4 * n; d++) {
        if (internal_splice[d] < 0 || vis[d]) continue;
        loops++;
        int cur = d;
        while (internal_splice[cur] >= 0 && !vis[cur]) {
            vis[cur] = 1;
            int nxt = internal_splice[cur];
            vis[nxt] = 1;
            cur = opp[nxt];
        }
    }

    memcpy(out_opp, new_opp, sizeof(int) * 4 * m);
    *extra_loops += loops;
    return m;
}

static int greedy_reduce_c(
    int n_in,
    const int *opp_in,
    const int *vt_in,
    int *opp_out,
    int *vt_out,
    int *extra_loops_out
) {
    int cur_n = n_in;
    int cur_opp[4 * MAX_N];
    int cur_vt[MAX_N];
    memcpy(cur_opp, opp_in, sizeof(int) * 4 * n_in);
    memcpy(cur_vt, vt_in, sizeof(int) * n_in);
    int loops = 0;
    int splice[4 * MAX_N];

    while (cur_n > 0) {
        /* 1. Omega_1- / Omega_6- */
        int found = 0;
        for (int v = 0; v < cur_n && !found; v++) {
            for (int s = 0; s < 4; s++) {
                int s_adj = (s + 1) & 3;
                if (cur_opp[4 * v + s] == 4 * v + s_adj) {
                    for (int i = 0; i < 4 * cur_n; i++) splice[i] = -1;
                    int da = 4 * v + ((s + 2) & 3);
                    int db = 4 * v + ((s + 3) & 3);
                    splice[da] = db;
                    splice[db] = da;
                    cur_n = remove_spliced_vertices_c(cur_n, cur_opp, cur_vt, 1u << v, splice, cur_opp, cur_vt, &loops);
                    found = 1;
                    break;
                }
            }
        }
        if (found) continue;

        /* 2. Omega_2- */
        for (int u = 0; u < cur_n && !found; u++) {
            if (cur_vt[u] >= 2) continue;
            for (int su = 0; su < 4; su++) {
                int od1 = cur_opp[4 * u + su];
                int w = od1 >> 2;
                if (w <= u || cur_vt[w] >= 2) continue;
                int sw = od1 & 3;
                int su2 = (su + 1) & 3;
                int sw2 = (sw - 1) & 3;
                if (cur_opp[4 * u + su2] != ((w << 2) | sw2)) continue;
                if (c_is_over(cur_vt[u], su) != c_is_over(cur_vt[w], sw)) continue;

                for (int i = 0; i < 4 * cur_n; i++) splice[i] = -1;
                int ue1 = (u << 2) | (su ^ 2), we1 = (w << 2) | (sw ^ 2);
                int ue2 = (u << 2) | (su2 ^ 2), we2 = (w << 2) | (sw2 ^ 2);
                splice[ue1] = we1; splice[we1] = ue1;
                splice[ue2] = we2; splice[we2] = ue2;
                cur_n = remove_spliced_vertices_c(cur_n, cur_opp, cur_vt, (1u << u) | (1u << w), splice, cur_opp, cur_vt, &loops);
                found = 1;
                break;
            }
        }
        if (!found) break;
    }

    memcpy(opp_out, cur_opp, sizeof(int) * 4 * cur_n);
    memcpy(vt_out, cur_vt, sizeof(int) * cur_n);
    *extra_loops_out = loops;
    return cur_n;
}

static int is_connected_c(int n, const int *opp) {
    if (n <= 1) return 1;
    uint32_t vis = 1u;
    int q[MAX_N], qh = 0, qt = 0;
    q[qt++] = 0;
    while (qh < qt) {
        int u = q[qh++];
        for (int s = 0; s < 4; s++) {
            int w = opp[4 * u + s] >> 2;
            if (!((vis >> w) & 1u)) {
                vis |= (1u << w);
                q[qt++] = w;
            }
        }
    }
    return (qt == n) ? 1 : 0;
}

/* Check if an (n+2)-vertex 1-gon-free connected diagram has a certifying 2-edge cut:
 * - Returns 1 if certified discard (3 <= |A| <= n-1, OR |A|==2 with nontrivial 2-vertex factor,
 *   OR |A|==2 with 0_1 factor and D_B reducing/cutting).
 * - If |A|==2 with 0_1 factor and valid n-vertex D_B, writes D_B to (bridge_opp, bridge_vt) and sets *has_bridge = 1!
 */
static int check_np2_cuts_c(
    int n,
    const int *opp,
    const int *vt,
    int *bridge_opp,
    int *bridge_vt,
    int *has_bridge
) {
    int m = n + 2;
    *has_bridge = 0;
    int e_a[2 * MAX_N], e_b[2 * MAX_N], ne = 0;
    for (int d = 0; d < 4 * m; d++) {
        int od = opp[d];
        if (d < od) {
            e_a[ne] = d;
            e_b[ne] = od;
            ne++;
        }
    }
    for (int i = 0; i < ne; i++) {
        for (int j = i + 1; j < ne; j++) {
            int b1 = e_a[i], b2 = e_b[i], b3 = e_a[j], b4 = e_b[j];
            uint32_t vis = 1u;
            int q[MAX_N], qh = 0, qt = 0;
            q[qt++] = 0;
            while (qh < qt) {
                int u = q[qh++];
                for (int s = 0; s < 4; s++) {
                    int dd = (u << 2) | s;
                    if (dd == b1 || dd == b2 || dd == b3 || dd == b4) continue;
                    int w = opp[dd] >> 2;
                    if (!((vis >> w) & 1u)) {
                        vis |= (1u << w);
                        q[qt++] = w;
                    }
                }
            }
            if (qt == m) continue;
            if (qt >= 3 && qt <= m - 3) {
                return 1; /* Both factors have <= n - 1 vertices! */
            }
            if (qt == 2 || qt == n) {
                /* Ensure both e_i and e_j actually cross between vis and !vis */
                int i_in1 = (vis >> (b1 >> 2)) & 1u, i_in2 = (vis >> (b2 >> 2)) & 1u;
                int j_in1 = (vis >> (b3 >> 2)) & 1u, j_in2 = (vis >> (b4 >> 2)) & 1u;
                if (i_in1 == i_in2 || j_in1 == j_in2) continue;

                int small_flag = (qt == 2) ? 1 : 0;
                int cap_s1 = (((vis >> (b1 >> 2)) & 1u) == (uint32_t)small_flag) ? b1 : b2;
                int cap_s2 = (((vis >> (b3 >> 2)) & 1u) == (uint32_t)small_flag) ? b3 : b4;
                int cap_l1 = (cap_s1 == b1) ? b2 : b1;
                int cap_l2 = (cap_s2 == b3) ? b4 : b3;

                /* Build 2-vertex factor D_2 */
                int o2n_s[MAX_N], idx_s = 0;
                int s_vt[2], s_opp[8];
                for (int v = 0; v < m; v++) {
                    if ((int)((vis >> v) & 1u) == small_flag) {
                        o2n_s[v] = idx_s;
                        s_vt[idx_s] = vt[v];
                        idx_s++;
                    } else {
                        o2n_s[v] = -1;
                    }
                }
                for (int v = 0; v < m; v++) {
                    int nv = o2n_s[v];
                    if (nv < 0) continue;
                    for (int s = 0; s < 4; s++) {
                        int od_orig = (v << 2) | s;
                        int tgt = (od_orig == cap_s1) ? cap_s2 : ((od_orig == cap_s2) ? cap_s1 : opp[od_orig]);
                        s_opp[4 * nv + s] = (o2n_s[tgt >> 2] << 2) | (tgt & 3);
                    }
                }
                int r2_opp[8], r2_vt[2], r2_loops = 0;
                int n2 = greedy_reduce_c(2, s_opp, s_vt, r2_opp, r2_vt, &r2_loops);
                int is_trivial_s2 = (n2 == 0 && r2_loops == 1) ? 1 : 0;
                if (!is_trivial_s2 && n2 == 2 && r2_loops == 0) {
                    int s2_nbr_o[16][4 * MAX_N];
                    int s2_nbr_v[16][MAX_N];
                    int s2_cnt = list_all_level_neighbors_c(2, r2_opp, r2_vt, s2_nbr_o, s2_nbr_v, 16);
                    for (int si = 0; si < s2_cnt; si++) {
                        int sub_o[8], sub_v[2], sub_loops = 0;
                        int sub_n = greedy_reduce_c(2, s2_nbr_o[si], s2_nbr_v[si], sub_o, sub_v, &sub_loops);
                        if (sub_n == 0 && sub_loops == 1) {
                            is_trivial_s2 = 1;
                            break;
                        }
                    }
                }
                if (!is_trivial_s2) {
                    /* 2-vertex factor is NOT 0_1 (e.g. 2_1^1 or 2_1^-1, or 2-sphere unlink)! */
                    return 1;
                }

                /* D_2 == 0_1: build n-vertex factor D_n! */
                int o2n_l[MAX_N], idx_l = 0;
                int l_vt[MAX_N], l_opp[4 * MAX_N];
                for (int v = 0; v < m; v++) {
                    if ((int)((vis >> v) & 1u) != small_flag) {
                        o2n_l[v] = idx_l;
                        l_vt[idx_l] = vt[v];
                        idx_l++;
                    } else {
                        o2n_l[v] = -1;
                    }
                }
                for (int v = 0; v < m; v++) {
                    int nv = o2n_l[v];
                    if (nv < 0) continue;
                    for (int s = 0; s < 4; s++) {
                        int od_orig = (v << 2) | s;
                        int tgt = (od_orig == cap_l1) ? cap_l2 : ((od_orig == cap_l2) ? cap_l1 : opp[od_orig]);
                        l_opp[4 * nv + s] = (o2n_l[tgt >> 2] << 2) | (tgt & 3);
                    }
                }
                int rn_opp[4 * MAX_N], rn_vt[MAX_N], rn_loops = 0;
                int nn = greedy_reduce_c(n, l_opp, l_vt, rn_opp, rn_vt, &rn_loops);
                if (nn < n || rn_loops > 0 || !is_connected_c(n, rn_opp) || (n >= 4 && has_0_or_2_edge_cut_c(n, rn_opp))) {
                    return 1;
                }
                if (!(*has_bridge)) {
                    memcpy(bridge_opp, rn_opp, sizeof(int) * 4 * n);
                    memcpy(bridge_vt, rn_vt, sizeof(int) * n);
                    *has_bridge = 1;
                }
            }
        }
    }
    return 0;
}

/* Evaluate one stabilized state `(m_stab, stab_opp, stab_vt)` (`n+1` or `n+2`) via a deduplicated
 * bounded BFS of up to `state_budget` states (max 480):
 * - Returns 1 if `D` is certified non-minimal (`< n` vertices), split, or composite!
 * - Appends any newly discovered `n`-vertex bridge neighbors to `bridges` (up to `max_bridges`).
 */
static void try_add_bridge_c(
    int n,
    const int *b_opp,
    const int *b_vt,
    CanonKey *bridge_keys,
    int bridge_opps[][4 * MAX_N],
    int bridge_vts[][MAX_N],
    int *num_bridges,
    int max_bridges
) {
    if (*num_bridges >= max_bridges) return;
    CanonKey bk;
    canonical_sym_key_c(n, b_opp, b_vt, &bk);
    for (int i = 0; i < *num_bridges; i++) {
        if (canon_key_eq(n, &bridge_keys[i], &bk)) return;
    }
    bridge_keys[*num_bridges] = bk;
    memcpy(bridge_opps[*num_bridges], b_opp, sizeof(int) * 4 * n);
    memcpy(bridge_vts[*num_bridges], b_vt, sizeof(int) * n);
    (*num_bridges)++;
}

static int probe_stabilized_state_c(
    int n,
    int m_stab,
    const int *stab_opp,
    const int *stab_vt,
    int depth_or_budget,
    int (*sq_opp)[4 * MAX_N],
    int (*sq_vt)[MAX_N],
    int *sq_n,
    CanonKey *bridge_keys,
    int bridge_opps[][4 * MAX_N],
    int bridge_vts[][MAX_N],
    int *num_bridges,
    int max_bridges
) {
    int budget;
    if (depth_or_budget <= 0) return 0;
    if (depth_or_budget == 1) budget = 24;
    else if (depth_or_budget == 2) budget = 64;
    else if (depth_or_budget < 16) budget = 160;
    else budget = (depth_or_budget > 480) ? 480 : depth_or_budget;

    uint8_t seen_used[1024];
    uint8_t seen_n[1024];
    CanonKey seen_keys[1024];
    memset(seen_used, 0, sizeof(seen_used));

    int nbr_opps[128][4 * MAX_N];
    int nbr_vts[128][MAX_N];
    int red_opp[4 * MAX_N], red_vt[MAX_N];
    int b_opp[4 * MAX_N], b_vt[MAX_N];

    memcpy(sq_opp[0], stab_opp, sizeof(int) * 4 * m_stab);
    memcpy(sq_vt[0], stab_vt, sizeof(int) * m_stab);
    sq_n[0] = m_stab;
    int sqh = 0, sqt = 1;

    while (sqh < sqt) {
        int cur_m = sq_n[sqh];
        const int *cur_opp = sq_opp[sqh];
        const int *cur_vt = sq_vt[sqh];
        sqh++;

        int cnt = list_all_level_neighbors_c(cur_m, cur_opp, cur_vt, nbr_opps, nbr_vts, 128);
        for (int i = 0; i < cnt; i++) {
            int loops = 0;
            int nr = greedy_reduce_c(cur_m, nbr_opps[i], nbr_vts[i], red_opp, red_vt, &loops);
            if (nr < n) return 1;
            if (loops > 0 || !is_connected_c(nr, red_opp)) return 1;

            if (nr == n) {
                if (n >= 4 && has_0_or_2_edge_cut_c(n, red_opp)) return 1;
                try_add_bridge_c(n, red_opp, red_vt, bridge_keys, bridge_opps, bridge_vts, num_bridges, max_bridges);
                continue;
            }

            /* nr is n + 1 or n + 2 (1-gon-free and connected): deduplicate before expensive cut checks! */
            CanonKey rk;
            canonical_sym_key_c(nr, red_opp, red_vt, &rk);
            size_t h_idx = ((size_t)canon_key_hash(nr, &rk) ^ ((size_t)nr * 1315423911ULL)) & 1023ULL;
            int dup = 0;
            while (seen_used[h_idx]) {
                if (seen_n[h_idx] == (uint8_t)nr && canon_key_eq(nr, &seen_keys[h_idx], &rk)) {
                    dup = 1;
                    break;
                }
                h_idx = (h_idx + 1) & 1023ULL;
            }
            if (dup) continue;

            if (nr == n + 1 && n >= 3 && has_0_or_2_edge_cut_c(nr, red_opp)) {
                return 1;
            }
            if (nr == n + 2 && n >= 3) {
                int hb = 0;
                if (check_np2_cuts_c(n, red_opp, red_vt, b_opp, b_vt, &hb)) return 1;
                if (hb) {
                    try_add_bridge_c(n, b_opp, b_vt, bridge_keys, bridge_opps, bridge_vts, num_bridges, max_bridges);
                }
            }

            if (sqt < budget) {
                seen_used[h_idx] = 1;
                seen_n[h_idx] = (uint8_t)nr;
                seen_keys[h_idx] = rk;
                memcpy(sq_opp[sqt], red_opp, sizeof(int) * 4 * nr);
                memcpy(sq_vt[sqt], red_vt, sizeof(int) * nr);
                sq_n[sqt] = nr;
                sqt++;
            }
        }
    }
    return 0;
}

static int run_single_probe_stabilizations_c(
    int n,
    const int *po,
    const int *pv,
    int depth_plus1,
    int depth_plus2,
    int (*sq_opp)[4 * MAX_N],
    int (*sq_vt)[MAX_N],
    int *sq_n,
    CanonKey *bridge_keys,
    int bridge_opps[][4 * MAX_N],
    int bridge_vts[][MAX_N],
    int *num_bridges,
    int max_bridges
) {
    if (depth_plus1 >= 1 && n + 1 <= MAX_N) {
        for (int d1 = 0; d1 < 4 * n; d1++) {
            int d2 = po[d1];
            if (d2 < d1) continue;
            for (int side = 0; side < 2; side++) {
                int s_opp[4 * MAX_N];
                int s_vt[MAX_N];
                memcpy(s_opp, po, sizeof(int) * 4 * n);
                memcpy(s_vt, pv, sizeof(int) * n);
                int u = n;
                if (side == 0) {
                    s_opp[4 * u + 0] = 4 * u + 1;
                    s_opp[4 * u + 1] = 4 * u + 0;
                    s_opp[4 * u + 2] = d1; s_opp[d1] = 4 * u + 2;
                    s_opp[4 * u + 3] = d2; s_opp[d2] = 4 * u + 3;
                } else {
                    s_opp[4 * u + 1] = 4 * u + 2;
                    s_opp[4 * u + 2] = 4 * u + 1;
                    s_opp[4 * u + 3] = d1; s_opp[d1] = 4 * u + 3;
                    s_opp[4 * u + 0] = d2; s_opp[d2] = 4 * u + 0;
                }
                for (int vt_u = 0; vt_u < 4; vt_u++) {
                    s_vt[u] = vt_u;
                    if (probe_stabilized_state_c(n, n + 1, s_opp, s_vt, depth_plus1, sq_opp, sq_vt, sq_n, bridge_keys, bridge_opps, bridge_vts, num_bridges, max_bridges)) {
                        return 1;
                    }
                }
            }
        }
    }

    if (depth_plus2 >= 1 && n + 2 <= MAX_N) {
        uint8_t vis_d[4 * MAX_N];
        memset(vis_d, 0, 4 * n);
        for (int sd = 0; sd < 4 * n; sd++) {
            if (vis_d[sd]) continue;
            int f_darts[4 * MAX_N], flen = 0;
            int cur = sd;
            while (!vis_d[cur]) {
                vis_d[cur] = 1;
                f_darts[flen++] = cur;
                int od = po[cur];
                cur = ((od >> 2) << 2) | (((od & 3) - 1) & 3);
            }
            if (flen < 2) continue;
            for (int a = 0; a < flen; a++) {
                int da = f_darts[a], oda = po[da];
                int ea = (da < oda) ? da : oda;
                for (int b = a + 1; b < flen; b++) {
                    int db = f_darts[b], odb = po[db];
                    int eb = (db < odb) ? db : odb;
                    if (ea == eb) continue;
                    int s_opp[4 * MAX_N];
                    int s_vt[MAX_N];
                    memcpy(s_opp, po, sizeof(int) * 4 * n);
                    memcpy(s_vt, pv, sizeof(int) * n);
                    int u = n, w = n + 1;
                    s_opp[4 * u + 0] = da;  s_opp[da] = 4 * u + 0;
                    s_opp[4 * u + 1] = 4 * w + 3;
                    s_opp[4 * u + 2] = 4 * w + 2;
                    s_opp[4 * u + 3] = odb; s_opp[odb] = 4 * u + 3;
                    s_opp[4 * w + 0] = oda; s_opp[oda] = 4 * w + 0;
                    s_opp[4 * w + 1] = db;  s_opp[db] = 4 * w + 1;
                    s_opp[4 * w + 2] = 4 * u + 2;
                    s_opp[4 * w + 3] = 4 * u + 1;
                    for (int vt_uw = 0; vt_uw < 2; vt_uw++) {
                        s_vt[u] = vt_uw;
                        s_vt[w] = vt_uw;
                        if (probe_stabilized_state_c(n, n + 2, s_opp, s_vt, depth_plus2, sq_opp, sq_vt, sq_n, bridge_keys, bridge_opps, bridge_vts, num_bridges, max_bridges)) {
                            return 1;
                        }
                    }
                }
            }
        }
    }
    return 0;
}

/* Run +1 (Omega_1+, Omega_6+) and +2 (Omega_2+) mountain-pass tunneling across `num_in` level-kept
 * representatives of size `n`:
 * - Discards any component that reduces to `< n` vertices, splits, or exposes a nontrivial 2-edge cut!
 * - Merges components connected via `+1/+2` mountain-pass bridges + level moves!
 */
int tunnel_and_merge_diagrams_c(
    int n,
    int num_in,
    const int *in_opp,
    const int *in_vt,
    int depth_plus1,   /* state budget or depth code (1->96, 2->160) after Omega_1+/6+ */
    int depth_plus2,   /* state budget or depth code (1->96, 2->160) after Omega_2+ */
    int max_bfs_states,
    int *out_opp,
    int *out_vt,
    int max_kept
) {
    if (num_in <= 0) return 0;
    if (n <= 2) {
        return merge_and_filter_diagrams_c(n, num_in, in_opp, in_vt, max_bfs_states, out_opp, out_vt, max_kept);
    }

    /* Step 1: Map every input representative `i` (`0..num_in-1`) and its level-BFS states to `i`! */
    CanonHashSet key_to_comp;
    hashset_init(&key_to_comp, 1ULL << 19);

    size_t q_alloc = (max_bfs_states < 512) ? 512 : (size_t)max_bfs_states;
    int (*q_opp)[4 * MAX_N] = malloc(q_alloc * sizeof(*q_opp));
    int (*q_vt)[MAX_N] = malloc(q_alloc * sizeof(*q_vt));
    int nbr_opps[128][4 * MAX_N];
    int nbr_vts[128][MAX_N];

    /* Disjoint-set union over `0..num_in-1`, plus `dead[root]` flag and best rep per root */
    int *parent = (int *)malloc((size_t)num_in * sizeof(int));
    int *dead = (int *)calloc((size_t)num_in, sizeof(int));
    CanonKey *best_keys = (CanonKey *)malloc((size_t)num_in * sizeof(CanonKey));
    int (*best_opps)[4 * MAX_N] = malloc((size_t)num_in * sizeof(*best_opps));
    int (*best_vts)[MAX_N] = malloc((size_t)num_in * sizeof(*best_vts));

    int *slot_owner = (int *)malloc(key_to_comp.cap * sizeof(int));
    for (size_t i = 0; i < key_to_comp.cap; i++) slot_owner[i] = -1;

    #define FIND_ROOT(x_var, out_r) do { \
        int _r = (x_var); \
        while (parent[_r] != _r) { parent[_r] = parent[parent[_r]]; _r = parent[_r]; } \
        (out_r) = _r; \
    } while (0)

    #define ENSURE_HASH_CAP() do { \
        if ((key_to_comp.count + 1) * 2 >= key_to_comp.cap) { \
            size_t old_cap = key_to_comp.cap; \
            HashSlot *old_s = key_to_comp.slots; \
            int *old_o = slot_owner; \
            key_to_comp.cap = old_cap * 2; \
            key_to_comp.slots = (HashSlot *)calloc(key_to_comp.cap, sizeof(HashSlot)); \
            slot_owner = (int *)malloc(key_to_comp.cap * sizeof(int)); \
            for (size_t _z = 0; _z < key_to_comp.cap; _z++) slot_owner[_z] = -1; \
            size_t _msk = key_to_comp.cap - 1; \
            for (size_t _z = 0; _z < old_cap; _z++) { \
                if (old_s[_z].used) { \
                    size_t _idx = (size_t)canon_key_hash(n, &old_s[_z].key) & _msk; \
                    while (key_to_comp.slots[_idx].used) _idx = (_idx + 1) & _msk; \
                    key_to_comp.slots[_idx] = old_s[_z]; \
                    slot_owner[_idx] = old_o[_z]; \
                } \
            } \
            free(old_s); \
            free(old_o); \
        } \
    } while (0)

    #define UNION_OWNERS(o1, o2) do { \
        int _r1, _r2; \
        FIND_ROOT((o1), _r1); \
        FIND_ROOT((o2), _r2); \
        if (_r1 != _r2) { \
            parent[_r1] = _r2; \
            dead[_r2] |= dead[_r1]; \
            if (canon_key_less(n, &best_keys[_r1], &best_keys[_r2])) { \
                best_keys[_r2] = best_keys[_r1]; \
                memcpy(best_opps[_r2], best_opps[_r1], sizeof(int) * 4 * n); \
                memcpy(best_vts[_r2], best_vts[_r1], sizeof(int) * n); \
            } \
        } \
    } while (0)

    for (int i = 0; i < num_in; i++) {
        parent[i] = i;
        const int *d_opp = in_opp + i * (4 * n);
        const int *d_vt = in_vt + i * n;
        canonical_sym_key_c(n, d_opp, d_vt, &best_keys[i]);
        memcpy(best_opps[i], d_opp, sizeof(int) * 4 * n);
        memcpy(best_vts[i], d_vt, sizeof(int) * n);
    }

    for (int i = 0; i < num_in; i++) {
        const int *d_opp = in_opp + i * (4 * n);
        const int *d_vt = in_vt + i * n;

        memcpy(q_opp[0], d_opp, sizeof(int) * 4 * n);
        memcpy(q_vt[0], d_vt, sizeof(int) * n);
        int qh = 0, qt = 1;

        {
            CanonKey k0 = best_keys[i];
            ENSURE_HASH_CAP();
            size_t msk = key_to_comp.cap - 1;
            size_t idx = (size_t)canon_key_hash(n, &k0) & msk;
            int already = 0;
            while (key_to_comp.slots[idx].used) {
                if (canon_key_eq(n, &key_to_comp.slots[idx].key, &k0)) {
                    UNION_OWNERS(i, slot_owner[idx]);
                    already = 1;
                    break;
                }
                idx = (idx + 1) & msk;
            }
            if (already) continue;
            key_to_comp.slots[idx].used = 1;
            key_to_comp.slots[idx].key = k0;
            slot_owner[idx] = i;
            key_to_comp.count++;
        }

        int bfs_cap = (max_bfs_states < 512) ? max_bfs_states : 512;
        while (qh < qt) {
            const int *co = q_opp[qh];
            const int *cv = q_vt[qh];
            qh++;
            int nnbr = list_all_level_neighbors_c(n, co, cv, nbr_opps, nbr_vts, 128);
            for (int ni = 0; ni < nnbr; ni++) {
                if (has_immediate_reduction(n, nbr_opps[ni], nbr_vts[ni]) || !is_connected_c(n, nbr_opps[ni])) {
                    int ri; FIND_ROOT(i, ri); dead[ri] = 1;
                    continue;
                }
                CanonKey nk;
                canonical_sym_key_c(n, nbr_opps[ni], nbr_vts[ni], &nk);
                ENSURE_HASH_CAP();
                size_t msk = key_to_comp.cap - 1;
                size_t idx = (size_t)canon_key_hash(n, &nk) & msk;
                int hit = 0;
                while (key_to_comp.slots[idx].used) {
                    if (canon_key_eq(n, &key_to_comp.slots[idx].key, &nk)) {
                        UNION_OWNERS(i, slot_owner[idx]);
                        hit = 1;
                        break;
                    }
                    idx = (idx + 1) & msk;
                }
                if (hit) continue;
                if (n >= 4 && has_0_or_2_edge_cut_c(n, nbr_opps[ni])) {
                    int ri; FIND_ROOT(i, ri); dead[ri] = 1;
                    continue;
                }
                key_to_comp.slots[idx].used = 1;
                key_to_comp.slots[idx].key = nk;
                slot_owner[idx] = i;
                key_to_comp.count++;

                int ri; FIND_ROOT(i, ri);
                if (canon_key_less(n, &nk, &best_keys[ri])) {
                    best_keys[ri] = nk;
                    memcpy(best_opps[ri], nbr_opps[ni], sizeof(int) * 4 * n);
                    memcpy(best_vts[ri], nbr_vts[ni], sizeof(int) * n);
                }
                if (qt < bfs_cap) {
                    memcpy(q_opp[qt], nbr_opps[ni], sizeof(int) * 4 * n);
                    memcpy(q_vt[qt], nbr_vts[ni], sizeof(int) * n);
                    qt++;
                }
            }
        }
    }

    /* Step 2: Tunnel from level-states of each surviving root `i`, plus 2nd-hop probes on new n-vertex bridges! */
    const int max_bridges = 1024;
    CanonKey *bridge_keys = (CanonKey *)malloc((size_t)max_bridges * sizeof(CanonKey));
    int (*bridge_opps)[4 * MAX_N] = malloc((size_t)max_bridges * sizeof(*bridge_opps));
    int (*bridge_vts)[MAX_N] = malloc((size_t)max_bridges * sizeof(*bridge_vts));
    int (*sq_opp)[4 * MAX_N] = malloc(512 * sizeof(*sq_opp));
    int (*sq_vt)[MAX_N] = malloc(512 * sizeof(*sq_vt));
    int *sq_n = (int *)malloc(512 * sizeof(int));
    CanonKey *b_seen = (CanonKey *)malloc(256 * sizeof(CanonKey));

    const int max_probes = (depth_plus1 <= 32 && depth_plus2 <= 32) ? 1 : 8;
    const int max_extra = (depth_plus1 <= 32 && depth_plus2 <= 32) ? 0 : 8;

    for (int i = 0; i < num_in; i++) {
        int ri;
        FIND_ROOT(i, ri);
        if (dead[ri]) continue;

        int probe_opps[8][4 * MAX_N];
        int probe_vts[8][MAX_N];
        memcpy(probe_opps[0], in_opp + i * (4 * n), sizeof(int) * 4 * n);
        memcpy(probe_vts[0], in_vt + i * n, sizeof(int) * n);
        int n_probes = 1;
        if (max_probes > 1) {
            int nn0 = list_all_level_neighbors_c(n, probe_opps[0], probe_vts[0], nbr_opps, nbr_vts, 64);
            for (int p_i = 0; p_i < nn0 && n_probes < max_probes; p_i++) {
                if (!has_immediate_reduction(n, nbr_opps[p_i], nbr_vts[p_i]) && is_connected_c(n, nbr_opps[p_i])) {
                    memcpy(probe_opps[n_probes], nbr_opps[p_i], sizeof(int) * 4 * n);
                    memcpy(probe_vts[n_probes], nbr_vts[p_i], sizeof(int) * n);
                    n_probes++;
                }
            }
        }

        int num_bridges = 0;
        int killed = 0;
        for (int pr = 0; pr < n_probes && !killed; pr++) {
            if (run_single_probe_stabilizations_c(
                n, probe_opps[pr], probe_vts[pr], depth_plus1, depth_plus2,
                sq_opp, sq_vt, sq_n, bridge_keys, bridge_opps, bridge_vts, &num_bridges, max_bridges
            )) {
                killed = 1;
            }
        }

        FIND_ROOT(i, ri);
        if (killed) {
            dead[ri] = 1;
            continue;
        }

        /* Process discovered n-vertex bridges; if a bridge lands in a brand-new n-vertex level component,
         * register that level component into key_to_comp and run up to max_extra second-hop probes! */
        int extra_probes_done = 0;
        for (int b_i = 0; b_i < num_bridges && !killed; b_i++) {
            CanonKey bk = bridge_keys[b_i];
            size_t msk = key_to_comp.cap - 1;
            size_t idx = (size_t)canon_key_hash(n, &bk) & msk;
            int found_owner = -1;
            while (key_to_comp.slots[idx].used) {
                if (canon_key_eq(n, &key_to_comp.slots[idx].key, &bk)) {
                    found_owner = slot_owner[idx];
                    break;
                }
                idx = (idx + 1) & msk;
            }
            if (found_owner >= 0) {
                UNION_OWNERS(i, found_owner);
                FIND_ROOT(i, ri);
                if (dead[ri]) { killed = 1; break; }
                continue;
            }

            /* Bridge missed `key_to_comp`: explore up to 256 level states from `bridge_opps[b_i]`! */
            memcpy(q_opp[0], bridge_opps[b_i], sizeof(int) * 4 * n);
            memcpy(q_vt[0], bridge_vts[b_i], sizeof(int) * n);
            b_seen[0] = bk;
            int bqh = 0, bqt = 1;
            while (bqh < bqt && !killed && found_owner < 0) {
                const int *bo = q_opp[bqh];
                const int *bv = q_vt[bqh];
                bqh++;
                int bnn = list_all_level_neighbors_c(n, bo, bv, nbr_opps, nbr_vts, 128);
                for (int ni = 0; ni < bnn; ni++) {
                    if (has_immediate_reduction(n, nbr_opps[ni], nbr_vts[ni]) || !is_connected_c(n, nbr_opps[ni])) {
                        killed = 1;
                        break;
                    }
                    if (n >= 4 && has_0_or_2_edge_cut_c(n, nbr_opps[ni])) {
                        killed = 1;
                        break;
                    }
                    CanonKey nk;
                    canonical_sym_key_c(n, nbr_opps[ni], nbr_vts[ni], &nk);
                    size_t nidx = (size_t)canon_key_hash(n, &nk) & msk;
                    while (key_to_comp.slots[nidx].used) {
                        if (canon_key_eq(n, &key_to_comp.slots[nidx].key, &nk)) {
                            found_owner = slot_owner[nidx];
                            break;
                        }
                        nidx = (nidx + 1) & msk;
                    }
                    if (found_owner >= 0) break;
                    int dup = 0;
                    for (int s = 0; s < bqt; s++) {
                        if (canon_key_eq(n, &b_seen[s], &nk)) { dup = 1; break; }
                    }
                    if (!dup && bqt < 256) {
                        memcpy(q_opp[bqt], nbr_opps[ni], sizeof(int) * 4 * n);
                        memcpy(q_vt[bqt], nbr_vts[ni], sizeof(int) * n);
                        b_seen[bqt] = nk;
                        bqt++;
                    }
                }
            }
            if (killed) {
                FIND_ROOT(i, ri);
                dead[ri] = 1;
                break;
            }
            if (found_owner >= 0) {
                UNION_OWNERS(i, found_owner);
                FIND_ROOT(i, ri);
                if (dead[ri]) { killed = 1; break; }
            } else {
                /* Register all `bqt` states in `key_to_comp` owned by `i`, and run a 2nd-hop probe! */
                for (int s = 0; s < bqt; s++) {
                    ENSURE_HASH_CAP();
                    msk = key_to_comp.cap - 1;
                    size_t sidx = (size_t)canon_key_hash(n, &b_seen[s]) & msk;
                    int ex = 0;
                    while (key_to_comp.slots[sidx].used) {
                        if (canon_key_eq(n, &key_to_comp.slots[sidx].key, &b_seen[s])) { ex = 1; break; }
                        sidx = (sidx + 1) & msk;
                    }
                    if (!ex) {
                        key_to_comp.slots[sidx].used = 1;
                        key_to_comp.slots[sidx].key = b_seen[s];
                        slot_owner[sidx] = i;
                        key_to_comp.count++;
                    }
                }
                if (extra_probes_done < max_extra) {
                    extra_probes_done++;
                    if (run_single_probe_stabilizations_c(
                        n, bridge_opps[b_i], bridge_vts[b_i], depth_plus1, depth_plus2,
                        sq_opp, sq_vt, sq_n, bridge_keys, bridge_opps, bridge_vts, &num_bridges, max_bridges
                    )) {
                        killed = 1;
                        FIND_ROOT(i, ri);
                        dead[ri] = 1;
                        break;
                    }
                }
            }
        }
    }

    int kept_cnt = 0;
    for (int i = 0; i < num_in; i++) {
        if (parent[i] == i && !dead[i]) {
            if (kept_cnt < max_kept) {
                memcpy(out_opp + kept_cnt * (4 * n), best_opps[i], sizeof(int) * 4 * n);
                memcpy(out_vt + kept_cnt * n, best_vts[i], sizeof(int) * n);
            }
            kept_cnt++;
        }
    }

    #undef UNION_OWNERS
    #undef ENSURE_HASH_CAP
    #undef FIND_ROOT
    free(b_seen);
    free(sq_n);
    free(sq_vt);
    free(sq_opp);
    free(bridge_vts);
    free(bridge_opps);
    free(bridge_keys);
    free(slot_owner);
    free(best_vts);
    free(best_opps);
    free(best_keys);
    free(dead);
    free(parent);
    free(q_vt);
    free(q_opp);
    hashset_free(&key_to_comp);
    return kept_cnt;
}


