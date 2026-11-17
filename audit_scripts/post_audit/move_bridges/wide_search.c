#include "/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/worktrees/1dadd8f0-f304-455a-b682-76a831212ce6/subagent-Deep-Soundness---Move-Bridge-Auditor-self-bc48a04a/code/src/surfaceknots/_fast_enum.c"

typedef struct {
    uint8_t nv;
    uint8_t owner;
    uint8_t expanded_lvl;
    uint8_t pad;
    int32_t parent;
    uint8_t opp[4 * MAX_N];
    uint8_t vt[MAX_N];
} WideNode;

typedef struct {
    int32_t node_idx; /* -1 if empty */
    uint8_t nv;
    CanonKey key;
} WideHSlot;

typedef struct {
    int32_t owner_a;
    int32_t parent_node_a;
    int32_t mid_nv;
    int32_t node_b;
    int32_t wave_h;
    int32_t mid_opp[4 * MAX_N];
    int32_t mid_vt[MAX_N];
} WideBridgeRecord;

typedef struct {
    int32_t collapsed; /* 0 = no, 1 = reduced (< n), 2 = split/loops, 3 = 2-cut at n */
    int32_t parent_node;
    int32_t term_nv;
    int32_t wave_h;
    int32_t term_opp[4 * MAX_N];
    int32_t term_vt[MAX_N];
} WideCollapseRecord;

static inline int dsu_find(int *dsu, int x) {
    while (dsu[x] != x) {
        dsu[x] = dsu[dsu[x]];
        x = dsu[x];
    }
    return x;
}

static inline void unpack_node(const WideNode *nd, int *opp_out, int *vt_out) {
    int m = (int)nd->nv;
    for (int i = 0; i < 4 * m; i++) opp_out[i] = (int)nd->opp[i];
    for (int i = 0; i < m; i++) vt_out[i] = (int)nd->vt[i];
}

static inline void pack_node(WideNode *nd, int nv, int owner, int parent, const int *opp_in, const int *vt_in) {
    nd->nv = (uint8_t)nv;
    nd->owner = (uint8_t)owner;
    nd->expanded_lvl = 0;
    nd->pad = 0;
    nd->parent = (int32_t)parent;
    for (int i = 0; i < 4 * nv; i++) nd->opp[i] = (uint8_t)opp_in[i];
    for (int i = 0; i < nv; i++) nd->vt[i] = (uint8_t)vt_in[i];
}

/* Insert or match state (nr, r_opp, r_vt) generated from `parent_idx` (or -1 for seed) by `owner`.
 * Returns:
 *   >= 0 : newly allocated node index
 *   -1   : duplicate within same DSU component (or table full)
 *   -2   : bridged two distinct DSU components!
 */
static int try_insert_or_bridge(
    int n_base,
    int wave_h,
    int nr,
    const int *r_opp,
    const int *r_vt,
    int owner,
    int parent_idx,
    WideNode *nodes,
    int *num_nodes,
    int max_nodes,
    WideHSlot *ht,
    size_t ht_mask,
    int *dsu,
    WideBridgeRecord *bridges,
    int *num_bridges,
    int max_bridges
) {
    CanonKey rk;
    canonical_sym_key_c(nr, r_opp, r_vt, &rk);
    uint64_t h = canon_key_hash(nr, &rk) ^ ((uint64_t)nr * 0x9e3779b97f4a7c15ULL);
    size_t idx = (size_t)h & ht_mask;

    while (ht[idx].node_idx >= 0) {
        if (ht[idx].nv == (uint8_t)nr && canon_key_eq(nr, &ht[idx].key, &rk)) {
            int other_node = ht[idx].node_idx;
            int other_owner = (int)nodes[other_node].owner;
            int ra = dsu_find(dsu, owner);
            int rb = dsu_find(dsu, other_owner);
            if (ra != rb) {
                if (*num_bridges < max_bridges) {
                    WideBridgeRecord *br = &bridges[*num_bridges];
                    br->owner_a = owner;
                    br->parent_node_a = parent_idx;
                    br->mid_nv = nr;
                    br->node_b = other_node;
                    br->wave_h = wave_h;
                    for (int k = 0; k < 4 * nr; k++) br->mid_opp[k] = r_opp[k];
                    for (int k = 0; k < nr; k++) br->mid_vt[k] = r_vt[k];
                    (*num_bridges)++;
                }
                if (ra < rb) dsu[rb] = ra;
                else dsu[ra] = rb;
                return -2;
            }
            return -1;
        }
        idx = (idx + 1) & ht_mask;
    }

    if (*num_nodes >= max_nodes) {
        return -1;
    }
    int new_id = (*num_nodes)++;
    pack_node(&nodes[new_id], nr, owner, parent_idx, r_opp, r_vt);
    ht[idx].node_idx = new_id;
    ht[idx].nv = (uint8_t)nr;
    ht[idx].key = rk;
    return new_id;
}

static void record_collapse(
    WideCollapseRecord *col,
    int kind,
    int parent_node,
    int term_nv,
    int wave_h,
    const int *term_opp,
    const int *term_vt
) {
    if (col->collapsed != 0) return;
    col->collapsed = kind;
    col->parent_node = parent_node;
    col->term_nv = term_nv;
    col->wave_h = wave_h;
    for (int k = 0; k < 4 * term_nv; k++) col->term_opp[k] = term_opp[k];
    for (int k = 0; k < term_nv; k++) col->term_vt[k] = term_vt[k];
}

static int count_distinct_alive_roots(int num_owners, int *dsu, const WideCollapseRecord *cols) {
    uint8_t root_seen[256];
    memset(root_seen, 0, sizeof(root_seen));
    int cnt = 0;
    for (int o = 0; o < num_owners; o++) {
        int r = dsu_find(dsu, o);
        if (!root_seen[r]) {
            root_seen[r] = 1;
            /* Check if any member in this root collapsed */
            int any_col = 0;
            for (int j = 0; j < num_owners; j++) {
                if (dsu_find(dsu, j) == r && cols[j].collapsed != 0) {
                    any_col = 1;
                    break;
                }
            }
            if (!any_col) cnt++;
        }
    }
    return cnt;
}

int wide_bidirectional_search_c(
    int n_base,
    int num_owners,
    int num_seeds,
    const int *seed_opp,
    const int *seed_vt,
    const int *seed_owner,
    int max_headroom,          /* 1..4 (n_base + max_headroom <= MAX_N) */
    int wave_budget,           /* e.g. 120000 states per wave */
    int up1_per_owner,         /* e.g. 512 states per owner climbed via Omega_1+/6+ */
    int up2_per_owner,         /* e.g. 256 states per owner climbed via Omega_2+ */
    int stop_when_merged,      /* 1 = return as soon as <= 1 uncollapsed DSU root remains (or single owner collapsed) */
    /* Outputs: */
    int *out_dsu_roots,        /* [num_owners] */
    WideCollapseRecord *out_collapses, /* [num_owners] */
    WideBridgeRecord *out_bridges,     /* [max_bridges] */
    int max_bridges,
    int *out_num_bridges,
    /* Exported node chain buffer so Python can reconstruct any parent chain: */
    int *out_total_nodes,
    /* Callback-free chain extraction: copy all nodes up to out_nodes_cap if requested, or reconstruct chains of bridges/collapses! */
    int *out_chains_buf,       /* packed chains: [num_chains, then for each chain: len, (nv, opp[4*nv], vt[nv]) * len] */
    int max_chains_ints
) {
    if (n_base < 1 || n_base + max_headroom > MAX_N || num_owners < 1 || num_owners > 128) {
        return -1;
    }

    int max_nodes = num_seeds + (max_headroom + 1) * wave_budget + 8192;
    if (max_nodes > 1600000) max_nodes = 1600000;

    size_t ht_cap = 1ULL << 18;
    while (ht_cap < (size_t)max_nodes * 2ULL) ht_cap <<= 1;
    size_t ht_mask = ht_cap - 1ULL;

    WideNode *nodes = (WideNode *)malloc((size_t)max_nodes * sizeof(WideNode));
    WideHSlot *ht = (WideHSlot *)malloc(ht_cap * sizeof(WideHSlot));
    if (!nodes || !ht) {
        free(nodes);
        free(ht);
        return -2;
    }
    for (size_t i = 0; i < ht_cap; i++) ht[i].node_idx = -1;

    int dsu[128];
    for (int o = 0; o < num_owners; o++) {
        dsu[o] = o;
        memset(&out_collapses[o], 0, sizeof(WideCollapseRecord));
    }

    int num_nodes = 0;
    int num_bridges = 0;

    /* Per-owner priority queues: prio_q holds reduced states (nv < m_ceil), reg_q holds m_ceil states */
    int q_cap_per_owner = (wave_budget / (num_owners > 0 ? num_owners : 1)) * (max_headroom + 2) + num_seeds + 4096;
    if (q_cap_per_owner > max_nodes) q_cap_per_owner = max_nodes;
    int32_t **prio_q = (int32_t **)malloc((size_t)num_owners * sizeof(int32_t *));
    int32_t **reg_q = (int32_t **)malloc((size_t)num_owners * sizeof(int32_t *));
    int ph[128], pt[128], rh[128], rt[128];
    for (int o = 0; o < num_owners; o++) {
        prio_q[o] = (int32_t *)malloc((size_t)q_cap_per_owner * sizeof(int32_t));
        reg_q[o] = (int32_t *)malloc((size_t)q_cap_per_owner * sizeof(int32_t));
        ph[o] = pt[o] = rh[o] = rt[o] = 0;
    }

    #define ENQ_NODE(owner_v, nid_v, m_ceil_v) do { \
        int _ow = (owner_v); \
        int _nid = (nid_v); \
        if (_nid >= 0) { \
            if ((int)nodes[_nid].nv < (m_ceil_v)) { \
                if (pt[_ow] < q_cap_per_owner) prio_q[_ow][pt[_ow]++] = _nid; \
            } else { \
                if (rt[_ow] < q_cap_per_owner) reg_q[_ow][rt[_ow]++] = _nid; \
            } \
        } \
    } while (0)

    int nbr_opps[128][4 * MAX_N];
    int nbr_vts[128][MAX_N];
    int cur_opp[4 * MAX_N], cur_vt[MAX_N];
    int red_opp[4 * MAX_N], red_vt[MAX_N];

    /* Wave 0: Insert all n_base seeds and complete n_base level closure */
    for (int s = 0; s < num_seeds; s++) {
        const int *so = seed_opp + s * (4 * n_base);
        const int *sv = seed_vt + s * n_base;
        int ow = seed_owner[s];
        int loops = 0;
        int nr = greedy_reduce_c(n_base, so, sv, red_opp, red_vt, &loops);
        if (loops > 0 || !is_connected_c(nr, red_opp)) {
            record_collapse(&out_collapses[ow], 2, -1, nr, 0, red_opp, red_vt);
            continue;
        }
        if (nr < n_base) {
            record_collapse(&out_collapses[ow], 1, -1, nr, 0, red_opp, red_vt);
            continue;
        }
        if (n_base >= 4 && has_0_or_2_edge_cut_c(n_base, red_opp)) {
            record_collapse(&out_collapses[ow], 3, -1, nr, 0, red_opp, red_vt);
            continue;
        }
        int nid = try_insert_or_bridge(
            n_base, 0, n_base, red_opp, red_vt, ow, -1,
            nodes, &num_nodes, max_nodes, ht, ht_mask, dsu,
            out_bridges, &num_bridges, max_bridges
        );
        ENQ_NODE(ow, nid, n_base);
    }

    /* Finish Wave 0 level BFS at n_base */
    {
        int any_q = 1;
        while (any_q && num_nodes < max_nodes) {
            any_q = 0;
            for (int o = 0; o < num_owners; o++) {
                int nid = -1;
                if (ph[o] < pt[o]) nid = prio_q[o][ph[o]++];
                else if (rh[o] < rt[o]) nid = reg_q[o][rh[o]++];
                if (nid < 0) continue;
                any_q = 1;
                if (nodes[nid].expanded_lvl) continue;
                nodes[nid].expanded_lvl = 1;
                unpack_node(&nodes[nid], cur_opp, cur_vt);
                int cnt = list_all_level_neighbors_c(n_base, cur_opp, cur_vt, nbr_opps, nbr_vts, 128);
                for (int i = 0; i < cnt; i++) {
                    int loops = 0;
                    int nr = greedy_reduce_c(n_base, nbr_opps[i], nbr_vts[i], red_opp, red_vt, &loops);
                    if (loops > 0 || !is_connected_c(nr, red_opp)) {
                        record_collapse(&out_collapses[o], 2, nid, nr, 0, red_opp, red_vt);
                        break;
                    }
                    if (nr < n_base) {
                        record_collapse(&out_collapses[o], 1, nid, nr, 0, red_opp, red_vt);
                        break;
                    }
                    if (n_base >= 4 && has_0_or_2_edge_cut_c(n_base, red_opp)) {
                        record_collapse(&out_collapses[o], 3, nid, nr, 0, red_opp, red_vt);
                        break;
                    }
                    int cid = try_insert_or_bridge(
                        n_base, 0, n_base, red_opp, red_vt, o, nid,
                        nodes, &num_nodes, max_nodes, ht, ht_mask, dsu,
                        out_bridges, &num_bridges, max_bridges
                    );
                    ENQ_NODE(o, cid, n_base);
                }
            }
        }
    }

    /* Waves h = 1 .. max_headroom */
    for (int h = 1; h <= max_headroom; h++) {
        if (stop_when_merged) {
            if (num_owners == 1 && out_collapses[0].collapsed != 0) break;
            if (num_owners > 1 && count_distinct_alive_roots(num_owners, dsu, out_collapses) <= 1) break;
        }
        if (num_nodes >= max_nodes) break;

        int m_ceil = n_base + h;
        int wave_start_nodes = num_nodes;
        int scan_limit = num_nodes;

        /* Reset queues for wave h */
        for (int o = 0; o < num_owners; o++) {
            ph[o] = pt[o] = rh[o] = rt[o] = 0;
        }

        /* Collect up to up1_per_owner states at nv == m_ceil - 1, and up2_per_owner states at nv == m_ceil - 2 */
        int *up1_lists[128], *up2_lists[128];
        int up1_cnt[128], up2_cnt[128];
        for (int o = 0; o < num_owners; o++) {
            up1_lists[o] = (int *)malloc((size_t)(up1_per_owner + 1) * sizeof(int));
            up2_lists[o] = (int *)malloc((size_t)(up2_per_owner + 1) * sizeof(int));
            up1_cnt[o] = 0;
            up2_cnt[o] = 0;
        }
        for (int idx = 0; idx < scan_limit; idx++) {
            int o = (int)nodes[idx].owner;
            int nv = (int)nodes[idx].nv;
            if (nv == m_ceil - 1 && up1_cnt[o] < up1_per_owner) {
                up1_lists[o][up1_cnt[o]++] = idx;
            }
            if (h >= 2 && nv == m_ceil - 2 && up2_cnt[o] < up2_per_owner) {
                up2_lists[o][up2_cnt[o]++] = idx;
            }
        }

        /* 1A. Round-robin Omega_1+ / Omega_6+ climb from (m_ceil - 1) -> m_ceil */
        int m1 = m_ceil - 1;
        for (int step = 0; step < up1_per_owner && num_nodes < max_nodes && (num_nodes - wave_start_nodes) < wave_budget; step++) {
            int any_step = 0;
            for (int o = 0; o < num_owners; o++) {
                if (step >= up1_cnt[o]) continue;
                if (num_owners == 1 && out_collapses[0].collapsed != 0) break;
                any_step = 1;
                int p_nid = up1_lists[o][step];
                int po[4 * MAX_N], pv[MAX_N];
                unpack_node(&nodes[p_nid], po, pv);
                int u = m1;
                for (int d1 = 0; d1 < 4 * m1; d1++) {
                    int d2 = po[d1];
                    if (d2 < d1) continue;
                    for (int side = 0; side < 2; side++) {
                        int s_opp[4 * MAX_N], s_vt[MAX_N];
                        memcpy(s_opp, po, sizeof(int) * 4 * m1);
                        memcpy(s_vt, pv, sizeof(int) * m1);
                        if (side == 0) {
                            s_opp[4 * u + 0] = 4 * u + 1; s_opp[4 * u + 1] = 4 * u + 0;
                            s_opp[4 * u + 2] = d1; s_opp[d1] = 4 * u + 2;
                            s_opp[4 * u + 3] = d2; s_opp[d2] = 4 * u + 3;
                        } else {
                            s_opp[4 * u + 1] = 4 * u + 2; s_opp[4 * u + 2] = 4 * u + 1;
                            s_opp[4 * u + 3] = d1; s_opp[d1] = 4 * u + 3;
                            s_opp[4 * u + 0] = d2; s_opp[d2] = 4 * u + 0;
                        }
                        for (int vt_u = 0; vt_u < 4; vt_u++) {
                            s_vt[u] = vt_u;
                            int cid = try_insert_or_bridge(
                                n_base, h, m_ceil, s_opp, s_vt, o, p_nid,
                                nodes, &num_nodes, max_nodes, ht, ht_mask, dsu,
                                out_bridges, &num_bridges, max_bridges
                            );
                            ENQ_NODE(o, cid, m_ceil);
                        }
                    }
                }
            }
            if (!any_step) break;
            if (stop_when_merged && num_owners > 1 && count_distinct_alive_roots(num_owners, dsu, out_collapses) <= 1) break;
        }

        /* 1B. Round-robin Omega_2+ poke climb from (m_ceil - 2) -> m_ceil */
        if (h >= 2) {
            int m2 = m_ceil - 2;
            for (int step = 0; step < up2_per_owner && num_nodes < max_nodes && (num_nodes - wave_start_nodes) < wave_budget; step++) {
                if (stop_when_merged && num_owners > 1 && count_distinct_alive_roots(num_owners, dsu, out_collapses) <= 1) break;
                int any_step = 0;
                for (int o = 0; o < num_owners; o++) {
                    if (step >= up2_cnt[o]) continue;
                    if (num_owners == 1 && out_collapses[0].collapsed != 0) break;
                    any_step = 1;
                    int p_nid = up2_lists[o][step];
                    int po[4 * MAX_N], pv[MAX_N];
                    unpack_node(&nodes[p_nid], po, pv);
                    uint8_t vis_d[4 * MAX_N];
                    memset(vis_d, 0, (size_t)(4 * m2));
                    for (int sd = 0; sd < 4 * m2; sd++) {
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
                                int s_opp[4 * MAX_N], s_vt[MAX_N];
                                memcpy(s_opp, po, sizeof(int) * 4 * m2);
                                memcpy(s_vt, pv, sizeof(int) * m2);
                                int u = m2, w = m2 + 1;
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
                                    int cid = try_insert_or_bridge(
                                        n_base, h, m_ceil, s_opp, s_vt, o, p_nid,
                                        nodes, &num_nodes, max_nodes, ht, ht_mask, dsu,
                                        out_bridges, &num_bridges, max_bridges
                                    );
                                    ENQ_NODE(o, cid, m_ceil);
                                }
                            }
                        }
                    }
                }
                if (!any_step) break;
            }
        }

        for (int o = 0; o < num_owners; o++) {
            free(up1_lists[o]);
            free(up2_lists[o]);
        }

        /* 2. Fair Round-Robin Level + Greedy-Reduction Priority BFS at <= m_ceil */
        int any_q = 1;
        while (any_q && num_nodes < max_nodes && (num_nodes - wave_start_nodes) < wave_budget) {
            if (stop_when_merged) {
                if (num_owners == 1 && out_collapses[0].collapsed != 0) break;
                if (num_owners > 1 && count_distinct_alive_roots(num_owners, dsu, out_collapses) <= 1) break;
            }
            any_q = 0;
            for (int o = 0; o < num_owners; o++) {
                if (num_owners == 1 && out_collapses[0].collapsed != 0) break;
                int nid = -1;
                if (ph[o] < pt[o]) nid = prio_q[o][ph[o]++];
                else if (rh[o] < rt[o]) nid = reg_q[o][rh[o]++];
                if (nid < 0) continue;
                any_q = 1;
                if (nodes[nid].expanded_lvl) continue;
                nodes[nid].expanded_lvl = 1;

                int cur_m = (int)nodes[nid].nv;
                unpack_node(&nodes[nid], cur_opp, cur_vt);
                int cnt = list_all_level_neighbors_c(cur_m, cur_opp, cur_vt, nbr_opps, nbr_vts, 128);
                for (int i = 0; i < cnt; i++) {
                    int loops = 0;
                    int nr = greedy_reduce_c(cur_m, nbr_opps[i], nbr_vts[i], red_opp, red_vt, &loops);
                    if (loops > 0 || !is_connected_c(nr, red_opp)) {
                        record_collapse(&out_collapses[o], 2, nid, nr, h, red_opp, red_vt);
                        break;
                    }
                    if (nr < n_base) {
                        record_collapse(&out_collapses[o], 1, nid, nr, h, red_opp, red_vt);
                        break;
                    }
                    if (nr == n_base && n_base >= 4 && has_0_or_2_edge_cut_c(n_base, red_opp)) {
                        record_collapse(&out_collapses[o], 3, nid, nr, h, red_opp, red_vt);
                        break;
                    }
                    if (nr == n_base + 1 && n_base >= 3 && has_0_or_2_edge_cut_c(nr, red_opp)) {
                        record_collapse(&out_collapses[o], 4, nid, nr, h, red_opp, red_vt);
                        break;
                    }
                    if (nr == n_base + 2 && n_base >= 3) {
                        int hb = 0;
                        int b_opp[4 * MAX_N], b_vt[MAX_N];
                        if (check_np2_cuts_c(n_base, red_opp, red_vt, b_opp, b_vt, &hb)) {
                            record_collapse(&out_collapses[o], 5, nid, nr, h, red_opp, red_vt);
                            break;
                        }
                        if (hb) {
                            /* Do not set parent=nid directly on b_opp (since b_opp pinched a trivial 2-vertex S^2),
                             * or only insert if we want bridge discovery; wait: for path verification, if b_opp is
                             * just an S^2 blister pinch, standard Omega moves from red_opp will reach b_opp anyway! */
                        }
                    }
                    int cid = try_insert_or_bridge(
                        n_base, h, nr, red_opp, red_vt, o, nid,
                        nodes, &num_nodes, max_nodes, ht, ht_mask, dsu,
                        out_bridges, &num_bridges, max_bridges
                    );
                    ENQ_NODE(o, cid, m_ceil);
                }
            }
        }
    }

    for (int o = 0; o < num_owners; o++) {
        out_dsu_roots[o] = dsu_find(dsu, o);
        free(prio_q[o]);
        free(reg_q[o]);
    }
    free(prio_q);
    free(reg_q);

    *out_num_bridges = num_bridges;
    *out_total_nodes = num_nodes;

    /* Pack parent chains into out_chains_buf:
     * First `num_owners` chains: for each owner `o`, if collapsed != 0 and parent_node >= 0, chain from root seed -> parent_node (else len=0).
     * Next `2 * num_bridges` chains: for each bridge `b`, chain A (seed -> parent_node_a) and chain B (seed -> node_b).
     */
    if (out_chains_buf && max_chains_ints > 16) {
        int pos = 0;
        int total_chains = num_owners + 2 * num_bridges;
        out_chains_buf[pos++] = total_chains;

        int stack_ids[4096];
        #define EMIT_CHAIN(end_nid) do { \
            int _len = 0; \
            int _cur = (end_nid); \
            while (_cur >= 0 && _len < 4096) { \
                stack_ids[_len++] = _cur; \
                _cur = nodes[_cur].parent; \
            } \
            int _needed = 1; \
            for (int _k = 0; _k < _len; _k++) { \
                int _m = (int)nodes[stack_ids[_k]].nv; \
                _needed += 1 + 5 * _m; \
            } \
            if (pos + _needed > max_chains_ints) { \
                if (pos < max_chains_ints) out_chains_buf[pos++] = -1; \
            } else { \
                out_chains_buf[pos++] = _len; \
                for (int _k = _len - 1; _k >= 0; _k--) { \
                    const WideNode *_nd = &nodes[stack_ids[_k]]; \
                    int _m = (int)_nd->nv; \
                    out_chains_buf[pos++] = _m; \
                    for (int _j = 0; _j < 4 * _m; _j++) out_chains_buf[pos++] = (int)_nd->opp[_j]; \
                    for (int _j = 0; _j < _m; _j++) out_chains_buf[pos++] = (int)_nd->vt[_j]; \
                } \
            } \
        } while (0)

        for (int o = 0; o < num_owners; o++) {
            if (out_collapses[o].collapsed != 0 && out_collapses[o].parent_node >= 0) {
                EMIT_CHAIN(out_collapses[o].parent_node);
            } else {
                if (pos < max_chains_ints) out_chains_buf[pos++] = 0;
            }
        }
        for (int b = 0; b < num_bridges; b++) {
            EMIT_CHAIN(out_bridges[b].parent_node_a);
            EMIT_CHAIN(out_bridges[b].node_b);
        }
        #undef EMIT_CHAIN
    }

    #undef ENQ_NODE
    free(ht);
    free(nodes);
    return 0;
}
