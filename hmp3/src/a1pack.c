/*
 * a1pack.c - lossless re-coding of MPEG-1 Layer III granules (roadmap stage A1).
 *
 * See pub/a1pack.h. Every quantized value and every band's effective step size
 * stays as the encoder chose it; only the way they are written changes.
 *
 * Effective step of a long band b with nonzero values (quarter steps):
 *     E[b] = global_gain - 2*(1+scalefac_scale)*(sf[b] + preflag*pretab[b])
 * and E[21] = global_gain for the band without a scalefactor. Any choice of
 * global_gain, scalefac_scale, preflag and sf that keeps E[b] for every band with
 * nonzero values decodes to the same PCM; bands with no nonzero value can take
 * any scalefactor. Short blocks keep their gains; only free values and
 * scalefac_compress change there.
 *
 * This file is a contribution to the Helix MP3 encoder and is distributed under
 * the same terms as the files it accompanies (RPSL/RCSL, see hmp3/LICENSE.txt).
 */

#include <string.h>
#include <stdlib.h>

#include "l3e.h"
#include "a1pack.h"

#include "htable.h"     /* Helix Huffman tables, as used by l3pack.c */

#define A1_INF      0x3fffffff
#define A1_NPAIR    288
#define A1_WINDOW   80          /* big_values boundary search: values past the minimum */
#define A1_MAXCAND  64

/* scalefactor band starts, MPEG-1, by header sample-rate index: 44.1, 48, 32 kHz */
static const int a1_sfb_l[3][23] = {
    {0, 4, 8, 12, 16, 20, 24, 30, 36, 44, 52, 62, 74, 90, 110, 134, 162, 196, 238, 288, 342, 418, 576},
    {0, 4, 8, 12, 16, 20, 24, 30, 36, 42, 50, 60, 72, 88, 106, 128, 156, 190, 230, 276, 330, 384, 576},
    {0, 4, 8, 12, 16, 20, 24, 30, 36, 44, 54, 66, 82, 102, 126, 156, 194, 240, 296, 364, 448, 550, 576}
};
static const int a1_sfb_s[3][14] = {
    {0, 4, 8, 12, 16, 22, 30, 40, 52, 66, 84, 106, 136, 192},
    {0, 4, 8, 12, 16, 22, 28, 38, 50, 64, 80, 100, 126, 192},
    {0, 4, 8, 12, 16, 22, 30, 42, 58, 78, 104, 138, 180, 192}
};
static const int a1_pretab[22] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 2, 0};
static const int a1_slen1[16] = {0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4};
static const int a1_slen2[16] = {0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3};
static const int a1_grp[5] = {0, 6, 11, 16, 21};       /* scfsi groups of long bands */
static const int a1_quada_code[16] = {1, 5, 4, 5, 6, 5, 4, 4, 7, 3, 6, 0, 7, 2, 3, 1};
static const int a1_quada_len[16] = {1, 4, 4, 5, 4, 6, 5, 6, 4, 5, 5, 6, 5, 6, 6, 6};
/* largest value of each base table (without linbits); -1: table does not exist */
static const int a1_base_max[32] = {0, 1, 2, 2, -1, 3, 3, 5, 5, 5, 7, 7, 7, 15, -1, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15};
static const int a1_small[13] = {1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15};


/* first entry of a1_small[] whose table holds a pair maximum of m (m <= 15) */
static const int a1_small_start[16] = {0, 0, 1, 3, 5, 5, 8, 8, 11, 11, 11, 11, 11, 11, 11, 11};

static int tmax[32];                    /* largest codable value, linbits included */
static unsigned char tlen[32][16][16];  /* code lengths, base values */
static unsigned int tcode[32][16][16];  /* codes */
static int lg2[A1_NPAIR + 1];           /* floor(log2(n)) */
static int lcost_tab[5][4][16], lcomp_tab[5][4][16];    /* [w1][w2][scfsi] */
static int scost_tab[5][4], scomp_tab[5][4];
static int a1_ready = 0;                /* tables above built (once per process) */

/*--------------------------------------------------------------------*/
static void
get_code ( int t, int x, int y, int *code, int *len )
{
    void *tb = huff_case_table[t].table;

    switch ( huff_case_table[t].icase )
    {
    case 1:
        *code = ( ( BYTE_HUFF_STRUCT ( * )[2] ) tb )[x][y].code;
        *len = ( ( BYTE_HUFF_STRUCT ( * )[2] ) tb )[x][y].len;
        break;
    case 2:
        *code = ( ( BYTE_HUFF_STRUCT ( * )[4] ) tb )[x][y].code;
        *len = ( ( BYTE_HUFF_STRUCT ( * )[4] ) tb )[x][y].len;
        break;
    case 3:
        *code = ( ( BYTE_HUFF_STRUCT ( * )[8] ) tb )[x][y].code;
        *len = ( ( BYTE_HUFF_STRUCT ( * )[8] ) tb )[x][y].len;
        break;
    case 4:
        *code = ( ( BYTE_HUFF_STRUCT ( * )[16] ) tb )[x][y].code;
        *len = ( ( BYTE_HUFF_STRUCT ( * )[16] ) tb )[x][y].len;
        break;
    case 5:
        *code = ( ( HUFF_STRUCT ( * )[16] ) tb )[x][y].code;
        *len = ( ( HUFF_STRUCT ( * )[16] ) tb )[x][y].len;
        break;
    default:
        *code = 0;
        *len = 0;
        break;
    }
}

static int
widthv ( int v )        /* bits needed for a scalefactor value 0..15 */
{
    if ( v <= 0 )
        return 0;
    if ( v < 2 )
        return 1;
    if ( v < 4 )
        return 2;
    if ( v < 8 )
        return 3;
    return 4;
}

static void huff_tables_init ( void );

/*--------------------------------------------------------------------*/
/* code tables and scalefactor cost tables, the same for every encoder */
static void
a1_tables_init ( void )
{
    int t, x, y, c, l, w1, w2, r, g;

    if ( a1_ready )
        return;

    memset ( tlen, 0, sizeof ( tlen ) );
    for ( t = 0; t < 32; t++ )
    {
        int m = a1_base_max[t];

        if ( m < 0 )
        {
            tmax[t] = -1;
            continue;
        }
        tmax[t] = ( t >= 16 ) ? 15 + ( 1 << huff_case_table[t].linbits ) - 1 : m;
        if ( t == 0 )
            continue;
        for ( x = 0; x <= m; x++ )
            for ( y = 0; y <= m; y++ )
            {
                get_code ( t, x, y, &c, &l );
                tlen[t][x][y] = ( unsigned char ) l;
                tcode[t][x][y] = ( unsigned int ) c;
            }
    }
    lg2[0] = lg2[1] = 0;
    for ( x = 2; x <= A1_NPAIR; x++ )
        lg2[x] = lg2[x >> 1] + 1;
    huff_tables_init (  );

    /* cheapest scalefac_compress for long blocks: needed widths, scfsi-reused groups */
    for ( w1 = 0; w1 < 5; w1++ )
        for ( w2 = 0; w2 < 4; w2++ )
            for ( r = 0; r < 16; r++ )
            {
                int best = A1_INF, bc = 0;

                for ( c = 0; c < 16; c++ )
                {
                    int s1 = a1_slen1[c], s2 = a1_slen2[c], bits = 0;

                    if ( s1 < w1 || s2 < w2 )
                        continue;
                    for ( g = 0; g < 4; g++ )
                        if ( !( r & ( 8 >> g ) ) )
                            bits += ( a1_grp[g + 1] - a1_grp[g] ) * ( g < 2 ? s1 : s2 );
                    if ( bits < best )
                    {
                        best = bits;
                        bc = c;
                    }
                }
                lcost_tab[w1][w2][r] = best;
                lcomp_tab[w1][w2][r] = bc;
            }
    /* short blocks, non-mixed: 18 values per slen */
    for ( w1 = 0; w1 < 5; w1++ )
        for ( w2 = 0; w2 < 4; w2++ )
        {
            int best = A1_INF, bc = 0;

            for ( c = 0; c < 16; c++ )
            {
                if ( a1_slen1[c] < w1 || a1_slen2[c] < w2 )
                    continue;
                if ( 18 * ( a1_slen1[c] + a1_slen2[c] ) < best )
                {
                    best = 18 * ( a1_slen1[c] + a1_slen2[c] );
                    bc = c;
                }
            }
            scost_tab[w1][w2] = best;
            scomp_tab[w1][w2] = bc;
        }
    a1_ready = 1;
}

/*====================================================================*/
/* Huffman: exhaustive coding of fixed quantized values                */
/*====================================================================*/
/*
 * Prefix sums over pairs, 16 lanes per position:
 *   lanes 0..12  bits in table a1_small[lane], code + signs (pairs up to 15)
 *   lane 13      base bits in table 16 (same codes as 17..23), values clamped to 15
 *   lane 14      base bits in table 24 (same codes as 25..31)
 *   lane 15      count of values >= 15 (each takes linbits in an escape table)
 */
#define A1_NL 16
typedef struct
{
    int ps[A1_NPAIR + 1][A1_NL];
    int lbp[A1_NPAIR + 1];      /* prefix of each pair's cheapest cost in any table */
    int mx[9][A1_NPAIR];        /* sparse table of pair maxima */
    int qa[2][146], qb[2][146]; /* count1 prefix costs, quads aligned at 0 or 2 mod 4 */
}
HCTX;

typedef struct
{
    int ok;             /* 0: values cannot be coded (never seen): keep stock */
    int big_values, region0_count, region1_count;
    int table_select[3];
    int count1table_select, nquads;
    int bits;
}
HRES;

/* per-encoder state (a1_create) */
struct A1CTX
{
    const int *sfbl, *sfbs;     /* band starts for the stream's sample rate */
    int keep_gain;              /* global_gain kept as coded */
    HCTX hc;                    /* Huffman search scratch */
    A1FRAME fz;                 /* values as far as the stock side info codes them */
    unsigned char tmp[4096];    /* stock granule re-encoded by stock_check */
};

static int pcost[16][16][A1_NL];        /* lane costs of one pair (values clamped to 15) */
static int lbsmall[16][16];             /* cheapest cost of a pair with both values <= 15 */
/*
 * For a segment whose largest value is m <= 15, the lanes worth trying: every lane
 * whose table holds m, minus those another such lane beats or ties on every pair
 * of values up to m (a dominated table can never be strictly cheaper).
 */
static int cand_n[16], cand_lane[16][15], cand_mul[16][15], cand_tab[16][15];
static int esc_t[2][14];                /* smallest escape table of each family for linbits n */
static unsigned char esc_nb[8192];      /* bits of v: the linbits a value of 15 + v needs */

/* smallest table of escape family f (0: 16..23, 1: 24..31) that holds m */
static int
esc_table ( int f, int m )
{
    int t;

    for ( t = 16 + 8 * f; t < 24 + 8 * f; t++ )
        if ( tmax[t] >= m )
            return t;
    return -1;
}

static void
huff_tables_init ( void )
{
    int x, y, i;

    for ( x = 0; x < 16; x++ )
        for ( y = 0; y < 16; y++ )
        {
            int s = ( x != 0 ) + ( y != 0 ), m = x > y ? x : y, lb;
            int *c = pcost[x][y];

            for ( i = 0; i < 13; i++ )
                c[i] = tlen[a1_small[i]][x][y] + s;     /* unused where the table cannot hold m */
            c[13] = tlen[16][x][y] + s;
            c[14] = tlen[24][x][y] + s;
            c[15] = ( x == 15 ) + ( y == 15 );
            lb = c[13] + c[15];
            if ( c[14] + 4 * c[15] < lb )
                lb = c[14] + 4 * c[15];
            for ( i = a1_small_start[m]; i < 13; i++ )
                if ( c[i] < lb )
                    lb = c[i];
            lbsmall[x][y] = ( m == 0 ) ? 0 : lb;        /* (0,0) may sit in a table-0 region */
        }

    {
        int m, a, b;

        for ( m = 1; m < 16; m++ )
        {
            int ok[15], n = 0;

            for ( a = 0; a < 15; a++ )
                ok[a] = a >= 13 || tmax[a1_small[a]] >= m;
            for ( a = 0; a < 15; a++ )
            {
                if ( !ok[a] )
                    continue;
                for ( b = 0; b < 15 && ok[a]; b++ )
                {
                    int dom = 1, same = 1;

                    if ( b == a || !ok[b] )
                        continue;
                    /* does lane b beat or tie lane a on every pair up to m? */
                    for ( x = 0; x <= m && dom; x++ )
                        for ( y = 0; y <= m && dom; y++ )
                        {
                            int ca = pcost[x][y][a] + ( a == 13 ? pcost[x][y][15] : a == 14 ? 4 * pcost[x][y][15] : 0 );
                            int cb = pcost[x][y][b] + ( b == 13 ? pcost[x][y][15] : b == 14 ? 4 * pcost[x][y][15] : 0 );

                            if ( cb > ca )
                                dom = 0;
                            if ( cb != ca )
                                same = 0;
                        }
                    if ( dom && ( !same || b < a ) )
                        ok[a] = 0;      /* drop a; of two equal lanes keep the first */
                }
            }
            for ( a = 0; a < 15; a++ )
                if ( ok[a] )
                {
                    cand_lane[m][n] = a;
                    cand_mul[m][n] = a == 13 ? 1 : a == 14 ? 4 : 0;
                    cand_tab[m][n] = a < 13 ? a1_small[a] : a == 13 ? 16 : 24;
                    n++;
                }
            cand_n[m] = n;
        }
        for ( a = 0; a < 2; a++ )
            for ( b = 0; b < 14; b++ )
                esc_t[a][b] = esc_table ( a, 15 + ( 1 << b ) - 1 );
        esc_nb[0] = 0;
        for ( a = 1; a < 8192; a++ )
            esc_nb[a] = ( unsigned char ) ( esc_nb[a >> 1] + 1 );
    }
}

static void
huff_prepare ( HCTX * h, const int *ix, int np )
{
    int k, i, lv;

    memset ( h->ps[0], 0, sizeof ( h->ps[0] ) );
    h->lbp[0] = 0;
    for ( k = 0; k < np; k++ )
    {
        int x = ix[2 * k], y = ix[2 * k + 1];
        int m = x > y ? x : y, lb;
        int xx = x > 15 ? 15 : x, yy = y > 15 ? 15 : y;
        const int *c = pcost[xx][yy];
        const int *p0 = h->ps[k];
        int *p1 = h->ps[k + 1];

        for ( i = 0; i < A1_NL; i++ )
            p1[i] = p0[i] + c[i];
        h->mx[0][k] = m;
        if ( m <= 15 )
            lb = lbsmall[x][y];
        else
        {
            int t16 = esc_t[0][esc_nb[m - 15]], t24 = esc_t[1][esc_nb[m - 15]];

            lb = c[13] + huff_case_table[t16].linbits * c[15];
            if ( c[14] + huff_case_table[t24].linbits * c[15] < lb )
                lb = c[14] + huff_case_table[t24].linbits * c[15];
        }
        h->lbp[k + 1] = h->lbp[k] + lb;
    }
    for ( lv = 1; ( 1 << lv ) <= np; lv++ )
        for ( k = 0; k + ( 1 << lv ) <= np; k++ )
        {
            int a = h->mx[lv - 1][k], b = h->mx[lv - 1][k + ( 1 << ( lv - 1 ) )];

            h->mx[lv][k] = a > b ? a : b;
        }
}

/* count1 costs of the quads from position nmin on, both alignments (prefix sums
   start at the first quad at or after nmin; quads that would pass 576 get none) */
static void
quad_prepare ( HCTX * h, const int *ix, int nmin, int last_nz )
{
    int a, q;

    for ( a = 0; a < 2; a++ )
    {
        int q0 = ( nmin - 2 * a + 3 ) >> 2;

        if ( q0 < 0 )
            q0 = 0;
        h->qa[a][q0] = h->qb[a][q0] = 0;
        for ( q = q0; 2 * a + 4 * q + 4 <= 576 && 2 * a + 4 * q <= last_nz; q++ )
        {
            const int *p = ix + 2 * a + 4 * q;
            int v = ( ( p[0] != 0 ) << 3 ) | ( ( p[1] != 0 ) << 2 ) | ( ( p[2] != 0 ) << 1 ) | ( p[3] != 0 );
            int s = ( p[0] != 0 ) + ( p[1] != 0 ) + ( p[2] != 0 ) + ( p[3] != 0 );

            h->qa[a][q + 1] = h->qa[a][q] + a1_quada_len[v] + s;
            h->qb[a][q + 1] = h->qb[a][q] + 4 + s;
        }
    }
}

/* cheapest table for pairs [a,b) whose largest value is m; returns bits, table in *tab */
static int
seg_cost ( const HCTX * h, int a, int b, int m, int *tab )
{
    const int *pa = h->ps[a], *pb = h->ps[b];
    int ne = pb[15] - pa[15], best, bt, i, c;

    if ( m == 0 || a >= b )
    {
        *tab = 0;
        return 0;
    }
    if ( m <= 15 )
    {
        const int *cl = cand_lane[m], *cm = cand_mul[m], *ct = cand_tab[m];
        int n = cand_n[m];

        best = pb[cl[0]] - pa[cl[0]] + cm[0] * ne;
        bt = ct[0];
        for ( i = 1; i < n; i++ )
        {
            c = pb[cl[i]] - pa[cl[i]] + cm[i] * ne;
            if ( c < best )
            {
                best = c;
                bt = ct[i];
            }
        }
        *tab = bt;
        return best;
    }
    {   /* escape tables only: the smallest linbits that holds m in each family */
        int lb = esc_nb[m - 15];        /* m <= 8206 was checked */
        int t16 = esc_t[0][lb], t24 = esc_t[1][lb];

        best = pb[13] - pa[13] + huff_case_table[t16].linbits * ne;
        bt = t16;
        c = pb[14] - pa[14] + huff_case_table[t24].linbits * ne;
        if ( c < best )
        {
            best = c;
            bt = t24;
        }
        *tab = bt;
        return best;
    }
}

/* largest value of pairs [a,b) */
static int
seg_max ( const HCTX * h, int a, int b )
{
    int k, m;

    if ( a >= b )
        return 0;
    k = lg2[b - a];
    m = h->mx[k][a];
    return h->mx[k][b - ( 1 << k )] > m ? h->mx[k][b - ( 1 << k )] : m;
}

static int
seg_best ( const HCTX * h, int a, int b, int *tab )
{
    return seg_cost ( h, a, b, seg_max ( h, a, b ), tab );
}

/*
 * Cheapest coding of one granule: every big_values in [nmin, nmin + A1_WINDOW],
 * every region split, the cheapest table of each region, count1 table A or B.
 * Branch and bound keeps it exact: big(np), the cheapest big_values part for np
 * pairs, never decreases with np by less than the cheapest cost of the added
 * pairs (lbp), so a big_values value or a split whose bound cannot beat the best
 * so far is skipped without changing the result.
 */
static void
huff_best ( A1CTX * ctx, const int *ix, int ws, HRES * r, A1STATS * st )
{
    HCTX *h = &ctx->hc;
    int last_nz = -1, last_big = -1, i, j, nmin, nend, nmax, nbig;
    int Bp[23];
    int M0[23], tM0[23];
    int W[23], Wj0[23], tW[23];      /* best region0 + region1 ending at band j */
    int BL[23];         /* prefix of each band's own cheapest cost: a bound for any region */
    int bmax[23], pmax[23];     /* largest value of each band, of bands [0, j) */
    int best = A1_INF;
    int chain_np = -1, chain_lb = 0;    /* big(np) >= chain_lb + lbp[np] - lbp[chain_np] */

    memset ( r, 0, sizeof ( *r ) );
    for ( i = 575; i >= 0; i-- )
    {
        int v = ix[i];

        if ( v < 0 || v > tmax[31] )
            return;     /* r->ok = 0 */
        if ( v && last_nz < 0 )
            last_nz = i;
        if ( v > 1 && last_big < 0 )
            last_big = i;
    }
    r->ok = 1;
    if ( last_nz < 0 )
        return;         /* all zero: nothing to code */
    nmin = last_big >= 0 ? ( ( last_big + 2 ) & ~1 ) : 0;
    nend = ( last_nz + 2 ) & ~1;
    nmax = nmin + A1_WINDOW;
    if ( nmax > nend )
        nmax = nend;
    huff_prepare ( h, ix, nmax / 2 );
    quad_prepare ( h, ix, nmin, last_nz );

    for ( j = 0; j < 23; j++ )
        Bp[j] = ctx->sfbl[j] >> 1;
    if ( !ws )
    {
        int np_max = nmax / 2, t;

        BL[0] = 0;
        pmax[0] = 0;
        for ( j = 0; j < 22; j++ )
        {
            bmax[j] = Bp[j + 1] <= np_max ? seg_max ( h, Bp[j], Bp[j + 1] ) : 0;
            pmax[j + 1] = pmax[j] > bmax[j] ? pmax[j] : bmax[j];
            BL[j + 1] = BL[j] + ( Bp[j + 1] <= np_max ? seg_cost ( h, Bp[j], Bp[j + 1], bmax[j], &t ) : 0 );
        }
        for ( j = 1; j <= 16; j++ )
            if ( Bp[j] < np_max )
                M0[j] = seg_cost ( h, 0, Bp[j], pmax[j], &tM0[j] );
        for ( j = 2; j <= 22; j++ )
        {
            int j0, rm = 0;     /* rm: largest value of bands [j0, j) */

            W[j] = A1_INF;
            Wj0[j] = 0;
            tW[j] = 0;
            if ( Bp[j] >= np_max )
                continue;
            for ( j0 = j - 1; j0 > 16; j0-- )
                rm = rm > bmax[j0] ? rm : bmax[j0];
            for ( j0 = j - 1 > 16 ? 16 : j - 1; j0 >= 1 && j0 >= j - 8; j0-- )
            {
                int t, c;

                rm = rm > bmax[j0] ? rm : bmax[j0];
                if ( M0[j0] + BL[j] - BL[j0] >= W[j] )
                    continue;
                c = M0[j0] + seg_cost ( h, Bp[j0], Bp[j], rm, &t );
                if ( c < W[j] )
                {
                    W[j] = c;
                    Wj0[j] = j0;
                    tW[j] = t;
                }
            }
        }
    }

    for ( nbig = nmin; nbig <= nmax; nbig += 2 )
    {
        int np = nbig >> 1, nq = 0, ca = 0, cb = 0, c1, c1sel;
        int big, bound, r0 = 0, r1 = 0, t0 = 0, t1 = 0, t2 = 0;

        if ( last_nz >= nbig )
        {
            int a = ( nbig >> 1 ) & 1, q0 = ( nbig - 2 * a ) >> 2, q1;

            nq = ( last_nz - nbig ) / 4 + 1;
            if ( nbig + 4 * nq > 576 )
                continue;
            q1 = q0 + nq;
            ca = h->qa[a][q1] - h->qa[a][q0];
            cb = h->qb[a][q1] - h->qb[a][q0];
        }
        c1sel = cb < ca;
        c1 = c1sel ? cb : ca;
        if ( c1 >= best )
            continue;
        bound = best - c1;      /* big part needed to improve */
        if ( chain_np >= 0 && chain_lb + h->lbp[np] - h->lbp[chain_np] >= bound )
            continue;
        big = bound;

        if ( ws )
        {
            int e = np < 18 ? np : 18, ta, tb, c;

            c = seg_best ( h, 0, e, &ta ) + seg_best ( h, e, np, &tb );
            if ( c < big )
            {
                big = c;
                t0 = ta;
                t1 = tb;
            }
        }
        else
        {
            int jn = 1;         /* first band boundary at or past np */
            int smax[23];       /* largest value of pairs [Bp[j], np) */

            /* region0 covers everything */
            for ( j = 1; j <= 16; j++ )
                if ( Bp[j] >= np )
                {
                    int ta, c = seg_best ( h, 0, np, &ta );

                    if ( c < big )
                    {
                        big = c;
                        r0 = j - 1;
                        r1 = 0;
                        t0 = ta;
                        t1 = t2 = 0;
                    }
                    break;
                }
            /* region1 (two regions) or region2 (three) runs from band j to the end */
            while ( Bp[jn] < np )
                jn++;
            smax[jn - 1] = seg_max ( h, Bp[jn - 1], np );
            for ( j = jn - 2; j >= 1; j-- )
                smax[j] = smax[j + 1] > bmax[j] ? smax[j + 1] : bmax[j];
            for ( j = 1; j < jn; j++ )
            {
                /* bound for region [Bp[j], np): whole bands, then the pairs of band jn-1 */
                int lbn = BL[jn - 1] - BL[j] + h->lbp[np] - h->lbp[Bp[jn - 1]];
                int k = jn - j - 1, okB, okC, N, tN;

                okB = j <= 16 && k <= 7;        /* region1_count = k reaches np */
                okC = j >= 2 && W[j] < A1_INF;
                if ( ( !okB || M0[j] + lbn >= big ) && ( !okC || W[j] + lbn >= big ) )
                    continue;
                N = seg_cost ( h, Bp[j], np, smax[j], &tN );
                if ( okB && M0[j] + N < big )
                {
                    big = M0[j] + N;
                    r0 = j - 1;
                    r1 = k;
                    t0 = tM0[j];
                    t1 = tN;
                    t2 = 0;
                }
                if ( okC && W[j] + N < big )
                {
                    int j0 = Wj0[j];

                    big = W[j] + N;
                    r0 = j0 - 1;
                    r1 = j - j0 - 1;
                    t0 = tM0[j0];
                    t1 = tW[j];
                    t2 = tN;
                }
            }
        }
        chain_np = np;
        chain_lb = big;         /* exact if below bound, else a lower bound */
        if ( big < bound )
        {
            best = big + c1;
            r->big_values = np;
            r->region0_count = r0;
            r->region1_count = r1;
            r->table_select[0] = t0;
            r->table_select[1] = t1;
            r->table_select[2] = t2;
            r->count1table_select = c1sel;
            r->nquads = nq;
            r->bits = best;
        }
    }
    if ( best >= A1_INF )
        r->ok = 0;
    if ( st && r->big_values * 2 == nmin + A1_WINDOW )
        st->nbig_window_hits++;
}

/*====================================================================*/
/* bit writer                                                          */
/*====================================================================*/
typedef struct
{
    unsigned char *p;
    unsigned long long acc;
    int n;
    long bits;
}
BW;

/* up to 56 bits at once */
static void
bw_put ( BW * w, unsigned long long v, int nb )
{
    if ( nb <= 0 )
        return;
    w->acc = ( w->acc << nb ) | ( v & ( ( 1ull << nb ) - 1 ) );
    w->n += nb;
    w->bits += nb;
    while ( w->n >= 8 )
    {
        w->n -= 8;
        *w->p++ = ( unsigned char ) ( w->acc >> w->n );
    }
}

/* 8 bits starting at bit pos (reads one byte past them) */
static unsigned int
get8 ( const unsigned char *p, long pos )
{
    const unsigned char *q = p + ( pos >> 3 );

    return ( ( ( ( unsigned int ) q[0] << 8 ) | q[1] ) >> ( 8 - ( pos & 7 ) ) ) & 0xff;
}

static void
bw_copy ( BW * w, const unsigned char *src, long bitpos, long nbits )
{
    for ( ; nbits >= 8; nbits -= 8, bitpos += 8 )
        bw_put ( w, get8 ( src, bitpos ), 8 );
    if ( nbits > 0 )
        bw_put ( w, get8 ( src, bitpos ) >> ( 8 - nbits ), ( int ) nbits );
}

static int
bw_flush ( BW * w, unsigned char *start )
{
    if ( w->n > 0 )
    {
        *w->p++ = ( unsigned char ) ( w->acc << ( 8 - w->n ) );
        w->n = 0;
    }
    return ( int ) ( w->p - start );
}

/* one pair: code, linbits and signs in a single write (at most 47 bits) */
static void
put_pair ( BW * w, int t, int x, int y, int sx, int sy )
{
    int lb = huff_case_table[t].linbits, nb;
    int xx = x > 15 ? 15 : x, yy = y > 15 ? 15 : y;
    unsigned long long v;

    v = tcode[t][xx][yy];
    nb = tlen[t][xx][yy];
    if ( t >= 16 && xx == 15 )
    {
        v = ( v << lb ) | ( unsigned ) ( x - 15 );
        nb += lb;
    }
    if ( x )
    {
        v = ( v << 1 ) | ( sx & 1 );
        nb++;
    }
    if ( t >= 16 && yy == 15 )
    {
        v = ( v << lb ) | ( unsigned ) ( y - 15 );
        nb += lb;
    }
    if ( y )
    {
        v = ( v << 1 ) | ( sy & 1 );
        nb++;
    }
    bw_put ( w, v, nb );
}

static void
put_huff ( BW * w, const int *sfbl, const int *ix, const unsigned char *sg, int ws, const HRES * r )
{
    int k, np = r->big_values, e0, e1, q;

    if ( ws )
    {
        e0 = 18;
        e1 = A1_NPAIR;
    }
    else
    {
        int j1 = r->region0_count + r->region1_count + 2;

        e0 = sfbl[r->region0_count + 1] >> 1;
        e1 = sfbl[j1 > 22 ? 22 : j1] >> 1;
    }
    for ( k = 0; k < np; k++ )
    {
        int t = k < e0 ? r->table_select[0] : ( k < e1 ? r->table_select[1] : r->table_select[2] );

        if ( t == 0 )
            continue;   /* all-zero region */
        put_pair ( w, t, ix[2 * k], ix[2 * k + 1], sg[2 * k], sg[2 * k + 1] );
    }
    for ( q = 0; q < r->nquads; q++ )
    {
        int p = 2 * np + 4 * q, i;
        int v = ( ix[p] << 3 ) | ( ix[p + 1] << 2 ) | ( ix[p + 2] << 1 ) | ix[p + 3];
        unsigned long long c;
        int nb;

        if ( r->count1table_select )
        {
            c = v ^ 15;
            nb = 4;
        }
        else
        {
            c = a1_quada_code[v];
            nb = a1_quada_len[v];
        }
        for ( i = 0; i < 4; i++ )
            if ( ix[p + i] )
            {
                c = ( c << 1 ) | ( sg[p + i] & 1 );
                nb++;
            }
        bw_put ( w, c, nb );
    }
}

/*====================================================================*/
/* scalefactors                                                        */
/*====================================================================*/
typedef struct
{
    int g, ss, pf;
    int v[21];
    int gm[4];          /* group maxima */
    int cost;           /* part2 bits on its own (no scfsi) */
}
LCAND;

static int
lcost ( const int gm[4], int reuse, int *comp )
{
    int m1 = 0, m2 = 0;

    if ( !( reuse & 8 ) )
        m1 = gm[0];
    if ( !( reuse & 4 ) && gm[1] > m1 )
        m1 = gm[1];
    if ( !( reuse & 2 ) )
        m2 = gm[2];
    if ( !( reuse & 1 ) && gm[3] > m2 )
        m2 = gm[3];
    if ( m2 > 7 )
    {   /* cannot be coded (candidates never produce this) */
        if ( comp )
            *comp = 0;
        return A1_INF;
    }
    if ( comp )
        *comp = lcomp_tab[widthv ( m1 )][widthv ( m2 )][reuse];
    return lcost_tab[widthv ( m1 )][widthv ( m2 )][reuse];
}

static void
group_max ( const int *v, int gm[4] )
{
    int g, b;

    for ( g = 0; g < 4; g++ )
    {
        gm[g] = 0;
        for ( b = a1_grp[g]; b < a1_grp[g + 1]; b++ )
            if ( v[b] > gm[g] )
                gm[g] = v[b];
    }
}

/* all representations of a long granule that keep E[b] on the active bands */
static int
long_cands ( const int *E, const int *act, int s21, int gg, LCAND * c )
{
    int n = 0, ss, pf, b, any = 0;

    for ( b = 0; b < 21; b++ )
        any |= act[b];
    if ( !any )
    {
        memset ( &c[0], 0, sizeof ( c[0] ) );
        c[0].g = gg;    /* keeps band 21's step */
        return 1;
    }
    for ( ss = 0; ss < 2; ss++ )
        for ( pf = 0; pf < 2; pf++ )
        {
            int m = 2 << ss, rr = -1, ok = 1, g0 = 0, first = 1, g;

            for ( b = 0; b < 21; b++ )
                if ( act[b] )
                {
                    int cls = ( ( E[b] % m ) + m ) % m;
                    int need = E[b] + m * pf * a1_pretab[b];

                    if ( rr < 0 )
                        rr = cls;
                    else if ( cls != rr )
                    {
                        ok = 0;
                        break;
                    }
                    if ( first || need > g0 )
                    {
                        g0 = need;
                        first = 0;
                    }
                }
            if ( !ok )
                continue;
            if ( g0 < 0 )
                g0 = 0;
            g0 += ( ( rr - g0 % m ) % m + m ) % m;
            if ( s21 )
            {
                if ( gg < g0 || ( ( gg - rr ) % m + m ) % m )
                    continue;
                g0 = gg;
            }
            for ( g = g0; g <= 255 && n < A1_MAXCAND; g += m )
            {
                LCAND *q = &c[n];
                int bad = 0;

                q->g = g;
                q->ss = ss;
                q->pf = pf;
                for ( b = 0; b < 21; b++ )
                {
                    int v;

                    if ( !act[b] )
                    {
                        q->v[b] = 0;
                        continue;
                    }
                    v = ( g - E[b] ) / m - pf * a1_pretab[b];
                    if ( v > ( b < 11 ? 15 : 7 ) )
                    {
                        bad = 1;
                        break;
                    }
                    q->v[b] = v;
                }
                if ( bad )
                    break;
                group_max ( q->v, q->gm );
                q->cost = lcost ( q->gm, 0, NULL );
                n++;
                if ( s21 )
                    break;
            }
        }
    return n;
}

/*====================================================================*/
/* frame                                                               */
/*====================================================================*/
typedef struct
{
    int empty;          /* no nonzero value: nothing written */
    int is_short;
    int act[22];        /* long: band holds a nonzero value (21 = band without sf) */
    int act_s[3][13];   /* short: [window][band] */
    int E[22];          /* long effective steps, active bands */
    int dec_l[22];      /* decoded long scalefactors of the stock stream */
    int dec_s[3][13];
    /* result */
    int g, ss, pf, comp;
    int v[21];          /* long values written (or reused) */
    int vs[3][13];
    int part2;
    HRES h;
}
GSTATE;

static void
activity ( const A1CTX * ctx, GSTATE * s, const int *ix )
{
    const int *sfbl = ctx->sfbl, *sfbs = ctx->sfbs;
    int b, w, i;

    for ( i = 0; i < 576 && !ix[i]; i++ );
    s->empty = ( i == 576 );
    if ( !s->is_short )
    {
        for ( b = 0; b < 22; b++ )
        {
            s->act[b] = 0;
            for ( i = sfbl[b]; i < sfbl[b + 1]; i++ )
                if ( ix[i] )
                {
                    s->act[b] = 1;
                    break;
                }
        }
    }
    else
    {
        for ( b = 0; b < 13; b++ )
        {
            int wd = sfbs[b + 1] - sfbs[b];

            for ( w = 0; w < 3; w++ )
            {
                int st = 3 * sfbs[b] + w * wd;

                s->act_s[w][b] = 0;
                for ( i = st; i < st + wd; i++ )
                    if ( ix[i] )
                    {
                        s->act_s[w][b] = 1;
                        break;
                    }
            }
        }
    }
}

/* long part2 bits of the stock stream, from its own scalefac_compress and scfsi */
static int
stock_part2 ( const GR * g, int reuse )
{
    int s1 = a1_slen1[g->scalefac_compress], s2 = a1_slen2[g->scalefac_compress], grp, bits = 0;

    if ( g->window_switching_flag && g->block_type == 2 )
        return 18 * ( s1 + s2 );
    for ( grp = 0; grp < 4; grp++ )
        if ( !( reuse & ( 8 >> grp ) ) )
            bits += ( a1_grp[grp + 1] - a1_grp[grp] ) * ( grp < 2 ? s1 : s2 );
    return bits;
}

/* the values fit the regions and tables the stock side info names */
static int
stock_codable ( const GR * g, const int *ix, const int *sfbl )
{
    int k, q, e0, e1, np = g->big_values, nq = g->aux_nquads;

    if ( np < 0 || np > A1_NPAIR || nq < 0 || 2 * np + 4 * nq > 576 )
        return 0;
    if ( g->window_switching_flag )
    {
        e0 = 18;
        e1 = A1_NPAIR;
    }
    else
    {
        int j1 = g->region0_count + g->region1_count + 2;

        if ( g->region0_count < 0 || g->region0_count > 15 || g->region1_count < 0 || g->region1_count > 7 )
            return 0;
        e0 = sfbl[g->region0_count + 1] >> 1;
        e1 = sfbl[j1 > 22 ? 22 : j1] >> 1;
    }
    for ( k = 0; k < np; k++ )
    {
        int t = k < e0 ? g->table_select[0] : ( k < e1 ? g->table_select[1] : g->table_select[2] );
        int m = ix[2 * k] > ix[2 * k + 1] ? ix[2 * k] : ix[2 * k + 1];

        if ( t < 0 || t > 31 || tmax[t] < m )
            return 0;   /* also tables 4 and 14, and table 0 holding a value */
    }
    for ( q = 2 * np; q < 2 * np + 4 * nq; q++ )
        if ( ix[q] > 1 )
            return 0;
    return 1;
}

/*
 * Re-encode a stock granule exactly as its side info describes it and compare with
 * the stock bits. A match proves that a decoder of the stock stream reads these
 * scalefactors and these quantized values, which is what the re-coding keeps.
 */
static int
stock_check ( A1CTX * ctx, const GR * g, const SCALEFACT * sfv, int igr, int scfsi, const int *ix,
              const unsigned char *sg, const unsigned char *stock_main, long off )
{
    unsigned char *tmp = ctx->tmp;
    BW w;
    HRES r;
    int n = g->part2_3_length, s1, s2, b, k, grp;
    long i;

    if ( n == 0 )
        return g->big_values == 0 && g->scalefac_compress == 0;
    if ( !stock_codable ( g, ix, ctx->sfbl ) )
        return 0;
    s1 = a1_slen1[g->scalefac_compress];
    s2 = a1_slen2[g->scalefac_compress];
    w.p = tmp;
    w.acc = 0;
    w.n = 0;
    w.bits = 0;
    if ( g->window_switching_flag && g->block_type == 2 )
    {
        for ( b = 0; b < 12; b++ )
            for ( k = 0; k < 3; k++ )
                bw_put ( &w, sfv->s[k][b], b < 6 ? s1 : s2 );
    }
    else
    {
        for ( grp = 0; grp < 4; grp++ )
        {
            if ( igr == 1 && ( scfsi & ( 8 >> grp ) ) )
                continue;
            for ( b = a1_grp[grp]; b < a1_grp[grp + 1]; b++ )
                bw_put ( &w, sfv->l[b], grp < 2 ? s1 : s2 );
        }
    }
    memset ( &r, 0, sizeof ( r ) );
    r.big_values = g->big_values;
    r.region0_count = g->region0_count;
    r.region1_count = g->region1_count;
    r.table_select[0] = g->table_select[0];
    r.table_select[1] = g->table_select[1];
    r.table_select[2] = g->table_select[2];
    r.count1table_select = g->count1table_select;
    r.nquads = g->aux_nquads;
    put_huff ( &w, ctx->sfbl, ix, sg, g->window_switching_flag, &r );
    if ( w.bits != n )
        return 0;
    bw_flush ( &w, tmp );
    for ( i = 0; i + 8 <= n; i += 8 )
        if ( tmp[i >> 3] != get8 ( stock_main, off + i ) )
            return 0;
    if ( i < n && ( tmp[i >> 3] ^ get8 ( stock_main, off + i ) ) >> ( 8 - ( n - i ) ) )
        return 0;
    return 1;
}

/* re-code both granules of one channel; returns 0 if it must be copied instead */
static int
plan_channel ( A1CTX * ctx, const SIDE_INFO * si, const SCALEFACT sf[2][2], int ch,
               const A1FRAME * f, const unsigned char *stock_main, long off[2][2],
               GSTATE gs[2], int *scfsi_out, A1STATS * st )
{
    int igr, b, w;
    LCAND c0[A1_MAXCAND], c1[A1_MAXCAND];
    int n0 = 0, n1 = 0;
    int reuse_ok;

    for ( igr = 0; igr < 2; igr++ )
    {
        const GR *g = &si->gr[igr][ch];
        GSTATE *s = &gs[igr];
        int null_stock = ( g->part2_3_length == 0 );
        int c = g->scalefac_compress, s1 = a1_slen1[c], s2 = a1_slen2[c];
        int m = 2 << g->scalefac_scale;

        memset ( s, 0, sizeof ( *s ) );
        if ( g->mixed_block_flag )
            return 0;   /* never produced by Helix */
        if ( c < 0 || c > 15 ||
             !stock_check ( ctx, g, &sf[igr][ch], igr, si->scfsi[ch], f->ix[igr][ch], f->sign[igr][ch],
                            stock_main, off[igr][ch] ) )
        {
            st->stock_mismatches++;
            return 0;
        }
        s->is_short = g->window_switching_flag && g->block_type == 2;
        activity ( ctx, s, f->ix[igr][ch] );
        if ( null_stock && !s->empty )
            return 0;   /* values the stock stream does not carry */

        /* scalefactors as a decoder reads them from the stock stream */
        if ( s->is_short )
        {
            for ( w = 0; w < 3; w++ )
                for ( b = 0; b < 12; b++ )
                {
                    int v = null_stock ? 0 : sf[igr][ch].s[w][b];

                    if ( v < 0 || v >= ( 1 << ( b < 6 ? s1 : s2 ) ) )
                        return 0;
                    s->dec_s[w][b] = v;
                }
        }
        else
        {
            int grp;

            for ( grp = 0; grp < 4; grp++ )
            {
                int reused = ( igr == 1 ) && ( si->scfsi[ch] & ( 8 >> grp ) );

                for ( b = a1_grp[grp]; b < a1_grp[grp + 1]; b++ )
                {
                    int v;

                    if ( reused )
                    {
                        if ( gs[0].is_short )
                            return 0;
                        v = gs[0].dec_l[b];
                        if ( null_stock == 0 && gs[0].dec_l[b] != sf[0][ch].l[b] && gs[0].empty )
                            st->stock_scfsi_quirks++;
                    }
                    else
                    {
                        v = null_stock ? 0 : sf[igr][ch].l[b];
                        if ( v < 0 || v >= ( 1 << ( b < 11 ? s1 : s2 ) ) )
                            return 0;
                    }
                    s->dec_l[b] = v;
                }
            }
            s->dec_l[21] = 0;
            for ( b = 0; b < 22; b++ )
                s->E[b] = g->global_gain - m * ( s->dec_l[b] + g->preflag * a1_pretab[b] );
        }

        /* Huffman */
        huff_best ( ctx, f->ix[igr][ch], g->window_switching_flag, &s->h, st );
        if ( !s->h.ok )
            return 0;
    }

    /* scalefactors */
    for ( igr = 0; igr < 2; igr++ )
    {
        const GR *g = &si->gr[igr][ch];
        GSTATE *s = &gs[igr];

        s->g = g->global_gain;
        s->ss = g->scalefac_scale;
        s->pf = g->preflag;
        if ( s->empty )
        {
            s->comp = 0;
            s->part2 = 0;
            memset ( s->v, 0, sizeof ( s->v ) );
            memset ( s->vs, 0, sizeof ( s->vs ) );
            continue;
        }
        if ( s->is_short )
        {
            int m1 = 0, m2 = 0;

            for ( w = 0; w < 3; w++ )
                for ( b = 0; b < 12; b++ )
                {
                    int v = s->act_s[w][b] ? s->dec_s[w][b] : 0;

                    s->vs[w][b] = v;
                    if ( b < 6 && v > m1 )
                        m1 = v;
                    if ( b >= 6 && v > m2 )
                        m2 = v;
                }
            s->comp = scomp_tab[widthv ( m1 )][widthv ( m2 )];
            s->part2 = scost_tab[widthv ( m1 )][widthv ( m2 )];
        }
    }

    reuse_ok = !gs[0].is_short && !gs[1].is_short;
    if ( !gs[0].is_short && !gs[0].empty )
        n0 = long_cands ( gs[0].E, gs[0].act, gs[0].act[21] || ctx->keep_gain, si->gr[0][ch].global_gain, c0 );
    if ( !gs[1].is_short && !gs[1].empty )
        n1 = long_cands ( gs[1].E, gs[1].act, gs[1].act[21] || ctx->keep_gain, si->gr[1][ch].global_gain, c1 );
    if ( ( !gs[0].is_short && !gs[0].empty && n0 == 0 ) ||
         ( !gs[1].is_short && !gs[1].empty && n1 == 0 ) )
        return 0;       /* cannot happen: the stock representation is always valid */

    *scfsi_out = 0;
    if ( reuse_ok && !gs[1].empty )
    {
        /* joint choice for both granules: gr1 may reuse gr0's values group by group */
        int a, ai, cc, best = A1_INF, ba = -1, bc = -1, bS = 0;
        static const LCAND zero_cand;
        const LCAND *p0, *p1;
        int null0 = gs[0].empty;
        int c1S[A1_MAXCAND][16];        /* gr1 part2 bits for every reuse pattern */
        int fmc[A1_MAXCAND][4];         /* largest value gr0 must take in a free band */
        int idx0[A1_MAXCAND];

        if ( null0 )
            n0 = 1;     /* gr0 writes nothing: its scalefactors decode as 0 */
        for ( cc = 0; cc < n1; cc++ )
        {
            int S, grp;

            for ( S = 0; S < 16; S++ )
                c1S[cc][S] = lcost ( c1[cc].gm, S, NULL );
            for ( grp = 0; grp < 4; grp++ )
            {
                fmc[cc][grp] = 0;
                if ( !null0 )
                    for ( b = a1_grp[grp]; b < a1_grp[grp + 1]; b++ )
                        if ( gs[1].act[b] && !gs[0].act[b] && c1[cc].v[b] > fmc[cc][grp] )
                            fmc[cc][grp] = c1[cc].v[b];
            }
        }
        /* gr0 candidates by increasing cost, so the search can stop early */
        for ( a = 0; a < n0; a++ )
        {
            int k = a;

            while ( k > 0 && !null0 && c0[idx0[k - 1]].cost > c0[a].cost )
            {
                idx0[k] = idx0[k - 1];
                k--;
            }
            idx0[k] = a;
        }
        for ( ai = 0; ai < n0; ai++ )
        {
            a = idx0[ai];
            p0 = null0 ? &zero_cand : &c0[a];
            if ( p0->cost >= best )
                break;  /* gr0's own bits only grow with reuse */
            for ( cc = 0; cc < n1; cc++ )
            {
                int R = 0, grp, S;

                p1 = &c1[cc];
                for ( grp = 0; grp < 4; grp++ )
                {
                    int ok = 1;

                    for ( b = a1_grp[grp]; b < a1_grp[grp + 1] && ok; b++ )
                        if ( gs[1].act[b] &&
                             ( null0 ? p1->v[b] != 0 : ( gs[0].act[b] && p1->v[b] != p0->v[b] ) ) )
                            ok = 0;
                    if ( ok )
                        R |= 8 >> grp;
                }
                if ( p0->cost + c1S[cc][R] >= best )
                    continue;
                /* every subset of the reusable groups */
                S = R;
                for ( ;; )
                {
                    int tot = c1S[cc][S], k, gm0[4], grow = 0;

                    if ( !null0 )
                    {
                        for ( k = 0; k < 4; k++ )
                        {
                            gm0[k] = p0->gm[k];
                            if ( ( S & ( 8 >> k ) ) && fmc[cc][k] > gm0[k] )
                            {
                                gm0[k] = fmc[cc][k];
                                grow = 1;
                            }
                        }
                        tot += grow ? lcost ( gm0, 0, NULL ) : p0->cost;
                    }
                    if ( tot < best )
                    {
                        best = tot;
                        ba = a;
                        bc = cc;
                        bS = S;
                    }
                    if ( S == 0 )
                        break;
                    S = ( S - 1 ) & R;
                }
            }
        }
        if ( bc < 0 )
            return 0;
        /* apply */
        if ( !null0 )
        {
            int k, gm0[4];

            p0 = &c0[ba];
            gs[0].g = p0->g;
            gs[0].ss = p0->ss;
            gs[0].pf = p0->pf;
            memcpy ( gs[0].v, p0->v, sizeof ( gs[0].v ) );
            for ( k = 0; k < 4; k++ )
                if ( bS & ( 8 >> k ) )
                    for ( b = a1_grp[k]; b < a1_grp[k + 1]; b++ )
                        if ( gs[1].act[b] && !gs[0].act[b] )
                            gs[0].v[b] = c1[bc].v[b];
            group_max ( gs[0].v, gm0 );
            gs[0].part2 = lcost ( gm0, 0, &gs[0].comp );
        }
        p1 = &c1[bc];
        gs[1].g = p1->g;
        gs[1].ss = p1->ss;
        gs[1].pf = p1->pf;
        memcpy ( gs[1].v, p1->v, sizeof ( gs[1].v ) );
        {
            int k;

            for ( k = 0; k < 4; k++ )
                if ( bS & ( 8 >> k ) )
                    for ( b = a1_grp[k]; b < a1_grp[k + 1]; b++ )
                        gs[1].v[b] = gs[0].empty ? 0 : gs[0].v[b];     /* what a decoder will hold */
        }
        gs[1].part2 = lcost ( p1->gm, bS, &gs[1].comp );
        *scfsi_out = bS;
    }
    else
    {
        /* each long granule on its own */
        for ( igr = 0; igr < 2; igr++ )
        {
            LCAND *cl = igr ? c1 : c0;
            int n = igr ? n1 : n0, k, bk = 0;

            if ( gs[igr].is_short || gs[igr].empty )
                continue;
            for ( k = 1; k < n; k++ )
                if ( cl[k].cost < cl[bk].cost )
                    bk = k;
            gs[igr].g = cl[bk].g;
            gs[igr].ss = cl[bk].ss;
            gs[igr].pf = cl[bk].pf;
            memcpy ( gs[igr].v, cl[bk].v, sizeof ( gs[igr].v ) );
            gs[igr].part2 = lcost ( cl[bk].gm, 0, &gs[igr].comp );
        }
    }

    /* self-check: every active band decodes to the stock step */
    for ( igr = 0; igr < 2; igr++ )
    {
        GSTATE *s = &gs[igr];
        int m = 2 << s->ss, s1 = a1_slen1[s->comp], s2 = a1_slen2[s->comp];

        if ( s->empty || s->is_short )
        {
            if ( s->is_short && !s->empty )
                for ( w = 0; w < 3; w++ )
                    for ( b = 0; b < 12; b++ )
                        if ( s->vs[w][b] >= ( 1 << ( b < 6 ? s1 : s2 ) ) )
                        {
                            st->check_failures++;
                            return 0;
                        }
            continue;
        }
        for ( b = 0; b < 21; b++ )
        {
            int reused = igr == 1 && ( *scfsi_out & ( 8 >> ( b < 6 ? 0 : b < 11 ? 1 : b < 16 ? 2 : 3 ) ) );

            if ( !reused && s->v[b] >= ( 1 << ( b < 11 ? s1 : s2 ) ) )
            {
                st->check_failures++;
                return 0;
            }
            if ( s->act[b] && s->g - m * ( s->v[b] + s->pf * a1_pretab[b] ) != s->E[b] )
            {
                st->check_failures++;
                return 0;
            }
        }
        if ( s->act[21] && s->g != s->E[21] )
        {
            st->check_failures++;
            return 0;
        }
        if ( s->g < 0 || s->g > 255 )
        {
            st->check_failures++;
            return 0;
        }
    }
    return 1;
}

/*--------------------------------------------------------------------*/
static void
write_channel_gr ( BW * w, const int *sfbl, const GSTATE * s, int igr, int scfsi, const int *ix, const unsigned char *sg,
                   int ws )
{
    int b, k;

    if ( s->empty )
        return;
    if ( s->is_short )
    {
        int s1 = a1_slen1[s->comp], s2 = a1_slen2[s->comp];

        for ( b = 0; b < 12; b++ )
            for ( k = 0; k < 3; k++ )
                bw_put ( w, s->vs[k][b], b < 6 ? s1 : s2 );
    }
    else
    {
        int s1 = a1_slen1[s->comp], s2 = a1_slen2[s->comp], grp;

        for ( grp = 0; grp < 4; grp++ )
        {
            if ( igr == 1 && ( scfsi & ( 8 >> grp ) ) )
                continue;
            for ( b = a1_grp[grp]; b < a1_grp[grp + 1]; b++ )
                bw_put ( w, s->v[b], grp < 2 ? s1 : s2 );
        }
    }
    put_huff ( w, sfbl, ix, sg, ws, &s->h );
}

/*====================================================================*/
int
a1_recode_frame ( A1CTX * ctx, const SIDE_INFO * si, const SCALEFACT sf[2][2], int nchan,
                  const A1FRAME * f, const unsigned char *stock_main,
                  unsigned char *out, SIDE_INFO * si_out, A1STATS * st )
{
    A1FRAME *fz = &ctx->fz;
    GSTATE gs[2][2];    /* [ch][gr] */
    int recode[2] = { 0, 0 }, scfsi[2] = { 0, 0 };
    long off[2][2];
    long pos = 0;
    int ch, igr, attempt;
    BW w;

    *si_out = *si;
    for ( igr = 0; igr < 2; igr++ )
        for ( ch = 0; ch < nchan; ch++ )
        {
            const GR *gi = &si->gr[igr][ch];
            int nc = gi->part2_3_length ? 2 * gi->big_values + 4 * gi->aux_nquads : 0;

            /* values past the coded region are not in the stock stream: zero */
            if ( nc < 0 )
                nc = 0;
            if ( nc > 576 )
                nc = 576;       /* stock_codable() rejects such a granule */
            memcpy ( fz->ix[igr][ch], f->ix[igr][ch], nc * sizeof ( int ) );
            memset ( fz->ix[igr][ch] + nc, 0, ( 576 - nc ) * sizeof ( int ) );
            memcpy ( fz->sign[igr][ch], f->sign[igr][ch], nc );
            memset ( fz->sign[igr][ch] + nc, 0, 576 - nc );
            off[igr][ch] = pos;
            pos += gi->part2_3_length;
        }
    f = fz;

    for ( ch = 0; ch < nchan; ch++ )
    {
        int sb = si->gr[0][ch].part2_3_length + si->gr[1][ch].part2_3_length, ab;

        recode[ch] = plan_channel ( ctx, si, sf, ch, f, stock_main, off, gs[ch], &scfsi[ch], st );
        if ( recode[ch] )
        {
            int b0 = gs[ch][0].part2 + gs[ch][0].h.bits, b1 = gs[ch][1].part2 + gs[ch][1].h.bits;

            ab = b0 + b1;
            if ( ab >= sb || b0 > 4095 || b1 > 4095 )
                recode[ch] = 0;     /* no gain, or part2_3_length would overflow */
        }
    }

    for ( attempt = 0; attempt < 2; attempt++ )
    {
        int ok = 1;

        w.p = out;
        w.acc = 0;
        w.n = 0;
        w.bits = 0;
        for ( igr = 0; igr < 2; igr++ )
            for ( ch = 0; ch < nchan; ch++ )
            {
                GR *go = &si_out->gr[igr][ch];
                const GR *gi = &si->gr[igr][ch];
                long b0 = w.bits;

                *go = *gi;
                if ( !recode[ch] )
                {
                    bw_copy ( &w, stock_main, off[igr][ch], gi->part2_3_length );
                    si_out->scfsi[ch] = si->scfsi[ch];
                    continue;
                }
                {
                    GSTATE *s = &gs[ch][igr];

                    write_channel_gr ( &w, ctx->sfbl, s, igr, scfsi[ch], f->ix[igr][ch], f->sign[igr][ch],
                                       gi->window_switching_flag );
                    go->part2_3_length = ( int ) ( w.bits - b0 );
                    if ( go->part2_3_length != s->part2 + s->h.bits && !s->empty )
                        ok = 0;
                    if ( s->empty && go->part2_3_length != 0 )
                        ok = 0;
                    go->global_gain = s->g;
                    go->scalefac_scale = s->ss;
                    go->preflag = s->pf;
                    go->scalefac_compress = s->comp;
                    go->big_values = s->h.big_values;
                    go->table_select[0] = s->h.table_select[0];
                    go->table_select[1] = s->h.table_select[1];
                    go->table_select[2] = s->h.table_select[2];
                    if ( !gi->window_switching_flag )
                    {
                        go->region0_count = s->h.region0_count;
                        go->region1_count = s->h.region1_count;
                    }
                    go->count1table_select = s->h.count1table_select;
                    si_out->scfsi[ch] = scfsi[ch];
                }
            }
        if ( ok )
            break;
        /* a bit count did not match its model: copy the whole frame */
        st->check_failures++;
        recode[0] = recode[1] = 0;
    }

    /* statistics */
    st->frames++;
    for ( ch = 0; ch < nchan; ch++ )
    {
        if ( !recode[ch] )
            st->copied_channels++;
        for ( igr = 0; igr < 2; igr++ )
        {
            const GR *gi = &si->gr[igr][ch];
            int p2s = gi->part2_3_length ? stock_part2 ( gi, igr ? si->scfsi[ch] : 0 ) : 0;

            st->stock_bits += gi->part2_3_length;
            st->stock_part2 += p2s;
            st->stock_part3 += gi->part2_3_length - p2s;
            st->a1_bits += si_out->gr[igr][ch].part2_3_length;
            if ( recode[ch] )
            {
                st->a1_part2 += gs[ch][igr].empty ? 0 : gs[ch][igr].part2;
                st->a1_part3 += gs[ch][igr].empty ? 0 : gs[ch][igr].h.bits;
            }
            else
            {
                st->a1_part2 += p2s;
                st->a1_part3 += gi->part2_3_length - p2s;
            }
        }
    }

    return bw_flush ( &w, out );
}

/*====================================================================*/
A1CTX *
a1_create ( int sr_index, int keep_global_gain )
{
    A1CTX *c;

    a1_tables_init (  );
    c = ( A1CTX * ) calloc ( 1, sizeof ( A1CTX ) );
    if ( c == NULL )
        return NULL;
    if ( sr_index < 0 || sr_index > 2 )
        sr_index = 0;
    c->sfbl = a1_sfb_l[sr_index];
    c->sfbs = a1_sfb_s[sr_index];
    c->keep_gain = keep_global_gain != 0;
    return c;
}

void
a1_destroy ( A1CTX * c )
{
    free ( c );
}
