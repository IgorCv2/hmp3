/*
 * a1pack.h - lossless re-coding of MPEG-1 Layer III granules (roadmap stage A1).
 *
 * The encoder decides the quantized values and step sizes as usual. a1pack then
 * re-codes every frame for the fewest bits while keeping every quantized value
 * and every band's effective step size, so decoders output the same PCM:
 *
 *   Huffman data   - big_values/count1 boundary (every position up to 80 values
 *                    past the smallest legal one), region0/region1 split, the
 *                    cheapest table per region, count1 table A or B;
 *   scalefactors   - scalefac_scale and preflag re-expressed (global_gain too,
 *                    unless kept, see a1_create), any value in bands that hold
 *                    no nonzero value, the cheapest scalefac_compress, and scfsi
 *                    chosen for both granules of a channel together.
 *
 * Both searches are exact: the result is the cheapest coding among these choices.
 * Each granule is first re-encoded from the stock side info and compared with
 * the stock bits; a channel that does not match, or whose re-coding would not be
 * smaller, is copied bit for bit.
 *
 * This file is a contribution to the Helix MP3 encoder and is distributed under
 * the same terms as the files it accompanies (RPSL/RCSL, see hmp3/LICENSE.txt).
 */
#ifndef _A1PACK_H_
#define _A1PACK_H_

#include "l3e.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct
{
    int ix[2][2][576];                  /* [gr][ch] quantized magnitudes, coded order */
    unsigned char sign[2][2][576];      /* [gr][ch] sign bits */
}
A1FRAME;

typedef struct
{
    double frames;
    double stock_bits;          /* part2+part3 bits as Helix coded them */
    double a1_bits;             /* part2+part3 bits after re-coding */
    double stock_part2, a1_part2;
    double stock_part3, a1_part3;
    double copied_channels;     /* channels kept bit for bit */
    double stock_mismatches;    /* stock bits not reproduced from its side info (kept stock) */
    double stock_scfsi_quirks;  /* scfsi reuse from an empty granule-0 (decoded as zeros) */
    double check_failures;      /* self-check rejected a re-coding (kept stock) */
    double nbig_window_hits;    /* best big_values at the edge of the search window */
    /* filled by the encoder */
    double stock_stream_bytes, stock_stream_frames;
    double a1_stream_bytes, a1_stream_frames;
}
A1STATS;

typedef struct A1CTX A1CTX;    /* per-encoder state, opaque */

/*
 * One per encoder, MPEG-1 only; NULL if out of memory. sr_index as in the frame
 * header (0 44.1k, 1 48k, 2 32k).
 *
 * keep_global_gain = 1 keeps every granule's global_gain as coded and re-expresses
 * only scalefac_scale, preflag and the scalefactors around it. Then each band's
 * step exponent, as global_gain plus the integer (1+scalefac_scale)*(sf+preflag*
 * pretab), is unchanged term by term, so even decoders that build the step from
 * two floating-point factors (such as minimp3) output identical PCM. With 0
 * global_gain may move too, which saves a little more: the steps are still the
 * same, so decoders that work from the integer step (such as ffmpeg and mpg123)
 * still output identical PCM, while two-factor float decoders differ by float
 * rounding (below 1e-6 of full scale).
 *
 * The code tables shared by all contexts are built by the first call; make that
 * call before starting encoder threads. A context serves one encoder at a time.
 */
A1CTX *a1_create ( int sr_index, int keep_global_gain );
void a1_destroy ( A1CTX * c );

/*
 * Re-code one frame. si/sf hold the encoder's side info and scalefactors for both
 * granules, f the quantized values, stock_main the encoder's own main data for this
 * frame. Each granule is first re-encoded exactly as the stock side info describes
 * it and compared with stock_main bit for bit, so the re-coding starts from what a
 * decoder reads; a channel that does not match, or would not get smaller, is copied
 * from stock_main. Writes the new main data to out (byte aligned) and the new side
 * info to si_out. Returns the number of bytes, never more than the stock frame's.
 */
int a1_recode_frame ( A1CTX * c, const SIDE_INFO * si, const SCALEFACT sf[2][2], int nchan,
                      const A1FRAME * f, const unsigned char *stock_main,
                      unsigned char *out, SIDE_INFO * si_out, A1STATS * st );

#ifdef __cplusplus
}
#endif

#endif /* _A1PACK_H_ */
