#include "sketch_rearrange.h"
#include <err.h>

uint64_t ctxmask, tupmask, ho_mask_len, hc_mask_len, io_mask_len, ho_mask_left, hc_mask_left, io_mask, hc_mask_right, ho_mask_right;
uint8_t iolen, klen, hclen, holen;
uint32_t gid_mask;
bitslen_t Bitslen;
void const_comask_init(minco_sketch_stat_t *minco_stat)
{
    //  init all public vars ;
    klen = minco_stat->klen;
    if (minco_stat->coden_len > 0)
    {
        Bitslen.ctx = 4 * minco_stat->coden_len;
        Bitslen.obj = 2 * klen - Bitslen.ctx;
        Bitslen.gid = GID_NBITS;
        gid_mask = (1U << GID_NBITS) - 1;
        if (minco_stat_needs_ctxobj96(minco_stat) ||
            minco_stat_needs_ctxgidobj128(minco_stat)) {
            ctxmask = 0;
            tupmask = 0;
            ho_mask_len = hc_mask_len = io_mask_len = 0;
            ho_mask_left = hc_mask_left = io_mask = hc_mask_right = ho_mask_right = 0;
            hclen = holen = iolen = 0;
            return;
        }
        ctxmask = generate_coden_pattern64();
    }
    else
    {
        holen = minco_stat->holen;
        hclen = minco_stat->hclen;
        iolen = klen - 2 * (hclen + holen);

        Bitslen.ctx = 4 * hclen;
        ho_mask_len = holen == 0 ? 0 : UINT64_MAX >> (64 - 2 * holen);
        hc_mask_len = hclen == 0 ? 0 : UINT64_MAX >> (64 - 2 * hclen);
        io_mask_len = iolen == 0 ? 0 : UINT64_MAX >> (64 - 2 * iolen); // UINT64_MAX >> 64 is undefined  not 0

        ho_mask_left = ho_mask_len << (2 * (klen - holen));          //(2*(holen + hclen + iolen + hclen));
        hc_mask_left = hc_mask_len << (2 * (holen + hclen + iolen)); // 2*(holen + hclen + iolen));
        io_mask = io_mask_len << (2 * (holen + hclen));
        hc_mask_right = hc_mask_len << (2 * holen);
        ho_mask_right = ho_mask_len;

        ctxmask = hc_mask_left | hc_mask_right;
    }
    Bitslen.obj = 2 * klen - Bitslen.ctx;
    Bitslen.gid = GID_NBITS;
    tupmask = klen == 0 ? 0 : UINT64_MAX >> (64 - 2 * klen);
    gid_mask = (1U << GID_NBITS) - 1;
}

void sketch64_2ctxobj64(uint64_t *sketch64, uint32_t arrlen)
{
#pragma omp parallel for num_threads(32) schedule(guided)
    for (uint32_t ri = 0; ri < arrlen; ri++)
    {
        sketch64[ri] = uint64kmer2generic_ctxobj(sketch64[ri]);
    }
}

uint96_t *sketch64_2uint96co(uint64_t *sketch_index, uint64_t *sketch64, int infile_num, uint32_t arrlen)
{
    //  int iolen = klen - 2*(hclen + holen);
    uint96_t *uint96co = malloc(arrlen * sizeof(uint96_t));
#pragma omp parallel for num_threads(32) schedule(guided)
    for (uint32_t rn = 0; rn < infile_num; rn++)
    {
        for (uint32_t ri = sketch_index[rn]; ri < sketch_index[rn + 1]; ri++)
        {
            uint96co[ri] = uint64_kmer2uint96(sketch64[ri], rn);
        }
    }
    return uint96co;
}

ctxgidobj_t *ctxobj64_2ctxgidobj(uint64_t *sketch_index, uint64_t *ctxobj64, int infile_num, uint32_t arrlen)
{
    minco_sketch_stat_t active_stat = {
        .coden_len = 0,
        .klen = klen,
        .hclen = hclen,
        .holen = holen,
    };
    active_stat.coden_len = Bitslen.ctx ? (int)(Bitslen.ctx / 4u) : 0;
    if (Bitslen.ctx + Bitslen.obj > 64u || Bitslen.ctx + GID_NBITS > 64u || Bitslen.obj > 32u)
        minco_reject_extended_ctxobj_if_needed(__func__, &active_stat,
                                               "packed ctxobj64 to ctxgid64obj32 conversion");
    ctxgidobj_t *ctxgidobj = malloc(arrlen * sizeof(ctxgidobj_t));
    uint64_t obj_len_mask = (1LU << Bitslen.obj) - 1;
#pragma omp parallel for num_threads(32) schedule(guided)
    for (uint32_t rn = 0; rn < infile_num; rn++)
    {
        for (uint32_t ri = sketch_index[rn]; ri < sketch_index[rn + 1]; ri++)
        {
            ctxgidobj[ri] = uint64_ctxobj2ctxgidobj96(ctxobj64[ri], rn, obj_len_mask);
        }
    }
    return ctxgidobj;
}

ctxgidobj128_t *ctxobj96_2ctxgidobj128(uint64_t *sketch_index, ctxobj96_t *ctxobj96,
                                       int infile_num, uint32_t arrlen)
{
    ctxgidobj128_t *ctxgidobj = malloc((size_t)arrlen * sizeof(ctxgidobj[0]));
    if (!ctxgidobj)
        err(EXIT_FAILURE, "%s(): OOM ctxgidobj128", __func__);
#pragma omp parallel for num_threads(32) schedule(guided)
    for (uint32_t rn = 0; rn < (uint32_t)infile_num; rn++)
    {
        for (uint32_t ri = sketch_index[rn]; ri < sketch_index[rn + 1]; ri++)
        {
            ctxgidobj[ri] = ctxgidobj128_make(ctxobj96[ri].ctx, rn, ctxobj96[ri].obj);
        }
    }
    return ctxgidobj;
}

uint64kmer2generic_ctxobj_fn uint64kmer2generic_ctxobj = NULL;
void set_uint64kmer2generic_ctxobj(bool is_coden_ctxobj_pattern)
{
    
    uint64kmer2generic_ctxobj = is_coden_ctxobj_pattern ? reorder_unituple_by_coden_pattern64 : uint64_kmer2ctxobj;
}
