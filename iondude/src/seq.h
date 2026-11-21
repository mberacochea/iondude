#ifndef SEQ
#define SEQ

#include "../lib/klib/kseq.h"
#include <stdint.h>
#include <zlib.h>

#define BUFFER_SIZE 409
__KS_TYPE(gzFile)

__KS_BASIC(gzFile, BUFFER_SIZE)

__KSEQ_TYPE(gzFile)

__KSEQ_BASIC(static, gzFile)

#define A_HOM 40
#define C_HOM 50
#define G_HOM 60
#define T_HOM 70

#define Q_OFFSET 0
#define Q_MIN 1
#define Q_MAX 2

static const char typenames[4][10] = {
    {"Phred"}, {"Sanger"}, {"Solexa"}, {"Illumina"}};

// quality score
static const int quality_constants[4][3] = {
    /* offset, min, max */
    {0, 4, 60},
    /* PHRED */    // Default for most technologies.
    {33, 33, 126}, /* SANGER */
    {64, 58, 112},
    /* SOLEXA; this is an approx; the transform is non-linear */ // Deprecated
    {64, 64, 110}                                                /* ILLUMINA */
};

// get quality values
// Taken from https://github.com/najoshi/sickle
int get_quality_num(char qualchar, int qualtype, int pos);

// ================================ //
//      Code and Decode related
// ================================ //
// 2-bit code of a base: A 0, C 1, G 2, T 3, -1 for any other byte
// decode with "ACGT"[code]
extern int8_t base_code[256];

// fill base_code, call it once at the start
void init_base_codes(void);

int is_acgt_seq(const char *seq);

// length is the homopolymer length - 1 (0 to HOM_MATRIX_DIM - 1)
char encode_h_base(char base, int length);

// 6-bit symbol of an encoded homopolymer: the base code (2 bits) and the
// length - 1 (4 bits). -1 for any other value
int hom_symbol(char encoded);

kseq_t *encode_homopolymers_seq(kseq_t *seq);

kseq_t *decode_homopolymers_seq(kseq_t *seq);

// copy a kseq_t
// it wont copy seq->f (the stream)
kseq_t *kseq_dup(kseq_t *seq);

void kseq_free(kseq_t *ks);

char *append_str(const char *str1, const char *str2);

#endif
