#include <stdio.h>

#include "dude.h"
#include "seq.h"

int get_quality_num(char qualchar, int qualtype, int pos) {
  /*
    Return the adjusted quality, depending on quality type.
    Note that this uses the array in sickle.h, which *approximates*
    the SOLEXA (pre-1.3 pipeline) qualities as linear. This is
    inaccurate with low-quality bases.
  */
  int qual_value = (int)qualchar;

  if (qual_value < quality_constants[qualtype][Q_MIN] ||
      qual_value > quality_constants[qualtype][Q_MAX]) {
    fprintf(stderr,
            "ERROR: Quality value (%d) does not fall within correct range for "
            "%s encoding.\n",
            qual_value, typenames[qualtype]);
    fprintf(stderr, "Range for %s encoding: %d-%d\n", typenames[qualtype],
            quality_constants[qualtype][Q_MIN],
            quality_constants[qualtype][Q_MAX]);
    fprintf(stderr, "Quality char: '%c'\n", qualchar);
    fprintf(stderr, "Quality position: %d\n", pos + 1);
    exit(1);
  }
  return (qual_value - quality_constants[qualtype][Q_OFFSET]);
}

// encode a substitutions base (A -> 0)
int encode_s_base(char base) {
  switch (base) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return -1;
  }
}

// decode a substitutions base (0 -> A)
char decode_s_base(int index) {
  switch (index) {
  case 0:
    return 'A';
  case 1:
    return 'C';
  case 2:
    return 'G';
  case 3:
    return 'T';
  default:
    return 'N';
  }
}

// check that the sequence has only A, C, G or T
int is_acgt_seq(const char *seq) { return seq[strspn(seq, "ACGT")] == '\0'; }

// Encode the Homopolymers Run-Length encoded
// the sequence must have only A, C, G or T (see is_acgt_seq)
kseq_t *encode_homopolymers_seq(kseq_t *seq) {

  // copy
  kseq_t *res = kseq_dup(seq);

  size_t i = 0, p = 0;

  while (i < seq->seq.l) {
    // runs longer than HOM_MATRIX_DIM are split, ie. 25 A -> 10 A, 10 A, 5 A
    size_t run = 1;
    while (i + run < seq->seq.l && seq->seq.s[i + run] == seq->seq.s[i] &&
           run < HOM_MATRIX_DIM) {
      run++;
    }
    res->seq.s[p] = encode_h_base(seq->seq.s[i], run - 1);
    // save the first position on the homopolymer
    res->qual.s[p] = seq->qual.s[i];
    p++;
    i += run;
  }

  res->seq.s[p] = '\0';
  res->qual.s[p] = '\0';
  res->seq.l = p;
  res->qual.l = p;

  return res;
}

// encode an homopolymer length
char encode_h_base(char base, int length) {

  int hom = A_HOM;
  switch (base) {
  case 'A':
    hom = A_HOM;
    break;
  case 'C':
    hom = C_HOM;
    break;
  case 'G':
    hom = G_HOM;
    break;
  case 'T':
    hom = T_HOM;
    break;
  }
  return (char)hom + length;
}

// decode the RLE
kseq_t *decode_homopolymers_seq(kseq_t *seq) {

  kseq_t *res = kseq_dup(seq);

  if (res == NULL) {
    fprintf(stderr, "decode_homopolymers_seq. Seq '%s'", seq->name.s);
    exit(1);
  }

  // get the sequence length on decode
  int seq_len = 0, offset = 0;

  for (int i = 0; seq->seq.l > i; i++) {
    int value = seq->seq.s[i];
    if (value >= T_HOM) {
      offset = T_HOM;
    } else if (value >= G_HOM) {
      offset = G_HOM;
    } else if (value >= C_HOM) {
      offset = C_HOM;
    } else if (value >= A_HOM) {
      offset = A_HOM;
    }
    seq_len += (value - offset) + 1;
  }

  // resize if necessary
  if (seq_len > res->seq.l) {
    char *ts = (char *)realloc(res->seq.s, seq_len + 1);
    if (ts == NULL) {
      puts("Error (re)allocating memory");
      exit(1);
    }
    res->seq.s = ts;
    char *tq = (char *)realloc(res->qual.s, seq_len + 1);
    if (tq == NULL) {
      puts("Error (re)allocating memory");
      exit(1);
    }
    res->qual.s = tq;
  }

  offset = 0;
  char base = '~';

  for (int j = 0, p = 0; seq->seq.l > j; j++, p++) {
    int value = seq->seq.s[j];
    if (value >= T_HOM) {
      base = 'T';
      offset = T_HOM;
    } else if (value >= G_HOM) {
      base = 'G';
      offset = G_HOM;
    } else if (value >= C_HOM) {
      base = 'C';
      offset = C_HOM;
    } else if (value >= A_HOM) {
      base = 'A';
      offset = A_HOM;
    }

    // get the homopolymer length
    int homopolymer_size = value - offset;

    if (homopolymer_size > 0) {
      // write as many chars as required
      for (int k = 0; k <= homopolymer_size; k++) {
        res->seq.s[p + k] = base;
        res->qual.s[p + k] = seq->qual.s[j];
      }
      // move forward
      p += homopolymer_size;
    } else {
      // print char decoded
      res->seq.s[p] = base;
      res->qual.s[p] = seq->qual.s[j];
    }
  }

  res->seq.l = seq_len;
  res->qual.l = seq_len;
  res->seq.s[seq_len] = '\0';
  res->qual.s[seq_len] = '\0';

  return res;
}

void copy_kstring(kstring_t *origin, kstring_t *dest) {
  if (origin->l > 0) {
    dest->s = strdup(origin->s);
    dest->l = origin->l;
    dest->m = origin->m;
  } else {
    // initialize
    dest->s = (char *)calloc(1, 1);
    dest->l = 0;
    dest->m = 1;
  }
}

// copy a kseq
kseq_t *kseq_dup(kseq_t *seq) {
  if (seq == NULL)
    return NULL;
  kseq_t *r = (kseq_t *)malloc(sizeof(kseq_t));
  if (r == NULL) {
    fprintf(stderr, "Out of memory.");
    exit(1);
  }
  copy_kstring(&seq->seq, &r->seq);
  copy_kstring(&seq->qual, &r->qual);
  copy_kstring(&seq->name, &r->name);
  copy_kstring(&seq->comment, &r->comment);
  r->last_char = seq->last_char;
  return r;
}

// free the memory of a kseq dup
void kseq_free(kseq_t *ks) {
  if (!ks)
    return;
  free(ks->name.s);
  free(ks->comment.s);
  free(ks->seq.s);
  free(ks->qual.s);
  free(ks);
}

char *append_str(const char *str1, const char *str2) {
  char *new_str = malloc(strlen(str1) + strlen(str2) + 1);
  strcpy(new_str, str1);
  strcat(new_str, str2);
  return new_str;
}