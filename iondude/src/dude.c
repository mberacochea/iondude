#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dude.h"

// max difference allowed on the matrices validation
#define MATRIX_TOLERANCE 1e-3

static void matrix_error(char *header, char *message) {
  fprintf(stderr, "Invalid matrix '%.*s': %s\n", (int)strcspn(header, "\r\n"),
          header, message);
  exit(1);
}

// print the range or the quality score of a slot
static void print_slot(FILE *file, dude_t *dude, int slot) {
  if (dude->type == QUAL) {
    fprintf(file, "quality %d", slot);
  } else {
    fprintf(file, "range %d-%d", slot * dude->range, (slot + 1) * dude->range);
  }
}

// read the rows of a dim x dim matrix, the column names line first
// each row starts with the row name: A\t0.98\t0.01\t...
static void read_matrix_rows(FILE *file, float *matrix, int dim, char *header) {
  char line[1024];

  if (fgets(line, sizeof line, file) == NULL) {
    matrix_error(header, "missing column names");
  }

  for (int i = 0; i < dim; i++) {
    if (fgets(line, sizeof line, file) == NULL) {
      matrix_error(header, "missing rows");
    }
    // skip the row name
    char *p = strchr(line, '\t');
    if (p == NULL) {
      matrix_error(header, "row without values");
    }
    for (int j = 0; j < dim; j++) {
      char *end;
      matrix[i * dim + j] = strtof(p, &end);
      if (end == p) {
        matrix_error(header, "row with missing or invalid values");
      }
      p = end;
    }
  }
}

// get the slot for the matrix header (#SF 0-20 ...)
// -1 if the matrix is not used, ALL_SLOTS if the header has no range
int get_matrix_slot(dude_t *dude, char *header) {
  int min, max;
  if (sscanf(header, "#%*s %d-%d", &min, &max) != 2) {
    // one matrix for all the positions or quality scores
    return ALL_SLOTS;
  }

  if (dude->type == QUAL) {
    // one matrix per quality score: "#SF 20-21" is the matrix for Phred 20
    return (min >= 0 && min < MAX_QUAL) ? min : -1;
  }

  if (max - min != dude->range || min % dude->range != 0) {
    matrix_error(header, "the range size is different from the -r value");
  }

  // ranges past the max read length are not used
  int slot = min / dude->range;
  return slot < dude->slots ? slot : -1;
}

void load_matrices(char *estimation, dude_t *dude) {
  // this method load the matrices obtained with the estimator
  // currently we have 2 sets of matrices
  // substitutions and homopolymers changes
  //
  // for substitutions the estimation file contains the following lines
  // #S   - Substitutions transition matrix counters (this is not used)
  // #SF  - Substitutions transition matrix frequencies
  // #SFI - Substitutions transition matrix frequencies inverted
  //
  // for homopolymers the patter is the same #I, #IF and #IFI
  //
  // each matrix has a range, the position in the read (POSITION type)
  // or the quality score (QUAL type):
  // #SF 0-20 - this means from 0 to 20.
  // a matrix without a range (#SF) is used for all the positions or quality
  // scores

  // open file
  FILE *estimation_file = fopen(estimation, "r");

  if (!estimation_file) {
    fprintf(stderr, "Could not open the Estimator matrices file '%s'.\n",
            estimation);
    exit(1);
  }

  char line[1024];

  printf("Loading matrices...\n");

  // headers (##) and the counters (#S, #I) are not used
  while (fgets(line, sizeof line, estimation_file) != NULL) {
    if (strncmp(line, "#SF", 3) == 0) {
      load_substitution_matrix(dude, estimation_file, line,
                               strncmp(line, "#SFI", 4) != 0);
    } else if (strncmp(line, "#IF", 3) == 0) {
      load_homopolymers_matrix(dude, estimation_file, line,
                               strncmp(line, "#IFI", 4) != 0);
    }
  }

  fclose(estimation_file);

  printf("Finished loading matrices\n");

  complete_matrices(dude);

  printf("Validating...\n");

  for (int z = 0; z < dude->slots; z++) {
    validate_matrices(dude, z, &dude->subs_channel[z][0][0],
                      &dude->subs_channel_i[z][0][0], SUBS_MATRIX_DIM);
    validate_matrices(dude, z, &dude->hom_channel[z][0][0],
                      &dude->hom_channel_i[z][0][0], HOM_MATRIX_DIM);
  }

  printf("Validation completed\n");
}

// find the closest slot with the matrices loaded from the file
static int nearest_loaded_slot(dude_t *dude, int slot, int flags) {
  for (int d = 1; d < dude->slots; d++) {
    if (slot - d >= 0 && (dude->loaded[slot - d] & flags) == flags)
      return slot - d;
    if (slot + d < dude->slots && (dude->loaded[slot + d] & flags) == flags)
      return slot + d;
  }
  return -1;
}

void complete_matrices(dude_t *dude) {
  // every range needs a matrix in the file
  // quality scores without a matrix use the matrix of the closest quality score
  for (int z = 0; z < dude->slots; z++) {
    int subs = LOADED_SUBS | LOADED_SUBS_I;
    int hom = LOADED_HOM | LOADED_HOM_I;

    if ((dude->loaded[z] & subs) != subs) {
      int from = (dude->type == QUAL) ? nearest_loaded_slot(dude, z, subs) : -1;
      if (from < 0) {
        fprintf(stderr, "Missing substitutions matrices (#SF and #SFI) for ");
        print_slot(stderr, dude, z);
        fprintf(stderr, ".\n");
        exit(1);
      }
      memcpy(dude->subs_channel[z], dude->subs_channel[from],
             sizeof(subs_matrix_type));
      memcpy(dude->subs_channel_i[z], dude->subs_channel_i[from],
             sizeof(subs_matrix_type));
      memcpy(dude->subs_schur_product[z], dude->subs_schur_product[from],
             sizeof(subs_schur_product));
    }

    if ((dude->loaded[z] & hom) != hom) {
      int from = (dude->type == QUAL) ? nearest_loaded_slot(dude, z, hom) : -1;
      if (from < 0) {
        fprintf(stderr, "Missing homopolymers matrices (#IF and #IFI) for ");
        print_slot(stderr, dude, z);
        fprintf(stderr, ".\n");
        exit(1);
      }
      memcpy(dude->hom_channel[z], dude->hom_channel[from],
             sizeof(hom_matrix_type));
      memcpy(dude->hom_channel_i[z], dude->hom_channel_i[from],
             sizeof(hom_matrix_type));
      memcpy(dude->hom_schur_product[z], dude->hom_schur_product[from],
             sizeof(hom_schur_product));
    }
  }
}

void validate_matrices(dude_t *dude, int slot, float *channel, float *inverse,
                       int dim) {
  // the channel rows must sum to 1 and channel x inverse must be the identity
  for (int i = 0; i < dim; i++) {
    float sum = 0;
    for (int j = 0; j < dim; j++) {
      sum += channel[i * dim + j];

      float product = 0;
      for (int n = 0; n < dim; n++) {
        product += channel[i * dim + n] * inverse[n * dim + j];
      }
      if (fabsf(product - (i == j)) > MATRIX_TOLERANCE) {
        fprintf(stderr, "The inverted matrix (%dx%d) for ", dim, dim);
        print_slot(stderr, dude, slot);
        fprintf(stderr, " is not the inverse of the matrix.\n");
        exit(1);
      }
    }
    if (fabsf(sum - 1) > MATRIX_TOLERANCE) {
      fprintf(stderr, "The matrix (%dx%d) for ", dim, dim);
      print_slot(stderr, dude, slot);
      fprintf(stderr, ": row %d doesn't sum to 1.\n", i);
      exit(1);
    }
  }
}

void load_substitution_matrix(dude_t *dude, FILE *file, char *header,
                              int isChannel) {

  subs_matrix_type matrix;
  read_matrix_rows(file, &matrix[0][0], SUBS_MATRIX_DIM, header);

  int first = get_matrix_slot(dude, header), last = first;
  if (first == ALL_SLOTS) {
    first = 0;
    last = dude->slots - 1;
  } else if (first < 0) {
    return;
  }

  for (int position = first; position <= last; position++) {
    if (isChannel) {
      memcpy(dude->subs_channel[position], matrix, sizeof matrix);
      dude->loaded[position] |= LOADED_SUBS;
    } else {
      memcpy(dude->subs_channel_i[position], matrix, sizeof matrix);
      dude->loaded[position] |= LOADED_SUBS_I;
      continue;
    }

    // Cache the Schur product (pairwise product)
    // for each possible base in X
    for (int observed_symb = 0; observed_symb < SUBS_MATRIX_DIM;
         observed_symb++) {
      // build the columns for that possibility
      for (int x = 0; x < SUBS_MATRIX_DIM; x++) {
        for (int y = 0; y < SUBS_MATRIX_DIM; y++) {
          dude->subs_schur_product[position][observed_symb][x][y] =
              dude->subs_channel[position][y][observed_symb] *
              dude->subs_loss[y][x];
        }
      }
    }
  }
}

void load_homopolymers_matrix(dude_t *dude, FILE *file, char *header,
                              int isChannel) {

  hom_matrix_type matrix;
  read_matrix_rows(file, &matrix[0][0], HOM_MATRIX_DIM, header);

  int first = get_matrix_slot(dude, header), last = first;
  if (first == ALL_SLOTS) {
    first = 0;
    last = dude->slots - 1;
  } else if (first < 0) {
    return;
  }

  for (int position = first; position <= last; position++) {
    if (isChannel) {
      memcpy(dude->hom_channel[position], matrix, sizeof matrix);
      dude->loaded[position] |= LOADED_HOM;
    } else {
      memcpy(dude->hom_channel_i[position], matrix, sizeof matrix);
      dude->loaded[position] |= LOADED_HOM_I;
      continue;
    }

    // Cache the Schur product
    for (int f = 0; f < HOM_MATRIX_DIM; f++) {
      for (int j = 0; j < HOM_MATRIX_DIM; j++) {
        for (int n = 0; n < HOM_MATRIX_DIM; n++) {
          dude->hom_schur_product[position][f][j][n] =
              dude->hom_channel[position][n][f] * dude->hom_loss[n][j];
        }
      }
    }
  }
}

// create a new context helper
context_helper_t *create_helper(void) {
  context_helper_t *helper =
      (context_helper_t *)malloc(sizeof(context_helper_t));
  if (helper == NULL) {
    error_exit();
  }
  return helper;
}

// Create a new DUDE params struct
dude_t *create_dude_params(int k, int h, int read_len, int range, char *type) {

  dude_t *dude = (dude_t *)malloc(sizeof(dude_t));

  if (dude == NULL) {
    error_exit();
  }

  // substitutions context length
  dude->k = k;

  // indels context length
  dude->h = h;

  // mode
  dude->type = POSITION;
  if (strncmp(type, "qual", 5) == 0) {
    dude->type = QUAL;
  }

  // The range doesn't matter for QUAL denoise type
  // the last range can be shorter, ie. read_len 150 and range 20 -> 8 ranges
  dude->range = range;
  dude->range_max = (read_len + range - 1) / range;

  // ======================================
  //   Substitutions Hamming loss function
  // ======================================
  for (int i = 0; i < SUBS_MATRIX_DIM; i++) {
    for (int j = 0; j < SUBS_MATRIX_DIM; j++) {
      dude->subs_loss[i][j] = i != j;
    }
  }

  // ======================================
  //   Homopolymers Hamming loss function
  // ======================================
  for (int i = 0; i < HOM_MATRIX_DIM; i++) {
    for (int j = 0; j < HOM_MATRIX_DIM; j++) {
      dude->hom_loss[i][j] = i != j;
    }
  }

  dude->slots = (dude->type == QUAL) ? MAX_QUAL : dude->range_max;
  int size = dude->slots;

  dude->loaded = calloc(size, sizeof(int));

  dude->subs_channel = calloc(size, sizeof(subs_matrix_type));
  dude->subs_channel_i = calloc(size, sizeof(subs_matrix_type));
  dude->subs_schur_product = calloc(size, sizeof(subs_schur_product));

  dude->hom_channel = calloc(size, sizeof(hom_matrix_type));
  dude->hom_channel_i = calloc(size, sizeof(hom_matrix_type));
  dude->hom_schur_product = calloc(size, sizeof(hom_schur_product));

  if (!dude->loaded || !dude->subs_channel || !dude->subs_channel_i ||
      !dude->subs_schur_product || !dude->hom_channel || !dude->hom_channel_i ||
      !dude->hom_schur_product) {
    error_exit();
  }

  return dude;
}

int get_index_for_matrix(int position_in_read, dude_t *dude) {
  int p = position_in_read / dude->range;
  if (p >= dude->range_max) {
    // positions past the max read length use the last range
    return dude->range_max - 1;
  } else {
    return p;
  }
}

void destroy_helper(context_helper_t *helper) {
  if (helper == NULL)
    return;
  // free(helper->counter);
  // free(helper->counters_inv);
  free(helper);
}

void error_exit() {
  puts("Error on malloc");
  exit(-1);
}
