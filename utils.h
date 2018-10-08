#ifndef UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <zlib.h>
#include <stdbool.h>
#include <stdint.h>

//#define GENOME_LEN 3100000000    // length of human genome
//#define MER_LEN 36    // must be multiples of 4
#define GENOME_LEN 1000    // length of human genome
#define MER_LEN 36    // must be multiples of 4

typedef struct 
{
    uint8_t bnt[MER_LEN / 4];    // 2-bit format of mer. bnt(bit nucleotide)
    uint8_t chr;    // chromosome number. i.e. 1, 2, 3, ..., 22, 23(X), 24(Y), 25(M)
    uint8_t gc;    // GC count of this 36-mer
    uint32_t pos;    // chromosome position start at 0
} meridx_t;    // meridx(mer index)


uint8_t *mer2bytes(const char *mer);
char *revcmpl(const char *mer);
bool check_nt(char nt);
bool check_mer(const char *mer);
char *get_fn(char *path);


#define UTILS_H
#endif