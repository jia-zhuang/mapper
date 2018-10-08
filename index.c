#include "index.h"
#include "utils.h"
#include "ksort.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define meridx_lt(a, b) (memcmp((a).bnt, (b).bnt, MER_LEN / 4) < 0)
KSORT_INIT(MERIDX_T, meridx_t, meridx_lt)

static uint8_t get_chr_num(const char *chr) {
    /*
     * Conver chromosome name to number
     *     chr1 -> 1
     *     chr2 -> 2
     *     ...
     *     chrX -> 23
     *     chrY -> 24
     *     chrM -> 25 
     */
    if (strcmp(chr, "chrX") == 0) 
        return 23;
    else if (strcmp(chr, "chrY") == 0)
        return 24;
    else if (strcmp(chr, "chrM") == 0)
        return 25;
    else {
        return atoi(chr + 3);
    }
}

int mapper_index(int argc, char *argv[])
{
    if (argc != 2) {
        printf("Usage: mapper index <in.fasta>\n\n");
        return 1;
    }

    // get file name: /path/to/ref.fa -> ref.fa
    char *fn = get_fn(argv[1]);

    // create/open output index file. i.e. ref.fa.idx and ref.fa.sz
    char *fn_idx = (char *)malloc(strlen(fn) + 5);
    strcpy(fn_idx, fn); strcat(fn_idx, ".idx");
    FILE *fidx = fopen(fn_idx, "wb");

    char *fn_sz = (char *)malloc(strlen(fn) + 4);
    strcpy(fn_sz, fn); strcat(fn_sz, ".sz");
    FILE *fsz = fopen(fn_sz, "wb");

    // open reference genome file
    gzFile fa = gzopen(argv[1], "r");
    kseq_t *record = kseq_init(fa);

    // start indexing
    clock_t start = clock();    // record time

    int rcd_num = 0;
    char mer[MER_LEN];
    uint8_t *bnt;
    meridx_t *meridx = (meridx_t *)malloc(GENOME_LEN * sizeof(meridx_t));
    //meridx_t *meridx = (meridx_t *)calloc(GENOME_LEN, sizeof(meridx_t));
    meridx_t *tmp;
    uint32_t idx_len = 0;
    while (kseq_read(record) >= 0) {
        printf("Record number %d, chromosome name %s, chromosome number %u, chromosome length %zu\n", rcd_num++, record->name.s, get_chr_num(record->name.s), record->seq.l);

        for (uint32_t i = 0; i < record->seq.l - MER_LEN + 1; ++i) {
        /*
         *  Can be extended for Masked Fasta File in here
         */

            uint8_t gc = 0;
            uint8_t j;
            // printf("***debug*** MER_LEN: %d\n", MER_LEN);
            for (j = 0; j < MER_LEN; ++j) {
                mer[j] = record->seq.s[i + j];

                // printf("**debug*** check_nt(mer[%d]): %d\n", j, check_nt(mer[j]));
                if (!check_nt(mer[j])) break;

                if (mer[j] == 'G' || mer[j] == 'g' || mer[j] == 'C' || mer[j] == 'c') gc++;
            }

            // if mer pass check_nt, then j == MER_LEN
            // printf("***debug*** j = %u\n", j);
            if (j == MER_LEN) {
                bnt = mer2bytes(mer);
                tmp = meridx + idx_len;
                memcpy(tmp->bnt, bnt, MER_LEN / 4);
                tmp->chr = get_chr_num(record->name.s);
                tmp->gc = gc;
                tmp->pos = i + 1;
                idx_len++;
            }
        }
    }

    ks_introsort(MERIDX_T, idx_len, meridx);

    // remove duplicated and save index
    meridx_t *p = meridx;
    meridx_t *q = meridx + 1;
    uint32_t idx_len_rmdp = 0;
    while (p < meridx + idx_len) {
        if (memcmp(p->bnt, q->bnt, MER_LEN / 4) == 0) {
            while (q < meridx + idx_len) {
                q++;
                if (memcmp(p->bnt, q->bnt, MER_LEN / 4) != 0) {
                    p = q++;
                    break;
                }
            }
        }
        else {
            idx_len_rmdp++;
            fwrite(p, sizeof(meridx_t), 1, fidx);
            p = q++;
        }
    }

    printf("Length of index: %u\n", idx_len);
    printf("Length of unique index: %u\n", idx_len_rmdp);
    fwrite(&idx_len_rmdp, sizeof(uint32_t), 1, fsz);

    printf("Finish indexing in %fs\n", (double)(clock() - start) / CLOCKS_PER_SEC);

    // close file and free memory
    kseq_destroy(record);
    gzclose(fa);
    fclose(fidx);
    fclose(fsz);
    free(meridx); meridx = NULL;
    free(fn_idx); fn_idx = NULL;
    free(fn_sz); fn_sz = NULL;

    return 0;
}