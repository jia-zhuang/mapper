#include "align.h"
#include "utils.h"
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int meridx_cmp(const void *p, const void *q) {
    return memcmp(((meridx_t *)p)->bnt, ((meridx_t *)q)->bnt, MER_LEN / 4);
}

int mapper_align(int argc, char *argv[])
{
    if (argc < 3 || argc > 4) {
        printf("Usage: \n");
        printf("For SE: mapper align <ref.fa> <read.fq>\n");
        printf("For PE: mapper align <ref.fa> <read1.fq> <read2.fq>\n\n");
        return 1;
    }

    char *fn_ref = get_fn(argv[1]);
    char *fn_read = get_fn(argv[2]);

    // check if the index if in shared memory
    char *fn_sz = (char *)malloc(strlen(fn_ref) + 4);
    strcpy(fn_sz, fn_ref); strcat(fn_sz, ".sz");
    char *fn_idx = (char *)malloc(strlen(fn_ref) + 5);
    strcpy(fn_idx, fn_ref); strcat(fn_idx, ".idx");

    int fd_idx = shm_open(fn_idx, O_RDONLY, 00444);
    if (fd_idx == -1) {
        printf("[align] Error: index of %s is not in shared memory!\n", fn_ref);
        printf("Use command 'mapper shm -l %s' to load first\n", fn_ref);
        return 2;
    }

    // load index from shared memory
    int fd_sz = shm_open(fn_sz, O_RDONLY, 00444);
    uint32_t *psz = (uint32_t *)mmap(NULL, sizeof(uint32_t), PROT_READ, MAP_SHARED, fd_sz, 0);
    meridx_t *pidx = (meridx_t *)mmap(NULL, *psz * sizeof(meridx_t), PROT_READ, MAP_SHARED, fd_idx, 0);

    // create/open output file
    char *fn_out = (char *)malloc(strlen(fn_read) + 5);
    strcpy(fn_out, fn_read); strcat(fn_out, ".out");
    FILE *fout = fopen(fn_out, "w");

    // SE
    if (argc == 3) {
        // read in reads file 
        clock_t start = clock();   // record time
        gzFile fp_read = gzopen(argv[2], "r");
        kseq_t *record = kseq_init(fp_read);
        meridx_t target, *pt;
        uint8_t *bnt;
        fprintf(fout, "chr\tpos\tGC\n");
        while (kseq_read(record) >= 0) {
            if (record->seq.l == MER_LEN && check_mer(record->seq.s)) {
                // +
                bnt = mer2bytes(record->seq.s);
                memcpy(target.bnt, bnt, MER_LEN / 4);
                pt = (meridx_t *)bsearch(&target, pidx, *psz, sizeof(meridx_t), meridx_cmp);
                if (pt != NULL) {
                    fprintf(fout, "%u\t%u\t%u\n", pt->chr, pt->pos, pt->gc);
                }

                // - reverse complement
                bnt = mer2bytes(revcmpl(record->seq.s));
                memcpy(target.bnt, bnt, MER_LEN / 4);
                pt = (meridx_t *)bsearch(&target, pidx, *psz, sizeof(meridx_t), meridx_cmp);
                if (pt != NULL) {
                    fprintf(fout, "%u\t%u\t%u\n", pt->chr, pt->pos, pt->gc);
                }
            }
        }

        printf("Finish alignment in %fs\n", (double)(clock() - start) / CLOCKS_PER_SEC);
        kseq_destroy(record);
        gzclose(fp_read);
    }

    // PE
    else if (argc == 4) {
        // read in reads file
        clock_t start = clock();   // record time
        gzFile fp_r1 = gzopen(argv[2], "r");
        gzFile fp_r2 = gzopen(argv[3], "r");
        kseq_t *record1 = kseq_init(fp_r1); 
        kseq_t *record2 = kseq_init(fp_r2); 
        meridx_t t1, t2;
        meridx_t *p1, *p2;
        uint8_t *bnt;
        fprintf(fout, "chr1\tpos1\tGC1\tchr2\tpos2\tGC2\n");
        while (kseq_read(record1) >= 0 && kseq_read(record2) >= 0) {
            if (record1->seq.l == MER_LEN && record2->seq.l == MER_LEN && check_mer(record1->seq.s) && check_mer(record2->seq.s)) {
                // +
                bnt = mer2bytes(record1->seq.s);
                memcpy(t1.bnt, bnt, MER_LEN / 4);
                bnt = mer2bytes(revcmpl(record2->seq.s));
                memcpy(t2.bnt, bnt, MER_LEN / 4);
                p1 = (meridx_t *)bsearch(&t1, pidx, *psz, sizeof(meridx_t), meridx_cmp);
                p2 = (meridx_t *)bsearch(&t2, pidx, *psz, sizeof(meridx_t), meridx_cmp);
                if (p1 != NULL && p2 != NULL) {
                    fprintf(fout, "%u\t%u\t%u\t%u\t%u\t%u\n", p1->chr, p1->pos, p1->gc, p2->chr, p2->pos, p2->gc);
                }

                // -
                bnt = mer2bytes(revcmpl(record1->seq.s));
                memcpy(t1.bnt, bnt, MER_LEN / 4);
                bnt = mer2bytes(record2->seq.s);
                memcpy(t2.bnt, bnt, MER_LEN / 4);
                p1 = (meridx_t *)bsearch(&t1, pidx, *psz, sizeof(meridx_t), meridx_cmp);
                p2 = (meridx_t *)bsearch(&t2, pidx, *psz, sizeof(meridx_t), meridx_cmp);
                if (p1 != NULL && p2 != NULL) {
                    fprintf(fout, "%u\t%u\t%u\t%u\t%u\t%u\n", p1->chr, p1->pos, p1->gc, p2->chr, p2->pos, p2->gc);
                }
            }
        }

        printf("Finish alignment in %fs\n", (double)(clock() - start) / CLOCKS_PER_SEC);
        kseq_destroy(record1);
        kseq_destroy(record2);
        gzclose(fp_r1);
        gzclose(fp_r2);
    }


    // close file and free memory
    free(fn_out); fn_out = NULL;
    fclose(fout);
    munmap(fn_idx, *psz * sizeof(meridx_t));
    munmap(fn_sz, sizeof(uint32_t));
    close(fd_idx);
    close(fd_sz);
    free(fn_sz); fn_sz = NULL;
    free(fn_idx); fn_idx = NULL;

    return 0;
}