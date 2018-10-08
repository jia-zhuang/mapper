#define _XOPEN_SOURCE 500
#include "shm.h"
#include "utils.h"
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>

int mapper_shm(int argc, char *argv[])
{
    if (argc == 1 || (argc == 2 && strcmp(argv[1], "-p") != 0) || argc > 3) {
        printf("Usage: mapper shm [-l|-c|-p|-d] <ref.fa>\n");
        printf("Options: -l      load index of ref.fa into shared memory\n");
        printf("         -c      check whether index of ref.fa is in shared memory\n");
        printf("         -p      print files in /dev/shm/\n");
        printf("         -d      destroy index of ref.fa in shared memeory\n\n");
        return 1;
    }

    if (strcmp(argv[1], "-l") == 0) {
        // get file name with relative path
        char *path_fn_sz = (char *)malloc(strlen(argv[2]) + 4);
        strcpy(path_fn_sz, argv[2]); strcat(path_fn_sz, ".sz");
        char *path_fn_idx = (char *)malloc(strlen(argv[2]) + 5);
        strcpy(path_fn_idx, argv[2]); strcat(path_fn_idx, ".idx");

        // get file name
        char *fn_sz = get_fn(path_fn_sz);
        char *fn_idx = get_fn(path_fn_idx);

        // open sz file
        FILE *fsz = fopen(path_fn_sz, "rb");
        uint32_t idx_len;
        fread(&idx_len, sizeof(uint32_t), 1, fsz);

        // open shared memory
        int fd_idx = shm_open(fn_idx, O_RDWR | O_CREAT | O_EXCL, 00666);
        int fd_sz = shm_open(fn_sz, O_RDWR | O_CREAT | O_EXCL, 00666);
        if (fd_idx == -1 || fd_sz == -1) {
            printf("Error loading index to shm!\n");
            printf("Use command 'mapper shm -c %s' to check if the index is already in shared memory\n", argv[2]);
            return 2;
        }
        ftruncate(fd_idx, idx_len * sizeof(meridx_t));
        ftruncate(fd_sz, sizeof(uint32_t));
        void *pidx = mmap(0, idx_len*sizeof(meridx_t), PROT_READ | PROT_WRITE, MAP_SHARED, fd_idx, 0);
        void *psz = mmap(0, sizeof(uint32_t), PROT_READ | PROT_WRITE, MAP_SHARED, fd_sz, 0);        

        // open idx file
        FILE *fidx = fopen(path_fn_idx, "rb");
        printf("Write %zu items to shared memory\n", fread(pidx, sizeof(meridx_t), idx_len, fidx));
        memcpy(psz, &idx_len, sizeof(uint32_t));

        // close file and free memory
        munmap(fn_idx, idx_len * sizeof(meridx_t));
        munmap(fn_sz, sizeof(uint32_t));
        close(fd_idx);
        close(fd_sz);
        fclose(fidx);
        fclose(fsz);
        free(path_fn_sz); path_fn_sz = NULL;
        free(path_fn_idx); path_fn_idx = NULL;

    } else if (strcmp(argv[1], "-c") == 0) {
        char *fn = get_fn(argv[2]);
        char *fn_idx = (char *)malloc(strlen(fn) + 5);
        strcpy(fn_idx, fn); strcat(fn_idx, ".idx");

        int fidx = shm_open(fn_idx, O_RDWR | O_CREAT | O_EXCL, 00666);
        if (fidx == -1) {
            printf("Index of %s is already in shared memory\n", fn);
        } else {
            printf("Index of %s is *not* in shared memory\n", fn);
            shm_unlink(fn_idx);
            close(fidx);
        }

        free(fn_idx); fn_idx = NULL;
    } else if (strcmp(argv[1], "-p") == 0) {
        system("ls /dev/shm/");

    } else if (strcmp(argv[1], "-d") == 0) {
        char *fn = get_fn(argv[2]);
        // file name
        char *fn_sz = (char *)malloc(strlen(fn) + 4);
        strcpy(fn_sz, fn); strcat(fn_sz, ".sz");
        char *fn_idx = (char *)malloc(strlen(fn) + 5);
        strcpy(fn_idx, fn); strcat(fn_idx, ".idx");

        // check whether the index is in shared memory
        int fidx = shm_open(fn_idx, O_RDWR | O_CREAT | O_EXCL, 00666);
        if (fidx != -1) {
            printf("Index of %s is not in shared memory. Nothing to be done!\n", fn);
            shm_unlink(fn_idx);
            close(fidx);
            free(fn_sz); fn_sz = NULL;
            free(fn_idx); fn_idx = NULL;
            return 2;
        }

        // remove index
        shm_unlink(fn_idx);
        shm_unlink(fn_sz);
        free(fn_sz); fn_sz = NULL;
        free(fn_idx); fn_idx = NULL;
        printf("Successfully destroyed index of %s in shared memory\n", fn);

    } else {
        printf("[shm] wrong options: %s\n", argv[1]);
        return 1;
    }

    return 0;
}
