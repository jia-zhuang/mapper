#include "utils.h"
#include "index.h"
#include "align.h"
#include "shm.h"

int main(int argc, char *argv[])
{
    int ret = 0;
    if (argc < 2) {
        printf("Usage: mapper <command> [options]\n\n");
        printf("command: index    index sequences in the FASTA format\n");
        printf("         align    reads alignment\n");
        printf("         shm      manage indices in shared memory\n\n");
    } 
    else if (strcmp(argv[1], "index") == 0) ret = mapper_index(argc - 1, argv + 1);
    else if (strcmp(argv[1], "align") == 0) ret = mapper_align(argc - 1, argv + 1);
    else if (strcmp(argv[1], "shm") == 0) ret = mapper_shm(argc - 1, argv + 1);
    else {
        printf("[main] unrecognized command '%s'\n", argv[1]);
        ret = -1;
    }

    return ret;
}