#include "utils.h"

static int8_t nt2code(char nt) {
    /*
     * Conver nucleotide to number
     *     A/a -> 0
     *     C/c -> 1
     *     G/g -> 2
     *     T/t -> 3
     */
    switch (nt) {
        case 'A':
            return 0;
        case 'a':
            return 0;
        case 'C':
            return 1;
        case 'c':
            return 1;
        case 'G':
            return 2;
        case 'g':
            return 2;
        case 'T':
            return 3;
        case 't':
            return 3;
        default:
            return -1;    // this should not happen!
    }
}

static uint8_t four_nt2one_byte(const char *nt) {
    /*
     * every nucleotide has 4 state which can be represent by 2 bit
     * Convert 4 nucleotide sequence to 1 byte
     * length(nt) must be 4
     */
    return (
        nt2code(nt[0]) << 6 | \
        nt2code(nt[1]) << 4 | \
        nt2code(nt[2]) << 2 | \
        nt2code(nt[3]) 
    );
}

uint8_t *mer2bytes(const char *mer) {
    /*
     * Convert mer to bytes. e.g. 36-mer to 9 bytes
     * return a static array which will be changed by next calling, be carefull!
     */
    static uint8_t bytes[MER_LEN / 4];
    for (int i = 0; i < MER_LEN / 4; ++i) {
        bytes[i] = four_nt2one_byte(mer + 4*i);
    }
    return bytes;
}

static char cmpl(char nt) {
    /*
     * Nucleotide complement
     *     A <=> T
     *     C <=> G
     */
    switch (nt) {
        case 'A':
            return 'T';
        case 'a':
            return 't';
        case 'C':
            return 'G';
        case 'c':
            return 'g';
        case 'g':
            return 'c';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        case 't':
            return 'a';
        default:
            return -1;    // this should not happen!
    }
}

char *revcmpl(const char *mer) {
    /*
     * The reverse complement of 36-mer
     * return a static array which will be changed by next calling, be carefull!
     */
    static char rcmpl[MER_LEN];
    for (int i = 0; i < MER_LEN; ++i) {
        rcmpl[i] = cmpl(mer[MER_LEN - 1 - i]);
    }
    return rcmpl;
}

bool check_nt(char nt) {
    /*
     * Check whether nucleotide is valid. i.e. equals 'A/a' or 'C/c' or 'G/g' or 'T/t'
     */
    if (
        nt != 'A' && \
        nt != 'a' && \
        nt != 'C' && \
        nt != 'c' && \
        nt != 'G' && \
        nt != 'g' && \
        nt != 'T' && \
        nt != 't'
    ) 
        return false;
    else 
        return true;
}

bool check_mer(const char *mer) {
    /*
     * Check whether the mer is valid. i.e. contains only 'A/a', 'C/c', 'G/g' and 'T/t'
     */
    for (int i = 0; i < MER_LEN; ++i) {
        if (!check_nt(mer[i])) 
            return false;
    }
    return true;
}

char *get_fn(char *path) {
    /*
     *  Get file name
     *  /path/to/ref.fa  =>  ref.fa
     */
    char *fn = strrchr(path, '/');
    return (fn == NULL? path: fn + 1);
}

/*
 * Self Test
 */ 

/*
int main(int argc, char const *argv[])
{
    // test nt2code
    for (int i = 0; i < strlen(argv[1]); ++i) {
        printf("%d\n", nt2code(argv[1][i]));
    }


    // test four_nt2one_byte
    //printf("%02x\n", (uint8_t)four_nt2one_byte(argv[1]));


    // test mer2bytes
    int8_t *bytes = mer2bytes(argv[1]);
    for (int i = 0; i < MER_LEN / 4; ++i) {
        printf("%02x\n", (uint8_t)bytes[i]);
    }


    // test revcmpl
    printf("%s\n", revcmpl(argv[1]));

    // test check_mer
    printf("%d\n", check_mer(argv[1]));

    return 0;

}

*/