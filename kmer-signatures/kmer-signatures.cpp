#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <chrono>
#include "uthash.hpp"

typedef unsigned char byte;

#define SIGNATURE_LEN 64

int DENSITY = 21;
int PARTITION_SIZE;

int inverse[256];
const char* alphabet = "CSTPAGNDEQHRKMILVFYW";


void seed_random(char* term, int length);
short random_num(short max);
int compare_files(const char* filename1, const char* filename2);

int doc_sig[SIGNATURE_LEN];

int WORDLEN;
FILE* sig_file;

typedef struct
{
    char term[100];
    short sig[SIGNATURE_LEN];
    UT_hash_handle hh;
} hash_term;

hash_term* vocab = NULL;


short* compute_new_term_sig(char* term, short* term_sig)
{
    seed_random(term, WORDLEN);

    int non_zero = SIGNATURE_LEN * DENSITY / 100;

    int positive = 0;
    while (positive < non_zero / 2)
    {
        short pos = random_num(SIGNATURE_LEN);
        if (term_sig[pos] == 0)
        {
            term_sig[pos] = 1;
            positive++;
        }
    }

    int negative = 0;
    while (negative < non_zero / 2)
    {
        short pos = random_num(SIGNATURE_LEN);
        if (term_sig[pos] == 0)
        {
            term_sig[pos] = -1;
            negative++;
        }
    }
    return term_sig;
}

short* find_sig(char* term)
{
    hash_term* entry;
    HASH_FIND(hh, vocab, term, WORDLEN, entry);
    if (entry == NULL)
    {
        entry = (hash_term*)malloc(sizeof(hash_term));
        strncpy_s(entry->term, sizeof(entry->term), term, WORDLEN);
        memset(entry->sig, 0, sizeof(entry->sig));
        compute_new_term_sig(term, entry->sig);
        HASH_ADD(hh, vocab, term, WORDLEN, entry);
    }

    return entry->sig;
}


void signature_add(char* term)
{
    short* term_sig = find_sig(term);
    for (int i = 0; i < SIGNATURE_LEN; i++)
        doc_sig[i] += term_sig[i];
}

int doc = 0;

void compute_signature(char* sequence, int length)
{
    memset(doc_sig, 0, sizeof(doc_sig));

    for (int i = 0; i < length - WORDLEN + 1; i++)
        signature_add(sequence + i);

    // save document number to sig file
    fwrite(&doc, sizeof(int), 1, sig_file);

    // flatten and output to sig file
    for (int i = 0; i < SIGNATURE_LEN; i += 8)
    {
        byte c = 0;
        for (int j = 0; j < 8; j++)
            c |= (doc_sig[i + j] > 0) << (7 - j);
        fwrite(&c, sizeof(byte), 1, sig_file);
    }
}

#define min(a,b) ((a) < (b) ? (a) : (b))

void partition(char* sequence, int length)
{
    int i = 0;
    do
    {
        compute_signature(sequence + i, min(PARTITION_SIZE, length - i));
        i += PARTITION_SIZE / 2;
    } while (i + PARTITION_SIZE / 2 < length);
    doc++;
}

int power(int n, int e)
{
    int p = 1;
    for (int j = 0; j < e; j++)
        p *= n;
    return p;
}

int main(int argc, char* argv[])
{
    const char* filename = "qut2.fasta";
    const char* test_file = "test_release_qut2.fasta.part16_sigs03_64";
    //const char* filename = "qut3.fasta";
    //const char* test_file = "test_qut3.fasta.part16_sigs03_64";

    WORDLEN = 3;
    PARTITION_SIZE = 16;
    int WORDS = power(20, WORDLEN);

    for (int i = 0; i < strlen(alphabet); i++)
        inverse[alphabet[i]] = i;

    auto start = std::chrono::high_resolution_clock::now();

    FILE* file;
    errno_t OK = fopen_s(&file, filename, "r");

    if (OK != 0)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        return 1;
    }

    char outfile[256];
    sprintf_s(outfile, 256, "%s.part%d_sigs%02d_%d", filename, PARTITION_SIZE, WORDLEN, SIGNATURE_LEN);
    fopen_s(&sig_file, outfile, "w");

    char buffer[10000];
    while (!feof(file))
    {
        fgets(buffer, 10000, file); // skip meta data line
        fgets(buffer, 10000, file);
        int n = (int)strlen(buffer) - 1;
        buffer[n] = 0;
        partition(buffer, n);
    }
    fclose(file);

    fclose(sig_file);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    printf("%s %f seconds\n", filename, duration.count());

    if (compare_files(outfile, test_file) != 0) {
        fprintf(stderr, "Error: output file does not match test file\n");
        return 1;
    }
    else {
        printf("Output file matches test file\n");
    }

    return 0;
}

int compare_files(const char *filename1, const char *filename2) {

    FILE* file1, *file2;

    fopen_s(&file1, filename1, "r");
    fopen_s(&file2, filename2, "r");

    unsigned long pos;
    int c1, c2;
    for (pos = 0;; pos++) {
        c1 = getc(file1);
        c2 = getc(file2);
        if (c1 != c2 || c1 == EOF)
            break;
    }
    if (c1 == c2) {
        return 0;  // files are identical
    }
    else
        if (c1 == EOF) {
            return 1;
        }
        else
            if (c2 == EOF) {
                return 2;
            }
            else {
                return 3;
            }
}
