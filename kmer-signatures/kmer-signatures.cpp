#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <unordered_map>
#include <ppl.h>
#include <fstream> 
#include <chrono>
#include <array>
#include "uthash.hpp"
#include "robin_hood.h"
#include "ISAAC-rand.hpp"

typedef unsigned char byte;

#define SIGNATURE_LEN 64

#define THREADS 8

int DENSITY = 21;
int PARTITION_SIZE;

int inverse[256];
const char* alphabet = "CSTPAGNDEQHRKMILVFYW";

int WORDLEN;
FILE* sig_file;

byte** sigs;

char buffer[THREADS][10000];
int lengths[THREADS];
int sizes[THREADS];

typedef std::array<short, SIGNATURE_LEN> computed_sig;

//typedef std::unordered_map<char*, computed_sig> computed_sig_map;
typedef robin_hood::unordered_map<char*, computed_sig> computed_sig_map;

thread_local computed_sig_map sig_map;

int compare_files(const char* filename1, const char* filename2);

int doc = 0;

computed_sig compute_new_term_sig(char* term, computed_sig &term_sig)
{
    randctx* R = seed_random(term, WORDLEN);

    int non_zero = SIGNATURE_LEN * DENSITY / 100;

    int positive = 0;
    while (positive < non_zero / 2)
    {
        short pos = random_num(SIGNATURE_LEN, R);
        if (term_sig[pos] == 0)
        {
            term_sig[pos] = 1;
            positive++;
        }
    }

    int negative = 0;
    while (negative < non_zero / 2)
    {
        short pos = random_num(SIGNATURE_LEN, R);
        if (term_sig[pos] == 0)
        {
            term_sig[pos] = -1;
            negative++;
        }
    }
    free(R);
    return term_sig;
}

computed_sig find_sig(char* term)
{
    auto item = sig_map.find(term);
    if (item != sig_map.end())
    {
        return item->second;
    }
 
    computed_sig new_sig = computed_sig();
    compute_new_term_sig(term, new_sig);
    auto new_item = sig_map.insert({ term, new_sig });

    return new_item.first->second;
}

void signature_add(char* term, int doc_sig[])
{
    computed_sig term_sig = find_sig(term);
    for (int i = 0; i < SIGNATURE_LEN; i++)
    {
        doc_sig[i] += term_sig[i];
    }
}

void compute_signature(char* sequence, int length, int x, int n)
{
    int doc_sig[SIGNATURE_LEN];
    memset(doc_sig, 0, sizeof(doc_sig));

    for (int i = 0; i < length - WORDLEN + 1; i++)
    {
        signature_add(sequence + i, doc_sig);
    }

    // flatten and output to sig file
    for (int i = 0; i < SIGNATURE_LEN; i += 8)
    {
        sigs[x][n] = 0;
        for (int j = 0; j < 8; j++) {
            sigs[x][n] |= (doc_sig[i + j] > 0) << (7 - j);
        }
        n++;
    }
}

#define min(a,b) ((a) < (b) ? (a) : (b))

void partition(int x)
{
    int size = ((lengths[x] - 1) / (PARTITION_SIZE / 2)) * SIGNATURE_LEN / 8;

    sigs[x] = (byte*)calloc(size, sizeof(byte));

    int si = 0;
    int i = 0;
    do
    {
        compute_signature(buffer[x] + i, min(PARTITION_SIZE, lengths[x] - i), x, si);
        i += PARTITION_SIZE / 2;
        si += (SIGNATURE_LEN/8);
    } while (i + PARTITION_SIZE / 2 < lengths[x]);

    sizes[x] = size;
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
    OK = fopen_s(&sig_file, outfile, "w");
    if (OK != 0)
    {
        fprintf(stderr, "Error: failed to open sig file\n");
        return 1;
    }

    sigs = (byte**)malloc(sizeof(byte*) * THREADS);
    
    //setvbuf(sig_file, NULL, _IOFBF, 16000);

    while (!feof(file))
    {
        int work = 0;
        while (work < THREADS && !feof(file)) {
            fgets(buffer[work], 10000, file); // skip meta data line
            fgets(buffer[work], 10000, file);
            int n = (int)strlen(buffer[work]) - 1;
            lengths[work] = n;
            buffer[work][n] = 0;
            work++;
        }

        concurrency::parallel_for(int(0), work, [&](int i) {
            partition(i);
        }, concurrency::static_partitioner());

        for (int i = 0; i < work; i++) {

            for (int j = 0; j < sizes[i];j+=8) {
                fwrite(&doc, sizeof(int), 1, sig_file);
                fwrite(&sigs[i][j], sizeof(byte), 8, sig_file);
            }
            doc++;
            free(sigs[i]);
        }
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
        //printf("files are identical and have %lu bytes\n", pos);
        return 0;  // files are identical
    }
    else
        if (c1 == EOF) {
            printf("file1 is included in file2, the first %lu bytes are identical\n", pos);
            return 1;
        }
        else
            if (c2 == EOF) {
                printf("file2 is included in file1, the first %lu bytes are identical\n", pos);
                return 2;
            }
            else {
                printf("file1 and file2 differ at position %lu: %u <>%u\n", pos, c1, c2);
                return 3;
            }
}
