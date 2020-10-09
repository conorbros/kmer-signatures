#include <stdio.h>
#include <stdlib.h>
#include <ppl.h>
#include <chrono>
#include <array>
#include "robin_hood.h"
#include "ISAAC-rand.hpp"

#define min(a,b) ((a) < (b) ? (a) : (b))

constexpr int SIGNATURE_LEN = 64;

constexpr int CHUNKS = 2000;

constexpr int DENSITY = 21;

constexpr int PARTITION_SIZE = 16;

constexpr int WORDLEN = 3;

typedef unsigned char byte;
byte** sigs;

char buffer[CHUNKS][10000];
int lengths[CHUNKS];
int sizes[CHUNKS];

typedef std::array<short, SIGNATURE_LEN> computed_sig;

// Each thread has their own hashmap to store previously computed signatures
thread_local robin_hood::unordered_map<char*, computed_sig> sig_map;

int doc = 0;

int compare_files(const char* filename1, const char* filename2);

/**
 * Computes a new signature for the supplied term.
 * 
 * \param term term to compute the new signature for
 * \param term_sig array to populate with computed signature
 */
void compute_new_term_sig(char* term, computed_sig& term_sig)
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
}

/**
 * Finds the signature for the supplied term by retrieing from the hashmap or computing it.
 * 
 * \param term to get signature for
 * \return the computed sig
 */
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

/**
 * Adds signatures to the doc sig array.
 * 
 * \param term to add signature for
 * \param doc_sig array to write the signatures to
 */
void signature_add(char* term, int doc_sig[])
{
	computed_sig term_sig = find_sig(term);
	for (int i = 0; i < SIGNATURE_LEN; i++)
	{
		doc_sig[i] += term_sig[i];
	}
}

/**
 * Computes the signature for a sequence partition and adds the results to the signature array.
 * 
 * \param sequence partition to compute signature for
 * \param length length of the partition
 * \param x index of the thread in global variables
 * \param n index to write from in the sig array
 */
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

/**
 * .
 * 
 * \param sequence sequence to partition and compute signatures for
 * \param length length of the sequence
 * \param x index of the thread in global variables
 * \return 
 */
int partition(char * sequence, int length, int x)
{
	int size = ((length - 1) / (PARTITION_SIZE / 2)) * SIGNATURE_LEN / 8;

	sigs[x] = (byte*)calloc(size, sizeof(byte));

	int si = 0;
	int i = 0;
	do
	{
		compute_signature(sequence + i, min(PARTITION_SIZE, length - i), x, si);
		i += PARTITION_SIZE / 2;
		si += (SIGNATURE_LEN / 8);
	} while (i + PARTITION_SIZE / 2 < length);

	return size;
}

int power(int n, int e)
{
	int p = 1;
	for (int j = 0; j < e; j++)
		p *= n;
	return p;
}

/**
 * Creates a sig file with the supplied .fasta file.
 *
 * \param threads max number of threads to allocate
 * \return duration in seconds that it took to complete the processing
 */
double kmer_signatures(int threads) {
	concurrency::Scheduler* qsScheduler = concurrency::Scheduler::Create(concurrency::SchedulerPolicy(2, concurrency::MinConcurrency, 1, concurrency::MaxConcurrency, threads));
	qsScheduler->Attach();

	const char* filename = "qut2.fasta";
	const char* test_file = "test_release_qut2.fasta.part16_sigs03_64";
	//const char* filename = "qut3.fasta";
	//const char* test_file = "test_qut3.fasta.part16_sigs03_64";

	auto start = std::chrono::high_resolution_clock::now();

	FILE* file;
	errno_t OK = fopen_s(&file, filename, "r");
	if (OK != 0)
	{
		fprintf(stderr, "Error: failed to open file %s\n", filename);
		exit(EXIT_FAILURE);
	}

	FILE* sig_file;
	char outfile[256];
	sprintf_s(outfile, 256, "%s.part%d_sigs%02d_%d", filename, PARTITION_SIZE, WORDLEN, SIGNATURE_LEN);
	OK = fopen_s(&sig_file, outfile, "w");
	if (OK != 0)
	{
		fprintf(stderr, "Error: failed to open sig file\n");
		exit(EXIT_FAILURE);
	}

	sigs = (byte**)malloc(sizeof(byte*) * CHUNKS);

	while (!feof(file))
	{
		int work = 0;
		while (work < CHUNKS && !feof(file)) {
			fgets(buffer[work], 10000, file); // skip meta data line
			fgets(buffer[work], 10000, file);
			int n = (int)strlen(buffer[work]) - 1;
			lengths[work] = n;
			buffer[work][n] = 0;
			work++;
		}

		concurrency::parallel_for(int(0), work, [&](int i) {
			sizes[i] = partition(buffer[i], lengths[i], i);
		}, concurrency::static_partitioner());

		for (int i = 0; i < work; i++) {
			for (int j = 0; j < sizes[i];j += 8) {
				fwrite(&doc, sizeof(int), 1, sig_file);
				fwrite(&sigs[i][j], sizeof(byte), 8, sig_file);
			}
			doc++;
			free(sigs[i]);
		}
	}

	fclose(file);
	fclose(sig_file);

	concurrency::CurrentScheduler::Detach();

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;

	if (compare_files(outfile, test_file) != 0) {
		fprintf(stderr, "Error: output file does not match test file\n");
		return 1;
	}
	else {
		printf("Output file matches test file\n");
	}
	return duration.count();
}

int main(int argc, char* argv[])
{
	printf("%f seconds \n", kmer_signatures(8));

	return 0;
}

/**
 * Opens and compares two binary files byte by byte.
 *
 * \param filename1 first file to compare
 * \param filename2 second file to compare
 * \return 0 if identical, >0 if not identical
 */
int compare_files(const char* filename1, const char* filename2) {

	FILE* file1, * file2;

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
		return 0;
	}
	else {
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
}
