#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bwtexk.h"
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"
#include "bwa.h"
#include "bwase.h"
#include "kstring.h"
#include "klcp.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.1"
#endif

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

exk_opt_t *exk_init_opt()
{
	exk_opt_t *o;
	o = (exk_opt_t*)calloc(1, sizeof(exk_opt_t));
	o->mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
	o->n_threads = 1;
	o->trim_qual = 0;
	o->kmer_length = 14;
	o->use_klcp = 0;
	o->output_rids = 0;
	return o;
}

int bwt_cal_sa_coord(const bwt_t *bwt, int len, const ubyte_t *str, uint64_t* k, uint64_t* l, int start_pos)
{
	bwtint_t ok, ol;
	int i;
	*k = 0; *l = bwt->seq_len;

	//fprintf(stderr, "start k = %d, l = %d\n", *k, *l);
	for (i = start_pos; i < start_pos + len; ++i) {
		ubyte_t c = str[i];
		if (c < 4) {
			bwt_2occ(bwt, *k - 1, *l, c, &ok, &ol);
			*k = bwt->L2[c] + ok + 1;
			*l = bwt->L2[c] + ol;
			//fprintf(stderr, "after i = %d character, cur k = %d, l = %d\n", i, *k, *l);
		}
		if (*k > *l || c > 3) { // then restart
			return i - start_pos;
		}
	}
	return len;
}

int bwt_cal_sa_coord_continue(const bwt_t *bwt, int len, const ubyte_t *str,
															uint64_t* k, uint64_t* l, int start_pos, klcp_t* klcp)
{
	bwtint_t ok, ol;
	int i;
	//int start_k = *k;
	//int start_l = *l;
	//fprintf(stderr, "initial k = %d, l = %d\n", *k, *l);

	// if (*k != klcp->prev[start_k]) {
	// 	fprintf(stderr, "init k = %d, found k = %d, prev[k] = %d\n", start_k, *k, klcp->prev[start_k]);
	// 	fprintf(stderr, "klcp[start_k] = %d, klcp[*k] = %d\n", is_member(klcp->klcp, start_k), is_member(klcp->klcp, *k));
	// 	int a;
	// 	scanf("%d", &a);
	// }

	// if (*l != klcp->next[start_l]) {
	// 	fprintf(stderr, "init l = %d, found l = %d, next[l] = %d\n", start_l, *l, klcp->next[start_l]);
	// 	fprintf(stderr, "init k = %d\n", start_k);
	// 	fprintf(stderr, "klcp[start_l] = %d, klcp[*l] = %d\n", is_member(klcp->klcp, start_l), is_member(klcp->klcp, *l));
	// 	int a;
	// 	scanf("%d", &a);
	// }

	*k = decrease_k(klcp, *k);
	*l = increase_l(klcp, *l);

	//fprintf(stderr, "increased k = %d, l = %d\n", *k, *l);

	//fprintf(stderr, "start k = %d, l = %d\n", *k, *l);
	for (i = start_pos; i < start_pos + len; ++i) {
		ubyte_t c = str[i];
		if (c < 4) {
			bwt_2occ(bwt, *k - 1, *l, c, &ok, &ol);
			*k = bwt->L2[c] + ok + 1;
			*l = bwt->L2[c] + ol;
			//fprintf(stderr, "after i = %d character, cur k = %d, l = %d\n", i, *k, *l);
		}
		if (*k > *l || c > 3) { // then restart
			return i - start_pos;
		}
	}
	return len;
}

void output_chromosomes(const bwaidx_t* idx, const int seq_len, const uint64_t k,
												const uint64_t l, int8_t* seen_rids_marks) {
	int* seen_rids = malloc((l - k + 1) * sizeof(int));
	int rids_cnt = 0;
	int t;
	for(t = k; t <= l; ++t) {
		int strand;
		int pos = bwa_sa2pos(idx->bns, idx->bwt, t, seq_len, &strand);//bwt_sa(bwt, t);
		int rid = bns_pos2rid(idx->bns, pos);
		int seen = 0;
		if (rid != -1) {
			seen = seen_rids_marks[rid];
		}
		if (!seen && rid != -1) {
			seen_rids[rids_cnt] = rid;
			++rids_cnt;
			seen_rids_marks[rid] = 1;
			//fprintf(stderr, "t = %d, pos = %d, rid = %d\n", t, pos, rid);
		}
	}
	fprintf(stdout, "%d ", rids_cnt);
	int r;
	for(r = 0; r < rids_cnt; ++r) {
		fprintf(stdout, "%d ", seen_rids[r]);
		seen_rids_marks[seen_rids[r]] = 0;
	}
	fprintf(stdout, "\n");
}

void bwa_cal_sa(int tid, bwaidx_t* idx, int n_seqs, bwa_seq_t *seqs,
								const exk_opt_t *opt, klcp_t* klcp)
{
	bwase_initialize();
	int i, j;

	bwt_t* bwt = idx->bwt;

	int8_t* seen_rids_marks = malloc(idx->bns->n_seqs * sizeof(int8_t));
	for(i = 0; i < idx->bns->n_seqs; ++i) {
		seen_rids_marks[i] = 0;
	}
	fprintf(stdout, "\n");
	for (i = 0; i != n_seqs; ++i) {
		if (i % 1000 == 1) {
			fprintf(stderr, "processed %d reads in chunk\n", i);
		}
		bwa_seq_t *p = seqs + i;
		p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;

		// NEED TO UNDERSTAND
		// core function
		// for (j = 0; j < p->len; ++j) // we need to complement
		// 	p->seq[j] = p->seq[j] > 3? 4 : 3 - p->seq[j];

		if (opt->output_rids) {
			fprintf(stdout, "#");
			for(j = 0; j < p->len; ++j) {
				fprintf(stdout, "%c", "ACGTN"[p->seq[j]]);
			}
			fprintf(stdout, "\n");
		}
		uint64_t k, l;
		int start_pos = 0;
		int zero_streak = 0;
		int was_one = 0;
		while (start_pos <= p->len - opt->kmer_length) {
			if (start_pos == 0) {
				k = 0;
				l = 0;
				bwt_cal_sa_coord(bwt, opt->kmer_length, p->seq, &k, &l, start_pos);
			} else {
				if (opt->use_klcp && k <= l) {
					bwt_cal_sa_coord_continue(bwt, 1, p->seq, &k, &l, start_pos + opt->kmer_length - 1, klcp);
				} else {
					k = 0;
					l = 0;
					bwt_cal_sa_coord(bwt, opt->kmer_length, p->seq, &k, &l, start_pos);
				}
			}
			//fprintf(stderr, "start_pos = %d\n", start_pos);
			//fprintf(stderr, "found k = %llu, l = %llu\n", k, l);
			if (opt->output_rids) {
				output_chromosomes(idx, opt->kmer_length, k, l, seen_rids_marks);
			}
			if (opt->skip_after_fail) {
				if (k <= l) {
					was_one = 1;
					zero_streak = 0;
				} else {
					if (was_one) {
						if (zero_streak == 0) {
							zero_streak += opt->kmer_length - 2;
							if (opt->output_rids) {
								int ind;
								for(ind = 0; ind < opt->kmer_length - 2 && start_pos + ind < p->len - opt->kmer_length; ++ind) {
									fprintf(stdout, "0 \n");
								}
							}
							start_pos += opt->kmer_length - 2;
						} else {
							zero_streak++;
						}
					}
				}
			}
			start_pos++;
		}
		//fprintf(stdout, "#\n");


		free(p->name); free(p->seq); free(p->rseq); free(p->qual);
		p->name = 0; p->seq = p->rseq = p->qual = 0;
	}
}

bwa_seqio_t *bwa_open_reads_new(int mode, const char *fn_fa)
{
	bwa_seqio_t *ks;
	if (mode & BWA_MODE_BAM) { // open BAM
		int which = 0;
		if (mode & BWA_MODE_BAM_SE) which |= 4;
		if (mode & BWA_MODE_BAM_READ1) which |= 1;
		if (mode & BWA_MODE_BAM_READ2) which |= 2;
		if (which == 0) which = 7; // then read all reads
		ks = bwa_bam_open(fn_fa, which);
	} else ks = bwa_seq_open(fn_fa);
	return ks;
}

void bwa_exk_core(const char *prefix, const char *fn_fa, const exk_opt_t *opt) {
	int n_seqs, tot_seqs = 0;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	clock_t t;
	bwt_t *bwt;
	bwaidx_t* idx;

	{ // load BWT
		//fprintf(stderr, "%s\n", prefix);
		if ((idx = bwa_idx_load(prefix, BWA_IDX_ALL)) == 0) {
			fprintf(stderr, "Couldn't load idx from %s\n", prefix);
			return;
		}
		bwt = idx->bwt;
	}

	klcp_t* klcp = malloc(sizeof(klcp_t));
	klcp->klcp = malloc(sizeof(bitarray_t));
	if (opt->use_klcp) {
		char* fn = malloc((strlen(prefix) + 10) * sizeof(char));
	  strcpy(fn, prefix);
		strcat(fn, ".");
	  char* kmer_length_str = malloc(5 * sizeof(char));
	  sprintf(kmer_length_str, "%d", opt->kmer_length);
	  strcat(fn, kmer_length_str);
	  strcat(fn, ".bit.klcp");
		klcp_restore(fn, klcp);
	}

	ks = bwa_open_reads_new(opt->mode, fn_fa);
	float total_time = 0;
	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual)) != 0) {
		tot_seqs += n_seqs;
		t = clock();
		bwa_cal_sa(0, idx, n_seqs, seqs, opt, klcp);
		total_time += (float)(clock() - t) / CLOCKS_PER_SEC;
		bwa_free_read_seq(n_seqs, seqs);
	}
	fprintf(stderr, "match time: %.2f sec\n", total_time);
	//fprintf(stderr, "tot_seqs = %d\n", tot_seqs);
	//fprintf(stderr, "overall_increase = %llu\n", overall_increase);
	//fprintf(stderr, "increase per k-mer = %lf\n", 1.0 * overall_increase / (tot_seqs * (seq_len - opt->kmer_length + 1)));
	// destroy
	bwt_destroy(bwt);
	bwa_seq_close(ks);
}

int exk_match(int argc, char *argv[])
{
	int c, opte = -1;
	exk_opt_t *opt;
	char *prefix;

	opt = exk_init_opt();
	while ((c = getopt(argc, argv, "suvk:n:o:e:i:d:l:LR:t:NM:O:E:q:f:b012IYB:")) >= 0) {
		switch (c) {
		case 'v': opt->output_rids = 1; break;
		case 'u': opt->use_klcp = 1; break;
		case 'k': opt->kmer_length = atoi(optarg); break;
		case 's': opt->skip_after_fail = 1; break;
		case 'e': opte = atoi(optarg); break;
		case 't': opt->n_threads = atoi(optarg); break;
		case 'L': opt->mode |= BWA_MODE_LOGGAP; break;
		case 'q': opt->trim_qual = atoi(optarg); break;
		case 'N': opt->mode |= BWA_MODE_NONSTOP; break;
		case 'f': xreopen(optarg, "wb", stdout); break;
		case 'b': opt->mode |= BWA_MODE_BAM; break;
		case '0': opt->mode |= BWA_MODE_BAM_SE; break;
		case '1': opt->mode |= BWA_MODE_BAM_READ1; break;
		case '2': opt->mode |= BWA_MODE_BAM_READ2; break;
		case 'I': opt->mode |= BWA_MODE_IL13; break;
		case 'Y': opt->mode |= BWA_MODE_CFY; break;
		case 'B': opt->mode |= atoi(optarg) << 24; break;
		default: return 1;
		}
	}
	if (opte > 0) {
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   exk match [options] <prefix> <in.fq>\n\n");
		fprintf(stderr, "Options: -k INT    length of k-mer\n");
		fprintf(stderr, "         -u        use klcp for matching\n");
		fprintf(stderr, "         -v        output set of chromosomes for every k-mer\n");
		fprintf(stderr, "         -s        skip k-1 k-mers after failing matching k-mer\n");

		// fprintf(stderr, "         -t INT    number of threads [%d]\n", opt->n_threads);
		// fprintf(stderr, "         -B INT    length of barcode\n");
		// fprintf(stderr, "         -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
    // fprintf(stderr, "         -f FILE   file to write output to instead of stdout\n");
		// fprintf(stderr, "         -B INT    length of barcode\n");
		// fprintf(stderr, "         -I        the input is in the Illumina 1.3+ FASTQ-like format\n");
		// fprintf(stderr, "         -b        the input read file is in the BAM format\n");
		// fprintf(stderr, "         -0        use single-end reads only (effective with -b)\n");
		// fprintf(stderr, "         -1        use the 1st read in a pair (effective with -b)\n");
		// fprintf(stderr, "         -2        use the 2nd read in a pair (effective with -b)\n");
		// fprintf(stderr, "         -Y        filter Casava-filtered sequences\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if ((prefix = bwa_idx_infer_prefix(argv[optind])) == 0) {
		fprintf(stderr, "[%s] fail to locate the index %s\n", __func__, argv[optind]);
		free(opt);
		return 1;
	}
	bwa_exk_core(prefix, argv[optind+1], opt);
	free(opt); free(prefix);
	return 0;
}

int exk_index(int argc, char *argv[])
{
	int c, opte = -1;	exk_opt_t *opt;
	char *prefix;
	opt = exk_init_opt();
	while ((c = getopt(argc, argv, "k:n:o:e:i:d:l:k:LR:m:t:NM:O:E:q:f:b012IYB:")) >= 0) {
		switch (c) {
		case 'k': opt->kmer_length = atoi(optarg); break;
		case 'e': opte = atoi(optarg); break;
		case 't': opt->n_threads = atoi(optarg); break;
		case 'L': opt->mode |= BWA_MODE_LOGGAP; break;
		case 'q': opt->trim_qual = atoi(optarg); break;
		case 'N': opt->mode |= BWA_MODE_NONSTOP; break;
		case 'f': xreopen(optarg, "wb", stdout); break;
		case 'b': opt->mode |= BWA_MODE_BAM; break;
		case '0': opt->mode |= BWA_MODE_BAM_SE; break;
		case '1': opt->mode |= BWA_MODE_BAM_READ1; break;
		case '2': opt->mode |= BWA_MODE_BAM_READ2; break;
		case 'I': opt->mode |= BWA_MODE_IL13; break;
		case 'Y': opt->mode |= BWA_MODE_CFY; break;
		case 'B': opt->mode |= atoi(optarg) << 24; break;
		default: return 1;
		}
	}
	if (opte > 0) {
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   exk index <prefix>\n\n");
		// fprintf(stderr, "Options: -t INT    number of threads [%d]\n", opt->n_threads);
		// fprintf(stderr, "         -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
    // fprintf(stderr, "         -f FILE   file to write output to instead of stdout\n");

		fprintf(stderr, "Options:  -k INT    length of k-mer\n");

		// fprintf(stderr, "         -B INT    length of barcode\n");
		// fprintf(stderr, "         -I        the input is in the Illumina 1.3+ FASTQ-like format\n");
		// fprintf(stderr, "         -b        the input read file is in the BAM format\n");
		// fprintf(stderr, "         -0        use single-end reads only (effective with -b)\n");
		// fprintf(stderr, "         -1        use the 1st read in a pair (effective with -b)\n");
		// fprintf(stderr, "         -2        use the 2nd read in a pair (effective with -b)\n");
		// fprintf(stderr, "         -Y        filter Casava-filtered sequences\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if ((prefix = bwa_idx_infer_prefix(argv[optind])) == 0) {
		fprintf(stderr, "[%s] fail to locate the index %s\n", __func__, argv[optind]);
		return 1;
	}
	exk_index_core(prefix, argv[optind+1], opt);
	free(prefix);
	return 0;
}

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: exk (alignment of k-mers)\n");
	fprintf(stderr, "Usage:   exk command [options]\n\n");
	fprintf(stderr, "Command: index         construct klcp array\n");
	fprintf(stderr, "Command: match         align k-mers\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	extern char *bwa_pg;
	int i, ret = 0;
	kstring_t pg = {0,0,0};
	ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
	for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
	bwa_pg = pg.s;
	if (argc < 2) return usage();
	if (strcmp(argv[1], "index") == 0) ret = exk_index(argc - 1, argv + 1);
	else if (strcmp(argv[1], "match") == 0) ret = exk_match(argc-1, argv+1);

	err_fflush(stdout);
	err_fclose(stdout);

	free(bwa_pg);
	return ret;
}
