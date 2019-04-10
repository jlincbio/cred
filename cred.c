// chemred.c, version Apr 4 2019
#include <stdarg.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <getopt.h>
#include <zlib.h>
#include <math.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/khash.h"

// DATA STRUCTURE CONFIGURATIONS
KHASH_MAP_INIT_INT(genome, int);
typedef struct {     // auxiliary data structure
    samFile *fp;     // the file handle
    bam_hdr_t *hdr;  // the file header
    hts_itr_t *iter; // NULL if a region not specified
    int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;
struct entry_interval { 
    int start, end;
	int tid;
	double score;
	double p;
};
typedef struct entry_interval regions;
struct intr { 
    int x, i;
};
typedef struct intr int_r;

// CONSTANTS
#define EPS2 0.00000001
#define EPS3 0.000000000001 // convergence cutoff for Welch's t-test
#define LOG_EPS2 -18.42068
#define SCORE_MAX 485165195 // maximum score for enrichment
#define SCORE_LOG 20.00000
#define SCORE_MIN -20.00000 // minimum log score
#define INV_SQRT_2PI 0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#define KS_P_INIT -1.233700550136169749038 // -pi/2 * pi/4
// #define FEATURE_LIMIT 134217713 // maximum array size
#define FEATURE_LIMIT 1048576 // maximum array size
#define CUTOFF_RFRAC 0.20 // screening cutoff for fraction of nonzero reads
#define CUTOFF_MAPQ 30

// MATH MACROS
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#define ABSDIFF(X,Y) (fabs(X - Y)) // t-test
#define BETWEEN(Z, X, Y) (Z < Y && Z > X) // testing for overlap

// COMPARISON FUNCTIONS
int cmpfun_ranks(const void *a, const void *b){
	// array with ranks; used with comparator function below
	// sort index, but compare with original values
    int_r *p1 = (int_r*) a;
    int_r *p2 = (int_r*) b;
    return (p1->x < p2->x) ? -1 : (p1->x > p2->x);
}
int cmpfun_region(const void *v1, const void *v2) {
	// sorts regions
	regions *p1 = (regions*) v1;
	regions *p2 = (regions*) v2;
	return (p1 -> start) - (p2 -> start);
}

// PARAMETER DECLARATIONS
static char *fg_bam, *bg_bam, *f_output, *ref_faidx;
static int nthread, wsize_min, wsize_max;
static int wsize_step = 10;
static int wsize_consect = 10; // number of consecutive windows
static int max_depth = 8000; // coverage limit
static int baseQ = 0; // base quality threshold
static int mapQ = CUTOFF_MAPQ;
static int min_len = 0; // minimum query length; ignore reads shorter
static aux_t **fg_data, **bg_data;
static FILE *out_fh;
static double p_cutoff = 0.0001; // default p-value
static int status = EXIT_SUCCESS; // exit status
const char read_strand[2] = {'+', '-'};

// FUNCTIONAL DECLARATIONS
double array_maxabs(const double *x, const int n);
double array_mean(double x[], int nx);
double array_var(double x[], int nx, double x_mean);
double log_limit(double p);
double log_max(double x);
double parse_double(const char *str);
double p_ks2(double D, double n);
double p_student (double x[], int nx, double y[], int ny);
double region_tests(int tid, int head, int tail, khash_t(genome) *fg, khash_t(genome) *bg, double *p, const double r_xy, const int stage);
double stat_ks2(int *x1, int n1, int *x2, int n2);
double test_ks2(int x[], int y[], int n, double *p);
double test_student(int x[], int y[], int n, double *p);
int parse_int(const char *str);
int process_bam(khash_t(genome) *f_coverage, double *f_mdp, void *data, char **chrom_names, int *len_chrom, long int *sum_chrom, int n_chrom, int *clr_chrom);
int region_print(regions *features, char **names, int counts, double p_cutoff);
static int initialize_bam(aux_t **data, char *f_bam, const int n);
static int read_bam(void *data, bam1_t *b);
static void usage(const char *str);
void print_time(char *s);
void region_merge(regions *sites, int size, int counts, int *region_counts);
int set_window(regions *features, khash_t(genome) *depths, int tid, double chr_average, int chr_size, int wsize, int *feature_count, int *f_size);
	
// BASIC OPERATIONS
static void usage(const char *str) {
	fprintf(stderr, "CRED: Chem-seq Read Enrichment Discovery (Version 0.1, Apr 2019 Initial Release)\n");
    fprintf(stderr, "Command  :  cred [options] -t TREATMENT.BAM -c CONTROL.BAM > OUTPUT.BED\n");
    fprintf(stderr, "Required :\n");
	fprintf(stderr, "  -t  TREATMENT.BAM  Path to the treatment (\"pulldown\") track [BAM]\n");
	fprintf(stderr, "  -c  CONTROL.BAM    Path to the control (\"input\") Chem-seq track [BAM]\n");
	fprintf(stderr, "Optional :\n");
	fprintf(stderr, "  -p  [P-VALUE]      Significance level [default 0.0001]\n");
	fprintf(stderr, "  -q  [SCORE]        Minimum MAPQ quality for reads to count [default 30]\n");
	fprintf(stderr, "  -w  [INTEGER]      Size of differential windows [default 1200 bp]\n");
	fprintf(stderr, "  -k                 Evaluate site significance with Kolmogorov-Smirnov\n");
	fprintf(stderr, "Reminders:\n");
	fprintf(stderr, "   1. BAM files must be sorted and indexed.\n");
	fprintf(stderr, "   2. Use a pipe (\">\") to capture CRED output.\n\n");
	
	if (str != NULL) {
		// error message specified: exit with usage
		fprintf(stderr, "ERROR: %s\n", str);
		exit(EXIT_FAILURE);
	}
}
int parse_int(const char *str) {
    errno = 0; 
    char *end;
	int num = strtol(str, &end, 0);
    if (end != str && *end != '\0') {
		fprintf(stderr, "# [CRED] ERROR: can't convert '%s' to integer (got '%s')\n", str, end);
	}
	return(num);
}
double parse_double(const char *str) {
    errno = 0;
    char *end; double num;
	if (strchr(str, '/')) {
		// possible fraction
		int n1, n2;
		sscanf(str, "%d/%d", &n1, &n2);
		num = (1.0 * n1)/(1.0 * n2);
	} else {
		num = strtod(str, &end);
	    if (end != str && *end != '\0') {
			fprintf(stderr, "# [CRED] ERROR: can't convert '%s' to decimal (got '%s')\n", str, end);
		}
	}
    return(num);
}
void print_time(char *s) {
	time_t t; struct tm *tf;
	time(&t); tf = localtime(&t);
	fprintf(stderr, "# [CRED %04d/%02d/%02d %02d:%02d:%02d] %s\n",
		tf->tm_year + 1900,
		tf->tm_mon + 1,
		tf->tm_mday,
		tf->tm_hour,
		tf->tm_min,
		tf->tm_sec, s);
}

// BAM OPERATIONS
static int read_bam(void *data, bam1_t *b) {
    aux_t *aux = (aux_t*) data; // data is a pointer to an auxiliary structure
    int ret;
    while (1) {
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if ( ret<0 ) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        if ( (int)b->core.qual < aux->min_mapQ ) continue;
        if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
        break;
    }
    return ret;
}
static int initialize_bam(aux_t **data, char *f_bam, const int n) {
	// the core multi-pileup loop
	int i;
	for (i = 0; i < n; ++i) {
        int rf;
        data[i] = calloc(1, sizeof(aux_t));
		data[i]->fp = sam_open(f_bam, "r");
		
        if (data[i]->fp == NULL) {
            fprintf(stderr, "# [CRED] ERROR: cannot open %s\n", f_bam);
            status = EXIT_FAILURE;
            goto depth_end;
        }
        rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ;
        if (baseQ) rf |= SAM_QUAL;
        data[i]->min_mapQ = mapQ;                    // set the mapQ filter
        data[i]->min_len  = min_len;                 // set the qlen filter
        data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
        if (data[i]->hdr == NULL) {
            fprintf(stderr, "# [CRED] ERROR: cannot retrieve header from %s\n", f_bam);
            status = EXIT_FAILURE;
            goto depth_end;
        } else {
        	status = EXIT_SUCCESS;
        }
    }
	return(status);
	
	depth_end:
    if (fclose(out_fh) != 0) {
        if (status == EXIT_SUCCESS) {
			fprintf(stderr, "# [CRED] ERROR when trying to close %s\n",
				(f_output && strcmp(f_output, "-") != 0 ? f_output : "stdout"));
            status = EXIT_FAILURE;
        }
    }
    for (i = 0; i < n && data[i]; ++i) {
        bam_hdr_destroy(data[i]->hdr);
        if (data[i]->fp) sam_close(data[i]->fp);
        hts_itr_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data);
    return(status);
}
int process_bam(khash_t(genome) *f_coverage, double *f_mdp, void *data, char **chrom_names, int *len_chrom, long int *sum_chrom, int n_chrom, int *clr_chrom) {
	int i, w, ret, tid, pos; // chromosome id
    int n = 1, ml = 0; int beg = 0;
	int end = INT_MAX;
	
	// initialize and set maximum coverage depth
    bam_mplp_t mplp; mplp = bam_mplp_init(n, read_bam, (void **) data);
	if (0 < max_depth)
        bam_mplp_set_maxcnt(mplp,max_depth);
    else if (!max_depth)
        bam_mplp_set_maxcnt(mplp,INT_MAX);
	
	// n_plp[i] is the number of covering reads from the i-th BAM
    // plp[i] points to the array of covering reads (internal in mplp)
    int *n_plp; const bam_pileup1_t **plp;
	n_plp = calloc(n, sizeof(int));
	plp = calloc(n, sizeof(bam_pileup1_t*));
	
	int cur_tid = -1, n_nonzero = 1, n_cur_nonzero = 1;
	int len_genome = 0;
	fprintf(stderr, "# [CRED] Progress:");
	int max_nonzero = -1, max_cur_nonzero = -1; // switches for overflow warning
	while ((ret = bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) {
		// come to the next covered position
		if (cur_tid != tid) {
			// chromosome switched
			fprintf(stderr, "\n# [CRED %9d bp] %s...", len_chrom[tid], chrom_names[tid]);
			if (clr_chrom[tid] < 1) {
				// only increment when switching for the first time
				len_genome = len_genome + len_chrom[tid];
			} else {
				fprintf(stderr, "\n# [CRED] ERROR: is BAM sorted?");
				
			}
			if (cur_tid > -1) {
				f_mdp[cur_tid] = ((double) sum_chrom[cur_tid]) / ((double) n_cur_nonzero);
			}
			// fprintf(stderr, "%d ", tid);
			cur_tid = tid;
			n_cur_nonzero = 1;
			++clr_chrom[tid];
			// when tid changes
		}
		// tid is the same; same chromosome -> fine if BAM is sorted
	
		if (pos < beg || pos >= end) continue; // out of range; skip
        if (tid >= n_chrom) continue;     // diff # of @SQ lines per file?
    
        for (i = 0; i < n; ++i) {
			// filter out base, calculate depth and estimate mean
            int j, m = 0;
            for (j = 0; j < n_plp[i]; ++j) {
                const bam_pileup1_t *p = plp[i] + j;
				// having dels or refskips at tid:pos
                if (p->is_del || p->is_refskip) ++m;
                else if (p->qpos < p->b->core.l_qseq && 
					bam_get_qual(p->b)[p->qpos] < baseQ) ++m;
					// low base quality
            }
		
			int cur_coverage = n_plp[i] - m;
			if (cur_coverage) {
				int ret_cur_pos; // need to initialize hash of hash
				khiter_t stat_cur_pos = kh_put(genome, &f_coverage[tid], pos, &ret_cur_pos);
				// fg genome; use chromosomal pos as key
			
				if (ret_cur_pos > 0) {
					// printf("saving into hash... "); // empty, save coverage
					kh_value(&f_coverage[tid], stat_cur_pos) = cur_coverage;
				} else {
					// key present, increment
					// printf("appending at position... ");
					int cur_val = kh_value(&f_coverage[tid], stat_cur_pos) + cur_coverage;
					kh_value(&f_coverage[tid], stat_cur_pos) = cur_val;
					// printf("In hash = %d\n", kh_value(fg_coverage, pos));
				}
				sum_chrom[tid] = sum_chrom[tid] + cur_coverage;
			}
		}
		// right here
		if (n_nonzero < 0) {
			// overflow
			if (max_nonzero < 0) {
				max_nonzero = 1;
				n_nonzero = len_genome;
			}
		} else {
			++n_nonzero;
		}
		if ((n_cur_nonzero + 1) > len_chrom[tid]) {
			if (max_cur_nonzero < 0) {
				max_cur_nonzero = 1;
				n_cur_nonzero = len_chrom[tid];
				// fprintf(stderr, "\n# [CRED] WARNING: %s counter overflow...", chrom_names[tid]);
			}
		} else {
			++n_cur_nonzero;
		}
	}
	// last contig finished, port again
	f_mdp[cur_tid] = ((double) sum_chrom[cur_tid]) / ((double) n_cur_nonzero);
	// fprintf(stderr, "CUR = %d NEXT = %d SUM = %lu NONZERO = %d (R = %f)\n", 
	//	cur_tid, tid, sum_chrom[cur_tid], n_cur_nonzero, f_mdp[cur_tid]);
	
    if (ret < 0) status = EXIT_FAILURE;
    free(n_plp); free(plp); free(sum_chrom);
    bam_mplp_destroy(mplp);
	
	n_nonzero = abs(n_nonzero);
	if (n_nonzero > len_genome) {
		// fprintf(stderr, "\n# [CRED] WARNING: candidate upper limit exceeded (%u)", n_nonzero);
		n_nonzero = len_genome;
	}
	return(n_nonzero);
}

// SLIDING WINDOW OPERATIONS 
int set_window(regions *features, khash_t(genome) *depths, int tid, double chr_average, int chr_size, int wsize, int *feature_count, int *f_size) {
	int f_final = *f_size; int f_incr = 1;
	regions *f_buffer = NULL;
	
	int i, j, reg_head = 0, reg_tail = 0;
	int r_avg = (int) (chr_average * wsize);
	int gapl = wsize/5;
	if (chr_average > 0.01) {
		// minimum 1% of the contig should be covered before proceeding
		int r_retain = 0; int w_sum_prev = 0;
		int f_count = *feature_count;
		// fprintf(stderr, "Start counting from %d (or %d) on chromosome %d...\n", *feature_count, f_count, tid);
	
		khiter_t chr_head = kh_begin(depths);
		khiter_t chr_tail = kh_end(depths);
	
		for (i = chr_head; i < chr_tail - wsize; i = i + wsize_step) {
			//fprintf(stderr, "# [CRED] Searching in positions %u...", i);
			int w_sum = 0; int n_zeros = 0;
			int w_tail = i + (wsize - 1);
			if (w_tail >= (chr_size - 1)) {
				w_tail = chr_size - 1;
			}
		
			//fprintf(stderr, "# [CRED] Iterators: %u to %u\n", j, k1);
			for (j = i; j <= w_tail; ++j) {			
				int cur_depth = -1;
				if (kh_exist(depths, j)) {
					cur_depth = kh_value(depths, j);
					// fprintf(stderr, "# [CRED] Running sum = %u\n", w_sum);
				}
				// fprintf(stderr, "# [CRED] [%u - %u] %u pos. coverage: %u\n", i, w_tail, j, cur_depth);
			
				if (cur_depth > ((int) chr_average)) {
					w_sum = w_sum + cur_depth;
				} else {
					++n_zeros;
					if (n_zeros > gapl) break;
				}
				// if (depths[j] > 0) {
					// sum of depth within window
					// w_sum = w_sum + depths[j];
				//}
			}
			if ((w_sum > r_avg) && (((double) n_zeros/ wsize) < CUTOFF_RFRAC)) {
				++r_retain;
			} else {
				// overlapping region is broken
				r_retain = 0;
			}
		
			if (r_retain >= wsize_consect) {
				// >num of preset consecutive windows present; start tracking 
				reg_tail = j;
				if (reg_head < (i - r_retain * wsize_step)) {
					// head of the region is too far away from # of overlapping windows
					// printf("head of the region is too far away from # of overlapping windows\n");
					reg_head = i; 
				}
			
				if (f_count > 0) {
					// check if two regions are too close together
					int f0 = f_count - 1;
					if (BETWEEN(reg_head, features[f0].start, features[f0].end)) {	
						if (reg_tail > features[f0].end) {
							// overlap w/ tail, lengthen if w_sum > w_sum_prev
							if (w_sum > features[f0].score) {
								features[f0].end = reg_tail;
								features[f0].score = w_sum;
							} else if (BETWEEN(w_sum, (features[f0].score)/4, features[f0].score)) {
								features[f0].end = reg_tail;
								features[f0].score = MAX(w_sum, features[f0].score);
							} else {
								// falls outside tolerance range, reset
								r_retain = 0;
							}	
						}
						continue;
					}
				}
				
				if (f_count >= f_final) {
					// realloc required
					int f_cur = f_final;
					f_final = (int) f_final * (0.5 + f_incr);
					// fprintf(stderr, "# [CRED] Reallocating site counts to %d (from %d)...\n", f_final, f_count);
					f_buffer = realloc(features, f_final * sizeof(regions));
					if (f_buffer == NULL) {
						fprintf(stderr, "# [CRED] ERROR: maximum candidate size exceeded...\n");
						return -1;
					}
					features = f_buffer;
					++f_incr;
				}
				w_sum_prev = w_sum;
				features[f_count].start = reg_head;
				features[f_count].end = reg_tail;
				features[f_count].tid = tid;
				features[f_count].score = (double) w_sum_prev;
				++f_count;
			}
			// printf("Depth sum from %d to %d = %d\n", i, w_tail, w_sum);
		}
		fprintf(stderr, ": %d candidates.\n", f_count);
		*feature_count = f_count;
	}
	*f_size = f_final;
	return 1;
}
int region_print(regions *features, char **names, int counts, double p_cutoff) {
	int k, n = 0;
	for (k = 0; k < counts; ++k) {
		if (features[k].p <= p_cutoff) {
			++n;
			printf("%s\t%d\t%d\tpeak_%s_%d\t%.4f\t.\t%.4f\n",
				names[features[k].tid],
				features[k].start,
				features[k].end,
				names[features[k].tid],
				n, log_max(features[k].score),
				log_limit(features[k].p));
		}
	}
	return(n);
}
void region_merge(regions *sites, int size, int counts, int *region_counts) {
	fprintf(stderr, "# [CRED] Filtering and merging %d candidates...\n", counts);
	if (size == 0) {
		*region_counts = 0;
	} else {
		int curr = 0, i = 1;
		for(i = 1; i < counts; ++i) {
			// fprintf(stderr, "TID = %d (%d to %d) with score %d\n", sites[i].tid, sites[i].start, sites[i].end, sites[i].score);
			if ((sites[curr].end >= sites[i].start) && (sites[curr].tid == sites[i].tid)) {
				// fprintf(stderr, "Merging feature %d (start: %d -> end: %d) in position %d\n", i, sites[i].start, sites[i].end, curr);
				sites[curr].end = MAX(sites[curr].end, sites[i].end);
				sites[curr].tid = sites[i].tid;
				sites[curr].score = MAX(sites[curr].score, sites[i].score);
			} else {
				++curr;
				// fprintf(stderr, "Skipping to next feature (index %d from %d)...", curr, i);
				sites[curr] = sites[i];
			}
		}
		// fprintf(stderr, "Current region counts: %d\n", curr);
		*region_counts = curr;
	}
}

// STATISTICAL OPERATIONS
double log_max(double x) {
	if (x < EPS2) {
		return(SCORE_MIN);
	}
	if (x > SCORE_MAX) {
		return(SCORE_LOG);
	}
	return(log(x));
}
double log_limit(double p) {
	if (p < EPS2) {
		return(LOG_EPS2);
	}
	return(log(p));
}
double array_mean(double x[], int nx) {
	if (nx > 0) {
		int i; double sum = 0.0;
		for (i = 0; i < nx; ++i) {
			sum = sum + x[i];
		}
		return (sum/(double) nx);
	}
	return 0.0; // better return value is ?
}
double array_var(double x[], int nx, double x_mean) {
	if (nx > 1) {
		int i; double x_var = 0.0;
		for (i = 0; i < nx; ++i) {
			//1st part of added unbiased_sample_variance
			x_var = x_var + (x[i] - x_mean) * (x[i] - x_mean);
		}
		x_var = x_var/(nx - 1);
		return (x_var);
	}
	return 0.0; // needs a better return value
}
double array_maxabs(const double *x, const int n) {
	double max = 0.0;
	int i;
	for (i = 0; i < n; ++i) {
		double z = fabs(x[i]);
		if (z > max) {
			max = z;
		}
	}
	return max;
}

// WELCH'S T-TEST
double p_student (double x[], int nx, double y[], int ny) {
	if ((nx > 1) && (ny > 1)) {
		// evaluate a p-value only if both arrays have n > 1
		double x_mean = array_mean(x, nx);
		double y_mean = array_mean(y, ny);
		if (ABSDIFF(x_mean, y_mean) <= EPS2) {
			return 1.0; // no difference in mean, p = 1
		}
	
		double x_var = array_var(x, nx, x_mean);
		double y_var = array_var(y, ny, y_mean);
	
		double t = (x_mean - y_mean)/sqrt(x_var/nx + y_var/ny);
		double df = pow((x_var/nx + y_var/ny), 2.0)/((x_var*x_var)/(nx*nx*(nx-1))+(y_var*y_var)/(ny*ny*(ny-1)));
		double a = df/2;
		double pval = df/(t*t+df);
		
		if (t < -3.0) {
			// mean of x is not greater than y, just reject
			return 1.0;
		}
		
		if ((isinf(pval) != 0) || (isnan(pval) != 0)) {
			return 1.0;
		}
		
		// BETAIN computes the incomplete Beta function ratio.
		double beta = lgammal(a) + 0.57236494292470009 - lgammal(a + 0.5);
		int indx;
		double cx, pp, psq, qq, rx, xx;
		double ai = 1.0, term = 1.0; // starting point
		
		if (pval <= 0.0) {
			return 0.0;
		} else if (pval >= 1.0) {
			return 1.0;
		}
		
		psq = a + 0.5; cx = 1.0 - pval;
		if (a < psq * pval) {
			xx = cx; cx = pval; pp = 0.5;
			qq = a; indx = 1;
		} else {
			xx = pval; pp = a;
			qq = 0.5; indx = 0;
		}
		
		pval = 1.0;
		int ns = (int) (qq + cx * psq);
		
		// Soper reduction
		rx = xx/cx;
		double temp = qq - ai;
		if (ns == 0) {
			rx = xx;
		}
		
		int i = 0; // iteration count, at most 200?
		while (i < 200) {
			term = term * temp * rx/(pp + ai);
			pval = pval + term;;
			temp = fabs(term);
			
			// warning message about iteration
			// if (i == 199) fprintf(stderr, "# [CRED] WARNING: max iteration reached: p = %g\n", pval);
			
			if (temp <= EPS3 && temp <= EPS3 * pval) {
				// convergence?
				pval = pval * exp(pp * log(xx) + (qq - 1.0) * log(cx) - beta)/ pp;
				if (indx) {
					pval = 1.0 - pval;
				}
				break;
			}
			
			ai = ai + 1.0;
			ns = ns - 1;
			if (0 <= ns) {
				temp = qq - ai;
				if (ns == 0) {
					rx = xx;
				}
			} else {
				temp = psq;
				psq = psq + 1.0;
			}
			// printf("Current iteration: %d -> temp = %g and pval = %g\n", mm, temp, pval);
			++i;
		}
		if (t < 0.0) {
			// reverse p-value for a negative t-statistic
			return (1 - pval/2.0);
		}
		return pval/2.0; // one-tailed
	}
	return 1.0;
}
double test_student(int x[], int y[], int n, double *p) {
	double score = 0.0, x_sum = 0.01, y_sum = 0.01;
	double xf[n], yf[n];
	int i;
	for (i = 0; i < n; ++i) {
		x_sum = x_sum + (double) x[i];
		y_sum = y_sum + (double) y[i];
		xf[i] = (double) x[i];
		yf[i] = (double) y[i];
	}	
	score = x_sum/y_sum;
	*p = p_student(xf, sizeof(xf)/sizeof(*xf), yf, sizeof(yf)/sizeof(*yf));
	return(score);
}

// KOLMOGOROV-SMIRNOV TEST
double p_ks2(double D, double n) {
	double p0 = sqrt(n) * D;
    double p1_new, p1_old, s, w, z;
    int k, k_max = (int) sqrt(2.0 - LOG_EPS2);

    if (p0 < 1) {
	    z = KS_P_INIT / (p0 * p0);
	    w = log(p0);
	    s = 0;
	    for(k = 1; k < k_max; k += 2) {
			s += exp(k * k * z - w);
	    }
	    p1_new = s / INV_SQRT_2PI;
	} else {
	    z = -2 * p0 * p0;
	    s = -1;
	    k = 1;
	    p1_old = 0;
	    p1_new = 1;
	    while(fabs(p1_old - p1_new) > EPS2) {
			// fprintf(stderr, "p' = %g\n", new);
			p1_old = p1_new;
			p1_new += 2 * s * exp(z * k * k);
			s *= -1;
			++k;
	    }
	}
	return (1.0 - p1_new);
}
double stat_ks2(int *x1, int n1, int *x2, int n2) {
	int i, n = n1 + n2, Ns = 0;
	int_r *z = malloc(n * sizeof(int_r)); 
	
	// combine the two input arrays
	for (i = 0; i < n1; ++i) {
        z[i].x = x1[i];
		z[i].i = i;
	}
	for (i = n1; i < n; ++i) {
        z[i].x = x2[i - n1];
		z[i].i = i;
	}
	
	double n_eff = (double) (n1 * n2)/(n1 + n2); // effect size calc 
	qsort(z, n, sizeof(int_r), cmpfun_ranks);
	// z is sorted; use z.i to access ranks
  
	double csum = 0.0;
	double *s0 = malloc(n * sizeof(double));
	double *s1 = malloc(n * sizeof(double)); // c. sum without ties
	double *d = malloc((n - 1) * sizeof(double));
	
	for (i = 0; i < n; i++) {
		double v = -1.0/n2;
		if (z[i].i < n1) {
			v = 1.0/n1;
		}
		csum = csum + v;
		s0[i] = csum;
		
		// calculate differences in sorted data to identify ties
		if (i < (n - 1)) {
			d[i] = z[i + 1].x - z[i].x;
			if (d[i] > EPS2) {
				s1[Ns] = s0[i]; // no tie present, save to s1
				++Ns;
			}
		}
    }
	
	s1[Ns] = s0[n - 1];
	// two-sided D-statistic: max(abs(s))
	double D = array_maxabs(s1, Ns + 1); 
	double p = p_ks2(D, n_eff);
	free(z); free(s0); free(s1); free(d);
	return p;
}
double test_ks2(int x[], int y[], int n, double *p) {
	double score = 0.0, x_sum = 0.01, y_sum = 0.01;
	int i;
	for (i = 0; i < n; ++i) {
		x_sum = x_sum + (double) x[i];
		y_sum = y_sum + (double) y[i];
	}	
	score = x_sum/y_sum;
	*p = stat_ks2(x, n, y, n);
	return(score);
}

// CALLER FUNCTION FOR HYPOTHESIS TESTING
double region_tests(int tid, int head, int tail, khash_t(genome) *fg, khash_t(genome) *bg, double *p, const double r_xy, const int stage) {
	int i, j, *x, *y;
	int f_length = tail - head + 1;
	double score = 0.0;
	
	x = calloc(f_length, sizeof(int));
	y = calloc(f_length, sizeof(int));
	
	for (i = 0; i < f_length; ++i) {
		j = i + head; // position
		int fg_depth = 0;
		if (kh_exist(fg, j)) {
			fg_depth = kh_value(fg, j);
			// fprintf(stderr, "# [CRED] Running sum = %u\n", w_sum);
		}
		int bg_depth = 0;
		if (kh_exist(bg, j)) {
			bg_depth = kh_value(bg, j);
			// fprintf(stderr, "# [CRED] Running sum = %u\n", w_sum);
		}
		x[i] = fg_depth / r_xy; // normalize with per-chrom coverage ratio 
		y[i] = bg_depth;
	}
	
	if (stage > 1) {
		score = test_ks2(x, y, f_length, p);
	} else {
		score = test_student(x, y, f_length, p);
	}
	return(score);
}
int main(int argc, char *argv[]) {
    // argument defaults
	fg_bam = NULL; bg_bam = NULL; f_output = NULL;
	wsize_min = 300; wsize_max = 1200;
	int i, n = 2, cur_stage = 1; // Welch's test as default
    bam_hdr_t *fg_header = NULL; // BAM header of the 1st input
	
	int args_opt = 0;
	while ((args_opt = getopt(argc, argv, "t:c:m:w:p:q:k")) != -1) {
		switch(args_opt) {
			case 't': fg_bam = optarg; break;
			case 'c': bg_bam = optarg; break;		
			case 'w': wsize_max = parse_int(optarg); break;
			case 'p': p_cutoff = parse_double(optarg); break;
			case 'q': mapQ = parse_int(optarg); break;
			case 'k': cur_stage = 1; break;
			default: usage("Required parameters not specified!");
		}
	}	
	if (optind >  argc) usage("Required parameters not specified!");
	if (fg_bam == NULL) usage("Foreground BAM filename not specified!");
	if (bg_bam == NULL) usage("Background BAM filename not specified!");
		
	print_time("Initializing:");
	fprintf(stderr, "# [CRED] MAPQ tag quality cutoff set to %d\n", mapQ);
    // initialize the auxiliary data structures
	int *len_chrom, *fc_chrom, *bc_chrom; // xc_chrom keeps track of sorting
	long int *fg_sum, *bg_sum;
	double *fg_mdp, *bg_mdp;
	fg_data = calloc(n, sizeof(aux_t*));
	bg_data = calloc(n, sizeof(aux_t*));
	initialize_bam(fg_data, fg_bam, 1);
	initialize_bam(bg_data, bg_bam, 1);
		
	fg_header = fg_data[0]->hdr; // foreground BAM header
    int n_chrom = fg_header->n_targets;
	len_chrom = (int*) calloc(n_chrom, sizeof(int));
	fc_chrom = (int*) calloc(n_chrom, sizeof(int));
	bc_chrom = (int*) calloc(n_chrom, sizeof(int));
	fg_sum = (long int*) calloc(n_chrom, sizeof(long int));
	bg_sum = (long int*) calloc(n_chrom, sizeof(long int));
	fg_mdp = (double*) calloc(n_chrom, sizeof(double));
	bg_mdp = (double*) calloc(n_chrom, sizeof(double));
			
	char **chrom_names = malloc(n_chrom * sizeof(*chrom_names));
	khash_t(genome) *fg_coverage = (khash_t(genome)*) malloc (sizeof(khash_t(genome)) * n_chrom);	
	khash_t(genome) *bg_coverage = (khash_t(genome)*) malloc (sizeof(khash_t(genome)) * n_chrom);	

	for (i = 0; i < fg_header->n_targets; ++i) {
		int l_chrom = (int) fg_header->target_len[i];
		chrom_names[i] = fg_header->target_name[i];
		len_chrom[i] = l_chrom;
		khash_t(genome) *fg_cur_chrom = kh_init(genome);
		kh_resize(genome, fg_cur_chrom, l_chrom);
		fg_coverage[i] = *fg_cur_chrom;
		
		khash_t(genome) *bg_cur_chrom = kh_init(genome);
		kh_resize(genome, bg_cur_chrom, l_chrom);
		bg_coverage[i] = *bg_cur_chrom;		
	}
			
	// load BAM's and calculate number of positive coverage sites
	print_time("Loading treatment BAM...");
	int n_fg_nonzero = process_bam(fg_coverage, fg_mdp, (void **) fg_data, chrom_names, len_chrom, fg_sum, n_chrom, fc_chrom);
	fprintf(stderr, "\n# [CRED] Treatment BAM loaded into memory (%u bases).\n", n_fg_nonzero);
	print_time("Loading control BAM...");
	int n_bg_nonzero = process_bam(bg_coverage, bg_mdp, (void **) bg_data, chrom_names, len_chrom, bg_sum, n_chrom, bc_chrom);
	fprintf(stderr,"\n# [CRED] Control BAM loaded into memory (%u bases).\n", n_bg_nonzero);
	
	// creating candidate regions		
	fprintf(stderr, "# [CRED] Initializing search windows...\n");
	regions *features_list = NULL;
	if (abs(n_fg_nonzero) > FEATURE_LIMIT) {
		features_list = malloc(FEATURE_LIMIT * sizeof(regions));
		n_fg_nonzero = FEATURE_LIMIT;
	} else {
		features_list = malloc(n_fg_nonzero * sizeof(regions));
	}
	if (features_list == NULL) {
		fprintf(stderr, "# [CRED] WARNING: memory allocation error encountered. Retrying...");
		n_fg_nonzero = FEATURE_LIMIT;
		while (features_list == NULL) {
			--n_fg_nonzero;
			features_list = malloc(n_fg_nonzero * sizeof(regions));
		}
		fprintf(stderr, "%d final\n", n_fg_nonzero);
	}
	// fprintf(stderr, "Number of non-zero sites = %d\n", n_nonzero);

	print_time("Evaluating candidate regions...");
	fprintf(stderr, "# [CRED] Reference sliding window size: %d bp\n", wsize_max);
	int b_nonzero = n_fg_nonzero; int f_counts = 0; 
	int n_fsize = n_fg_nonzero;
	
	// recount chromosomes here
	int np_chrom = 0;
	for (i = 0; i < n_chrom; ++i) {
		if (fc_chrom[i] == 1) {
			np_chrom++;
		} else {
			fc_chrom[i] = 0;
		}
	}
	double *f_dpratio = calloc(n_chrom, sizeof(double));
	for (i = 0; i < n_chrom; ++i) {
		f_dpratio[i] = fg_mdp[i]/bg_mdp[i];
		if (fc_chrom[i] > 0) {
			fprintf(stderr, "# [CRED] %s (rel. T:C coverage %.3f:%.3f = %.3f)", chrom_names[i], fg_mdp[i], bg_mdp[i], f_dpratio[i]);
			int cur_fsize = set_window(features_list, &fg_coverage[i], i, fg_mdp[i], len_chrom[i], wsize_max, &f_counts, &n_fsize);
			if (cur_fsize < 0) {
				// memory allocation error occured
				fprintf(stderr, "# [CRED] Memory allocation error encountered; stopping discovery...\n");
				break;
			}
		}
	}
	
	// merge sites
	int region_counts = 0;
	// if (wsize_min < 300) wsize_min = 300;
	region_merge(features_list, b_nonzero, f_counts, &region_counts);
	regions *sites_step0 = malloc(region_counts * sizeof(regions));
	fprintf(stderr, "# [CRED] %d candidate regions identified.\n", region_counts);
	
	// hypothesis testing
	for (i = 0; i < region_counts; ++i) {
		int tid = features_list[i].tid;
		double f_p;
		if (isfinite(f_dpratio[tid])) {
			double f_score = region_tests(tid, features_list[i].start, features_list[i].end, &fg_coverage[tid], &bg_coverage[tid], &f_p, f_dpratio[tid], cur_stage);
			sites_step0[i] = features_list[i];
			sites_step0[i].score = f_score;
			sites_step0[i].p = f_p;
		}
	}
	free(features_list);
	
	// site assessment
	fprintf(stderr, "# [CRED] Assessing candidate region significance...\n");
	fprintf(stderr, "# [CRED] Method: %s\n",
		((cur_stage < 2) ? "Welch's t-test" : "Kolmogorov-Smirnov test"));
	fprintf(stderr, "# [CRED] Significance level set to p = %f\n", p_cutoff);
	if (p_cutoff < EPS2) {
		p_cutoff = EPS2;
		fprintf(stderr, "# [CRED] WARNING: significance level coerced to %.8f\n", p_cutoff);
	}
	int filter_counts = region_print(sites_step0, chrom_names, region_counts, p_cutoff);
	fprintf(stderr, "# [CRED] %d final enrichment sites identified.\n", filter_counts);
	
	free(sites_step0);
	print_time("Process completed.");
	return(status);
}

/*
References for Hypothesis Testing:
1. Welch's t-test:
	a. Incomplete Beta function ratio:
	Majumder KL & Bhattacharjee GP, Appl. Stat. 22(3): 409-411, 1973.
	C version by J. Burkardt, modified here exclusively for one-tailed case of x > y
	b. https://rosettacode.org/wiki/Welch%27s_t-test#C

2. Kolmogorov-Smirnov test:
	J. Durbin, Distribution Theory for Tests Based on the Sample Distribution Function (1973)
*/