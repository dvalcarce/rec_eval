/*
 Copyright (c) 2017 - Daniel Valcarce.

 Permission is granted for use and modification of this file for
 research, non-commercial purposes.
 */

#include "common.h"
#include "sysfunc.h"
#include "trec_eval.h"
#include "functions.h"
#include "trec_format.h"

static int
te_calc_q_measure(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval);
static PARAMS default_q_measure_gains = { NULL, 0, NULL };

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_q_measure =
		{ "q_measure",
				"    Q-measure.\n\
    TODO.\n\
    Gain values are set to the appropriate relevance level by default.  \n\
    The default gain can be overridden on the command line by having \n\
    comma separated parameters 'rel_level=gain'.\n\
    Eg, 'trec_eval -m ndcg.1=3.5,2=9.0,4=7.0 ...'\n\
    will give gains 3.5, 9.0, 3.0, 7.0 for relevance levels 1,2,3,4\n\
    respectively (level 3 remains at the default).\n\
    Gains are allowed to be 0 or negative, and relevance level 0\n\
    can be given a gain.\n\
	Cite:  Tetsuya Sakai and Noriko Kando: On information retrieval metrics\n\
	designed for evaluation with incomplete relevance assessments. In\n\
	Information	Retrieval 11, 5 (2008), 447-470.\n\
	DOI=http://dx.doi.org/10.1007/s10791-008-9059-7\n",
				te_init_meas_s_float_p_pair, te_calc_q_measure, te_acc_meas_s,
				te_calc_avg_meas_s, te_print_single_meas_s_float,
				te_print_final_meas_s_float_p, &default_q_measure_gains, -1 };

/* Keep track of valid rel_levels and associated gains */
/* Initialized in setup_gains */
typedef struct {
	long rel_level;
	long num_at_level;
	double gain;
} REL_GAIN;

typedef struct {
	REL_GAIN *rel_gains;
	long num_gains;
	long total_num_at_levels;
} GAINS;

static int setup_gains(const TREC_MEAS *tm, const RES_RELS *res_rels,
		GAINS *gains);
static double get_gain(const long rel_level, const GAINS *gains);
static int comp_rel_gain();

static int te_calc_q_measure(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval) {

	RES_RELS res_rels;
	long i;
	long rel_so_far = 0;
	long cur_lvl, lvl_count = 0;
	double q_measure = 0.0;
	double cum_gain = 0.0;
	GAINS gains;
	double *cgi;

	if (UNDEF == te_form_res_rels(epi, rel_info, results, &res_rels)) {
		return (UNDEF);
	}

	if (UNDEF == setup_gains(tm, &res_rels, &gains)) {
		return (UNDEF);
	}

	/* Precompute ideal cumulative gains (cgi) at each position */
	cgi = (double *) malloc(res_rels.num_rel * sizeof(double));
	lvl_count = 0;
	cur_lvl = res_rels.num_rel_levels - 1;
	for (i = 0; 1; i++) {
		lvl_count++;
		while (cur_lvl > 0 && lvl_count > res_rels.rel_levels[cur_lvl]) {
			cur_lvl--;
			lvl_count = 1;
		}
		if (cur_lvl == 0) {
			break;
		}
		cgi[i] = get_gain(cur_lvl, &gains);
		if (i > 0) {
			cgi[i] += cgi[i - 1];
		}
	}

	/* Compute Q-measure */
	for (i = 0; i < res_rels.num_ret; i++) {
		if (res_rels.results_rel_list[i] >= epi->relevance_level) {
			rel_so_far++;
			cum_gain += get_gain(res_rels.results_rel_list[i], &gains);
			q_measure += (double) (cum_gain + rel_so_far)
					/ (double) (i + 1 + cgi[i]);
		}
	}

	if (res_rels.num_rel) {
		q_measure /= res_rels.num_rel;
	}

	free(cgi);

	eval->values[tm->eval_index].value = q_measure;
	return (1);
}

static int setup_gains(const TREC_MEAS *tm, const RES_RELS *res_rels,
		GAINS *gains) {
	FLOAT_PARAM_PAIR *pairs = NULL;
	long num_pairs = 0;
	long i, j;
	long num_gains;

	if (tm->meas_params) {
		pairs = (FLOAT_PARAM_PAIR *) tm->meas_params->param_values;
		num_pairs = tm->meas_params->num_params;
	}

	if (NULL
			== (gains->rel_gains = Malloc(res_rels->num_rel_levels + num_pairs,
					REL_GAIN)))
		return (UNDEF);
	num_gains = 0;
	for (i = 0; i < num_pairs; i++) {
		gains->rel_gains[num_gains].rel_level = atol(pairs[i].name);
		gains->rel_gains[num_gains].gain = (double) pairs[i].value;
		gains->rel_gains[num_gains].num_at_level = 0;
		num_gains++;
	}

	for (i = 0; i < res_rels->num_rel_levels; i++) {
		for (j = 0; j < num_gains && gains->rel_gains[j].rel_level != i; j++)
			;
		if (j < num_gains)
			/* Was included in list of parameters. Update occurrence info */
			gains->rel_gains[j].num_at_level = res_rels->rel_levels[i];
		else {
			/* Not included in list of parameters. New gain level */
			gains->rel_gains[num_gains].rel_level = i;
			gains->rel_gains[num_gains].gain = (double) i;
			gains->rel_gains[num_gains].num_at_level = res_rels->rel_levels[i];
			num_gains++;
		}
	}

	/* Sort gains by increasing gain value */
	qsort((char *) gains->rel_gains, (int) num_gains, sizeof(REL_GAIN),
			comp_rel_gain);

	gains->total_num_at_levels = 0;
	for (i = 0; i < num_gains; i++)
		gains->total_num_at_levels += gains->rel_gains[i].num_at_level;

	gains->num_gains = num_gains;
	return (1);
}

static int comp_rel_gain(REL_GAIN *ptr1, REL_GAIN *ptr2) {
	return (ptr1->gain - ptr2->gain);
}

static double get_gain(const long rel_level, const GAINS *gains) {
	long i;
	for (i = 0; i < gains->num_gains; i++)
		if (rel_level == gains->rel_gains[i].rel_level)
			return (gains->rel_gains[i].gain);
	return (0.0); /* Print Error ?? */
}

