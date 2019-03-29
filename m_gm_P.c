/*
 Copyright (c) 2008 - Chris Buckley.
 Copyright (c) 2019 - Daniel Valcarce.

 Permission is granted for use and modification of this file for
 research, non-commercial purposes.
 */
#include "common.h"
#include "sysfunc.h"
#include "trec_eval.h"
#include "functions.h"
#include "trec_format.h"

static int
te_calc_gm_P(const EPI *epi, const REL_INFO *rel_info, const RESULTS *results,
		const TREC_MEAS *tm, TREC_EVAL *eval);
static long long_cutoff_array[] = { 5, 10, 15, 20, 30, 100 };
static PARAMS default_P_cutoffs = {
NULL, sizeof(long_cutoff_array) / sizeof(long_cutoff_array[0]),
		&long_cutoff_array[0] };

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_gm_P = { "gm_P",
		"    Precision using geometric mean over the topics.\n",
		te_init_meas_a_float_cut_long, te_calc_gm_P, te_acc_meas_a_cut,
		te_calc_avg_meas_a_cut, te_print_single_meas_empty,
		te_print_final_meas_a_cut, (void *) &default_P_cutoffs, -1 };

static int te_calc_gm_P(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval) {
	long *cutoffs = (long *) tm->meas_params->param_values;
	long cutoff_index = 0;
	long i;
	RES_RELS res_rels;
	long rel_so_far = 0;
	double prec;

	if (UNDEF == te_form_res_rels(epi, rel_info, results, &res_rels))
		return (UNDEF);

	for (i = 0; i < res_rels.num_ret; i++) {
		if (i == cutoffs[cutoff_index]) {
			/* Calculate previous cutoff threshold.
			 Note all guaranteed to be positive by init_meas */
			prec = (double) rel_so_far / (double) i;
			eval->values[tm->eval_index + cutoff_index].value = (double) log(
					(double) (MAX(prec, MIN_GEO_MEAN)));
			if (++cutoff_index == tm->meas_params->num_params)
				break;
		}
		if (res_rels.results_rel_list[i] >= epi->relevance_level)
			rel_so_far++;
	}
	/* calculate values for those cutoffs not achieved */
	while (cutoff_index < tm->meas_params->num_params) {
		prec = (double) rel_so_far / (double) cutoffs[cutoff_index];
		eval->values[tm->eval_index + cutoff_index].value = (double) log(
				(double) (MAX(prec, MIN_GEO_MEAN)));
		cutoff_index++;
	}
	return (1);
}
