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
te_calc_gm_recip_rank(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval);

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_gm_recip_rank = { "gm_recip_rank",
		"    Reciprocal Rank using geometric mean over the topics..\n",
		te_init_meas_s_float, te_calc_gm_recip_rank, te_acc_meas_s,
		te_calc_avg_meas_s, te_print_single_meas_empty,
		te_print_final_meas_s_float,
		NULL, -1 };

static int te_calc_gm_recip_rank(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval) {
	RES_RELS res_rels;
	long i;

	if (UNDEF == te_form_res_rels(epi, rel_info, results, &res_rels))
		return (UNDEF);

	for (i = 0; i < res_rels.num_ret; i++) {
		if (res_rels.results_rel_list[i] >= epi->relevance_level)
			break;
	}
	if (i < res_rels.num_ret) {
		eval->values[tm->eval_index].value =
				(double) log(
						(double) (MAX(((double ) 1.0 / (double ) (i + 1)),
								MIN_GEO_MEAN)));
		;
	}
	return (1);
}
