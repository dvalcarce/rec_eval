/* 
 Copyright (c) 2008 - Chris Buckley. 

 Permission is granted for use and modification of this file for
 research, non-commercial purposes. 
 */

#include "common.h"
#include "sysfunc.h"
#include "trec_eval.h"
#include "functions.h"
#include "trec_format.h"

static int
te_calc_gm_set_recall(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval);

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_gm_set_recall = { "gm_set_recall",
		"    Recall using geometric mean over the topics.\n",
		te_init_meas_s_float, te_calc_gm_set_recall, te_acc_meas_s,
		te_calc_avg_meas_s_gm, te_print_single_meas_empty,
		te_print_final_meas_s_float,
		NULL, -1 };

static int te_calc_gm_set_recall(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval) {
	RES_RELS res_rels;
	double recall;

	if (UNDEF == te_form_res_rels(epi, rel_info, results, &res_rels))
		return (UNDEF);

	if (res_rels.num_rel) {
		recall = (double) res_rels.num_rel_ret / (double) res_rels.num_rel;
		eval->values[tm->eval_index].value = (double) log(
				(double) (MAX(recall, MIN_GEO_MEAN)));
	}

	return (1);
}
