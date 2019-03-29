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
te_calc_gm_set_P(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval);

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_gm_set_P = { "gm_set_P",
		"    Precision using geometric mean over the topics.\n",
		te_init_meas_s_float, te_calc_gm_set_P, te_acc_meas_s,
		te_calc_avg_meas_s_gm, te_print_single_meas_empty,
		te_print_final_meas_s_float, NULL, -1 };

static int te_calc_gm_set_P(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval) {
	RES_RELS res_rels;
	double prec;

	if (UNDEF == te_form_res_rels(epi, rel_info, results, &res_rels))
		return (UNDEF);

	if (res_rels.num_ret) {
		prec = (double) res_rels.num_rel_ret / (double) res_rels.num_ret;
		eval->values[tm->eval_index].value = (double) log(
				(double) (MAX(prec, MIN_GEO_MEAN)));
	}
	return (1);
}
