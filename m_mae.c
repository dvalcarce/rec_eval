/*
 Copyright (c) 2016 - Daniel Valcarce.

 Permission is granted for use and modification of this file for
 research, non-commercial purposes.
 */

#include "common.h"
#include "sysfunc.h"
#include "trec_eval.h"
#include "functions.h"
#include "trec_format.h"

static int
te_calc_mae(const EPI *epi, const REL_INFO *rel_info, const RESULTS *results,
		const TREC_MEAS *tm, TREC_EVAL *eval);

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_mae =
		{ "mae",
				"    Mean Absolute Error\n\
		MAE is a standard error metric.\n\
		It measures the difference between of the predicted rating\n",
				te_init_meas_s_float, te_calc_mae, te_acc_meas_s,
				te_calc_avg_meas_s, te_print_single_meas_s_float,
				te_print_final_meas_s_float, NULL, -1 };

static int te_calc_mae(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval) {
	RES_RELS res_rels;
	double mae;
	double actual_score;
	double predicted_score;
	int i;
	int n;

	if (UNDEF == te_form_res_rels(epi, rel_info, results, &res_rels))
		return (UNDEF);

	mae = 0.0;
	n = 0;
	for (i = 0; i < res_rels.num_ret; i++) {
		predicted_score =
				((TEXT_RESULTS_INFO *) results->q_results)->text_results->sim;
		actual_score = res_rels.results_rel_list[i];
		if (actual_score >= 0.0) {
			n++;
			mae += fabs(actual_score - predicted_score);
		}
	}

	if (res_rels.num_ret)
		eval->values[tm->eval_index].value = n > 0 ? mae / n : 0.0;

	return (1);
}
