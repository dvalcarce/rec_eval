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
te_calc_gm_set_F(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval);
static double set_F_param_array[] = { 1.0 };
static PARAMS default_set_F_params = {
NULL, sizeof(set_F_param_array) / sizeof(set_F_param_array[0]),
		&set_F_param_array[0] };

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_gm_set_F = { "gm_set_F",
		"      F-measure using geometric mean over the topics.\n",
		te_init_meas_s_float_p_float, te_calc_gm_set_F, te_acc_meas_s,
		te_calc_avg_meas_s_gm, te_print_single_meas_empty,
		te_print_final_meas_s_float_p, (void *) &default_set_F_params, -1 };

static int te_calc_gm_set_F(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval) {
	double *params = (double *) tm->meas_params->param_values;
	RES_RELS res_rels;
	double x, P, R, f1;

	if (UNDEF == te_form_res_rels(epi, rel_info, results, &res_rels))
		return (UNDEF);

	if (res_rels.num_rel_ret) {
		P = (double) res_rels.num_rel_ret / (double) res_rels.num_ret;
		R = (double) res_rels.num_rel_ret / (double) res_rels.num_rel;
		x = params[0];
		f1 = (x + 1.0) * P * R / (x * P + R);
		eval->values[tm->eval_index].value = (double) log(
				(double) (MAX(f1, MIN_GEO_MEAN)));
	}
	return (1);
}
