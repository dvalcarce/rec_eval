/*
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
te_calc_gm_err(const EPI *epi, const REL_INFO *rel_info, const RESULTS *results,
		const TREC_MEAS *tm, TREC_EVAL *eval);

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_gm_err = { "gm_err",
		"    ERR using geometric mean over the topics.\n", te_init_meas_s_float,
		te_calc_gm_err, te_acc_meas_s, te_calc_avg_meas_s,
		te_print_single_meas_empty, te_print_final_meas_s_float, NULL, -1 };

static double compute_gain(const long rel_level, const long max_rel);

static int te_calc_gm_err(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval) {

	RES_RELS res_rels;
	double p;
	double err;
	double gain;
	long r;

	if (UNDEF == te_form_res_rels(epi, rel_info, results, &res_rels))
		return (UNDEF);

	p = 1.0;
	err = 0.0;
	for (r = 0; r < res_rels.num_ret; r++) {
		gain = compute_gain(res_rels.results_rel_list[r],
				res_rels.num_rel_levels);
		err += p * gain / (double) (r + 1);
		p *= 1.0 - gain;
	}
	eval->values[tm->eval_index].value = (double) log(
			(double) (MAX(err, MIN_GEO_MEAN)));
	;

	return (1);

}

static inline double compute_gain(const long rel_level, const long max_rel) {
	return rel_level > 0 ?
			(pow(2.0, rel_level) - 1.0) / pow(2.0, (max_rel - 1)) : 0.0;
}
