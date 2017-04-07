/*
 Copyright (c) 2008 - Chris Buckley.
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
te_calc_infap2(const EPI *epi, const REL_INFO *rel_info, const RESULTS *results,
		const TREC_MEAS *tm, TREC_EVAL *eval);

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_infAP2 =
		{ "infAP2",
				"    Inferred AP 2 (modified) \n\
    A measure that allows sampling of judgements: Qrels/results divided\n\
    into unjudged, judged_rel and pooled_judged_nonrel.\n\
    Cite:    'Estimating Average Precision with Incomplete and Imperfect\n\
    Judgments', Emine Yilmaz and Javed A. Aslam. CIKM \n",
				te_init_meas_s_float, te_calc_infap2, te_acc_meas_s,
				te_calc_avg_meas_s, te_print_single_meas_s_float,
				te_print_final_meas_s_float,
				NULL, -1 };

static int te_calc_infap2(const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval) {
	RES_RELS res_rels;
	long j;
	long nonrel_so_far, rel_so_far, pool_unjudged_so_far;
	double inf_ap = 0.0;

	if (UNDEF == te_form_res_rels(epi, rel_info, results, &res_rels))
		return (UNDEF);

	nonrel_so_far = 0;
	rel_so_far = 0;
	pool_unjudged_so_far = 0;
	for (j = 0; j < res_rels.num_ret; j++) {
		if (res_rels.results_rel_list[j] == RELVALUE_NONPOOL) {
			/* document not in pool -> unjudged */
			pool_unjudged_so_far++;
			continue;
		}
		if (res_rels.results_rel_list[j] == RELVALUE_UNJUDGED) {
			/* document in pool but unjudged. */
			pool_unjudged_so_far++;
			continue;
		}

		if (res_rels.results_rel_list[j] >= 0
				&& res_rels.results_rel_list[j] < epi->relevance_level) {
			nonrel_so_far++;
		} else {
			/* Judged Rel doc */
			rel_so_far++;
			/* Special case nonrel_so_far == 0 to avoid division by 0 */
			/* inf_ap */
			if (0 == j)
				inf_ap += 1.0;
			else {
				double fj = (double) j;
				inf_ap += 1.0 / (fj + 1.0)
						+ (fj / (fj + 1.0))
								* ((rel_so_far - 1 + nonrel_so_far
										+ pool_unjudged_so_far) / fj)
								* ((rel_so_far - 1 + INFAP_EPSILON)
										/ (rel_so_far - 1 + nonrel_so_far
												+ 2 * INFAP_EPSILON));
			}
		}
	}
	if (res_rels.num_rel) {
		inf_ap /= res_rels.num_rel;
	}
	eval->values[tm->eval_index].value = inf_ap;

	return (1);
}
