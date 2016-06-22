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
te_calc_err (const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm,
		TREC_EVAL *eval);
static double err_param_array[] = {5.0};
static PARAMS default_err_param = {
	NULL, sizeof (err_param_array) / sizeof (err_param_array[0]),
	&err_param_array[0]};

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_err =
{"err",
	"    Expected Reciprocal Rank of the first relevant retrieved doc.\n\
		Measure is most useful for tasks in which there is only one relevant\n\
		doc, or the user only wants one relevant doc using graded relevance.\n\
		This measure accepts as a param the maximum rating in the collection;\n\
		otherwise, 5.0 is assumed as maximum value.\n\
		Use:  trec_eval -m err.4.0\n",
	te_init_meas_s_float_p_float,
	te_calc_err,
	te_acc_meas_s,
	te_calc_avg_meas_s,
	te_print_single_meas_s_float,
	te_print_final_meas_s_float_p,
	(void *) &default_err_param, -1};

static double setup (const TREC_MEAS *tm);
static double compute_gain (const long rel_level, const double max_rating);

	static int
te_calc_err (const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm,
		TREC_EVAL *eval)
{
	RES_RELS res_rels;
	double max_rating;
	double p;
	double err;
	double gain;
	long r;

	if (UNDEF == te_form_res_rels (epi, rel_info, results, &res_rels))
		return (UNDEF);

	max_rating = setup (tm);

	p = 1.0;
	err = 0.0;
	for (r = 0; r < res_rels.num_ret; r++) {
		gain = compute_gain(res_rels.results_rel_list[r], max_rating);
		err = err + p * gain / (double) (r + 1.0);
		p = p * (1.0 - gain);
	}
	eval->values[tm->eval_index].value = err;

	return (1);
}


	static double
setup (const TREC_MEAS *tm)
{
	if (tm->meas_params)
		return *((double *) tm->meas_params->param_values);
	else
		return (5.0); // Default value
}

	static double
compute_gain (const long rel_level, const double max_rating)
{
	return rel_level <= 0 ? 0.0 :
		(pow(2.0, (double) rel_level) - 1.0) / pow(2.0, max_rating);
}
