// Copyright (C) 2017 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

int ZeroEG_d2(int pathNum, point_data_d *Final, point_d last_approx, point_data_d *lastSample, int *cycle_num, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

int ZeroEG_mp2(int pathNum, point_data_mp *Final, point_mp last_approx, point_data_mp *lastSample, int *cycle_num, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

int ZeroEG_amp2(int pathNum, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d *lastSample_d, point_data_mp *lastSample_mp, int *cycle_num, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));


int ZeroEG_d(int pathNum, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks directly to endpoint in fixed double precision  *
\***************************************************************/
{
  int retVal, cycle_num;
  point_data_d lastSample;

  init_point_data_d(&lastSample, 0);

  // setup last_approx in case of failure
  point_cp_d(last_approx, Start->point);

  retVal = ZeroEG_d2(pathNum, Final, last_approx, &lastSample, &cycle_num, Start, T, OUT, midOUT, ED, eval_func, find_dehom);

  clear_point_data_d(&lastSample);

  return retVal;
}

int ZeroEG_rank_d(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks directly to endpoint in fixed double precision  *
\***************************************************************/
{
  int i, retVal, cycle_num, samples_per_loop = T->num_PSEG_sample_points + 1;
  comp_d finalTime;
  point_data_d lastSample;

  init_point_data_d(&lastSample, 0);

  // setup finalTime
  set_double_d(finalTime, T->targetT, 0);

  // setup last_approx in case of failure
  point_cp_d(last_approx, Start->point);

  // perform the endgame - returns the last sample point so that we can use it for rank determination
  retVal = ZeroEG_d2(pathNum, Final, last_approx, &lastSample, &cycle_num, Start, T, OUT, midOUT, ED, eval_func, find_dehom);

  if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
  { // continue on with rank determination
    retVal = Cauchy_rank_main_d(condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, retVal, Final, last_approx, finalTime, &lastSample, cycle_num, samples_per_loop, T, OUT, ED, eval_func);
  }
  else
  { // setup condNum since we have failure
    *condNum = -1;
    *smallest_nonzero_SV = *largest_zero_SV = 0;

    // setup last_approx
    point_cp_d(last_approx, Final->point);
    for (i = 0; i < last_approx->size; i++)
    {
      get_comp_rand_d(finalTime);
      mul_rdouble_d(finalTime, finalTime, T->final_tolerance);
      add_d(&last_approx->coord[i], &last_approx->coord[i], finalTime);
    }

    // assume failure since we could not converge
    if (rankType == 0)
      *rankDef = -1;
    else
      *corank = 1;
  }

  clear_point_data_d(&lastSample);

  return retVal;
}

int ZeroEG_d2(int pathNum, point_data_d *Final, point_d last_approx, point_data_d *lastSample, int *cycle_num, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks directly to endpoint in fixed double precision  *
\***************************************************************/
{
  int i, retVal;
  comp_d endTime;
  point_data_d Pt;

  // copy Start to Pt
  init_point_data_d(&Pt, Start->point->size); 
  point_data_cp_d(&Pt, Start);

  // setup last_approx in case of failure
  point_cp_d(last_approx, Start->point);

  // first, we need to track to the endgameBoundary - it is assumed the T->endgameBoundary >= 0 == T->targetT
  set_double_d(endTime, T->endgameBoundary, 0.0);
  T->currentNewtonTol = T->basicNewtonTol;
  T->minStepSize = T->minStepSizeBeforeEndGame;
  T->endgameSwitch = 0;
  T->error_at_latest_sample_point = 0.0;

  // track to the endgame boundary
  retVal = track_d(&Pt, &Pt, endTime, T, OUT, ED, eval_func);

  // print the ending point to midOUT
  fprintf(midOUT, "%d\n", pathNum);
  for (i = 0; i < Pt.point->size; i++)
    fprintf(midOUT, "%.15e %.15e\n", Pt.point->coord[i].r, Pt.point->coord[i].i);
  fprintf(midOUT, "\n");

  // since we have printed to midOUT, this path is considered to be "in the endgame"
  T->endgameSwitch = 1;

  // check to see if we made it to the endgame boundary successfully
  if (retVal)
  { // failure
    fprintf(OUT, "NOTE: Zeroth endgame never started!\n");

    // copy the last sample to Final
    point_data_cp_d(Final, &Pt);
    Final->cycle_num = *cycle_num = 0;
    T->t_val_at_latest_sample_point = Pt.time->r;
  }
  else
  { // success - so run the actual endgame now that we are at the endgame boundary

    // set the endgame tolerances
    set_double_d(endTime, T->targetT, 0.0);
    T->currentNewtonTol = T->endgameNewtonTol;
    T->minStepSize = T->minStepSizeDuringEndGame;

    // track to the end 
    retVal = track_d(Final, &Pt, endTime, T, OUT, ED, eval_func);
    
    // setup lastSample & cycle_num
    *cycle_num = Final->cycle_num = 1; 
    point_data_cp_d(lastSample, Final);
    point_cp_d(last_approx, Final->point);
    T->t_val_at_latest_sample_point = Final->time->r;

    if (retVal == 0)
    { // refine to precision requested in FinalTol 
      int sharpenDigits = ceil(-log10(T->final_tolerance));
      retVal = sharpen_d(T->outputLevel, sharpenDigits, &T->latest_newton_residual_d, Final, Final, OUT, ED, eval_func);
    }
  }

  // free the memory
  clear_point_data_d(&Pt);

  return retVal;
}

///////////// MAIN MP VERSIONS //////////////////

int ZeroEG_mp(int pathNum, point_data_mp *Final, point_mp last_approx, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks directly to endpoint in fixed multi precision   *
\***************************************************************/
{
  int retVal, cycle_num;
  point_data_mp lastSample;

  init_point_data_mp(&lastSample, 0);

  // setup last_approx in case of failure
  point_cp_mp(last_approx, Start->point);

  retVal = ZeroEG_mp2(pathNum, Final, last_approx, &lastSample, &cycle_num, Start, T, OUT, midOUT, ED, eval_func, find_dehom);

  // clear
  clear_point_data_mp(&lastSample);

  return retVal;
}

int ZeroEG_rank_mp(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_mp *Final, point_mp last_approx, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks directly to endpoint in fixed multi precision   *
\***************************************************************/
{
  int i, retVal, cycle_num, samples_per_loop = T->num_PSEG_sample_points + 1;
  comp_mp finalTime;
  point_data_mp lastSample;

  init_mp(finalTime);
  init_point_data_mp(&lastSample, 0);

  // setup finalTime
  set_double_mp(finalTime, T->targetT, 0);

  // setup last_approx in case of failure
  point_cp_mp(last_approx, Start->point);

  // perform the endgame - returns that last sample so that it can be used for rank determination
  retVal = ZeroEG_mp2(pathNum, Final, last_approx, &lastSample, &cycle_num, Start, T, OUT, midOUT, ED, eval_func, find_dehom);

  if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
  { // continue on with rank determination
    retVal = Cauchy_rank_main_mp(condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, retVal, Final, last_approx, finalTime, &lastSample, cycle_num, samples_per_loop, T, OUT, ED, eval_func);
  }
  else
  { // setup condNum since we have failure
    *condNum = -1;
    *smallest_nonzero_SV = *largest_zero_SV = 0;

    // setup last_approx
    point_cp_mp(last_approx, Final->point);
    for (i = 0; i < last_approx->size; i++)
    {
      get_comp_rand_mp(finalTime);
      mul_rdouble_mp(finalTime, finalTime, T->final_tolerance);
      add_mp(&last_approx->coord[i], &last_approx->coord[i], finalTime);
    }
  }

  // clear
  clear_mp(finalTime);
  clear_point_data_mp(&lastSample);

  return retVal;
}

int ZeroEG_mp2(int pathNum, point_data_mp *Final, point_mp last_approx, point_data_mp *lastSample, int *cycle_num, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks directly to endpoint using fixed multi precision*
\***************************************************************/
{
  int i, retVal;
  comp_mp endTime;
  point_data_mp Pt; 

  // initialize endTime
  init_mp(endTime);

  // initialize Pt 
  init_point_data_mp(&Pt, Start->point->size);
  point_data_cp_mp(&Pt, Start);

  // setup last_approx in case of failure
  point_cp_mp(last_approx, Start->point);

  // first, we need to track to the endgameBoundary - it is assumed the T->endgameBoundary >= 0 == T->targetT
  set_double_mp(endTime, T->endgameBoundary, 0.0);
  T->currentNewtonTol = T->basicNewtonTol;
  T->minStepSize = T->minStepSizeBeforeEndGame;
  T->endgameSwitch = 0;
  T->error_at_latest_sample_point = 0.0;

  // track
  retVal = track_mp(&Pt, &Pt, endTime, T, OUT, ED, eval_func);

  // print the ending point to midOUT
  fprintf(midOUT, "%d\n", pathNum);
  for (i = 0; i < Pt.point->size; i++)
  {
    print_mp(midOUT, 0, &Pt.point->coord[i]);
    fprintf(midOUT, "\n");
  }
  fprintf(midOUT, "\n");

  // since we have printed to midOUT, this path is considered to be "in the endgame"
  T->endgameSwitch = 1;

  // check to see if we made it to the endgame boundary successfully
  if (retVal)
  { // failure
    fprintf(OUT, "NOTE: Zeroth endgame never started!\n");

    // copy the last sample to Final
    point_data_cp_mp(Final, &Pt);
    Final->cycle_num = *cycle_num = 0;
    T->t_val_at_latest_sample_point = mpf_get_d(Pt.time->r);
  }
  else
  { // success - so run the actual endgame now that we are at the endgame boundary

    // set the endgame tolerances
    set_double_mp(endTime, T->targetT, 0.0);
    T->currentNewtonTol = T->endgameNewtonTol;
    T->minStepSize = T->minStepSizeDuringEndGame;

    // track to the end 
    retVal = track_mp(Final, &Pt, endTime, T, OUT, ED, eval_func);

    // setup lastSample & cycle_num
    *cycle_num = Final->cycle_num = 1; 
    point_data_cp_mp(lastSample, Final);
    point_cp_mp(last_approx, Final->point);
    T->t_val_at_latest_sample_point = mpf_get_d(Final->time->r);

    if (retVal == 0)
    { // refine to precision requested in FinalTol
      int sharpenDigits = ceil(-log10(T->final_tolerance));
      retVal = sharpen_mp(T->outputLevel, sharpenDigits, &T->latest_newton_residual_d, T->latest_newton_residual_mp, Final, Final, OUT, ED, eval_func);
    }
  }

  // free the memory
  clear_mp(endTime);
  clear_point_data_mp(&Pt);

  return retVal;
}

///////////// MAIN AMP VERSIONS //////////////////

int ZeroEG_amp(int pathNum, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, int *last_approx_prec, point_d last_approx_d, point_mp last_approx_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: prec - the ending precision - < 64 - Final_d   *
*                                             >= 64 - Final_mp  *
* time_first_increase - the time of the last increase to MP     *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal, cycle_num;
  point_data_d lastSample_d;
  point_data_mp lastSample_mp;

  init_point_data_d(&lastSample_d, 0);
  init_point_data_mp(&lastSample_mp, 0);

  // setup last_approx in case of failure
  *last_approx_prec = prec_in;
  if (prec_in < 64)
  {
    point_cp_d(last_approx_d, Start_d->point);
  }
  else
  {
    setprec_point_mp(last_approx_mp, prec_in);
    point_cp_mp(last_approx_mp, Start_mp->point);
  }

  retVal = ZeroEG_amp2(pathNum, prec, time_first_increase, Final_d, Final_mp, last_approx_d, last_approx_mp, last_approx_prec, &lastSample_d, &lastSample_mp, &cycle_num, Start_d, Start_mp, prec_in, T, OUT, midOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

  // clear
  clear_point_data_d(&lastSample_d);
  clear_point_data_mp(&lastSample_mp);

  return retVal;
}

int ZeroEG_rank_amp(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks directly to endpoint using adaptive multi prec  *
\***************************************************************/
{
  int i, retVal, cycle_num, samples_per_loop = T->num_PSEG_sample_points + 1;
  comp_d finalTime_d;
  comp_mp finalTime_mp;
  point_data_d lastSample_d;
  point_data_mp lastSample_mp;

  init_mp(finalTime_mp);
  init_point_data_d(&lastSample_d, 0);
  init_point_data_mp(&lastSample_mp, 0);

  // setup finalTime
  set_double_d(finalTime_d, T->targetT, 0);
  d_to_mp(finalTime_mp, finalTime_d);

  // setup last_approx in case of failure
  *last_approx_prec = prec_in;
  if (prec_in < 64)
  {
    point_cp_d(last_approx_d, Start_d->point);
  }
  else
  {
    setprec_point_mp(last_approx_mp, prec_in);
    point_cp_mp(last_approx_mp, Start_mp->point);
  }

  // perform the endgame - return that last sample so that it can be used for rank determination
  retVal = ZeroEG_amp2(pathNum, prec, time_first_increase, Final_d, Final_mp, last_approx_d, last_approx_mp, last_approx_prec, &lastSample_d, &lastSample_mp, &cycle_num, Start_d, Start_mp, prec_in, T, OUT, midOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

  if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
  { // continue on with rank determination
    retVal = Cauchy_rank_main_amp(condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, retVal, prec, time_first_increase, Final_d, Final_mp, last_approx_d, last_approx_mp, *last_approx_prec, finalTime_d, finalTime_mp, &lastSample_d, &lastSample_mp, *prec, cycle_num, samples_per_loop, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
  }
  else
  { // setup condNum since we have failure
    *condNum = -1;
    *smallest_nonzero_SV = *largest_zero_SV = 0;

    // setup last_approx
    *last_approx_prec = *prec;
    if (*prec < 64)
    { // setup _d
      point_cp_d(last_approx_d, Final_d->point);
      for (i = 0; i < last_approx_d->size; i++)
      {
        get_comp_rand_d(finalTime_d);
        mul_rdouble_d(finalTime_d, finalTime_d, T->final_tolerance);
        add_d(&last_approx_d->coord[i], &last_approx_d->coord[i], finalTime_d);
      }
    }
    else
    { // setup _mp
      setprec_point_mp(last_approx_mp, *last_approx_prec);
      point_cp_mp(last_approx_mp, Final_mp->point);
      for (i = 0; i < last_approx_mp->size; i++)
      {
        get_comp_rand_mp(finalTime_mp);
        mul_rdouble_mp(finalTime_mp, finalTime_mp, T->final_tolerance);
        add_mp(&last_approx_mp->coord[i], &last_approx_mp->coord[i], finalTime_mp);
      }
    }
  }

  // clear
  clear_mp(finalTime_mp);
  clear_point_data_d(&lastSample_d);
  clear_point_data_mp(&lastSample_mp);

  return retVal;
}

int ZeroEG_amp2(int pathNum, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d *lastSample_d, point_data_mp *lastSample_mp, int *cycle_num, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: prec - the ending precision - < 64 - Final_d   *
*                                             >= 64 - Final_mp  *
* time_first_increase - the time of the last increase to MP     *
* NOTES: returns the last sample point                          *
\***************************************************************/
{
  int i, retVal;
  comp_d endTime_d;
  comp_mp endTime_mp;

  point_data_d Pt_d;
  point_data_mp Pt_mp;

  // start off with atleast 64 bit precision
  T->Precision = MAX(64, prec_in);
  initMP(T->Precision);
  change_prec(ED_mp, T->Precision);

  // initialize
  i = prec_in < 64 ? Start_d->point->size : Start_mp->point->size; 
  init_point_data_d(&Pt_d, i);
  init_point_data_mp(&Pt_mp, i); 
  init_mp(endTime_mp);

  // initialze prec & time_first_increase
  *prec = prec_in; // start off at the intial precision
  *time_first_increase = 0; // no increase so far!
  Final_d->cycle_num = Final_mp->cycle_num = 0;

  // setup to track to the endgame boundary
  set_double_d(endTime_d, T->endgameBoundary, 0.0);
  set_double_mp(endTime_mp, T->endgameBoundary, 0.0);
  T->currentStepSize = T->maxStepSize;
  T->currentNewtonTol = T->basicNewtonTol;
  T->minStepSize = T->minStepSizeBeforeEndGame;
  T->endgameSwitch = 0;
  T->error_at_latest_sample_point = 0.0;

  // setup last_approx in case of failure
  *last_approx_prec = prec_in;
  if (prec_in < 64)
  {
    point_cp_d(last_approx_d, Start_d->point);
  }
  else
  {
    setprec_point_mp(last_approx_mp, prec_in);
    point_cp_mp(last_approx_mp, Start_mp->point);
  }

  // try to track to the endgame boundary - this will handle everything
  retVal = AMPtrack(&Pt_d, &Pt_mp, prec, time_first_increase, Start_d, Start_mp, prec_in, endTime_d, endTime_mp, prec_in, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  // so, either made it to the endgame boundary or path failure - first print to midOUT
  fprintf(midOUT, "%d\n", pathNum);
  if (*prec == 52)
  { // print to midOUT using samples_d
    for (i = 0; i < Pt_d.point->size; i++)
      fprintf(midOUT, "%.15e %.15e\n", Pt_d.point->coord[i].r, Pt_d.point->coord[i].i);
    fprintf(midOUT, "\n");
  }
  else
  { // print to midOUT using samples_mp
    for (i = 0; i < Pt_mp.point->size; i++)
    {
      print_mp(midOUT, 0, &Pt_mp.point->coord[i]);
      fprintf(midOUT, "\n");
    }
    fprintf(midOUT, "\n");
  }
  // since we have printed to midOUT, this path is considered to be "in the endgame"
  T->endgameSwitch = 1;

  // check to see if we made it to the endgame boundary successfully
  if (retVal)
  { // had path failure
    fprintf(OUT, "NOTE: Zeroth endgame never started!\n");

    // copy the last sample to Final
    if (*prec == 52)
    {
      point_data_cp_d(Final_d, &Pt_d);
      Final_d->cycle_num = *cycle_num = 0;
      T->t_val_at_latest_sample_point = Pt_d.time->r;
    }
    else
    {
      setprec_point_mp(Final_mp->point, *prec);
      setprec_mp(Final_mp->time, *prec);
      point_data_cp_mp(Final_mp, &Pt_mp);
      Final_mp->cycle_num = *cycle_num = 0;
      T->t_val_at_latest_sample_point = mpf_get_d(Pt_mp.time->r);
    }
  }
  else
  { // success - setup to run the acutal endgame!

    // set the endgame tolerances
    set_double_d(endTime_d, T->targetT, 0.0);
    set_double_mp(endTime_mp, T->targetT, 0.0);
    T->currentNewtonTol = T->endgameNewtonTol;
    T->minStepSize = T->minStepSizeDuringEndGame;

    // track  
    retVal = AMPtrack(Final_d, Final_mp, prec, time_first_increase, &Pt_d, &Pt_mp, *prec, endTime_d, endTime_mp, *prec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // setup last_approx, lastSample & cycle_num
    *last_approx_prec = *prec; 
    if (*prec < 64)
    { // use _d
      *cycle_num = Final_d->cycle_num = 1;
      point_cp_d(last_approx_d, Final_d->point);
      point_data_cp_d(lastSample_d, Final_d);
      T->t_val_at_latest_sample_point = Final_d->time->r;
    }
    else
    { // use _mp
      *cycle_num = Final_mp->cycle_num = 1;
      setprec_point_mp(last_approx_mp, *last_approx_prec);
      point_cp_mp(last_approx_mp, Final_mp->point);

      // setup lastSample_mp
      setprec_point_mp(lastSample_mp->point, *prec);
      setprec_mp(lastSample_mp->time, *prec);
      point_data_cp_mp(lastSample_mp, Final_mp); 
      T->t_val_at_latest_sample_point = mpf_get_d(Final_mp->time->r);
    }

    if (retVal == 0)
    { // refine to precision requested in FinalTol
      int sharpenDigits = ceil(-log10(T->final_tolerance));
      retVal = sharpen_amp(T->outputLevel, sharpenDigits, &T->latest_newton_residual_d, T->latest_newton_residual_mp, T->currentNewtonTol, NULL, 52, Final_d, Final_mp, prec, Final_d, Final_mp, *prec, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
    }
  }

  // clear the MP data types
  clear_mp(endTime_mp);
  clear_point_data_d(&Pt_d);
  clear_point_data_mp(&Pt_mp);

  return retVal;
}

