// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "parallel.h"

#ifdef _HAVE_MPI

// creates the MPI datatype mpi_worker_info
void create_worker_info(MPI_Datatype *mpi_worker_info)
{
  // arrays for length, displacement and datatypes in mpi_worker_info
  int worker_info_length[1] = {1};
  MPI_Datatype worker_info_datatypes[1] = {MPI_INT};
  MPI_Aint worker_info_displacements[1] = {0};

  // build the mpi datatype mpi_comd_d
  MPI_Type_create_struct(1, worker_info_length, worker_info_displacements, worker_info_datatypes, mpi_worker_info);
  MPI_Type_commit(mpi_worker_info);

  return;
}

// creates the MPI datatype mpi_comp_d
void create_comp_d(MPI_Datatype *mpi_comp_d)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  comp_d tempComplexNum;

  // arrays for length, displacement and datatypes in mpi_comp_d
  int comp_d_length[2] = {1, 1};
  MPI_Datatype comp_d_datatypes[2] = {MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint comp_d_displacements[2];

  // calculate displacements
  comp_d_displacements[0] = 0;
  MPI_Get_address(&tempComplexNum->r, &startaddress);
  MPI_Get_address(&tempComplexNum->i, &address);
  comp_d_displacements[1] = address - startaddress;

  // build the mpi datatype mpi_comd_d
  MPI_Type_create_struct(2, comp_d_length, comp_d_displacements, comp_d_datatypes, mpi_comp_d);
  MPI_Type_commit(mpi_comp_d);

  return;
}

// creates the MPI datatype mpi_comp_rat_int
void create_comp_rat_int(MPI_Datatype *mpi_comp_rat_int)
{
  // arrays for length, displacement and datatypes in mpi_comp_rat_int
  int comp_rat_length[1] = {2};
  MPI_Datatype comp_rat_datatypes[1] = {MPI_INT};
  MPI_Aint comp_rat_displacements[1] = {0};

  // build the mpi datatype mpi_comd_rat_int
  MPI_Type_create_struct(1, comp_rat_length, comp_rat_displacements, comp_rat_datatypes, mpi_comp_rat_int);
  MPI_Type_commit(mpi_comp_rat_int);

  return;
}

// creates the MPI datatype mpi_point_d_int
void create_point_d_int(MPI_Datatype *mpi_point_d_int)
{
  // arrays for length, displacement and datatypes in mpi_point_d
  int point_d_length[1] = {1};
  MPI_Datatype point_d_datatypes[1] = {MPI_INT};
  MPI_Aint point_d_displacements[1] = {0};

  // build the mpi datatype mpi_point_d
  MPI_Type_create_struct(1, point_d_length, point_d_displacements, point_d_datatypes, mpi_point_d_int);
  MPI_Type_commit(mpi_point_d_int);

  return;
}

// creates the MPI datatype mpi_comp_d & mpi_point_d
void create_comp_point_d_int(MPI_Datatype *mpi_comp_d, MPI_Datatype *mpi_point_d_int)
{
  // create mpi_comp_d datatype
  create_comp_d(mpi_comp_d);
 
  // create mpi_point_d datatype
  create_point_d_int(mpi_point_d_int);

  return;
}

// creates the MPI datatype mpi_point_data_d_int
void create_point_data_d_int(MPI_Datatype *mpi_point_data_d_int)
{
  MPI_Datatype mpi_point_d_int, mpi_comp_d;

  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  point_data_d_int tempPointData;

  // create mpi_comp_d & mpi_point_d datatype
  create_comp_point_d_int(&mpi_comp_d, &mpi_point_d_int);

  // arrays for length, displacement and datatypes in mpi_point_data_d
  int point_data_d_length[3] = {1, 1, 1};
  MPI_Datatype point_data_d_datatypes[3] = {mpi_point_d_int, mpi_comp_d, MPI_INT};
  MPI_Aint point_data_d_displacements[3];

  // calculate displacements
  point_data_d_displacements[0] = 0;
  MPI_Get_address(&tempPointData.point_int, &startaddress);
  MPI_Get_address(&tempPointData.time, &address);
  point_data_d_displacements[1] = address - startaddress;
  MPI_Get_address(&tempPointData.cycle_num, &address);
  point_data_d_displacements[2] = address - startaddress;

  // build the mpi datatype mpi_point_data_d
  MPI_Type_create_struct(3, point_data_d_length, point_data_d_displacements, point_data_d_datatypes, mpi_point_data_d_int);
  MPI_Type_commit(mpi_point_data_d_int);

  // clear the datatypes
  MPI_Type_free(&mpi_comp_d);
  MPI_Type_free(&mpi_point_d_int);

  return;
}

// creates the MPI datatype mpi_mat_d_int
void create_mat_d_int(MPI_Datatype *mpi_mat_d_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  mat_d_int tempMat;

  // arrays for length, displacement and datatypes in mpi_mat_d_int
  int mat_d_int_length[2] = {1, 1};
  MPI_Datatype mat_d_int_datatypes[2] = {MPI_INT, MPI_INT};
  MPI_Aint mat_d_int_displacements[2];

  // calculate displacements
  mat_d_int_displacements[0] = 0;
  MPI_Get_address(&tempMat.rows, &startaddress);
  MPI_Get_address(&tempMat.cols, &address);
  mat_d_int_displacements[1] = address - startaddress;

  // build the mpi datatype mpi_mat_d_int
  MPI_Type_create_struct(2, mat_d_int_length, mat_d_int_displacements, mat_d_int_datatypes, mpi_mat_d_int);
  MPI_Type_commit(mpi_mat_d_int);

  return;
}

//////////////// MP VERSIONS OF MPI DATATYPES - integers signifying the lengths of the string /////////

// creates the MPI datatype mpi_mpf_int
void create_mpf_int(MPI_Datatype *mpi_mpf_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  mpf_int tempMPF;

  // arrays for length, displacement and datatypes in mpi_comp_mp_int
  int mpf_int_length[2] = {1, 1};
  MPI_Datatype mpf_int_datatypes[2] = {MPI_INT, MPI_INT};
  MPI_Aint mpf_int_displacements[2];

  // calculate displacements
  mpf_int_displacements[0] = 0;
  MPI_Get_address(&tempMPF.prec, &startaddress);
  MPI_Get_address(&tempMPF.totalLength, &address);
  mpf_int_displacements[1] = address - startaddress;

  // build the mpi datatype mpi_comp_mp_int
  MPI_Type_create_struct(2, mpf_int_length, mpf_int_displacements, mpf_int_datatypes, mpi_mpf_int);
  MPI_Type_commit(mpi_mpf_int);

  return;
}

// creates the MPI datatype mpi_comp_mp_int
void create_comp_mp_int(MPI_Datatype *mpi_comp_mp_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  comp_mp_int tempComp;

  // arrays for length, displacement and datatypes in mpi_comp_mp_int
  int comp_mp_int_length[2] = {1, 1};
  MPI_Datatype comp_mp_int_datatypes[2] = {MPI_INT, MPI_INT};
  MPI_Aint comp_mp_int_displacements[2];

  // calculate displacements
  comp_mp_int_displacements[0] = 0;
  MPI_Get_address(&tempComp.prec, &startaddress);
  MPI_Get_address(&tempComp.totalLength, &address);
  comp_mp_int_displacements[1] = address - startaddress;

  // build the mpi datatype mpi_comp_mp_int
  MPI_Type_create_struct(2, comp_mp_int_length, comp_mp_int_displacements, comp_mp_int_datatypes, mpi_comp_mp_int);
  MPI_Type_commit(mpi_comp_mp_int);

  return;
}

// creates the MPI datatype mpi_point_mp_int
void create_point_mp_int(MPI_Datatype *mpi_point_mp_int)
{
  MPI_Datatype mpi_comp;

  create_comp_point_mp_int(&mpi_comp, mpi_point_mp_int);

  MPI_Type_free(&mpi_comp);

  return;
}

// creates the MPI datatypes mpi_comp_mp_int & mpi_point_mp_int
void create_comp_point_mp_int(MPI_Datatype *mpi_comp_mp_int, MPI_Datatype *mpi_point_mp_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  point_mp_int tempPoint;

  // create mpi_comp_mp_int datatype
  create_comp_mp_int(mpi_comp_mp_int);

  // arrays for length, displacement and datatypes in mpi_point_mp_int
  int point_mp_int_length[3] = {1, 1, 1};
  MPI_Datatype point_mp_int_datatypes[3] = {MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint point_mp_int_displacements[3];

  // calculate displacements
  point_mp_int_displacements[0] = 0;
  MPI_Get_address(&tempPoint.prec, &startaddress);
  MPI_Get_address(&tempPoint.size, &address);
  point_mp_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempPoint.totalLength, &address);
  point_mp_int_displacements[2] = address - startaddress;

  // build the mpi datatype mpi_point_mp_int
  MPI_Type_create_struct(3, point_mp_int_length, point_mp_int_displacements, point_mp_int_datatypes, mpi_point_mp_int);
  MPI_Type_commit(mpi_point_mp_int);

  return;
}

// creates the MPI datatype mpi_point_data_mp_int
void create_point_data_mp_int(MPI_Datatype *mpi_point_data_mp_int)
{
  MPI_Datatype mpi_comp_mp_int, mpi_point_mp_int;

  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  point_data_mp_int tempPoint;

  // create datatypes mpi_comp_mp_int & mpi_point_mp_int
  create_comp_point_mp_int(&mpi_comp_mp_int, &mpi_point_mp_int);

  // arrays for length, displacement and datatypes in mpi_point_data_mp_int
  int point_data_mp_int_length[4] = {1, 1, 1, 1};
  MPI_Datatype point_data_mp_int_datatypes[4] = {mpi_point_mp_int, mpi_comp_mp_int, MPI_INT, MPI_INT};
  MPI_Aint point_data_mp_int_displacements[4];

  // calculate displacements
  point_data_mp_int_displacements[0] = 0;
  MPI_Get_address(&tempPoint.point_int, &startaddress);
  MPI_Get_address(&tempPoint.time_int, &address);
  point_data_mp_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempPoint.cycle_num, &address);
  point_data_mp_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempPoint.totalLength, &address);
  point_data_mp_int_displacements[3] = address - startaddress;

  // build the mpi datatype mpi_point_data_mp_int
  MPI_Type_create_struct(4, point_data_mp_int_length, point_data_mp_int_displacements, point_data_mp_int_datatypes, mpi_point_data_mp_int);
  MPI_Type_commit(mpi_point_data_mp_int);

  // clear the datatypes
  MPI_Type_free(&mpi_comp_mp_int);
  MPI_Type_free(&mpi_point_mp_int);

  return;
}

// creates the MPI datatype mpi_mat_mp_int
void create_mat_mp_int(MPI_Datatype *mpi_mat_mp_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  mat_mp_int tempMat;

  // arrays for length, displacement and datatypes in mpi_mat_mp_int
  int mat_mp_int_length[4] = {1, 1, 1, 1};
  MPI_Datatype mat_mp_int_datatypes[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint mat_mp_int_displacements[4];

  // calculate displacements
  mat_mp_int_displacements[0] = 0;
  MPI_Get_address(&tempMat.prec, &startaddress);
  MPI_Get_address(&tempMat.rows, &address);
  mat_mp_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempMat.cols, &address);
  mat_mp_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempMat.totalLength, &address);
  mat_mp_int_displacements[3] = address - startaddress;

  // build the mpi datatype mpi_mat_mp_int
  MPI_Type_create_struct(4, mat_mp_int_length, mat_mp_int_displacements, mat_mp_int_datatypes, mpi_mat_mp_int);
  MPI_Type_commit(mpi_mat_mp_int);

  return;
}

// creates the MPI datatype mpi_point_rat_int
void create_point_rat_int(MPI_Datatype *mpi_point_rat_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  point_rat_int tempVec;

  // arrays for length, displacement and datatypes in mpi_point_rat_int
  int point_rat_int_length[2] = {1, 1};
  MPI_Datatype point_rat_int_datatypes[2] = {MPI_INT, MPI_INT};
  MPI_Aint point_rat_int_displacements[2];

  // calculate displacements
  point_rat_int_displacements[0] = 0;
  MPI_Get_address(&tempVec.size, &startaddress);
  MPI_Get_address(&tempVec.totalLength, &address);
  point_rat_int_displacements[1] = address - startaddress;

  // build the mpi datatype mpi_point_rat_int
  MPI_Type_create_struct(2, point_rat_int_length, point_rat_int_displacements, point_rat_int_datatypes, mpi_point_rat_int);
  MPI_Type_commit(mpi_point_rat_int);

  return;
}

// creates the MPI datatype mpi_mat_rat_int
void create_mat_rat_int(MPI_Datatype *mpi_mat_rat_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  mat_rat_int tempMat;

  // arrays for length, displacement and datatypes in mpi_mat_rat_int
  int mat_rat_int_length[3] = {1, 1, 1};
  MPI_Datatype mat_rat_int_datatypes[3] = {MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint mat_rat_int_displacements[3];

  // calculate displacements
  mat_rat_int_displacements[0] = 0;
  MPI_Get_address(&tempMat.rows, &startaddress);
  MPI_Get_address(&tempMat.cols, &address);
  mat_rat_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempMat.totalLength, &address);
  mat_rat_int_displacements[2] = address - startaddress;

  // build the mpi datatype mpi_mat_rat_int
  MPI_Type_create_struct(3, mat_rat_int_length, mat_rat_int_displacements, mat_rat_int_datatypes, mpi_mat_rat_int);
  MPI_Type_commit(mpi_mat_rat_int);

  return;
}

////////////////// OTHER DATATYPES ///////////////////////

// creates the MPI datatype mpi_endgame_data_t
void create_endgame_data_t_int(MPI_Datatype *mpi_endgame_data_t)
{
  MPI_Datatype mpi_point_data_d_int, mpi_point_data_mp_int, mpi_point_d_int, mpi_point_mp_int, mpi_mpf_int;

  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  endgame_data_t_int tempEndgame;

  // create datatypes mpi_point_data_d_int, mpi_point_data_mp_int, mpi_point_d_int, mpi_point_mp_int, mpi_mpf_int
  create_point_data_d_int(&mpi_point_data_d_int);
  create_point_data_mp_int(&mpi_point_data_mp_int);
  create_point_d_int(&mpi_point_d_int);
  create_point_mp_int(&mpi_point_mp_int);
  create_mpf_int(&mpi_mpf_int);

  // arrays for length, displacement and datatypes in mpi_endgame_data_t
  int endgame_data_t_length[20] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype endgame_data_t_datatypes[20] = {MPI_INT, mpi_point_data_d_int, mpi_point_data_mp_int, MPI_INT, mpi_point_d_int, mpi_point_mp_int, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, mpi_mpf_int, MPI_DOUBLE, mpi_mpf_int, MPI_DOUBLE, mpi_mpf_int, MPI_DOUBLE, mpi_mpf_int, MPI_INT};
  MPI_Aint endgame_data_t_displacements[20];

  // calculate displacements
  endgame_data_t_displacements[0] = 0;
  MPI_Get_address(&tempEndgame.prec, &startaddress);
  MPI_Get_address(&tempEndgame.PD_d_int, &address);
  endgame_data_t_displacements[1] = address - startaddress;
  MPI_Get_address(&tempEndgame.PD_mp_int, &address);
  endgame_data_t_displacements[2] = address - startaddress;
  MPI_Get_address(&tempEndgame.last_approx_prec, &address);
  endgame_data_t_displacements[3] = address - startaddress;
  MPI_Get_address(&tempEndgame.last_approx_d_int, &address);
  endgame_data_t_displacements[4] = address - startaddress;
  MPI_Get_address(&tempEndgame.last_approx_mp_int, &address);
  endgame_data_t_displacements[5] = address - startaddress;
  MPI_Get_address(&tempEndgame.retVal, &address);
  endgame_data_t_displacements[6] = address - startaddress;
  MPI_Get_address(&tempEndgame.pathNum, &address);
  endgame_data_t_displacements[7] = address - startaddress;
  MPI_Get_address(&tempEndgame.codim, &address);
  endgame_data_t_displacements[8] = address - startaddress;
  MPI_Get_address(&tempEndgame.first_increase, &address);
  endgame_data_t_displacements[9] = address - startaddress;
  MPI_Get_address(&tempEndgame.condition_number, &address);
  endgame_data_t_displacements[10] = address - startaddress;
  MPI_Get_address(&tempEndgame.function_residual_d, &address);
  endgame_data_t_displacements[11] = address - startaddress;
  MPI_Get_address(&tempEndgame.function_residual_mp_int, &address);
  endgame_data_t_displacements[12] = address - startaddress;
  MPI_Get_address(&tempEndgame.latest_newton_residual_d, &address);
  endgame_data_t_displacements[13] = address - startaddress;
  MPI_Get_address(&tempEndgame.latest_newton_residual_mp_int, &address);
  endgame_data_t_displacements[14] = address - startaddress;
  MPI_Get_address(&tempEndgame.t_val_at_latest_sample_point_d, &address);
  endgame_data_t_displacements[15] = address - startaddress;
  MPI_Get_address(&tempEndgame.t_val_at_latest_sample_point_mp_int, &address);
  endgame_data_t_displacements[16] = address - startaddress;
  MPI_Get_address(&tempEndgame.error_at_latest_sample_point_d, &address);
  endgame_data_t_displacements[17] = address - startaddress;
  MPI_Get_address(&tempEndgame.error_at_latest_sample_point_mp_int, &address);
  endgame_data_t_displacements[18] = address - startaddress;
  MPI_Get_address(&tempEndgame.totalLength, &address);
  endgame_data_t_displacements[19] = address - startaddress;

  // build the mpi datatype mpi_endgame_data_t
  MPI_Type_create_struct(20, endgame_data_t_length, endgame_data_t_displacements, endgame_data_t_datatypes, mpi_endgame_data_t);
  MPI_Type_commit(mpi_endgame_data_t);

  // clear the datatypes
  MPI_Type_free(&mpi_point_data_d_int);
  MPI_Type_free(&mpi_point_data_mp_int);
  MPI_Type_free(&mpi_point_d_int);
  MPI_Type_free(&mpi_point_mp_int);
  MPI_Type_free(&mpi_mpf_int);

  return;
}

// creates the MPI datatype mpi_trackBack_samples_t_int
void create_trackBack_samples_t_int(MPI_Datatype *mpi_trackBack_samples_t)
{
  MPI_Datatype mpi_endgame_data_t_int;

  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  trackBack_samples_t_int tempTrack;

  // create datatype endgame_data_t_int
  create_endgame_data_t_int(&mpi_endgame_data_t_int);

  // arrays for length, displacement and datatypes in mpi_endgame_data_t
  int trackBack_data_t_length[8] = {1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype trackBack_data_t_datatypes[8] = {mpi_endgame_data_t_int, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint trackBack_data_t_displacements[8];

  // calculate displacements
  trackBack_data_t_displacements[0] = 0;
  MPI_Get_address(&tempTrack.EG_int, &startaddress);
  MPI_Get_address(&tempTrack.numSamples, &address);
  trackBack_data_t_displacements[1] = address - startaddress;
  MPI_Get_address(&tempTrack.samplePts_prec, &address);
  trackBack_data_t_displacements[2] = address - startaddress;
  MPI_Get_address(&tempTrack.midPt_prec, &address);
  trackBack_data_t_displacements[3] = address - startaddress;
  MPI_Get_address(&tempTrack.pointSize, &address);
  trackBack_data_t_displacements[4] = address - startaddress;
  MPI_Get_address(&tempTrack.num_double, &address);
  trackBack_data_t_displacements[5] = address - startaddress;
  MPI_Get_address(&tempTrack.num_comp_d, &address);
  trackBack_data_t_displacements[6] = address - startaddress;
  MPI_Get_address(&tempTrack.totalLength, &address);
  trackBack_data_t_displacements[7] = address - startaddress;

  // build the mpi datatype mpi_trackBack_samples_t
  MPI_Type_create_struct(8, trackBack_data_t_length, trackBack_data_t_displacements, trackBack_data_t_datatypes, mpi_trackBack_samples_t);
  MPI_Type_commit(mpi_trackBack_samples_t);

  // clear the datatypes
  MPI_Type_free(&mpi_endgame_data_t_int);

  return;
}

// creates the MPI datatype mpi_tracker_config_t_relevant - send over MPI the relevant data
void create_tracker_config_t(MPI_Datatype *mpi_tracker_config_t_relevant)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  tracker_config_t_relevant tempTrackerData;

  int counter;

  // arrays for length, displacement and datatypes in mpi_tracker_config_t_relevant
  int tracker_config_t_length[63];
  MPI_Datatype tracker_config_t_datatypes[63] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
  MPI_Aint tracker_config_t_displacements[63];

  // set array of lengths
  for (counter = 0; counter < 63; counter++)
    tracker_config_t_length[counter] = 1;

  // calculate displacements
  tracker_config_t_displacements[0] = 0;
  MPI_Get_address(&tempTrackerData.numVars, &startaddress);
  MPI_Get_address(&tempTrackerData.numPathVars, &address);
  tracker_config_t_displacements[1] = address - startaddress;
  MPI_Get_address(&tempTrackerData.numParams, &address);
  tracker_config_t_displacements[2] = address - startaddress;
  MPI_Get_address(&tempTrackerData.numFuncs, &address);
  tracker_config_t_displacements[3] = address - startaddress;

  MPI_Get_address(&tempTrackerData.maxStepSize, &address);
  tracker_config_t_displacements[4] = address - startaddress;
  MPI_Get_address(&tempTrackerData.minStepSizeBeforeEndGame, &address);
  tracker_config_t_displacements[5] = address - startaddress;
  MPI_Get_address(&tempTrackerData.minStepSizeDuringEndGame, &address);
  tracker_config_t_displacements[6] = address - startaddress;

  MPI_Get_address(&tempTrackerData.minTrackT, &address);
  tracker_config_t_displacements[7] = address - startaddress;
  MPI_Get_address(&tempTrackerData.basicNewtonTol, &address);
  tracker_config_t_displacements[8] = address - startaddress;
  MPI_Get_address(&tempTrackerData.endgameNewtonTol, &address);
  tracker_config_t_displacements[9] = address - startaddress;

  MPI_Get_address(&tempTrackerData.cSecInc, &address);
  tracker_config_t_displacements[10] = address - startaddress;
  MPI_Get_address(&tempTrackerData.maxNewtonIts, &address);
  tracker_config_t_displacements[11] = address - startaddress;
  MPI_Get_address(&tempTrackerData.MPType, &address);
  tracker_config_t_displacements[12] = address - startaddress;
  MPI_Get_address(&tempTrackerData.Precision, &address);
  tracker_config_t_displacements[13] = address - startaddress;
  MPI_Get_address(&tempTrackerData.outputLevel, &address);
  tracker_config_t_displacements[14] = address - startaddress;
  MPI_Get_address(&tempTrackerData.screenOut, &address);
  tracker_config_t_displacements[15] = address - startaddress;
  MPI_Get_address(&tempTrackerData.targetT, &address);
  tracker_config_t_displacements[16] = address - startaddress;
  MPI_Get_address(&tempTrackerData.endgameBoundary, &address);
  tracker_config_t_displacements[17] = address - startaddress;

  MPI_Get_address(&tempTrackerData.goingToInfinity, &address);
  tracker_config_t_displacements[18] = address - startaddress;
  MPI_Get_address(&tempTrackerData.maxNumSteps, &address);
  tracker_config_t_displacements[19] = address - startaddress;
  MPI_Get_address(&tempTrackerData.endgameNumber, &address);
  tracker_config_t_displacements[20] = address - startaddress;

  MPI_Get_address(&tempTrackerData.power_series_sample_factor, &address);
  tracker_config_t_displacements[21] = address - startaddress;
  MPI_Get_address(&tempTrackerData.cycle_num_max, &address);
  tracker_config_t_displacements[22] = address - startaddress;
  MPI_Get_address(&tempTrackerData.num_PSEG_sample_points, &address);
  tracker_config_t_displacements[23] = address - startaddress;

  MPI_Get_address(&tempTrackerData.final_tolerance, &address);
  tracker_config_t_displacements[24] = address - startaddress;
  MPI_Get_address(&tempTrackerData.real_threshold, &address);
  tracker_config_t_displacements[25] = address - startaddress;
  MPI_Get_address(&tempTrackerData.endgameOnly, &address);
  tracker_config_t_displacements[26] = address - startaddress;

  MPI_Get_address(&tempTrackerData.AMP_bound_on_abs_vals_of_coeffs, &address);
  tracker_config_t_displacements[27] = address - startaddress;
  MPI_Get_address(&tempTrackerData.AMP_bound_on_degree, &address);
  tracker_config_t_displacements[28] = address - startaddress;
  MPI_Get_address(&tempTrackerData.AMP_eps, &address);
  tracker_config_t_displacements[29] = address - startaddress;
  MPI_Get_address(&tempTrackerData.AMP_Phi, &address);
  tracker_config_t_displacements[30] = address - startaddress;
  MPI_Get_address(&tempTrackerData.AMP_Psi, &address);
  tracker_config_t_displacements[31] = address - startaddress;
  MPI_Get_address(&tempTrackerData.AMP_safety_digits_1, &address);
  tracker_config_t_displacements[32] = address - startaddress;
  MPI_Get_address(&tempTrackerData.AMP_safety_digits_2, &address);
  tracker_config_t_displacements[33] = address - startaddress;
  MPI_Get_address(&tempTrackerData.AMP_max_prec, &address);
  tracker_config_t_displacements[34] = address - startaddress;

  MPI_Get_address(&tempTrackerData.sing_val_zero_tol, &address);
  tracker_config_t_displacements[35] = address - startaddress;
  MPI_Get_address(&tempTrackerData.cond_num_threshold, &address);
  tracker_config_t_displacements[36] = address - startaddress;

  MPI_Get_address(&tempTrackerData.step_fail_factor, &address);
  tracker_config_t_displacements[37] = address - startaddress;   
  MPI_Get_address(&tempTrackerData.step_success_factor, &address);
  tracker_config_t_displacements[38] = address - startaddress; 

  MPI_Get_address(&tempTrackerData.max_num_pts_for_trace, &address);
  tracker_config_t_displacements[39] = address - startaddress;
  MPI_Get_address(&tempTrackerData.max_num_mon_linears, &address);
  tracker_config_t_displacements[40] = address - startaddress;
  MPI_Get_address(&tempTrackerData.max_num_bad_loops_in_mon, &address);
  tracker_config_t_displacements[41] = address - startaddress;

  MPI_Get_address(&tempTrackerData.final_tol_multiplier, &address);
  tracker_config_t_displacements[42] = address - startaddress; 
  MPI_Get_address(&tempTrackerData.final_tol_times_mult, &address);
  tracker_config_t_displacements[43] = address - startaddress;

  MPI_Get_address(&tempTrackerData.sharpenOnly, &address);
  tracker_config_t_displacements[44] = address - startaddress; 
  MPI_Get_address(&tempTrackerData.sharpenDigits, &address);
  tracker_config_t_displacements[45] = address - startaddress;

  MPI_Get_address(&tempTrackerData.regen_remove_inf, &address);
  tracker_config_t_displacements[46] = address - startaddress;
  MPI_Get_address(&tempTrackerData.regen_higher_dim_check, &address);
  tracker_config_t_displacements[47] = address - startaddress;
  MPI_Get_address(&tempTrackerData.sliceBasicNewtonTol, &address);
  tracker_config_t_displacements[48] = address - startaddress;
  MPI_Get_address(&tempTrackerData.sliceEndgameNewtonTol, &address);
  tracker_config_t_displacements[49] = address - startaddress;
  MPI_Get_address(&tempTrackerData.sliceFinalTol, &address);
  tracker_config_t_displacements[50] = address - startaddress;

  MPI_Get_address(&tempTrackerData.minCycleTrackBack, &address);
  tracker_config_t_displacements[51] = address - startaddress;
  MPI_Get_address(&tempTrackerData.junkRemovalTest, &address);
  tracker_config_t_displacements[52] = address - startaddress;
  MPI_Get_address(&tempTrackerData.maxDepthLDT, &address);
  tracker_config_t_displacements[53] = address - startaddress;
  MPI_Get_address(&tempTrackerData.odePredictor, &address);
  tracker_config_t_displacements[54] = address - startaddress;
  MPI_Get_address(&tempTrackerData.securityLevel, &address);
  tracker_config_t_displacements[55] = address - startaddress;
  MPI_Get_address(&tempTrackerData.securityMaxNorm, &address);
  tracker_config_t_displacements[56] = address - startaddress;
  MPI_Get_address(&tempTrackerData.cutoffCycleTime, &address);
  tracker_config_t_displacements[57] = address - startaddress;
  MPI_Get_address(&tempTrackerData.cutoffRatioTime, &address);
  tracker_config_t_displacements[58] = address - startaddress;
  MPI_Get_address(&tempTrackerData.finiteThreshold, &address);
  tracker_config_t_displacements[59] = address - startaddress;
  MPI_Get_address(&tempTrackerData.funcResTol, &address);
  tracker_config_t_displacements[60] = address - startaddress;
  MPI_Get_address(&tempTrackerData.ratioTol, &address);
  tracker_config_t_displacements[61] = address - startaddress;
  MPI_Get_address(&tempTrackerData.maxStepsBeforeNewton, &address);
  tracker_config_t_displacements[62] = address - startaddress;

  // build the mpi datatype mpi_tracker_config_t_relevant
  MPI_Type_create_struct(63, tracker_config_t_length, tracker_config_t_displacements, tracker_config_t_datatypes, mpi_tracker_config_t_relevant);
  MPI_Type_commit(mpi_tracker_config_t_relevant);

  return;
}

// creates the MPI datatype mpi_prog_t_int
void create_prog_t_int(MPI_Datatype *mpi_prog_t_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  prog_t_int tempProg;

  int counter;

  // arrays for length, displacement and datatypes in mpi_prog_t_int
  int prog_t_int_length[31];
  MPI_Datatype prog_t_int_datatypes[31];
  MPI_Aint prog_t_int_displacements[31];

  for (counter = 0; counter < 31; counter++)
  {
    prog_t_int_length[counter] = 1;
    prog_t_int_datatypes[counter] = MPI_INT;
  }

  // calculate displacements
  prog_t_int_displacements[0] = 0;
  prog_t_int_datatypes[0] = MPI_LONG;
  MPI_Get_address(&tempProg.size, &startaddress);
  MPI_Get_address(&tempProg.memSize, &address);
  prog_t_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempProg.precision, &address);
  prog_t_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempProg.num_var_gps, &address);
  prog_t_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempProg.index_of_first_number_for_proj_trans, &address);
  prog_t_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempProg.numInstAtEndUpdate, &address);
  prog_t_int_displacements[5] = address - startaddress;
  prog_t_int_datatypes[5] = MPI_LONG;
  MPI_Get_address(&tempProg.numInstAtEndParams, &address);
  prog_t_int_displacements[6] = address - startaddress;
  prog_t_int_datatypes[6] = MPI_LONG;
  MPI_Get_address(&tempProg.numInstAtEndPDeriv, &address);
  prog_t_int_displacements[7] = address - startaddress;
  prog_t_int_datatypes[7] = MPI_LONG;
  MPI_Get_address(&tempProg.numInstAtEndFnEval, &address);
  prog_t_int_displacements[8] = address - startaddress;
  prog_t_int_datatypes[8] = MPI_LONG;
  MPI_Get_address(&tempProg.numInstAtEndJvEval, &address);
  prog_t_int_displacements[9] = address - startaddress;
  prog_t_int_datatypes[9] = MPI_LONG;
  MPI_Get_address(&tempProg.numVars, &address);
  prog_t_int_displacements[10] = address - startaddress;
  MPI_Get_address(&tempProg.numPathVars, &address);
  prog_t_int_displacements[11] = address - startaddress;
  MPI_Get_address(&tempProg.numNums, &address);
  prog_t_int_displacements[12] = address - startaddress;
  MPI_Get_address(&tempProg.numConsts, &address);
  prog_t_int_displacements[13] = address - startaddress;
  MPI_Get_address(&tempProg.numPars, &address);
  prog_t_int_displacements[14] = address - startaddress;
  MPI_Get_address(&tempProg.numFuncs, &address);
  prog_t_int_displacements[15] = address - startaddress;
  MPI_Get_address(&tempProg.numSubfuncs, &address);
  prog_t_int_displacements[16] = address - startaddress;
  MPI_Get_address(&tempProg.inpVars, &address);
  prog_t_int_displacements[17] = address - startaddress;
  MPI_Get_address(&tempProg.inpPathVars, &address);
  prog_t_int_displacements[18] = address - startaddress;
  MPI_Get_address(&tempProg.IAddr, &address);
  prog_t_int_displacements[19] = address - startaddress;
  MPI_Get_address(&tempProg.numAddr, &address);
  prog_t_int_displacements[20] = address - startaddress;
  MPI_Get_address(&tempProg.constAddr, &address);
  prog_t_int_displacements[21] = address - startaddress;
  MPI_Get_address(&tempProg.evalPars, &address);
  prog_t_int_displacements[22] = address - startaddress;
  MPI_Get_address(&tempProg.evalDPars, &address);
  prog_t_int_displacements[23] = address - startaddress;
  MPI_Get_address(&tempProg.evalFuncs, &address);
  prog_t_int_displacements[24] = address - startaddress;
  MPI_Get_address(&tempProg.evalJVars, &address);
  prog_t_int_displacements[25] = address - startaddress;
  MPI_Get_address(&tempProg.evalJPars, &address);
  prog_t_int_displacements[26] = address - startaddress;
  MPI_Get_address(&tempProg.evalSubs, &address);
  prog_t_int_displacements[27] = address - startaddress;
  MPI_Get_address(&tempProg.evalJSubsV, &address);
  prog_t_int_displacements[28] = address - startaddress;
  MPI_Get_address(&tempProg.evalJSubsP, &address);
  prog_t_int_displacements[29] = address - startaddress;
  MPI_Get_address(&tempProg.totalLength, &address);
  prog_t_int_displacements[30] = address - startaddress;

  // build the mpi datatype mpi_prog_t_int
  MPI_Type_create_struct(31, prog_t_int_length, prog_t_int_displacements, prog_t_int_datatypes, mpi_prog_t_int);
  MPI_Type_commit(mpi_prog_t_int);

  return;
}

///////////////////// EVALUATION DATA STRUCTURES ////////////////////////

// creates the MPI datatype mpi_preproc_data_int
void create_preproc_data_int(MPI_Datatype *mpi_preproc_data_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  preproc_data_int tempPPD;

  // arrays for length, displacement and datatypes in mpi_preproc_data_int
  int preproc_data_int_length[3] = {1, 1, 1};
  MPI_Datatype preproc_data_int_datatypes[3] = {MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint preproc_data_int_displacements[3];

  // calculate displacements
  preproc_data_int_displacements[0] = 0;
  MPI_Get_address(&tempPPD.num_funcs, &startaddress);
  MPI_Get_address(&tempPPD.num_hom_var_gp, &address);
  preproc_data_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempPPD.num_var_gp, &address);
  preproc_data_int_displacements[2] = address - startaddress;

  // build the mpi datatype mpi_preproc_data_int
  MPI_Type_create_struct(3, preproc_data_int_length, preproc_data_int_displacements, preproc_data_int_datatypes, mpi_preproc_data_int);
  MPI_Type_commit(mpi_preproc_data_int);

  return;
}

// creates the MPI datatype mpi_patch_eval_data_d_int
void create_patch_eval_data_d_int(MPI_Datatype *mpi_patch_eval_data_d_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  patch_eval_data_d_int tempPED;

  // arrays for length, displacement and datatypes in mpi_patch_eval_data_d
  int patch_eval_data_d_length[3] = {1, 1, 1};
  MPI_Datatype patch_eval_data_d_datatypes[3] = {MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint patch_eval_data_d_displacements[3];

  // calculate displacements
  patch_eval_data_d_displacements[0] = 0;
  MPI_Get_address(&tempPED.num_patches, &startaddress);
  MPI_Get_address(&tempPED.patchCoeff_rows, &address);
  patch_eval_data_d_displacements[1] = address - startaddress;
  MPI_Get_address(&tempPED.patchCoeff_cols, &address);
  patch_eval_data_d_displacements[2] = address - startaddress;

  // build the mpi datatype mpi_patch_eval_data_d
  MPI_Type_create_struct(3, patch_eval_data_d_length, patch_eval_data_d_displacements, patch_eval_data_d_datatypes, mpi_patch_eval_data_d_int);
  MPI_Type_commit(mpi_patch_eval_data_d_int);

  return;
}

// creates the MPI datatype mpi_patch_eval_data_mp_int
void create_patch_eval_data_mp_int(MPI_Datatype *mpi_patch_eval_data_mp_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  patch_eval_data_mp_int tempPED;

  // arrays for length, displacement and datatypes in mpi_patch_eval_data_mp_int
  int patch_eval_data_mp_int_length[5] = {1, 1, 1, 1, 1};
  MPI_Datatype patch_eval_data_mp_int_datatypes[5] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint patch_eval_data_mp_int_displacements[5];

  // calculate displacements
  patch_eval_data_mp_int_displacements[0] = 0;
  MPI_Get_address(&tempPED.prec, &startaddress);
  MPI_Get_address(&tempPED.num_patches, &address);
  patch_eval_data_mp_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempPED.patchCoeff_rows, &address);
  patch_eval_data_mp_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempPED.patchCoeff_cols, &address);
  patch_eval_data_mp_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempPED.totalLength, &address);
  patch_eval_data_mp_int_displacements[4] = address - startaddress;

  // build the mpi datatype mpi_patch_eval_data_mp_int
  MPI_Type_create_struct(5, patch_eval_data_mp_int_length, patch_eval_data_mp_int_displacements, patch_eval_data_mp_int_datatypes, mpi_patch_eval_data_mp_int);
  MPI_Type_commit(mpi_patch_eval_data_mp_int);

  return;
}

// creates the MPI datatype mpi_start_system_eval_data_d_int
void create_start_system_eval_data_d_int(MPI_Datatype *mpi_ssed_d_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  MPI_Datatype mpi_comp_d;

  start_system_eval_data_d_int tempSSED;

  // create mpi_comp_d
  create_comp_d(&mpi_comp_d);

  // arrays for length, displacement and datatypes in mpi_ssed_d_int
  int ssed_d_int_length[6] = {1, 1, 1, 1, 1, 1};
  MPI_Datatype ssed_d_int_datatypes[6] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, mpi_comp_d};
  MPI_Aint ssed_d_int_displacements[6];

  // calculate displacements
  ssed_d_int_displacements[0] = 0;
  MPI_Get_address(&tempSSED.startSystemType, &startaddress);
  MPI_Get_address(&tempSSED.size_r, &address);
  ssed_d_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempSSED.max_degree, &address);
  ssed_d_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempSSED.coeff_cols, &address);
  ssed_d_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempSSED.num_coeff, &address);
  ssed_d_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempSSED.gamma, &address);
  ssed_d_int_displacements[5] = address - startaddress;

  // build the mpi datatype mpi_ssed_d_int
  MPI_Type_create_struct(6, ssed_d_int_length, ssed_d_int_displacements, ssed_d_int_datatypes, mpi_ssed_d_int);
  MPI_Type_commit(mpi_ssed_d_int);

  // free the datatypes
  MPI_Type_free(&mpi_comp_d);

  return;
}

// creates the MPI datatype mpi_start_system_eval_data_mp_int
void create_start_system_eval_data_mp_int(MPI_Datatype *mpi_ssed_mp_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  start_system_eval_data_mp_int tempSSED;

  // arrays for length, displacement and datatypes in mpi_ssed_mp_int
  int ssed_mp_int_length[7] = {1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype ssed_mp_int_datatypes[7] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint ssed_mp_int_displacements[7];

  // calculate displacements
  ssed_mp_int_displacements[0] = 0;
  MPI_Get_address(&tempSSED.startSystemType, &startaddress);
  MPI_Get_address(&tempSSED.size_r, &address);
  ssed_mp_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempSSED.max_degree, &address);
  ssed_mp_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempSSED.coeff_cols, &address);
  ssed_mp_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempSSED.coeff_strLength, &address);
  ssed_mp_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempSSED.prec, &address);
  ssed_mp_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempSSED.totalLength, &address);
  ssed_mp_int_displacements[6] = address - startaddress;

  // build the mpi datatype mpi_ssed_mp_int
  MPI_Type_create_struct(7, ssed_mp_int_length, ssed_mp_int_displacements, ssed_mp_int_datatypes, mpi_ssed_mp_int);
  MPI_Type_commit(mpi_ssed_mp_int);

  return;
}

// creates the MPI datatype mpi_square_system_eval_data_d_int
void create_square_system_eval_data_d_int(MPI_Datatype *mpi_ssed_d_int)
{
  MPI_Datatype mpi_prog_t_int;

  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  square_system_eval_data_d_int tempSSED;

  // create datatype mpi_prog_t_int
  create_prog_t_int(&mpi_prog_t_int);

  // arrays for length, displacement and datatypes in mpi_ssed_d_int
  int ssed_d_int_length[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype ssed_d_int_datatypes[12] = {mpi_prog_t_int, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint ssed_d_int_displacements[12];

  // calculate displacements
  ssed_d_int_displacements[0] = 0;
  MPI_Get_address(&tempSSED.Prog_int, &startaddress);
  MPI_Get_address(&tempSSED.size_f, &address);
  ssed_d_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempSSED.B_rows, &address);
  ssed_d_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempSSED.B_cols, &address);
  ssed_d_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempSSED.B_perp_rows, &address);
  ssed_d_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempSSED.B_perp_cols, &address);
  ssed_d_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempSSED.noChanges, &address);
  ssed_d_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempSSED.max_of_W, &address);
  ssed_d_int_displacements[7] = address - startaddress;
  MPI_Get_address(&tempSSED.A_rows, &address);
  ssed_d_int_displacements[8] = address - startaddress;
  MPI_Get_address(&tempSSED.A_cols, &address);
  ssed_d_int_displacements[9] = address - startaddress;
  MPI_Get_address(&tempSSED.size_r, &address);
  ssed_d_int_displacements[10] = address - startaddress;
  MPI_Get_address(&tempSSED.num_comp_d, &address);
  ssed_d_int_displacements[11] = address - startaddress;

  // build the mpi datatype mpi_ssed_d_int
  MPI_Type_create_struct(12, ssed_d_int_length, ssed_d_int_displacements, ssed_d_int_datatypes, mpi_ssed_d_int);
  MPI_Type_commit(mpi_ssed_d_int);

  // free datatypes
  MPI_Type_free(&mpi_prog_t_int);

  return;
}

// creates the MPI datatype mpi_square_system_eval_data_mp_int
void create_square_system_eval_data_mp_int(MPI_Datatype *mpi_ssed_mp_int)
{
  MPI_Datatype mpi_prog_t_int;

  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  square_system_eval_data_mp_int tempSSED;

  // create datatypes mpi_prog_t_int
  create_prog_t_int(&mpi_prog_t_int);

  // arrays for length, displacement and datatypes in mpi_ssed_mp_int
  int ssed_mp_int_length[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype ssed_mp_int_datatypes[16] = {mpi_prog_t_int, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint ssed_mp_int_displacements[16];

  // calculate displacements
  ssed_mp_int_displacements[0] = 0;
  MPI_Get_address(&tempSSED.Prog_int, &startaddress);
  MPI_Get_address(&tempSSED.size_f, &address);
  ssed_mp_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempSSED.B_rows, &address);
  ssed_mp_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempSSED.B_cols, &address);
  ssed_mp_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempSSED.B_strLength, &address);
  ssed_mp_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempSSED.B_perp_rows, &address);
  ssed_mp_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempSSED.B_perp_cols, &address); 
  ssed_mp_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempSSED.B_perp_strLength, &address);
  ssed_mp_int_displacements[7] = address - startaddress;
  MPI_Get_address(&tempSSED.noChanges, &address);
  ssed_mp_int_displacements[8] = address - startaddress;
  MPI_Get_address(&tempSSED.max_of_W, &address);
  ssed_mp_int_displacements[9] = address - startaddress;
  MPI_Get_address(&tempSSED.A_rows, &address);
  ssed_mp_int_displacements[10] = address - startaddress;
  MPI_Get_address(&tempSSED.A_cols, &address); 
  ssed_mp_int_displacements[11] = address - startaddress;
  MPI_Get_address(&tempSSED.A_strLength, &address);
  ssed_mp_int_displacements[12] = address - startaddress;
  MPI_Get_address(&tempSSED.size_r, &address);
  ssed_mp_int_displacements[13] = address - startaddress;
  MPI_Get_address(&tempSSED.prec, &address);
  ssed_mp_int_displacements[14] = address - startaddress;
  MPI_Get_address(&tempSSED.totalLength, &address);
  ssed_mp_int_displacements[15] = address - startaddress;

  // build the mpi datatype mpi_ssed_mp_int
  MPI_Type_create_struct(16, ssed_mp_int_length, ssed_mp_int_displacements, ssed_mp_int_datatypes, mpi_ssed_mp_int);
  MPI_Type_commit(mpi_ssed_mp_int);

  // free mpi_prog_t_int
  MPI_Type_free(&mpi_prog_t_int);

  return;
}

// creates the MPI datatype mpi_basic_eval_data_d_int
void create_basic_eval_data_d_int(MPI_Datatype *mpi_bed_d_int)
{
  MPI_Datatype mpi_sq_d_int, mpi_patch_d_int, mpi_st_d_int, mpi_preproc_int;

  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  basic_eval_data_d_int tempBED;

  // create datatypes needed
  create_square_system_eval_data_d_int(&mpi_sq_d_int);
  create_patch_eval_data_d_int(&mpi_patch_d_int);
  create_start_system_eval_data_d_int(&mpi_st_d_int);
  create_preproc_data_int(&mpi_preproc_int);

  // arrays for length, displacement and datatypes in mpi_bed_mp_int
  int bed_d_int_length[5] = {1, 1, 1, 1, 1};
  MPI_Datatype bed_d_int_datatypes[5] = {mpi_sq_d_int, mpi_patch_d_int, mpi_st_d_int, mpi_preproc_int, MPI_INT};
  MPI_Aint bed_d_int_displacements[5];

  // calculate displacements
  bed_d_int_displacements[0] = 0;
  MPI_Get_address(&tempBED.squareSystem_int, &startaddress);
  MPI_Get_address(&tempBED.patch, &address);
  bed_d_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempBED.startSystem_int, &address);
  bed_d_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempBED.preProcData_int, &address);
  bed_d_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempBED.MPType, &address);
  bed_d_int_displacements[4] = address - startaddress;

  // build the mpi datatype mpi_bed_mp_int
  MPI_Type_create_struct(5, bed_d_int_length, bed_d_int_displacements, bed_d_int_datatypes, mpi_bed_d_int);
  MPI_Type_commit(mpi_bed_d_int);

  // free the datatypes
  MPI_Type_free(&mpi_sq_d_int);
  MPI_Type_free(&mpi_patch_d_int);
  MPI_Type_free(&mpi_st_d_int);
  MPI_Type_free(&mpi_preproc_int);

  return;
}

// creates the MPI datatype mpi_basic_eval_data_mp_int
void create_basic_eval_data_mp_int(MPI_Datatype *mpi_bed_mp_int)
{
  MPI_Datatype mpi_sq_mp_int, mpi_patch_mp_int, mpi_st_mp_int, mpi_preproc_int;

  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  basic_eval_data_mp_int tempBED;

  // create datatypes needed
  create_square_system_eval_data_mp_int(&mpi_sq_mp_int);
  create_patch_eval_data_mp_int(&mpi_patch_mp_int);
  create_start_system_eval_data_mp_int(&mpi_st_mp_int);
  create_preproc_data_int(&mpi_preproc_int);

  // arrays for length, displacement and datatypes in mpi_bed_mp_int
  int bed_mp_int_length[4] = {1, 1, 1, 1};
  MPI_Datatype bed_mp_int_datatypes[4] = {mpi_sq_mp_int, mpi_patch_mp_int, mpi_st_mp_int, mpi_preproc_int};
  MPI_Aint bed_mp_int_displacements[4];

  // calculate displacements
  bed_mp_int_displacements[0] = 0;
  MPI_Get_address(&tempBED.squareSystem_int, &startaddress);
  MPI_Get_address(&tempBED.patch_int, &address);
  bed_mp_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempBED.startSystem_int, &address);
  bed_mp_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempBED.preProcData_int, &address);
  bed_mp_int_displacements[3] = address - startaddress;

  // build the mpi datatype mpi_bed_mp_int
  MPI_Type_create_struct(4, bed_mp_int_length, bed_mp_int_displacements, bed_mp_int_datatypes, mpi_bed_mp_int);
  MPI_Type_commit(mpi_bed_mp_int);

  // free the datatypes
  MPI_Type_free(&mpi_sq_mp_int);
  MPI_Type_free(&mpi_patch_mp_int);
  MPI_Type_free(&mpi_st_mp_int);
  MPI_Type_free(&mpi_preproc_int);

  return;
}

// creates the MPI datatype mpi_eqData_t_int
void create_eqData_t_int(MPI_Datatype *mpi_eqd_t_int)
{
  MPI_Datatype mpi_comp_d;

  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  eqData_t_int tempEqD;

  // create datatypes mpi_comp_d
  create_comp_d(&mpi_comp_d);

  // arrays for length, displacement and datatypes in mpi_eqd_d_int
  int eqd_t_int_length[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype eqd_t_int_datatypes[11] = {mpi_comp_d, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint eqd_t_int_displacements[11];

  // calculate displacements
  eqd_t_int_displacements[0] = 0;
  MPI_Get_address(&tempEqD.gamma_d, &startaddress);
  MPI_Get_address(&tempEqD.curr_precision, &address);
  eqd_t_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempEqD.num_funcs, &address);
  eqd_t_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempEqD.num_subsystems, &address);
  eqd_t_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempEqD.num_var_gps, &address);
  eqd_t_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempEqD.num_vars, &address);
  eqd_t_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempEqD.num_coeff, &address);
  eqd_t_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempEqD.totalLength, &address);
  eqd_t_int_displacements[7] = address - startaddress;
  MPI_Get_address(&tempEqD.noChanges, &address);
  eqd_t_int_displacements[8] = address - startaddress;
  MPI_Get_address(&tempEqD.numSubFuncs, &address);
  eqd_t_int_displacements[9] = address - startaddress;
  MPI_Get_address(&tempEqD.numInts, &address);
  eqd_t_int_displacements[10] = address - startaddress;

  // build the mpi datatype mpi_eqd_t_int
  MPI_Type_create_struct(11, eqd_t_int_length, eqd_t_int_displacements, eqd_t_int_datatypes, mpi_eqd_t_int);
  MPI_Type_commit(mpi_eqd_t_int);

  // free datatypes
  MPI_Type_free(&mpi_comp_d);

  return;
}

// creates the MPI datatype mpi_witnessData_t_int
void create_witnessData_t_int(MPI_Datatype *mpi_wd_t_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  witnessData_t_int tempWD;

  // arrays for length, displacement and datatypes in mpi_wd_d_int
  int wd_t_int_length[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype wd_t_int_datatypes[10] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint wd_t_int_displacements[10];

  // calculate displacements
  wd_t_int_displacements[0] = 0;
  MPI_Get_address(&tempWD.curr_precision, &startaddress);
  MPI_Get_address(&tempWD.startFunction, &address);
  wd_t_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempWD.depth, &address);
  wd_t_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempWD.num_paths, &address);
  wd_t_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempWD.num_linears, &address);
  wd_t_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempWD.B_rows, &address);
  wd_t_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempWD.B_cols, &address);
  wd_t_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempWD.p_size, &address);
  wd_t_int_displacements[7] = address - startaddress;
  MPI_Get_address(&tempWD.num_comp_d, &address);
  wd_t_int_displacements[8] = address - startaddress;
  MPI_Get_address(&tempWD.totalLength, &address);
  wd_t_int_displacements[9] = address - startaddress;

  // build the mpi datatype mpi_wd_t_int
  MPI_Type_create_struct(10, wd_t_int_length, wd_t_int_displacements, wd_t_int_datatypes, mpi_wd_t_int);
  MPI_Type_commit(mpi_wd_t_int);

  return;
}

// creates the MPI datatype mpi_stageData_t_int
void create_stageData_t_int(MPI_Datatype *mpi_sd_t_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  stageData_t_int tempSD;

  // arrays for length, displacement and datatypes in mpi_sd_d_int
  int sd_t_int_length[13] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype sd_t_int_datatypes[13] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint sd_t_int_displacements[13];

  // calculate displacements
  sd_t_int_displacements[0] = 0;
  MPI_Get_address(&tempSD.curr_precision, &startaddress);
  MPI_Get_address(&tempSD.depth_x, &address);
  sd_t_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempSD.depth_y, &address);
  sd_t_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempSD.num_paths, &address);
  sd_t_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempSD.useIntrinsicSlice, &address);
  sd_t_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempSD.B_rows, &address);
  sd_t_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempSD.B_cols, &address);
  sd_t_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempSD.p_size, &address);
  sd_t_int_displacements[7] = address - startaddress;
  MPI_Get_address(&tempSD.Bt_rows, &address);
  sd_t_int_displacements[8] = address - startaddress;
  MPI_Get_address(&tempSD.Bt_cols, &address);
  sd_t_int_displacements[9] = address - startaddress;
  MPI_Get_address(&tempSD.pt_size, &address);
  sd_t_int_displacements[10] = address - startaddress;
  MPI_Get_address(&tempSD.num_comp_d, &address);
  sd_t_int_displacements[11] = address - startaddress;
  MPI_Get_address(&tempSD.totalLength, &address);
  sd_t_int_displacements[12] = address - startaddress;

  // build the mpi datatype mpi_sd_t_int
  MPI_Type_create_struct(13, sd_t_int_length, sd_t_int_displacements, sd_t_int_datatypes, mpi_sd_t_int);
  MPI_Type_commit(mpi_sd_t_int);

  return;
}

// creates the MPI datatype mpi_codim_t_int
void create_codim_t_int(MPI_Datatype *mpi_codim_t_int)
{
  MPI_Datatype mpi_prog_t, mpi_preproc;

  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  codim_t_int tempCD;

  // create datatypes mpi_prog_t & mpi_preproc
  create_prog_t_int(&mpi_prog_t);
  create_preproc_data_int(&mpi_preproc);

  // arrays for length, displacement and datatypes in mpi_cd_t_int
  int cd_t_int_length[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype cd_t_int_datatypes[10] = {MPI_INT, mpi_prog_t, mpi_preproc, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint cd_t_int_displacements[10];

  // calculate displacements
  cd_t_int_displacements[0] = 0;
  MPI_Get_address(&tempCD.curr_precision, &startaddress);
  MPI_Get_address(&tempCD.Prog_int, &address);
  cd_t_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempCD.PPD_int, &address);
  cd_t_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempCD.num_funcs, &address);
  cd_t_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempCD.system_rank, &address);
  cd_t_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempCD.num_codim, &address);
  cd_t_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempCD.orig_variables, &address);
  cd_t_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempCD.new_variables, &address);
  cd_t_int_displacements[7] = address - startaddress;
  MPI_Get_address(&tempCD.num_comp_d, &address);
  cd_t_int_displacements[8] = address - startaddress;
  MPI_Get_address(&tempCD.totalLength, &address);
  cd_t_int_displacements[9] = address - startaddress;

  // build the mpi datatype mpi_codim_t_int
  MPI_Type_create_struct(10, cd_t_int_length, cd_t_int_displacements, cd_t_int_datatypes, mpi_codim_t_int);
  MPI_Type_commit(mpi_codim_t_int);

  return;
}

// creates the MPI datatype mpi_codimData_t_int
void create_codimData_t_int(MPI_Datatype *mpi_codimData_t_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  codimData_t_int tempCD;

  // arrays for length, displacement and datatypes in mpi_codimData_t_int
  int cd_t_int_length[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype cd_t_int_datatypes[11] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint cd_t_int_displacements[11];

  // calculate displacements
  cd_t_int_displacements[0] = 0;
  MPI_Get_address(&tempCD.curr_precision, &startaddress);
  MPI_Get_address(&tempCD.codim, &address);
  cd_t_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempCD.A_W_rows, &address);
  cd_t_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempCD.A_W_cols, &address);
  cd_t_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempCD.H_size, &address);
  cd_t_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempCD.useIntrinsicSlice, &address);
  cd_t_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempCD.B_rows, &address);
  cd_t_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempCD.B_cols, &address);
  cd_t_int_displacements[7] = address - startaddress;
  MPI_Get_address(&tempCD.p_size, &address);
  cd_t_int_displacements[8] = address - startaddress;
  MPI_Get_address(&tempCD.num_comp_d, &address);
  cd_t_int_displacements[9] = address - startaddress;
  MPI_Get_address(&tempCD.totalLength, &address);
  cd_t_int_displacements[10] = address - startaddress;

  // build the mpi datatype mpi_codimData_t_int
  MPI_Type_create_struct(11, cd_t_int_length, cd_t_int_displacements, cd_t_int_datatypes, mpi_codimData_t_int);
  MPI_Type_commit(mpi_codimData_t_int);

  return;
}

// creates the MPI datatype mpi_cascade_t_int
void create_cascade_t_int(MPI_Datatype *mpi_cascade_t_int)
{
  MPI_Datatype mpi_prog_t, mpi_preproc;

  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  cascade_t_int tempCD;

  // create datatypes mpi_prog_t & mpi_preproc
  create_prog_t_int(&mpi_prog_t);
  create_preproc_data_int(&mpi_preproc);

  // arrays for length, displacement and datatypes in mpi_cd_t_int
  int cd_t_int_length[22] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype cd_t_int_datatypes[22] = {MPI_INT, mpi_prog_t, mpi_preproc, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint cd_t_int_displacements[22];

  // calculate displacements
  cd_t_int_displacements[0] = 0;
  MPI_Get_address(&tempCD.curr_precision, &startaddress);
  MPI_Get_address(&tempCD.Prog_int, &address);
  cd_t_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempCD.PPD_int, &address);
  cd_t_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempCD.num_funcs, &address);
  cd_t_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempCD.system_rank, &address);
  cd_t_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempCD.num_codim, &address);
  cd_t_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempCD.orig_variables, &address);
  cd_t_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempCD.new_variables, &address);
  cd_t_int_displacements[7] = address - startaddress;
  MPI_Get_address(&tempCD.A_W_rows, &address);
  cd_t_int_displacements[8] = address - startaddress;
  MPI_Get_address(&tempCD.A_W_cols, &address);
  cd_t_int_displacements[9] = address - startaddress;
  MPI_Get_address(&tempCD.H_size, &address);
  cd_t_int_displacements[10] = address - startaddress;
  MPI_Get_address(&tempCD.C_rows, &address);
  cd_t_int_displacements[11] = address - startaddress;
  MPI_Get_address(&tempCD.C_cols, &address);
  cd_t_int_displacements[12] = address - startaddress;
  MPI_Get_address(&tempCD.R_rows, &address);
  cd_t_int_displacements[13] = address - startaddress;
  MPI_Get_address(&tempCD.R_cols, &address);
  cd_t_int_displacements[14] = address - startaddress;
  MPI_Get_address(&tempCD.T_size, &address);
  cd_t_int_displacements[15] = address - startaddress;
  MPI_Get_address(&tempCD.B_rows, &address);
  cd_t_int_displacements[16] = address - startaddress;
  MPI_Get_address(&tempCD.B_cols, &address);
  cd_t_int_displacements[17] = address - startaddress;
  MPI_Get_address(&tempCD.p_size, &address);
  cd_t_int_displacements[18] = address - startaddress;
  MPI_Get_address(&tempCD.num_int, &address);
  cd_t_int_displacements[19] = address - startaddress;
  MPI_Get_address(&tempCD.num_comp_d, &address);
  cd_t_int_displacements[20] = address - startaddress;
  MPI_Get_address(&tempCD.totalLength, &address);
  cd_t_int_displacements[21] = address - startaddress;

  // build the mpi datatype mpi_cascade_t_int
  MPI_Type_create_struct(22, cd_t_int_length, cd_t_int_displacements, cd_t_int_datatypes, mpi_cascade_t_int);
  MPI_Type_commit(mpi_cascade_t_int);

  return;
}

// creates the MPI datatype mpi_cascadeCodim_t_int
void create_cascadeCodim_t_int(MPI_Datatype *mpi_cascadeCodim_t_int)
{
  // arrays for length, displacement and datatypes in mpi_cascadeCodim_t_int
  int cd_t_int_length[1] = {1};
  MPI_Datatype cd_t_int_datatypes[1] = {MPI_INT};
  MPI_Aint cd_t_int_displacements[1] = {0};

  // build the mpi datatype mpi_cascadeCodim_t_int
  MPI_Type_create_struct(1, cd_t_int_length, cd_t_int_displacements, cd_t_int_datatypes, mpi_cascadeCodim_t_int);
  MPI_Type_commit(mpi_cascadeCodim_t_int);

  return;
}

// creates the MPI datatype mpi_slice_moving_int
void create_membership_slice_moving_t_int(MPI_Datatype *mpi_slice_moving_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  membership_slice_moving_t_int tempSlice;

  // create mpi_prog_t_int
  MPI_Datatype mpi_prog_t_int;
  create_prog_t_int(&mpi_prog_t_int);

  // arrays for length, displacement and datatypes in mpi_slice_moving_int
  int slice_t_int_length[17] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype slice_t_int_datatypes[17] = {mpi_prog_t_int, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint slice_t_int_displacements[17];

  // calculate displacements
  slice_t_int_displacements[0] = 0;
  MPI_Get_address(&tempSlice.Prog_int, &startaddress);
  MPI_Get_address(&tempSlice.curr_codim, &address);
  slice_t_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempSlice.orig_variables, &address);
  slice_t_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempSlice.curr_precision, &address);
  slice_t_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempSlice.A_rows, &address);
  slice_t_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempSlice.A_cols, &address);
  slice_t_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempSlice.B_rows, &address);
  slice_t_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempSlice.B_cols, &address);
  slice_t_int_displacements[7] = address - startaddress;
  MPI_Get_address(&tempSlice.p_size, &address);
  slice_t_int_displacements[8] = address - startaddress;
  MPI_Get_address(&tempSlice.startSliceVec_size, &address);
  slice_t_int_displacements[9] = address - startaddress;
  MPI_Get_address(&tempSlice.startSliceVec_init, &address);
  slice_t_int_displacements[10] = address - startaddress;
  MPI_Get_address(&tempSlice.targetSliceVec_size, &address);
  slice_t_int_displacements[11] = address - startaddress;
  MPI_Get_address(&tempSlice.targetSliceVec_init, &address);
  slice_t_int_displacements[12] = address - startaddress;
  MPI_Get_address(&tempSlice.K_rows, &address);
  slice_t_int_displacements[13] = address - startaddress;
  MPI_Get_address(&tempSlice.K_cols, &address);
  slice_t_int_displacements[14] = address - startaddress;
  MPI_Get_address(&tempSlice.num_comp_d, &address);
  slice_t_int_displacements[15] = address - startaddress;
  MPI_Get_address(&tempSlice.totalLength, &address);
  slice_t_int_displacements[16] = address - startaddress;

  // build the mpi datatype mpi_slice_moving_int
  MPI_Type_create_struct(17, slice_t_int_length, slice_t_int_displacements, slice_t_int_datatypes, mpi_slice_moving_int);
  MPI_Type_commit(mpi_slice_moving_int);

  // free the datatypes
  MPI_Type_free(&mpi_prog_t_int);

  return;
}

// creates the MPI datatype mpi_endpoint_data_d_int
void create_endpoint_data_d_int(MPI_Datatype *mpi_endpoint_d_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  endpoint_data_d_int tempEP;

  // arrays for length, displacement and datatypes in mpi_bed_mp_int
  int ep_d_int_length[8] = {1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype ep_d_int_datatypes[8] = {MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint ep_d_int_displacements[8];

  // calculate displacements
  ep_d_int_displacements[0] = 0;
  MPI_Get_address(&tempEP.cond_num, &startaddress);
  MPI_Get_address(&tempEP.corank, &address);
  ep_d_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempEP.smallest_nonzero_SV, &address);
  ep_d_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempEP.largest_zero_SV, &address);
  ep_d_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempEP.retVal, &address);
  ep_d_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempEP.pt_size, &address);
  ep_d_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempEP.last_approx_size, &address);
  ep_d_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempEP.num_comp_d, &address);
  ep_d_int_displacements[7] = address - startaddress;

  // build the mpi datatype mpi_endpoint_d_int
  MPI_Type_create_struct(8, ep_d_int_length, ep_d_int_displacements, ep_d_int_datatypes, mpi_endpoint_d_int);
  MPI_Type_commit(mpi_endpoint_d_int);

  return;
}

// creates the MPI datatype mpi_endpoint_data_mp_int
void create_endpoint_data_mp_int(MPI_Datatype *mpi_endpoint_mp_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  endpoint_data_mp_int tempEP;

  // arrays for length, displacement and datatypes in mpi_bed_mp_int
  int ep_mp_int_length[8] = {1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype ep_mp_int_datatypes[8] = {MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint ep_mp_int_displacements[8];

  // calculate displacements
  ep_mp_int_displacements[0] = 0;
  MPI_Get_address(&tempEP.cond_num, &startaddress);
  MPI_Get_address(&tempEP.corank, &address);
  ep_mp_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempEP.smallest_nonzero_SV, &address);
  ep_mp_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempEP.largest_zero_SV, &address);
  ep_mp_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempEP.retVal, &address);
  ep_mp_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempEP.pt_size, &address);
  ep_mp_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempEP.last_approx_size, &address);
  ep_mp_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempEP.totalLength, &address);
  ep_mp_int_displacements[7] = address - startaddress;

  // build the mpi datatype mpi_endpoint_mp_int
  MPI_Type_create_struct(8, ep_mp_int_length, ep_mp_int_displacements, ep_mp_int_datatypes, mpi_endpoint_mp_int);
  MPI_Type_commit(mpi_endpoint_mp_int);

  return;
}

// creates the MPI datatype mpi_endpoint_data_amp_int
void create_endpoint_data_amp_int(MPI_Datatype *mpi_endpoint_amp_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  endpoint_data_amp_int tempEP;

  // arrays for length, displacement and datatypes in mpi_bed_mp_int
  int ep_amp_int_length[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype ep_amp_int_datatypes[11] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint ep_amp_int_displacements[11];

  // calculate displacements
  ep_amp_int_displacements[0] = 0;
  MPI_Get_address(&tempEP.curr_prec, &startaddress);
  MPI_Get_address(&tempEP.last_approx_prec, &address);
  ep_amp_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempEP.cond_num, &address);
  ep_amp_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempEP.corank, &address);
  ep_amp_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempEP.smallest_nonzero_SV, &address);
  ep_amp_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempEP.largest_zero_SV, &address);
  ep_amp_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempEP.retVal, &address);
  ep_amp_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempEP.pt_size, &address);
  ep_amp_int_displacements[7] = address - startaddress;
  MPI_Get_address(&tempEP.last_approx_size, &address);
  ep_amp_int_displacements[8] = address - startaddress;
  MPI_Get_address(&tempEP.num_comp_d, &address);
  ep_amp_int_displacements[9] = address - startaddress;
  MPI_Get_address(&tempEP.totalLength, &address);
  ep_amp_int_displacements[10] = address - startaddress;

  // build the mpi datatype mpi_endpoint_amp_int
  MPI_Type_create_struct(11, ep_amp_int_length, ep_amp_int_displacements, ep_amp_int_datatypes, mpi_endpoint_amp_int);
  MPI_Type_commit(mpi_endpoint_amp_int);

  return;
}

// creates the MPI datatype mpi_witnessCodim_int
void create_witnessCodim_t_int(MPI_Datatype *mpi_witnessCodim_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  witnessCodim_t_int tempWit;

  // arrays for length, displacement and datatypes in mpi_witnessCodim_int
  int wc_int_length[13] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype wc_int_datatypes[13] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint wc_int_displacements[13];

  // calculate displacements
  wc_int_displacements[0] = 0;
  MPI_Get_address(&tempWit.codim, &startaddress);
  MPI_Get_address(&tempWit.num_set, &address);
  wc_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempWit.num_nonsing, &address);
  wc_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempWit.num_sing, &address);
  wc_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempWit.A_rows, &address);
  wc_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempWit.A_cols, &address);
  wc_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempWit.H_size, &address);
  wc_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempWit.B_rows, &address);
  wc_int_displacements[7] = address - startaddress;
  MPI_Get_address(&tempWit.B_cols, &address);
  wc_int_displacements[8] = address - startaddress;
  MPI_Get_address(&tempWit.p_size, &address);
  wc_int_displacements[9] = address - startaddress;
  MPI_Get_address(&tempWit.curr_prec, &address);
  wc_int_displacements[10] = address - startaddress;
  MPI_Get_address(&tempWit.num_comp_d, &address);
  wc_int_displacements[11] = address - startaddress;
  MPI_Get_address(&tempWit.totalLength, &address);
  wc_int_displacements[12] = address - startaddress;

  // build the mpi datatype mpi_witnessCodim_int
  MPI_Type_create_struct(13, wc_int_length, wc_int_displacements, wc_int_datatypes, mpi_witnessCodim_int);
  MPI_Type_commit(mpi_witnessCodim_int);

  return;
}

// creates the MPI datatype mpi_witness_int
void create_witness_t_int(MPI_Datatype *mpi_witness_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  witness_t_int tempWit;

  MPI_Datatype mpi_prog_int, mpi_PPD_int;

  // create the MPI datatypes
  create_prog_t_int(&mpi_prog_int);
  create_preproc_data_int(&mpi_PPD_int);

  // arrays for length, displacement and datatypes in mpi_witnessCodim_int
  int wc_int_length[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype wc_int_datatypes[10] = {mpi_prog_int, mpi_PPD_int, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint wc_int_displacements[10];

  // calculate displacements
  wc_int_displacements[0] = 0;
  MPI_Get_address(&tempWit.Prog_int, &startaddress);
  MPI_Get_address(&tempWit.PPD_int, &address);
  wc_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempWit.num_funcs, &address);
  wc_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempWit.num_codim, &address);
  wc_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempWit.curr_precision, &address);
  wc_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempWit.system_rank, &address);
  wc_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempWit.orig_variables, &address);
  wc_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempWit.new_variables, &address);
  wc_int_displacements[7] = address - startaddress;
  MPI_Get_address(&tempWit.num_comp_d, &address);
  wc_int_displacements[8] = address - startaddress;
  MPI_Get_address(&tempWit.totalLength, &address);
  wc_int_displacements[9] = address - startaddress;

  // build the mpi datatype mpi_witnessCodim_int
  MPI_Type_create_struct(10, wc_int_length, wc_int_displacements, wc_int_datatypes, mpi_witness_int);
  MPI_Type_commit(mpi_witness_int);

  return;
}

// creates the MPI datatype mpi_regen_pos_dim_t_int
void create_regen_pos_dim_t_int(MPI_Datatype *mpi_regen_pos_dim_t_int)
{
  MPI_Datatype mpi_prog_t, mpi_preproc;

  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  regen_pos_dim_t_int tempRPD;

  // create datatypes mpi_prog_t & mpi_preproc
  create_prog_t_int(&mpi_prog_t);
  create_preproc_data_int(&mpi_preproc);

  // arrays for length, displacement and datatypes in mpi_rpd_t_int
  int rpd_t_int_length[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype rpd_t_int_datatypes[16] = {mpi_prog_t, mpi_preproc, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint rpd_t_int_displacements[16];

  // calculate displacements
  rpd_t_int_displacements[0] = 0;
  MPI_Get_address(&tempRPD.Prog_int, &startaddress);
  MPI_Get_address(&tempRPD.PPD_int, &address);
  rpd_t_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempRPD.system_rank, &address);
  rpd_t_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempRPD.orig_variables, &address);
  rpd_t_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempRPD.new_variables, &address);
  rpd_t_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempRPD.C_rows, &address);
  rpd_t_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempRPD.C_cols, &address);
  rpd_t_int_displacements[6] = address - startaddress;
  MPI_Get_address(&tempRPD.H_size, &address);
  rpd_t_int_displacements[7] = address - startaddress;
  MPI_Get_address(&tempRPD.patchCoeff_size, &address);
  rpd_t_int_displacements[8] = address - startaddress;
  MPI_Get_address(&tempRPD.num_funcs, &address);
  rpd_t_int_displacements[9] = address - startaddress;
  MPI_Get_address(&tempRPD.num_codim, &address);
  rpd_t_int_displacements[10] = address - startaddress;
  MPI_Get_address(&tempRPD.curr_precision, &address);
  rpd_t_int_displacements[11] = address - startaddress;
  MPI_Get_address(&tempRPD.sameA, &address);
  rpd_t_int_displacements[12] = address - startaddress;
  MPI_Get_address(&tempRPD.num_int, &address);
  rpd_t_int_displacements[13] = address - startaddress;
  MPI_Get_address(&tempRPD.num_comp_d, &address);
  rpd_t_int_displacements[14] = address - startaddress;
  MPI_Get_address(&tempRPD.totalLength, &address);
  rpd_t_int_displacements[15] = address - startaddress;

  // build the mpi datatype mpi_regen_pos_dim_t_int
  MPI_Type_create_struct(16, rpd_t_int_length, rpd_t_int_displacements, rpd_t_int_datatypes, mpi_regen_pos_dim_t_int);
  MPI_Type_commit(mpi_regen_pos_dim_t_int);

  // free the datatypes
  MPI_Type_free(&mpi_prog_t);
  MPI_Type_free(&mpi_preproc);
 
  return;
}

// creates the MPI datatype mpi_regenCodim_t_int
void create_regenCodim_t_int(MPI_Datatype *mpi_regenCodim_t_int)
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  regenCodim_t_int tempRPD;

  // arrays for length, displacement and datatypes in mpi_rpd_t_int
  int rpd_t_int_length[7] = {1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype rpd_t_int_datatypes[7] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint rpd_t_int_displacements[7];

  // calculate displacements
  rpd_t_int_displacements[0] = 0;
  MPI_Get_address(&tempRPD.codim, &startaddress);
  MPI_Get_address(&tempRPD.useIntrinsicSlice, &address);
  rpd_t_int_displacements[1] = address - startaddress;
  MPI_Get_address(&tempRPD.B_rows, &address);
  rpd_t_int_displacements[2] = address - startaddress;
  MPI_Get_address(&tempRPD.B_cols, &address);
  rpd_t_int_displacements[3] = address - startaddress;
  MPI_Get_address(&tempRPD.p_size, &address);
  rpd_t_int_displacements[4] = address - startaddress;
  MPI_Get_address(&tempRPD.num_comp_d, &address);
  rpd_t_int_displacements[5] = address - startaddress;
  MPI_Get_address(&tempRPD.totalLength, &address);
  rpd_t_int_displacements[6] = address - startaddress;

  // build the mpi datatype mpi_regenCodim_t_int
  MPI_Type_create_struct(7, rpd_t_int_length, rpd_t_int_displacements, rpd_t_int_datatypes, mpi_regenCodim_t_int);
  MPI_Type_commit(mpi_regenCodim_t_int);

  return;
}


#endif

