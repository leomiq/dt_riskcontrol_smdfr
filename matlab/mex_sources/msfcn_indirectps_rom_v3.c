#define S_FUNCTION_NAME  msfcn_indirectps_rom_v3
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

static char_T dworkNames[9][16];

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 0);

    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return;
    }

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    // Input ports
    if (!ssSetNumInputPorts(S, 10)) return;
    ssSetInputPortWidth(S, 0, 1);     // t
    ssSetInputPortWidth(S, 1, 1);     // Ncomps
    ssSetInputPortWidth(S, 2, 5);     // Nstates
    ssSetInputPortMatrixDimensions(S, 3, 5, 2);  // mat_states_index
    ssSetInputPortWidth(S, 4, 20);    // mat_trans
    ssSetInputPortMatrixDimensions(S, 5, 5, 5);  // failed_states
    ssSetInputPortWidth(S, 6, 5);     // initial_state
    ssSetInputPortWidth(S, 7, 2);     // initial_branch_state
    ssSetInputPortWidth(S, 8, 1);     // complementary_branch_flag
    ssSetInputPortWidth(S, 9, 1);     // parent_last_time

    for (int i = 0; i < 10; i++)
        ssSetInputPortRequiredContiguous(S, i, 1);

    for (int i = 0; i < 10; i++)
        ssSetInputPortDirectFeedThrough(S, i, 1);

    // Output ports (1-to-1 with DWorks)
    if (!ssSetNumOutputPorts(S, 9)) return;
    ssSetOutputPortWidth(S, 0, 5);   // current_sys_state
    ssSetOutputPortWidth(S, 1, 1);   // current_state_comp
    ssSetOutputPortWidth(S, 2, 5);   // trans_sys
    ssSetOutputPortWidth(S, 3, 1);   // new_trans
    ssSetOutputPortWidth(S, 4, 1);   // new_ttf
    ssSetOutputPortWidth(S, 5, 1);   // prev_ttf
    ssSetOutputPortWidth(S, 6, 1);   // failure_flag
    ssSetOutputPortWidth(S, 7, 10);  // mat_ttf
    ssSetOutputPortWidth(S, 8, 1);   // trans_c

    // DWork vectors
    ssSetNumDWork(S, 9);
    for (int i = 0; i < 9; i++) {
        int width = (i == 0 || i == 2) ? 5 : (i == 7 ? 10 : 1);
        ssSetDWorkWidth(S, i, width);
        ssSetDWorkDataType(S, i, SS_DOUBLE);
        snprintf(dworkNames[i], 16, "xS%d", i + 1);
        ssSetDWorkName(S, i, dworkNames[i]);
        ssSetDWorkUsedAsDState(S, i, 1);
    }

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 1);  // step_counter
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);
    ssSetOptions(S, 0);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S)
{
    for (int i = 0; i < 9; ++i) {
        real_T *x = ssGetDWork(S, i);
        int width = ssGetDWorkWidth(S, i);
        for (int j = 0; j < width; ++j)
            x[j] = 0.0;
    }

    real_T *x0 = ssGetDWork(S, 0);  // current_sys_state
    for (int i = 0; i < 5; ++i)
        x0[i] = 1.0;

    ssGetIWork(S)[0] = 0; // step_counter
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    for (int i = 0; i < 9; i++) {
        real_T *out = ssGetOutputPortRealSignal(S, i);
        real_T *dwork = ssGetDWork(S, i);
        int width = ssGetOutputPortWidth(S, i);
        for (int j = 0; j < width; j++)
            out[j] = dwork[j];
    }
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{
    real_T *current_sys_state    = ssGetDWork(S, 0);
    real_T *current_state_comp   = ssGetDWork(S, 1);
    real_T *trans_sys            = ssGetDWork(S, 2);
    real_T *new_trans            = ssGetDWork(S, 3);
    real_T *new_ttf              = ssGetDWork(S, 4);
    real_T *prev_ttf             = ssGetDWork(S, 5);
    real_T *failure_flag         = ssGetDWork(S, 6);
    real_T *mat_ttf              = ssGetDWork(S, 7);
    real_T *trans_c              = ssGetDWork(S, 8);
    int    *step_counter         = ssGetIWork(S);

    // Inputs
    real_T t = *ssGetInputPortRealSignal(S, 0);
    int Ncomps = (int)*ssGetInputPortRealSignal(S, 1);
    int Nstates = (int)*ssGetInputPortRealSignal(S, 2);
    const real_T *mat_states_index = ssGetInputPortRealSignal(S, 3);
    const real_T *mat_trans = ssGetInputPortRealSignal(S, 4);
    const real_T *failed_states = ssGetInputPortRealSignal(S, 5);
    real_T complementary_branch_flag = *ssGetInputPortRealSignal(S, 8);

    const real_T tsim = 1400.0;
    const real_T thres_sampled_time = 2.5;

    (*step_counter)++;
    int comp_branch_event = (complementary_branch_flag == 1 && *step_counter <= 1) ? 1 : 0;

    if (new_ttf[0] == 0) {
        for (int c = 0; c < Ncomps; c++) {
            int current_state = (int)current_sys_state[c] - 1;
            for (int j = 0; j < Nstates; j++) {
                if ((j != current_state) && (mat_states_index[c * Nstates + j] != 0)) {
                    real_T lambda = mat_trans[c * 4 + current_state * 2 + j];
                    real_T ttf = t - log((real_T)rand() / (real_T)RAND_MAX) / lambda;
                    mat_ttf[c * Nstates + j] = ttf;
                } else {
                    mat_ttf[c * Nstates + j] = 0.0;
                }
            }
        }

        real_T min_ttf = 1e10;
        int min_c = 0, min_j = 0;
        for (int c = 0; c < Ncomps; c++) {
            for (int j = 0; j < Nstates; j++) {
                real_T val = mat_ttf[c * Nstates + j];
                if (val > 0 && val < min_ttf) {
                    min_ttf = val;
                    min_c = c;
                    min_j = j;
                }
            }
        }

        prev_ttf[0] = t;
        new_ttf[0] = min_ttf;

        if (min_ttf > (tsim - thres_sampled_time) && min_ttf < tsim) {
            for (int k = 0; k < Ncomps * Nstates; k++) {
                if (mat_ttf[k] > 0)
                    mat_ttf[k] += thres_sampled_time;
            }
            new_ttf[0] += thres_sampled_time;
        }

        for (int i = 0; i < Ncomps; i++)
            trans_sys[i] = current_sys_state[i];

        trans_sys[min_c] = min_j + 1;
        new_trans[0] = min_j + 1;
        trans_c[0] = min_c + 1;
    } else if (comp_branch_event == 0 && t >= new_ttf[0]) {
        for (int i = 0; i < Ncomps; i++)
            current_sys_state[i] = trans_sys[i];

        current_state_comp[0] = new_trans[0];

        failure_flag[0] = 0;
        for (int i = 0; i < Ncomps; i++) {
            int state = (int)current_sys_state[i];
            if (failed_states[i * Nstates + state - 1] == 1) {
                failure_flag[0] = 1;
                break;
            }
        }

        prev_ttf[0] = t;
        new_ttf[0] = 0.0;
    }
}

static void mdlTerminate(SimStruct *S) {}

#ifdef MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
