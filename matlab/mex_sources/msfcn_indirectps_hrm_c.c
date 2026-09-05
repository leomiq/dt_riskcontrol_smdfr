#define S_FUNCTION_NAME  msfcn_indirectps_hrm_c
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>
#include <stdlib.h>

/* Function: mdlInitializeSizes ==============================================
 * Abstract:
 *   Setup sizes of input/output/DWork vectors
 */

static char_T dworkNames[11][16];  /* MUST be static so it survives beyond this call */

/*
 * D-Work State Variables (xS1…xS11):
 *
 *  S1  (DWork 0, width=5)   – current_sys_state[1..5]
 *      • Current system state vector (one entry per component)
 *
 *  S2  (DWork 1, width=1)   – current_comp
 *      • Index of the component currently undergoing transition
 *
 *  S3  (DWork 2, width=5)   – future_sys_state[1..5]
 *      • Next (future) system state vector after the upcoming transition
 *
 *  S4  (DWork 3, width=1)   – future_comp
 *      • Index of the component selected for the next transition
 *
 *  S5  (DWork 4, width=1)   – future_ttf
 *      • Time-to-failure sampled for that next transition
 *
 *  S6  (DWork 5, width=1)   – prev_ttf
 *      • Previous transition timestamp (used for Δ-time calculations)
 *
 *  S7  (DWork 6, width=1)   – failure_flag
 *      • 0 = system healthy, 1 = system has failed
 *
 *  S8  (DWork 7, width=10)  – mat_ttf
 *      • Flattened matrix of all sampled TTFs (size = Ncomps × maxStates)
 *
 *  S9  (DWork 8, width=1)   – trans_c
 *      • Index of the component whose TTF was minimal (next to transition)
 *
 *  S10 (DWork 9, width=1)   – comp_branch_event
 *      • Complementary-branch event flag (1 = branch event active, else 0)
 *
 *  S11 (DWork 10, width=1)  – step_counter
 *      • Incremental counter of how many Update calls have occurred
 */

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 0); // Adjust based on parameters

    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch reported by Simulink */
    }

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    /* Define input ports */
    if (!ssSetNumInputPorts(S, 10)) return;
    ssSetInputPortWidth(S, 0, 1);     // Current simulation time
    ssSetInputPortRequiredContiguous(S,0,1);
    ssSetInputPortWidth(S, 1, 1);     // Ncomps
    ssSetInputPortRequiredContiguous(S,1,1);
    ssSetInputPortWidth(S, 2, 5);     // Nstates
    ssSetInputPortRequiredContiguous(S,2,1);
    ssSetInputPortMatrixDimensions(S, 3, 5, 2);  // mat_states_index
    ssSetInputPortRequiredContiguous(S,3,1);
    ssSetInputPortWidth(S, 4, 20);    // mat_trans (flattened 2x2x5)
    ssSetInputPortRequiredContiguous(S,4,1);
    ssSetInputPortMatrixDimensions(S, 5, 5, 5);  // failed_states
    ssSetInputPortRequiredContiguous(S,5,1);
    ssSetInputPortWidth(S, 6, 5);     // initial_state
    ssSetInputPortRequiredContiguous(S,6,1);
    ssSetInputPortWidth(S, 7, 2);     // initial BRANCH state
    ssSetInputPortRequiredContiguous(S,7,1);
    ssSetInputPortWidth(S, 8, 1);     // complementary_branch_flag
    ssSetInputPortRequiredContiguous(S,8,1);
    ssSetInputPortWidth(S, 9, 1);     // parent_last_time
    ssSetInputPortRequiredContiguous(S,9,1);

    for (int i=0; i<10; i++)
        ssSetInputPortDirectFeedThrough(S, i, 1);

    /* Define output ports */
    if (!ssSetNumOutputPorts(S,9)) return;
    ssSetOutputPortWidth(S, 0, 5);
    ssSetOutputPortWidth(S, 1, 1);
    ssSetOutputPortWidth(S, 2, 5);
    ssSetOutputPortWidth(S, 3, 1);
    ssSetOutputPortWidth(S, 4, 1);
    ssSetOutputPortWidth(S, 5, 1);
    ssSetOutputPortWidth(S, 6, 10);
    ssSetOutputPortWidth(S, 7, 1);
    ssSetOutputPortWidth(S, 8, 1);

    /* DWork vectors */
    ssSetNumDWork(S, 11);
    for (int i = 0; i < 11; i++) {
        /* widths & types as before */
        ssSetDWorkWidth   (S, i, (i==0||i==2) ? 5 : (i==7 ? 10 : 1));
        ssSetDWorkDataType(S, i, SS_DOUBLE);

        /* fill the persistent buffer */
        snprintf(dworkNames[i], 16, "xS%d", i+1);

        /* now give it to Simulink */
        ssSetDWorkName(S, i, dworkNames[i]);
        ssSetDWorkUsedAsDState(S, i, 1);
    }

    /* One sample time */
    ssSetNumSampleTimes(S, 1);

    /* Explicitly zero out any extra work vectors */
    ssSetNumRWork       (S, 0);
    ssSetNumIWork       (S, 0);
    ssSetNumPWork       (S, 0);
    ssSetNumModes       (S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);
    ssSetOptions(S, 0);
}

/* Function: mdlInitializeSampleTimes ========================================
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

/* Function: mdlInitializeConditions ========================================
 */
#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S)
{
    real_T *xS_1 = ssGetDWork(S,0);
    real_T *xS_2 = ssGetDWork(S,1);
    real_T *xS_3 = ssGetDWork(S,2);
    real_T *xS_4 = ssGetDWork(S,3);
    real_T *xS_5 = ssGetDWork(S,4);
    real_T *xS_6 = ssGetDWork(S,5);   /* prev_ttf */
    real_T *xS_7 = ssGetDWork(S,6);   /* failure_flag */

    for (int i=0; i<5; i++){
        xS_1[i] = 1.0;
        xS_3[i] = 0.0;
    }
    *xS_2 = 1.0;
    *xS_4 = 1.0;
    *xS_5 = 0.0;
    xS_6[0] = 0.0;   /* initialize previous‐time‐to‐failure */
    xS_7[0] = 0.0;
}

/* Function: mdlOutputs ======================================================
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    for (int i=0; i<9; i++){
        real_T *out = ssGetOutputPortRealSignal(S,i);
        real_T *dwork = ssGetDWork(S,i);
        int width = ssGetOutputPortWidth(S,i);
        for (int j=0; j<width; j++)
            out[j] = dwork[j];
    }
}

/* Function: mdlUpdate =======================================================
 */
#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{
    const real_T t = *ssGetInputPortRealSignal(S,0);
    const int Ncomps = (int)*ssGetInputPortRealSignal(S,1);
    const int Nstates = (int)*ssGetInputPortRealSignal(S,2);
    const real_T *mat_states_index = ssGetInputPortRealSignal(S,3);
    const real_T *mat_trans = ssGetInputPortRealSignal(S,4);
    const real_T *failed_states = ssGetInputPortRealSignal(S,5);
    const real_T complementary_branch_flag = *ssGetInputPortRealSignal(S,8);

    real_T *current_sys_state = ssGetDWork(S,0);
    real_T *current_state_comp = ssGetDWork(S,1);
    real_T *trans_sys = ssGetDWork(S,2);
    real_T *new_trans = ssGetDWork(S,3);
    real_T *new_ttf = ssGetDWork(S,4);
    real_T *prev_ttf = ssGetDWork(S,5);
    real_T *failure_flag = ssGetDWork(S,6);
    real_T *mat_ttf = ssGetDWork(S,7);
    real_T *trans_c = ssGetDWork(S,8);
    real_T *comp_branch_event = ssGetDWork(S,9);
    real_T *step_counter = ssGetDWork(S,10);

    step_counter[0]++;
    comp_branch_event[0] = (complementary_branch_flag == 1 && step_counter[0] <= 1) ? 1 : 0;

    if (new_ttf[0] == 0) {
        for (int c = 0; c < Ncomps; c++) {
            int curr_state = (int)current_sys_state[c];
            for (int j = 0; j < Nstates; j++) {
                if ((j != curr_state) && (mat_states_index[c * Nstates + j] != 0)) {
                    int idx = c + 2 * (curr_state + 2 * j); // 2x2x5 flattened index
                    double lambda = mat_trans[idx];
                    double rand_num = (double)rand() / RAND_MAX;
                    mat_ttf[c * Nstates + j] = -log(rand_num) / lambda;
                }
            }
        }
        double min_ttf = 1e10;
        int min_c = 0, min_s = 0;
        for (int c = 0; c < Ncomps; c++) {
            for (int j = 0; j < Nstates; j++) {
                double val = mat_ttf[c * Nstates + j];
                if (val > 0 && val < min_ttf) {
                    min_ttf = val;
                    min_c = c;
                    min_s = j;
                }
            }
        }
        trans_c[0] = min_c;
        new_trans[0] = min_s;
        prev_ttf[0] = t;
        new_ttf[0] = min_ttf;

        for (int i = 0; i < Ncomps; i++)
            trans_sys[i] = current_sys_state[i];
        trans_sys[min_c] = min_s;
    }
    else if (comp_branch_event[0] == 0 && t >= new_ttf[0]) {
        for (int i = 0; i < Ncomps; i++)
            current_sys_state[i] = trans_sys[i];
        current_state_comp[0] = new_trans[0];

        failure_flag[0] = 0;
        for (int i = 0; i < Ncomps; i++) {
            if (failed_states[i * Ncomps + (int)current_sys_state[i] - 1] == 1) {
                failure_flag[0] = 1;
                break;
            }
        }

        prev_ttf[0] = t;
        new_ttf[0] = 0;
    }
}

static void mdlTerminate(SimStruct *S) {}

#ifdef MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
