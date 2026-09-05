#define S_FUNCTION_NAME  msfcn_schedule_hrm_c
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>
#include <string.h>

/*
 * Level-2 S-Function: msfcn_schedule_hrm_c
 *
 * D-Work State Variables (S#):
 *  S1 (DWork 0, width=1) – act_PS1: activation flag for the protection system
 *  S2 (DWork 1, width=2) – state_PS1: [demand, availability] current state vector
 *  S3 (DWork 2, width=1) – step_counter: counts Update calls
 */

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return;

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 7)) return;
    /* Ports 0,1,2,4,5,6 are width=1; port 3 is width=2 */
    for (int i=0; i<7; i++) {
        ssSetInputPortWidth           (S, i, (i==3 ? 2 : 1));
        ssSetInputPortDirectFeedThrough(S, i, 1);
        ssSetInputPortRequiredContiguous(S, i, 1);
    }

    if (!ssSetNumOutputPorts(S,2)) return;
    ssSetOutputPortWidth(S, 0, 1);
    ssSetOutputPortWidth(S, 1, 2);

    ssSetNumDWork(S, 3);

    ssSetDWorkWidth(S, 0, 1);  // act_PS1
    ssSetDWorkDataType(S, 0, SS_DOUBLE);
    ssSetDWorkUsageType(S, 0, SS_DWORK_USED_AS_DSTATE);
    ssSetDWorkName(S, 0, "Ssch_1");

    ssSetDWorkWidth(S, 1, 2);  // state_PS1
    ssSetDWorkDataType(S, 1, SS_DOUBLE);
    ssSetDWorkUsageType(S, 1, SS_DWORK_USED_AS_DSTATE);
    ssSetDWorkName(S, 1, "Ssch_2");

    ssSetDWorkWidth(S, 2, 1);  // step_counter2
    ssSetDWorkDataType(S, 2, SS_DOUBLE);
    ssSetDWorkUsageType(S, 2, SS_DWORK_USED_AS_DSTATE);
    ssSetDWorkName(S, 2, "Ssch_3");

     /* Register that you have one sample-time */
    ssSetNumSampleTimes(S, 1);

    /* Explicitly zero out the extra work vectors */
    ssSetNumRWork       (S, 0);
    ssSetNumIWork       (S, 0);
    ssSetNumPWork       (S, 0);
    ssSetNumModes       (S, 0);
    ssSetNumNonsampledZCs(S, 0);

    /* Optional: make your block behave like a built-in */
    ssSetOperatingPointCompliance     (S, USE_DEFAULT_OPERATING_POINT);
    ssSetRuntimeThreadSafetyCompliance(S, RUNTIME_THREAD_SAFETY_COMPLIANCE_TRUE);
    ssSetOptions                      (S, SS_OPTION_EXCEPTION_FREE_CODE);

    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S)
{
    real_T *act_PS1 = ssGetDWork(S,0);
    real_T *state_PS1 = ssGetDWork(S,1);
    real_T *step_counter2 = ssGetDWork(S,2);

    *act_PS1 = 0;
    state_PS1[0] = 0;
    state_PS1[1] = 1;
    *step_counter2 = 0;
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T *y0 = ssGetOutputPortRealSignal(S,0);
    real_T *y1 = ssGetOutputPortRealSignal(S,1);

    real_T *act_PS1 = ssGetDWork(S,0);
    real_T *state_PS1 = ssGetDWork(S,1);

    y0[0] = *act_PS1;
    memcpy(y1, state_PS1, 2 * sizeof(real_T));
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{
    /* bump your counter */
    real_T *step_counter = ssGetDWork(S,2);
    (*step_counter)++;

    /* read inputs */
    real_T  state_var_ps             = *ssGetInputPortRealSignal(S,0);
    real_T  failure_flag             = *ssGetInputPortRealSignal(S,1);
    real_T  tres_demand_ps           = *ssGetInputPortRealSignal(S,2);
    const real_T *prev_state         = ssGetInputPortRealSignal(S,3);
    real_T  parent_last_time         = *ssGetInputPortRealSignal(S,4);
    real_T  t                        = *ssGetInputPortRealSignal(S,5);
    real_T  complementary_branch_flag = *ssGetInputPortRealSignal(S,6);

    /* compute ava_ps and demand_ps exactly as MATLAB */
    int ava_ps    = (failure_flag == 0) ? 1 : 0;
    int demand_ps = (state_var_ps == tres_demand_ps) ? 1 : 0;

    /* pack current_state */
    int current_state[2] = { demand_ps, ava_ps };
    bool   all_same = (current_state[0]==prev_state[0])
                   && (current_state[1]==prev_state[1]);

    /* pointers to D-Works */
    real_T *act_PS1   = ssGetDWork(S,0);
    real_T *state_PS1 = ssGetDWork(S,1);

    if (complementary_branch_flag == 1) {
        /* ——— exactly MATLAB’s “branch==1” code ——— */
        if (demand_ps == 1) {
            if (ava_ps == 1) {
                *act_PS1 = 1.0;
                if (!all_same && *step_counter >= 2) {
                    ssSetStopRequested(S,1);
                }
            } else {
                *act_PS1 = 0.0;
                if (!all_same && *step_counter >= 3) {
                    ssSetStopRequested(S,1);
                }
            }
        } else {
            *act_PS1 = 0.0;
            if (prev_state[0] == 1) {
                ssSetStopRequested(S,1);
            }
        }
    }
    else if (complementary_branch_flag == 0) {
        /* ——— exactly MATLAB’s “branch==0” code ——— */
        if (demand_ps == 1) {
            if (ava_ps == 1) {
                *act_PS1 = 1.0;
                /* no ssSetStopRequested here (commented out in MATLAB) */
            } else {
                *act_PS1 = 0.0;
                /* no ssSetStopRequested here either */
            }
        } else {
            *act_PS1 = 0.0;
            if (prev_state[0] == 1) {
                ssSetStopRequested(S,1);
            }
        }
    }

    /* finally, store the new state vector */
    state_PS1[0] = (real_T)current_state[0];  /* demand */
    state_PS1[1] = (real_T)current_state[1];  /* availability */
}

static void mdlTerminate(SimStruct *S) {}

#ifdef MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
