#define S_FUNCTION_NAME  msfcn_schedule_rom_v2
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>
#include <string.h>

/*
 * Level-2 S-Function: msfcn_schedule_rom_v2
 *
 * D-Work State Variables (S#):
 *  S1 (DWork 0, width=1) – act_PS1: activation flag for the protection system
 *  S2 (DWork 1, width=2) – state_PS1: [demand, availability] current state vector
 *  S3 (DWork 2, width=1) – step_counter2: counts Update calls
 */

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return;

    /* No continuous or discrete states */
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    /* Inputs: 7 ports, port 3 is width=2 */
    if (!ssSetNumInputPorts(S, 7)) return;
    for (int i = 0; i < 7; i++) {
        ssSetInputPortWidth(S, i, (i == 3 ? 2 : 1));
        ssSetInputPortDirectFeedThrough(S, i, 1);
        ssSetInputPortRequiredContiguous(S, i, 1);
    }

    /* Outputs: 2 ports, port 1 width=1, port 2 width=2 */
    if (!ssSetNumOutputPorts(S, 2)) return;
    ssSetOutputPortWidth(S, 0, 1);
    ssSetOutputPortWidth(S, 1, 2);

    /* D-Work: 3 vectors, all used as discrete states */
    ssSetNumDWork(S, 3);
    ssSetDWorkWidth(S, 0, 1);  /* S1: act_PS1 */
    ssSetDWorkDataType(S, 0, SS_DOUBLE);
    ssSetDWorkUsageType(S, 0, SS_DWORK_USED_AS_DSTATE);
    ssSetDWorkName(S, 0, "Ssch_1");

    ssSetDWorkWidth(S, 1, 2);  /* S2: state_PS1 */
    ssSetDWorkDataType(S, 1, SS_DOUBLE);
    ssSetDWorkUsageType(S, 1, SS_DWORK_USED_AS_DSTATE);
    ssSetDWorkName(S, 1, "Ssch_2");

    ssSetDWorkWidth(S, 2, 1);  /* S3: step_counter2 */
    ssSetDWorkDataType(S, 2, SS_DOUBLE);
    ssSetDWorkUsageType(S, 2, SS_DWORK_USED_AS_DSTATE);
    ssSetDWorkName(S, 2, "Ssch_3");

    /* One sample time */
    ssSetNumSampleTimes(S, 1);

    /* No other work vectors */
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    /* Compliance and thread-safety */
    ssSetOperatingPointCompliance(S, USE_DEFAULT_OPERATING_POINT);
    ssSetRuntimeThreadSafetyCompliance(S, RUNTIME_THREAD_SAFETY_COMPLIANCE_TRUE);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
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
    real_T *act_PS1      = ssGetDWork(S, 0);
    real_T *state_PS1    = ssGetDWork(S, 1);
    real_T *step_counter = ssGetDWork(S, 2);

    /* Initialize discrete states */
    act_PS1[0]       = 0.0;
    state_PS1[0]     = 0.0;
    state_PS1[1]     = 1.0;
    step_counter[0]  = 0.0;
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T *act_PS1   = ssGetDWork(S, 0);
    real_T *state_PS1 = ssGetDWork(S, 1);

    real_T *y0 = ssGetOutputPortRealSignal(S, 0);
    real_T *y1 = ssGetOutputPortRealSignal(S, 1);

    /* Output activation and state vector */
    y0[0]     = act_PS1[0];
    y1[0]     = state_PS1[0];
    y1[1]     = state_PS1[1];
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{
    real_T *step_counter       = ssGetDWork(S, 2);
    real_T *act_PS1            = ssGetDWork(S, 0);
    real_T *state_PS1          = ssGetDWork(S, 1);

    /* Update counter */
    (*step_counter)++;

    /* Read inputs */
    real_T state_var_ps              = *ssGetInputPortRealSignal(S, 0);
    real_T failure_flag              = *ssGetInputPortRealSignal(S, 1);
    real_T tres_demand_ps            = *ssGetInputPortRealSignal(S, 2);
    const real_T *prev_state         = ssGetInputPortRealSignal(S, 3);
    real_T parent_last_time          = *ssGetInputPortRealSignal(S, 4);
    real_T t                         = *ssGetInputPortRealSignal(S, 5);
    real_T complementary_branch_flag = *ssGetInputPortRealSignal(S, 6);

    /* Compute flags */
    int ava_ps    = (failure_flag == 0) ? 1 : 0;
    int demand_ps = (state_var_ps == tres_demand_ps) ? 1 : 0;

    /* Activation logic */
    if (complementary_branch_flag == 1) {
        if (demand_ps && ava_ps) {
            *act_PS1 = 1.0;
            if ((demand_ps != prev_state[0] || ava_ps != prev_state[1]) && *step_counter >= 2) {
                ssSetStopRequested(S, 1);
            }
        } else {
            *act_PS1 = 0.0;
            if (demand_ps) {
                if ((demand_ps != prev_state[0] || ava_ps != prev_state[1]) && *step_counter >= 3) {
                    ssSetStopRequested(S, 1);
                }
            } else if (prev_state[0] == 1) {
                ssSetStopRequested(S, 1);
            }
        }
    } else {
        if (demand_ps) {
            *act_PS1 = ava_ps ? 1.0 : 0.0;
        } else {
            *act_PS1 = 0.0;
            if (prev_state[0] == 1) {
                ssSetStopRequested(S, 1);
            }
        }
    }

    /* Update state vector */
    state_PS1[0] = (real_T)demand_ps;
    state_PS1[1] = (real_T)ava_ps;
}

static void mdlTerminate(SimStruct *S) {}

#ifdef MATLAB_MEX_FILE
# include "simulink.c"
#else
# include "cg_sfun.h"
#endif
