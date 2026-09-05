#define S_FUNCTION_NAME  msfcn_limintm_hrm_c
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>
#include <string.h>

/*
 * Level-2 S-Function: msfcn_limintm_hrm_c
 *  - Echoes input port 0 to output port 0
 *  - Maintains two discrete D-Work states for internal logic
 *  - Outputs:
 *      y0: echo of input u
 *      y1: xD (current state)
 *      y2: xD_2 (max state reached)
 *      y3: flag based on 'reloj'
 *      y4[2]: mirrors 2-element input 'a'
 */

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))
        return;

    /* No continuous or discrete Simulink states */
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    /* Inputs: 4 ports */
    if (!ssSetNumInputPorts(S, 4)) return;
    ssSetInputPortWidth(S, 0, 1);
    ssSetInputPortWidth(S, 1, 2);
    ssSetInputPortWidth(S, 2, 1);
    ssSetInputPortWidth(S, 3, 1);
    for (int i = 0; i < 4; i++) {
        ssSetInputPortDirectFeedThrough(S, i, 1);
        ssSetInputPortRequiredContiguous(S, i, 1);
    }

    /* Outputs: 5 ports */
    if (!ssSetNumOutputPorts(S, 5)) return;
    ssSetOutputPortWidth(S, 0, 1);
    ssSetOutputPortWidth(S, 1, 1);
    ssSetOutputPortWidth(S, 2, 1);
    ssSetOutputPortWidth(S, 3, 1);
    ssSetOutputPortWidth(S, 4, 2);

    /* Two D-Work vectors */
    ssSetNumDWork(S, 2);
    for (int i = 0; i < 2; i++) {
        ssSetDWorkWidth(S, i, 1);
        ssSetDWorkDataType(S, i, SS_DOUBLE);
        ssSetDWorkUsageType(S, i, SS_DWORK_USED_AS_DSTATE);
    }

    /* Name the D-Works so they appear in loggedStates */
    ssSetDWorkName(S, 0, "Slim_1");  /* last input value (u) */
    ssSetDWorkName(S, 1, "Slim_2");  /* maximum state reached */

    /* Sample times */
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S)
{
    real_T *xD   = ssGetDWork(S, 0);
    real_T *xD_2 = ssGetDWork(S, 1);
    xD[0]   = 0.0;
    xD_2[0] = 0.0;
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    /* Read inputs (contiguous) */
    const real_T *u      = ssGetInputPortRealSignal(S, 0);
    const real_T *a      = ssGetInputPortRealSignal(S, 1);
    const real_T *reloj  = ssGetInputPortRealSignal(S, 2);

    /* D-Work states */
    real_T *xD   = ssGetDWork(S, 0);
    real_T *xD_2 = ssGetDWork(S, 1);

    /* Outputs */
    real_T *y0 = ssGetOutputPortRealSignal(S, 0);
    real_T *y1 = ssGetOutputPortRealSignal(S, 1);
    real_T *y2 = ssGetOutputPortRealSignal(S, 2);
    real_T *y3 = ssGetOutputPortRealSignal(S, 3);
    real_T *y4 = ssGetOutputPortRealSignal(S, 4);

    /* Echo and expose states */
    y0[0]     = u[0];
    y1[0]     = xD[0];
    y2[0]     = xD_2[0];
    y3[0]     = (reloj[0] > 5.0) ? 1.0 : 0.0;
    y4[0]     = a[0];
    y4[1]     = a[1];

    /* Debug print */
    //ssPrintf("mdlOutputs: u=%f, state=%f, maxstate=%f, flag=%f\n",
    //         u[0], xD[0], xD_2[0], y3[0]);
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{
    /* Read input signals (all marked contiguous) */
    const real_T  u     = *ssGetInputPortRealSignal(S, 0);
    const real_T *a     = ssGetInputPortRealSignal(S, 1);
    const real_T  lb_ps = *ssGetInputPortRealSignal(S, 3);

    /* Pointers to D-Work storage */
    real_T *xD   = ssGetDWork(S, 0);  /* holds current state index */
    real_T *xD_2 = ssGetDWork(S, 1);  /* holds max state reached */

    /* Load previous max-state */
    int_T states_max = (int_T)xD_2[0];
    int_T new_state  = 0;

    /* Compute bounds from the 2-element vector `a` */
    real_T ub = (a[0] < a[1]) ? a[0] : a[1];
    real_T lb = 0.9995;

    /*— Transition logic —*/
    if (u <= a[0]) {
        new_state = 0;
    }
    else if (u > a[0] && u <= a[1]) {
        new_state = 1;
    }
    else if (u > a[1]) {
        new_state = 2;
        /* Terminal condition: activation of end‐of‐protection */
        ssSetStopRequested(S, 1);
    }
    else {
        new_state = 0;
    }

    /* Store the new state index */
    xD[0] = (real_T)new_state;

    /* Update max‐state seen so far */
    states_max = (new_state > states_max) ? new_state : states_max;

    /* Reset condition: if `u` falls below lb_ps after entering any state */
    if (u < lb_ps && states_max > 0) {
        states_max = 0;
        xD[0]      = 0.0;
        /* Terminal condition: deactivation of protection */
        //ssSetStopRequested(S, 1);
    }

    /* Write back the updated max‐state */
    xD_2[0] = (real_T)states_max;
}


static void mdlTerminate(SimStruct *S) {}

#ifdef MATLAB_MEX_FILE
# include "simulink.c"
#else
# include "cg_sfun.h"
#endif
