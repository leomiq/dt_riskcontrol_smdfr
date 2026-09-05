#define S_FUNCTION_NAME  msfcn_limintm_rom_v2
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

/* Number of ports and DWorks */
#define NINPUTS   4
#define NOUTPUTS  5
#define NDWORKS   2

/* Function: mdlInitializeSizes ===============================================*/
static void mdlInitializeSizes(SimStruct *S)
{
    /* Parameters */
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))
        return;

    /* States: no continuous, NDWORKS discrete */
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, NDWORKS);

    /* Input ports */
    if (!ssSetNumInputPorts(S, NINPUTS)) return;
    ssSetInputPortWidth(S, 0, 1);
    ssSetInputPortWidth(S, 1, 2);
    ssSetInputPortWidth(S, 2, 1);
    ssSetInputPortWidth(S, 3, 1);
    for (int i = 0; i < NINPUTS; i++) {
        ssSetInputPortDirectFeedThrough(S, i, 1);
        ssSetInputPortRequiredContiguous(S, i, 1);
    }

    /* Output ports */
    if (!ssSetNumOutputPorts(S, NOUTPUTS)) return;
    ssSetOutputPortWidth(S, 0, 1);
    ssSetOutputPortWidth(S, 1, 1);
    ssSetOutputPortWidth(S, 2, 1);
    ssSetOutputPortWidth(S, 3, 1);
    ssSetOutputPortWidth(S, 4, 2);

    /* DWork vectors */
    ssSetNumDWork(S, NDWORKS);
    for (int i = 0; i < NDWORKS; i++) {
        ssSetDWorkWidth(S, i, 1);
        ssSetDWorkDataType(S, i, SS_DOUBLE);
        ssSetDWorkUsageType(S, i, SS_DWORK_USED_AS_DSTATE);
    }
    ssSetDWorkName(S, 0, "Slim_1");
    ssSetDWorkName(S, 1, "Slim_2");

    /* Sample times and work vectors */
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    /* SimState and options */
    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

/* Function: mdlInitializeSampleTimes =========================================*/
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

/* Function: mdlInitializeConditions ==========================================*/
static void mdlInitializeConditions(SimStruct *S)
{
    /* Zero all discrete states (DWorks) */
    for (int i = 0; i < NDWORKS; i++) {
        real_T *d = ssGetDWork(S, i);
        *d = 0.0;
    }
}

/* Function: mdlOutputs =======================================================*/
static void mdlOutputs(SimStruct *S, int_T tid)
{
    /* Input pointers */
    InputRealPtrsType u0 = ssGetInputPortRealSignalPtrs(S, 0);
    InputRealPtrsType u1 = ssGetInputPortRealSignalPtrs(S, 1);
    InputRealPtrsType u2 = ssGetInputPortRealSignalPtrs(S, 2);
    InputRealPtrsType u3 = ssGetInputPortRealSignalPtrs(S, 3);

    /* Output pointers */
    real_T *y0 = ssGetOutputPortRealSignal(S, 0);
    real_T *y1 = ssGetOutputPortRealSignal(S, 1);
    real_T *y2 = ssGetOutputPortRealSignal(S, 2);
    real_T *y3 = ssGetOutputPortRealSignal(S, 3);
    real_T *y4 = ssGetOutputPortRealSignal(S, 4);

    /* Discrete states */
    real_T *dwork0 = ssGetDWork(S, 0);
    real_T *dwork1 = ssGetDWork(S, 1);

    /* Simple pass-throughs and logic */
    y0[0]    = **u0;
    y1[0]    = *dwork0;
    y2[0]    = *dwork1;
    y3[0]    = (**u2 <= 5.0) ? 0.0 : 1.0;

    /* Lower-bound forward */
    y4[0]    = **u3;
    y4[1]    = 0.0;
}

/* Function: mdlUpdate ========================================================*/
static void mdlUpdate(SimStruct *S, int_T tid)
{
    InputRealPtrsType u0 = ssGetInputPortRealSignalPtrs(S, 0);
    InputRealPtrsType u1 = ssGetInputPortRealSignalPtrs(S, 1);
    InputRealPtrsType u3 = ssGetInputPortRealSignalPtrs(S, 3);

    real_T *dwork0 = ssGetDWork(S, 0);
    real_T *dwork1 = ssGetDWork(S, 1);

    real_T u      = **u0;
    real_T lb_ps  = **u3;
    int_T  nThresh = ssGetInputPortWidth(S, 1);
    real_T states_max = *dwork1;

    for (int i = 0; i < nThresh; i++) {
        real_T lb   = (i == 0) ? 0.9995 : (*u1[i-1]);
        real_T ub   = *u1[i];
        if (u >= lb && u <= ub) {
            *dwork0 = i;
            if (states_max < i) states_max = i;
        }
        if (i == nThresh-1 && u > ub) {
            *dwork0 = i+1;
            if (states_max < i+1) states_max = i+1;
            ssSetStopRequested(S, 1);
        }
    }

    if (u < lb_ps && states_max > 0) {
        states_max = 0;
    }
    *dwork1 = states_max;
}

/* Function: mdlTerminate =====================================================*/
static void mdlTerminate(SimStruct *S)
{
    /* Nothing to do */
}

#ifdef MATLAB_MEX_FILE
# include "simulink.c"
#else
# include "cg_sfun.h"
#endif
