#define S_FUNCTION_NAME  msfcn_indirectps_rom_v2
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

static char_T dworkNames[11][16];

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

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S)
{
    real_T *x;
    x=ssGetDWork(S,0); for(int j=0;j<5;j++) x[j]=1.0;
    x=ssGetDWork(S,1); x[0]=1.0;
    x=ssGetDWork(S,2); for(int j=0;j<5;j++) x[j]=0.0;
    x=ssGetDWork(S,3); x[0]=1.0;
    x=ssGetDWork(S,4); x[0]=0.0;
    x=ssGetDWork(S,5); x[0]=0.0;
    x=ssGetDWork(S,6); x[0]=0.0;
    x=ssGetDWork(S,7); for(int j=0;j<10;j++) x[j]=0.0;
    x=ssGetDWork(S,8); x[0]=0.0;
    x=ssGetDWork(S,9); x[0]=0.0;
    x=ssGetDWork(S,10);x[0]=0.0;

}

static void mdlOutputs(SimStruct *S,int_T tid)
{
    memcpy(ssGetOutputPortRealSignal(S, 0), ssGetDWork(S, 0), 5*sizeof(real_T));
    memcpy(ssGetOutputPortRealSignal(S, 1), ssGetDWork(S, 1), 1*sizeof(real_T));
    memcpy(ssGetOutputPortRealSignal(S, 2), ssGetDWork(S, 2), 5*sizeof(real_T));
    memcpy(ssGetOutputPortRealSignal(S, 3), ssGetDWork(S, 3), 1*sizeof(real_T));
    memcpy(ssGetOutputPortRealSignal(S, 4), ssGetDWork(S, 4), 1*sizeof(real_T));
    memcpy(ssGetOutputPortRealSignal(S, 5), ssGetDWork(S, 5), 1*sizeof(real_T));
    memcpy(ssGetOutputPortRealSignal(S, 6), ssGetDWork(S, 7),10*sizeof(real_T));  // FIXED
    memcpy(ssGetOutputPortRealSignal(S, 7), ssGetDWork(S, 8), 1*sizeof(real_T));  // FIXED
    memcpy(ssGetOutputPortRealSignal(S, 8), ssGetDWork(S, 9), 1*sizeof(real_T));  // FIXED
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{
    /* D-Work pointers */
    real_T *current_sys_state = ssGetDWork(S, 0);      /* current_sys_state[5] */
    real_T *current_state_comp = ssGetDWork(S, 1);     /* current state of transitioning component */
    real_T *trans_sys = ssGetDWork(S, 2);              /* future_sys_state[5] */
    real_T *new_trans = ssGetDWork(S, 3);              /* future state of transitioning component */
    real_T *new_ttf = ssGetDWork(S, 4);                /* future_ttf */
    real_T *prev_ttf = ssGetDWork(S, 5);               /* prev_ttf */
    real_T *failure_flag = ssGetDWork(S, 6);           /* failure_flag */
    real_T *mat_ttf = ssGetDWork(S, 7);                /* mat_ttf[10] */
    real_T *trans_c = ssGetDWork(S, 8);                /* transitioning component */
    real_T *comp_branch_event = ssGetDWork(S, 9);      /* comp_branch_event */
    real_T *step_counter = ssGetDWork(S,10);           /* step_counter */

    /* Read inputs */
    real_T t = *ssGetInputPortRealSignal(S,0);
    int Ncomps = (int)*ssGetInputPortRealSignal(S,1);
    int Nstates = (int)*ssGetInputPortRealSignal(S,2);
    const real_T *mat_states_index = ssGetInputPortRealSignal(S,3);
    const real_T *mat_trans = ssGetInputPortRealSignal(S,4);
    const real_T *failed_states = ssGetInputPortRealSignal(S,5);
    real_T complementary_branch_flag = *ssGetInputPortRealSignal(S,8);

    const real_T tsim = 1400.0;
    const real_T thres_sampled_time = 2.5;

    /* increment step counter and check complementary branch event */
    step_counter[0] += 1;
    comp_branch_event[0] = (complementary_branch_flag == 1 && step_counter[0] <= 1) ? 1 : 0;
    failure_flag[0] = 0;  // default: healthy

    if (new_ttf[0] == 0) {
        for (int c = 0; c < Ncomps; c++) {
            int current_state = (int)current_sys_state[c] - 1;
            for (int j = 0; j < Nstates; j++) {
                if ((j != current_state) && (mat_states_index[c*Nstates+j] != 0)) {
                    real_T lambda = mat_trans[c*4 + current_state*2 + j];
                    real_T ttf = t - log((real_T)rand() / (real_T)RAND_MAX) / lambda;
                    mat_ttf[c*Nstates+j] = ttf;
                } else {
                    mat_ttf[c*Nstates+j] = 0.0;
                }
            }
        }

        /* find min TTF */
        real_T min_ttf = 1e10;
        int min_c = 0, min_j = 0;
        for (int c = 0; c < Ncomps; c++) {
            for (int j = 0; j < Nstates; j++) {
                real_T val = mat_ttf[c*Nstates+j];
                if (val > 0 && val < min_ttf) {
                    min_ttf = val;
                    min_c = c;
                    min_j = j;
                }
            }
        }

        prev_ttf[0] = t;
        new_ttf[0] = min_ttf;

        /* threshold adjustment */
        if (min_ttf > (tsim - thres_sampled_time) && min_ttf < tsim) {
            for (int k=0; k<Ncomps*Nstates; k++)
                if (mat_ttf[k] > 0) mat_ttf[k] += thres_sampled_time;
            new_ttf[0] += thres_sampled_time;
        }

        for (int i = 0; i < Ncomps; i++)
            trans_sys[i] = current_sys_state[i];

        trans_sys[min_c] = min_j+1;
        new_trans[0] = min_j+1;
        trans_c[0] = min_c+1;
    }
    else if (comp_branch_event[0] == 0 && t >= new_ttf[0]) {
        for (int i = 0; i < Ncomps; i++)
            current_sys_state[i] = trans_sys[i];

        current_state_comp[0] = new_trans[0];

        failure_flag[0] = 0;
        for (int i = 0; i < Ncomps; i++) {
            int state = (int)current_sys_state[i];
            if (failed_states[i*Nstates+state-1] == 1) {
                failure_flag[0] = 1;
                break;
            }
        }

        prev_ttf[0] = t;
        new_ttf[0] = 0.0;
    }
}

static void mdlTerminate(SimStruct *S) {}

/* Required for Simulink to handle operating points and logged states */
#define MDL_SIM_STATE
static mxArray* mdlGetSimState(SimStruct* S)
{
    mxArray *simState = mxCreateDoubleMatrix(1,0,mxREAL);
    return simState;
}

static void mdlSetSimState(SimStruct* S, const mxArray* simState){}

#ifdef MATLAB_MEX_FILE
# include "simulink.c"
#else
# include "cg_sfun.h"
#endif
