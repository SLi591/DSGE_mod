// addpath ("C:\Users\u477193\OneDrive - Tilburg University\dynare-6.3\matlab")

/* 
This file is used to use the Dynare's perfect solver to study the transition dynamics of nonstationary
Ramsey model. We assume there is only capital.  We study the effects of discount factor change: beta increases from 0.95 to 0.99.
It starts from the old steady state and transit to the new steady state.
*/

@#define simulation_periods=40
//****************************************************************************
//Define variables
//****************************************************************************

var C   ${c}$ (long_name='consumption')
    K   ${k}$ (long_name='capital stock')
    r
    ;

varexo beta  ${\beta}$ (long_name='discount factor');

parameters delta alpha;

//****************************************************************************
//Set parameter values
//****************************************************************************

alpha = 0.33;
delta = 0.1;

//****************************************************************************
//enter the model equations (model-block)
//****************************************************************************

model;
1 / C = beta * (1/C(+1) )* (1 + r(+1));
C + K - (1-delta) * K(-1) = K(-1)^alpha;
delta + r = alpha * K(-1)^(alpha-1);
end;

//****************************************************************************
// Initial Value: Old steady state
//****************************************************************************

initval;
beta = 0.95;
r = 1/beta-1;
K = ((delta + 1 / beta - 1) / alpha) ^ (1 / (alpha-1));
C = K ^ alpha - delta * K;
end;
steady;

resid;

//****************************************************************************
// End Value: New steady state
//****************************************************************************
endval;
beta = 0.99;
r = 1/beta - 1;
K = ((delta + 1 / beta - 1) / alpha) ^ (1 / (alpha-1));
C = K ^ alpha - delta * K;
end;
steady;

resid;

/* Use `check` after a `steady;` command (or after `initval;` if you are manually providing steady state values)
to check the local stability and determinacy conditions of the model around the current steady state. . 
*/
check;

//****************************************************************************
//perfect_foresight_solver: compute the solution
//****************************************************************************
perfect_foresight_setup(periods=@{simulation_periods});
perfect_foresight_solver;


rplot K;
rplot C;




