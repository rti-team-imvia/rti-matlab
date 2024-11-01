/*
 * File: main.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 29-Apr-2020 12:32:32
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "rt_nonfinite.h"
#include "LP_xyz2phitheta.h"
#include "main.h"
#include "LP_xyz2phitheta_terminate.h"
#include "LP_xyz2phitheta_initialize.h"

/* Function Declarations */
static void argInit_50x3_real_T(double result[150]);
static double argInit_real_T(void);
static void main_LP_xyz2phitheta(void);

/* Function Definitions */

/*
 * Arguments    : double result[150]
 * Return Type  : void
 */
static void argInit_50x3_real_T(double result[150])
{
  int idx0;
  double result_tmp;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 50; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_tmp = argInit_real_T();
    result[idx0] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 50] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 100] = argInit_real_T();
  }
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_LP_xyz2phitheta(void)
{
  double dv0[150];
  double dirs[100];

  /* Initialize function 'LP_xyz2phitheta' input arguments. */
  /* Initialize function input argument 'LP'. */
  /* Call the entry-point 'LP_xyz2phitheta'. */
  argInit_50x3_real_T(dv0);
  LP_xyz2phitheta(dv0, dirs);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  LP_xyz2phitheta_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_LP_xyz2phitheta();

  /* Terminate the application.
     You do not need to do this more than one time. */
  LP_xyz2phitheta_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
