/* vectorisation quiz
 * for practical course
 * ----------------------
 * autor: Markus Brenk
 * date: 2003-10-21
 * ================================================== */


#include <stdio.h>
#include <math.h>

#define dim1 100
#define dim2 200

static int bsp1_(double *, double *, double *, int);
static int bsp2_(double *, double *, double *, int);
static int bsp3_(double *, int, int);
static int bsp4_(double *, int , int );
static int bsp5_(double *, double *, int, int, int, int);
static int bsp6_(double *, double *, double *, int);
static int bsp7_(double *, double *, double *, int *, int);
static int bsp8_(double *, double *, double *, double *, int *, int, int);
static int bsp9_(double *, double *, double *, int *, int *, int);
static int bsp10_(double *, double *, int);
static int bsp11_(double a[][dim1], int);
static int bsp12_(double a[][dim1], double b[][dim1], int);
static int bsp13_(double*,double*,double*,double a[][dim1],double b[][dim1],double*,double*,double*,int);
static int bsp14_(double*,double*,double*,double a[][dim1],double b[][dim1],double*,double*,double*,int);
static int bsp15_(double a[][dim1], double b[][dim1], double c[][dim1], int);
static int bsp16_(double *, double *, double *, int *, int);
static int bsp17_(double *, double *, int);
static int bsp18_(double *, double *, int);
static int bsp19_(double *, double *, int);

static double f_(double);

/* ---------------------------------------------------------------------- *\
 *  The functions bsp1_() up to bsp19_() contain loops.                   * 
 *  Is it possible to optimize the different loops with SSE2 pipelining?  * 
 *  What do you expect?                                                   * 
 *  What does the compiler do?                                            * 
 *  Compile on an intel architecture with the intel compiler              *
 *  Use the following commands:                                           * 
 *  > icc -std=c99 -restrict -O3 -ip -ip_no_inlining\     * 
 *        -vec_report3 -opt_report_filevector.opt -opt_report_phaseall\   *
 *        vector.c                                                        *
 *  What is the difference?                                               *
 *  Can you help the compiler with pragmas or restructuring the code?     *
\* ---------------------------------------------------------------------- */


/* Main program */ 
int main()
{
    int n1 = dim1;
    int n2 = dim2;
    int i, j, l, m;
    int ix[dim1], iy[dim1];
    double a[dim1][dim1], b[dim1][dim1],c[dim1][dim1],e[dim1];
    double x[dim2], y[dim2], z[dim2], s1, s2, s3;

/* initialisation */

	// Loop was vectorized.
	for (i = 0; i < n2; i++)
      {
	x[i] = 1.;
	y[i] = 2.;
	z[i] = 3.;
      }
   
	
	// The compiler says this loop was not vectorized, because it is not an inner loop.
	// I found out that the compiler is merging this loop with the next loop, and
	// vectorizing only the inner loop of the next loop. Change the number of iterations of
	// either of the loop, and you will see that this loop is vectorized by the compiler.
    for (i=0; i < n1; i++) 
      {
	ix[i] = n1 - i;
	iy[i] = (int)(i - i * .9);
      }
    
    for (i=0; i < n1; i++) 
      {
		 // Loop was vectorized.
         for (j=0; j < n1; j++) 
           {
	     a[i][j] = 1.;
	     b[i][j] = 2.;
	     c[i][j] = 3.;
           }
      }
    s1 = 1.;
    s2 = 2.;
    s3 = 3.;



    bsp1_(x, y, e, n1);
    bsp2_(x, y, e, n1);
    bsp3_(x, n1, n2);
    bsp4_(x, n1, n2);
    l = 1;
    m = 2;
    /* -- or -- */ 
    /* l = 2; 
       m = 1; */
    bsp5_(y, z, n1, n2, l, m);

	// Loop was vectorized.
    for (i=0; i < n2; i++) 
      {
	z[i] = 3.;
      }
    bsp6_(x, y, z, n1);

	// Loop was vectorized.
    for (i=0; i<n2; i++) 
      {
	x[i] = 1.;
	y[i] = 2.;
      }
    bsp7_(x, y, e, ix, n1);
    bsp8_(x, y, z, e, iy, n1, n2);

	// Loop was vectorized.
    for (i=0; i<n2; i++) 
      {
	x[i] = 1.;
      }
    bsp9_(x, y, z, ix, iy, n1);
	
	// Loop was vectorized.
    for (i=0; i<n2; i++) 
      {
	x[i] = 1.;
	y[i] = 2.;
	z[i] = 2.;
      }
    bsp10_(x, e, n1);
    bsp11_(a, n1);

    for (i=0; i<n1; i++) 
      {
		 // Loop was vectorized.
         for (j=0; j<n1; j++) 
           {
	      a[i][j] = 1.;
           }
      }
    bsp12_(a, b, n1);
      
    for (i=0; i<n1; i++) 
      {
      	 // Loop was vectorized.
         for (j=0; j<n1; j++) 
           {
              a[i][j] = 1.;
           }
      }
    bsp13_(x, y, z, a, b, &s1, &s2, &s3, n1);

	// Loop was vectorized.
    for (i=0; i <n1; i++) 
      {
         x[i] = 1.;
         z[i] = 3.;
         for (j=0; j<n1; j++) 
           {
              a[i][j] = 1.;
           }
      }
    bsp14_(x, y, z, a, b, &s1, &s2, &s3, n1);

	// Loop was vectorized.
    for (i=0; i<n1; i++) 
      {
         x[i] = 1.;
         z[i] = 3.;
         for (j=0; j <n1; j++) 
           {
              a[i][j] = 1.;
           }
      }
    bsp15_(a, b, c, n1);
    bsp16_(x, y, z, ix, n1);
    bsp17_(x, &s1, n1);

	// Loop was vectorized.
    for (i=0; i<50; i++) 
      {
         x[i] = 1.;
      }
    bsp18_(x, y, n1);
    bsp19_(x, y, n1);


} /* end of main */



int bsp1_(double * restrict x, double * restrict y, double* restrict e, int n)
{
  int i;

  // Not vectorized due to if condition
  // TODO: Some magic could possibly be done to remove the if.
  for (i=0; i< n; i++) { 
    if (x[i] > 0.) { 
       e[i] = x[i] + y[i];
    } else {
       e[i] = x[i] - y[i];
    }
  }
  return 0;
} /* bsp1_ */




int bsp2_(double * restrict x, double * restrict y, double * restrict e, int n)
{
  int i;
  double r1;

// This section can be converted into a loop, as follows.
    // i = 0;
// L10: 
    // /* Computing 3rd power */
    // r1 = x[i] - y[i];
    // e[i] = r1 * (r1 * r1);
    // i++;
    // if (i < n) {
	// goto L10;
    // }

	// Apart from writing the code as a loop, the use of r1 has been removed as
	// use of a variable which is shared between iterations will not allow vectorization.
	for(i=0; i<n; i++)
	{
		e[i] = x[i] - y[i];
		e[i] *= e[i] * e[i];
	}
	
  return 0;
} /* bsp2_ */



int bsp3_(double * restrict x, int n, int nn)
{
  int i, l;

  l = 65;
  // Seems like not vectorized because of inefficient memory access pattern
  for (i=0; i < n; i++) { 
     x[i] = x[i + l];
  }
  return 0;
} /* bsp3_ */



int bsp4_(double * restrict x, int n, int nn)
{
  int i, l;

  l = -100;
  // Seems like not vectorized because of inefficient memory access pattern
  for (i = n ; i < nn; i++) {
     x[i] = x[i + l];
  }
  return 0;
} /* bsp4_ */



int bsp5_(double * restrict y, double * restrict z, int n, int nn, int l, int m)
{
  int i;

  // Loop was vectorized.
  // Worth noting, that the compiler checks the values of l and m.
  // While the loop is vectorized with l = 1 and m = 2, the compiler
  // does not vectorize with l = 2 and m = 1, since that introduces
  // vector dependence.
  for (i=0; i< n; i++) { 
     z[i + l] = z[i + m] + y[i];
  }
  return 0;
} /* bsp5_ */



int bsp6_(double * restrict x, double * restrict y, double * restrict z, int n)
{
  int i;

  // PARTIAL LOOP WAS VECTORIZED. This loop can be rewritten to improve vectorization.
  // for (i = 1; i < n; i++) { 
     // x[i] = y[i - 1];
     // y[i] = z[i];
  // }
  
  // Writing the loop as follows, improves the vectorization of the code.
  x[1] = y[0];
  for (i=1; i < n-1; i++)
  {
  	y[i] = z[i];
  	x[i + 1] = y[i];
  }
  y[n] = z[n];
  
  return 0;
} /* bsp6_ */



int bsp7_(double * restrict x, double * restrict y, double * restrict e, int * restrict ix, int n)
{
  int i;

  for (i=0; i < n; i++) {
       // This statement has no effect, as e[i] will be over written.
  	   // After this change, the loop is vectorized by the compiler.
       // e[i] = x[ix[i]];
       e[i] = x[i] * y[i];
  }
  return 0;
} /* bsp7_ */



int bsp8_(double * restrict x, double * restrict y, double * restrict z,
          double * restrict e, int * restrict iy, int n, int nn)
{
  int i;

// The loop was not vectorized possibly because the stride length for accessing y array, is not uniform
// therefore, the compiler can not efficiently prefetch the next value. The compiler estimates, that
// vectorization will be inefficient.

// If we change the subscript of y, to something that can be analysed the loop gets vectorized.

  for (i=0; i < n; i++) { 
     x[i] = x[i * 2] + z[i + i];
     e[i] = y[i] + y[iy[i] * iy[i]];
  }
  return 0;
} /* bsp8_ */



int bsp9_(double * restrict x, double * restrict y, double * restrict z,
          int * restrict ix, int * restrict iy, int n)
{
  int i;

// Loop was not vectorized because it has vector dependency.
// The dependency can be removed by splitting the loop into 3 loops.
// However, the compiler does not vectorize it, because the stride length is not uniform.
  for (i=0; i < n; i++) {
     x[i] = y[iy[i]] * 2.5;
     z[i] = y[ix[i]] * 3.;
     y[ix[i]] = 0.;
  }
  return 0;
} /* bsp9_ */



int bsp10_(double * restrict x, double * restrict e, int n)
{
  int i;

  for (i=0; i < n; i++) { 
       e[i] = f_(x[i]);
  }
  return 0;
} /* bsp10_ */

// Changed the function to inline, so as to vectorize the loop in bsp10_
inline double f_(double x) 
{
  double ret_val, r1;

  r1 = sin(x);
  ret_val = x - r1 * r1;
  return ret_val;
} /* f_ */



int bsp11_(double a[restrict][dim1], int n)
{
  int i, j, ip1;

  for (i=0; i < n; i++) { 
     ip1 = i + 1;
     // Loop was vectorized.
     for (j = ip1; j < n; j++) { 
        a[i][j] /= a[i][i];
     }
  }

  /* for (j=0; j < n; j++) { 
     jp1 = j + 1;
     for (i = jp1; i < n; i++) {
        a[i][j] /= a[j][j];
     }
     }*/
  return 0;
} /* bsp11_ */



int bsp12_(double a[restrict][dim1], double b[restrict][dim1], int n)
{
  int i, j;

  // loop interchanged and vectorized.
  for (j=0; j < n; j++) {
     for (i=0; i < n-1; i++) { 
        a[i+1][j] = a[i][j] * b[i][j]; 
     }
  }
  return 0;
} /* bsp12_ */




int bsp13_(double * restrict x, double * restrict y, double * restrict z, 
           double a[restrict dim1][dim1], double b[restrict][dim1], 
           double * restrict s1, double * restrict s2, double * restrict s3, int n)
{
  int i, j, k;

  for (i=0; i < n; i++) {
     z[i] = (x[i] - y[i]) / *s1;
     x[i] = z[i] * *s2;
     // Loops were fused and vectorized.
     for (k=0; k < n; k++) { 
        a[i][k] = b[i][k] + x[i];
     }
     for (j=0; j < n; j++) {
        a[i][j] = a[i][j] * x[i] - *s3;
     }
  }
  return 0;
} /* bsp13_ */



int bsp14_(double * restrict x, double * restrict y, double * restrict z, 
           double a[restrict][dim1], double b[restrict][dim1],
           double * restrict s1, double * restrict s2, double * restrict s3, int n)
{
  int i, j, k;

  for (i=0; i < n; i++) {
     z[i] = (x[i] - y[i]) / *s1;
     x[i] = z[i] * *s2;
     for (k=0; k < n; k++) { // Loop was vectorized.
        a[i][k] = b[i][k] + x[i];
     }
     // This loop was not fused with the above loop, because we access x[i] in above loop and x[j] in the following loop.
     // This may lead to dependency. If both had x[i] or above loop had x[k], the compiler would fuse the two loops together and vectorize.
     for (j=0; j < n; j++) {
        a[i][j] = a[i][j] * x[j] - *s3;
     }
  }
  return 0;
} /* bsp14_ */



int bsp15_(double a[restrict][dim1], double b[restrict][dim1], 
           double c[restrict][dim1], int n)
{
  int i, j;
  int l=3, m=2;

/* ------------------------------------------ *
 * In which case is LOOP COLLAPSING possible? *
 * case 1: l=i and m=j   case 2: l=j and m=i  *           I 
 * ------------------------------------------ */

  // Loops were interchanged and vectorized.
  for (j=0; j < n; j++) {
     /* m = j; */ 
     /* m = i; */
     for (i=0; i < n; i++) { 
        /* l = i; */
        /* l = j; */
        a[i][j] = b[i][j] + c[i][j] * b[l][m];
     }
  }
  return 0;
} /* bsp15_ */



int bsp16_(double * restrict x, double * restrict y, double * restrict z,
           int * restrict ix, int n)
{

  int i;

  // Output dependency assumed. We can force the compiler to ignore dependency and vectorize,
  // but it will be inefficient because of the non-uniform stride length for accessing array x.
  // we can also split the loop into two, but it will be less efficient as second loop will not
  // be vectorized efficiently.
  for (i = 0; i < n; i++) { 
     x[i] = y[i] + z[i];
     x[ix[i]] = 1.;
  }
  return 0;
} /* bsp16_ */



int bsp17_(double * restrict x, double * restrict s, int n)
{

  int i;


  *s = 0.;
  // Loop was vectorized.
  for (i = 0; i < n; i++) { 
     *s += x[i];
  }
  return 0;
} /* bsp17_ */



int bsp18_(double * restrict x, double * restrict y, int n)
{
  int i;

  for (i = 0; i < n; i++) { // Not vectorized due to branch
     y[i] = 1.;
     if (x[i] < 0.) break;
  }
  return 0;
} /* bsp18_ */


int bsp19_(double * restrict x, double * restrict y, int n)
{
  int i;

  for (i = 0; i < n; i++) { // Not vectorized due to branch
     if (x[i] < 0.) break;
     x[i+1] -=  y[i];
  }
  return 0;
} /* bsp19_ */
