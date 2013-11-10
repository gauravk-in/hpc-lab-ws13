/////////////////////////////////////////////
// CHANGE 1
/////////////////////////////////////////////

int bsp2_(double * restrict x, double * restrict y, double * restrict e, int n)
{
  int i;
  double r1;

    i = 0;
L10:
    /* Computing 3rd power */
    r1 = x[i] - y[i];
    e[i] = r1 * (r1 * r1);
    i++;
    if (i < n) {
	goto L10;
    }
  return 0;
} /* bsp2_ */

/////////////////////////////////////////////
// Rewritten as 
/////////////////////////////////////////////

// Converted to simple for loop, and removed use of r1. Loop vectorized.
int bsp2_(double * restrict x, double * restrict y, double * restrict e, int n)
{
  int i;

  for(i=0; i<n; i++)
  {
    e[i] = x[i] - y[i];
    e[i] *= e[i] * e[i];
  }

  return 0;
} /* bsp2_ */



/////////////////////////////////////////////
// CHANGE 2
/////////////////////////////////////////////

// The following loops were not vectorized. Possibly because of inefficient memory access pattern??
// TODO

int bsp3_(double * restrict x, int n, int nn)
{
  int i, l;

  l = 65;
  for (i=0; i < n; i++) { 
     x[i] = x[i + l];
  }
  return 0;
} /* bsp3_ */



int bsp4_(double * restrict x, int n, int nn)
{
  int i, l;

  l = -100;
  for (i = n ; i < nn; i++) {
     x[i] = x[i + l];
  }
  return 0;
} /* bsp4_ */




/////////////////////////////////////////////
// CHANGE 3
/////////////////////////////////////////////

// This loop was vectorized, but it is worth noting, that the compiler saw the values of l and m, before vectorizing the loop.
// While for l = 1 and m = 2, there are no dependencies, the compiler does not vectorize the loop when l = 2 and m = 1, 
// because that introduces vector dependency.
// Also note that for l = 3 and m = 1, there are no vector dependency and loop is vectorized.

int bsp5_(double * restrict y, double * restrict z, int n, int nn, int l, int m)
{
  int i;

  for (i=0; i< n; i++) {
     z[i + l] = z[i + m] + y[i];
  }
  return 0;
} /* bsp5_ */



/////////////////////////////////////////////
// CHANGE 4
/////////////////////////////////////////////

// PARTIAL LOOP WAS VECTORIZED.
int bsp6_(double * restrict x, double * restrict y, double * restrict z, int n)
{
  int i;

  for (i = 1; i < n; i++) {
     x[i] = y[i - 1];
     y[i] = z[i];
  }
  return 0;
} /* bsp6_ */

/////////////////////////////////////////////
// Rewritten as 
/////////////////////////////////////////////

int bsp6_(double * restrict x, double * restrict y, double * restrict z, int n)
{
  int i;

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




/////////////////////////////////////////////
// CHANGE 5
/////////////////////////////////////////////

// remark: loop was not vectorized: vectorization possible but seems inefficient.

int bsp7_(double * restrict x, double * restrict y, double * restrict e, int * restrict ix, int n)
{
  int i;

  for (i=0; i < n; i++) { 
       e[i] = x[ix[i]];
       e[i] = x[i] * y[i];
  }
  return 0;
} /* bsp7_ */

// first statement in the loop does not have any effect as e[i] is overwritten. removing it leads to vectorization.

int bsp7_(double * restrict x, double * restrict y, double * restrict e, int * restrict ix, int n)
{
  int i;

  for (i=0; i < n; i++) { 
       e[i] = x[i] * y[i];
  }
  return 0;
} /* bsp7_ */




/////////////////////////////////////////////
// CHANGE 6
/////////////////////////////////////////////

// remark: loop was not vectorized: vectorization possible but seems inefficient.
int bsp8_(double * restrict x, double * restrict y, double * restrict z,
          double * restrict e, int * restrict iy, int n, int nn)
{
  int i;

  #pragma vector always
  for (i=0; i < n; i++) { 
     x[i] = x[i * 2] + z[i + i];
     e[i] = y[i] + y[iy[i] * iy[i]];
  }
  return 0;
} /* bsp8_ */

// The loop was not vectorized possibly because the stride length for accessing y array, is not uniform
// therefore, the compiler can not efficiently prefetch the next value. The compiler estimates, that
// vectorization will be inefficient.

// We can force the compiler to vectorize using #pragma vector always.

// If we change the subscript of y, to something that can be analysed the loop gets vectorized.


/////////////////////////////////////////////
// CHANGE 7
/////////////////////////////////////////////

// Loop was not vectorized because it has vector dependency.
// The dependency can be removed by splitting the loop into 2 loops.
// However, the compiler does not vectorize it, because the stride length is not uniform.
// We can use #pragma vector always to make force it to vectorize.

int bsp9_(double * restrict x, double * restrict y, double * restrict z,
          int * restrict ix, int * restrict iy, int n)
{
  int i;

  for (i=0; i < n; i++) {
     x[i] = y[iy[i]] * 2.5;
     z[i] = y[ix[i]] * 3.;
     y[ix[i]] = 0.;
  }
  return 0;
} /* bsp9_ */

/////////////////////////////////////////////
// Rewritten as 
/////////////////////////////////////////////

int bsp9_(double * restrict x, double * restrict y, double * restrict z,
          int * restrict ix, int * restrict iy, int n)
{
  int i;

  #pragma vector always
  for (i=0; i < n; i++) {
     x[i] = y[iy[i]] * 2.5;
     z[i] = y[ix[i]] * 3.;
  }

  #pragma vector always
  for (i=0; i < n; i++) {
     y[ix[i]] = 0.;
  }

  return 0;
} /* bsp9_ */


/////////////////////////////////////////////
// CHANGE 8
/////////////////////////////////////////////

// This loop was not vectorized, because of a function call.
// If we change the function f_ to inline, the function is vectorized.
int bsp10_(double * restrict x, double * restrict e, int n)
{
  int i;

  for (i=0; i < n; i++) {
       e[i] = f_(x[i]);
  }
  return 0;
} /* bsp10_ */

//double f_(double x)
inline double f_(double x)
{
  double ret_val, r1;

  r1 = sin(x);
  ret_val = x - r1 * r1;
  return ret_val;
} /* f_ */

/////////////////////////////////////////////
// CHANGE 9
/////////////////////////////////////////////

// Compiler interchanged the loops for efficient memory access pattern, and vectorized the code automatically.

int bsp12_(double a[restrict][dim1], double b[restrict][dim1], int n)
{
  int i, j;

  for (j=0; j < n; j++) {
     for (i=0; i < n-1; i++) { 
        a[i+1][j] = a[i][j] * b[i][j]; 
     }
  }
  return 0;
} /* bsp12_ */



/////////////////////////////////////////////
// CHANGE 10
/////////////////////////////////////////////

// Compiler fused the two inner loops and vectorized the fused loop automatically.

int bsp13_(double * restrict x, double * restrict y, double * restrict z, 
           double a[restrict dim1][dim1], double b[restrict][dim1], 
           double * restrict s1, double * restrict s2, double * restrict s3, int n)
{
  int i, j, k;

  for (i=0; i < n; i++) {
     z[i] = (x[i] - y[i]) / *s1;
     x[i] = z[i] * *s2;
     for (k=0; k < n; k++) { 
        a[i][k] = b[i][k] + x[i];
     }
     for (j=0; j < n; j++) {
        a[i][j] = a[i][j] * x[i] - *s3;
     }
  }
  return 0;
} /* bsp13_ */




/////////////////////////////////////////////
// CHANGE 11
/////////////////////////////////////////////

int bsp14_(double * restrict x, double * restrict y, double * restrict z, 
           double a[restrict][dim1], double b[restrict][dim1],
           double * restrict s1, double * restrict s2, double * restrict s3, int n)
{
  int i, j, k;

  for (i=0; i < n; i++) {
     z[i] = (x[i] - y[i]) / *s1;
     x[i] = z[i] * *s2;
     // Loop was vectorized.
     for (k=0; k < n; k++) { 
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




/////////////////////////////////////////////
// CHANGE 12
/////////////////////////////////////////////

// Compiler interchanged the loops and vectorized the inner loop.

int bsp15_(double a[restrict][dim1], double b[restrict][dim1], 
           double c[restrict][dim1], int n)
{
  int i, j;
  int l=3, m=2;

/* ------------------------------------------ *
 * In which case is LOOP COLLAPSING possible? *
 * case 1: l=i and m=j   case 2: l=j and m=i  *           I 
 * ------------------------------------------ */

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

/////////////////////////////////////////////
// CHANGE 13
/////////////////////////////////////////////

// Loop was not vectorized, because the compiler assumed existence of output dependency.
// Based on the  logic of the application, we can force the compiler to ignore the possibility and vectorize the code.
// However, it will be inefficient because of the non uniform stride length for accessing x.
// We can also split the loop into two, but again the second loop will not get vectorized, and it will be less efficient.

int bsp16_(double * restrict x, double * restrict y, double * restrict z,
           int * restrict ix, int n)
{

  int i;

  for (i = 0; i < n; i++) {
     x[i] = y[i] + z[i];
     x[ix[i]] = 1.;
  }
  return 0;
} /* bsp16_ */

/////////////////////////////////////////////
// Rewritten as 
/////////////////////////////////////////////

int bsp16_(double * restrict x, double * restrict y, double * restrict z,
           int * restrict ix, int n)
{

  int i;

  for (i = 0; i < n; i++) {
     x[i] = y[i] + z[i];
  }

  #pragma vector always
  for (i = 0; i < n; i++) {
     x[ix[i]] = 1.;
  }
  return 0;
} /* bsp16_ */
