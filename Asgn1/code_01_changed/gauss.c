/* row wise gaus algorithm
 * pattern for practical course
 * -------------------------
 * autor: Markus Brenk
 * date: 2002-09-25
 * =================================================== */

//#define NGLS 32768
//NGLS changed to 2000 for 2000 different right-hand sides
#define NGLS 2000

#define NUM_LOOPS 4

#include <stdio.h>
#include <math.h>
#include "timer.h"

/*
	Version containing optimizations for vectorization.
*/

/* print a 3x3 matrix */
void print_matrix(char* name, double matrix[3][3]);


/* print a 3d vector */
void print_vector(char* name, double vec[3]);


/**
 *  initialisation: generates the following LGS:
 *  ( 3  1  1)       (5)                         (1)
 *  ( 1  4  1) * X = (6)        => solution  X = (1)
 *  ( 1  1  5)       (7)                         (1)
 */
void init(double a[3][3], double b[3], double x[3], int os);

/** performs gauss elimination */
void gauss_elimination(double a[3][3], double b[3], double x[3]);


double a[NGLS][3][3];
double b[NGLS][3];
double x[NGLS][3];

int main() {

	int i, j;
	int n=3;
	char * name;
	char nm;
	
	name = &nm;

	time_marker_t time;

	for (i = 0; i < NGLS; i++) {
		init(a[i], b[i], x[i], i);
	}

	int flops;
	//flops = NGLS*0.5*n*(n^4+3*n^2-4*n+2);
	flops = NGLS*(0.5*(n*n*n*n) + 2*n*n - 1.5*n);
	flops = flops/1000; //KFLOPS
	
	printf("Input data\n");
	nm = 'A';
	print_matrix(name, a[NGLS-1]);
	nm = 'b';
	print_vector(name, b[NGLS-1]);
	printf("Variable flops/1000: %d \n", flops);
	
	time = get_time();
	for (i = 0; i < NGLS; i++) {
		gauss_elimination(a[i], b[i], x[i]);
	}

	print_flops(flops, time);
	printf("NAIV: Time elapsed. time: %f   ticks: %f\n", get_ToD_diff_time(time), get_ticks_diff_time(time));
	//printf("KFLOPS calc: %e \n", (double) flops/get_ToD_diff_time(time));
	
	//new output
	printf("Output data\n");
	nm = 'A';
	print_matrix(name, a[NGLS-1]);
	nm = 'x';
	print_vector(name, x[NGLS-1]);
		
	return(0);
}


//parameter os helps create different right hand sides
void init(double a[3][3], double b[3], double x[3], int os) {
	int i,j;
	int n = 3;

	for (j=0;j<n;j++)	//outer loop - cannot be vectorized 
	{
		b[j]=(float)(2*n-2)+(float)(j+1)+(float)(os)/10;
		a[j][j]=(float)(n-1)+(float)(j+1); //a[0][0]=3; a[1][1]=4; a[2][2]=5;
		x[j]=0.;

		for (i=j+1;i<n;i++) {	//loop was unrolled, vectorized
			a[i][j]=1.;	
		}
		for (i=0;i<j;i++) {		//loop was unrolled, vectorized
			a[i][j]=1.;
		}
	}
}

/*	
	The algorithm brings the matrix of coefficients to reduced row echelon form.
	It starts by first dividing the first row of the matrix by the value
	of the leading (non-zero) coefficient, so that the leading coefficient's value is set to equal one.
	Then, it adds a scalar multiple of the first row to all subsequent rows,
	so that the column of the leading coefficient contains only zeroes in all subsequent rows.
	It then proceeds by doing the same thing for all other rows in the matrix.
	The leading coefficient that gets selected for each row is an element on the main diagonal of the matrix.
	In the end, the matrix of coefficients A is transformed into a right triangular matrix 
	with the additional property of its main diagonal containing only the value of 1 (reduced row echelon form).
*/
void gauss_elimination(double a[3][3], double b[3], double x[3]) {
	int n = 3;
	int i,j,k;

	//outer for-loop selects current row 
	for (i = 0; i < n; i++) {	//outer loop - no vectorization

		//temporary variable to store row divisor
		double div = a[i][i];
	
		//current row is divided by the leading coefficient
		//first column then row - here inefficient memory access
		for (j = i+1; j < n; j++) {	//loop now vectorized!
		  	a[j][i] = a[j][i] / div;
		}
		b[i] = b[i] / div;	//right hand side must also be divided by the leading coefficient
		//a[i][i] never reached - does not change --> error?
		//correction: a[i][i] = 1.

		//perform row operations on all subsequent rows (j=i+1)
		//first columns, then rows - ensures efficient memory access
		
		//tmp storage of factors necessary
		double factors[n-i-1];
		
		//split for loop in 2 loops
		//dependency checking since arrays are used
		#pragma ivdep
		for (j = i+1; j < n; j++) {	
			/*
				choose factor in such a way, so that subtracting the current row from the row j
				eliminates all values in the same column as the leading coefficient in all rows underneath row i (current row).
			*/
			//double factor = a[j][i]; //old code
			
			factors[j-i-1] = a[i][j];
			b[j] = b[j] - a[i][j]* b[i];	//perform subtraction operation on right-hand side for row j
		}
		
		//second part of the for loop - first columns then rows ensures more efficient mem access
		for(k=i; k < n; k++){
					
			//run through all columns and perform the subtraction operation: row[j] = row[j] - multiple of row[i]
			#pragma ivdep
			for (j = i+1; j < n; j++) { //loop was vectorized
				a[k][j] = a[k][j] - a[k][i] * factors[j-i-1];
			}
			
		}
	}

	//after matrix A has been brought to reduced row echelon form
	//backwards-substitution is used to compute the vector of solutions x
	//loop interchange helps vectorize inner loop
	//columns, then rows are iterated
	for (j = n-1; j >= 0; j--) {	//outer loop - not vectorized
		x[j] = b[j];	//assignment of final solutions
		//no dependency checking necessary
		#pragma ivdep 
		for(i = 0; i < j; i++) {	//not vectorized due to vector dependencies
			b[i] = b[i] - a[j][i] * x[j];	//temporary storage of values in vector b
		}
	}
}

//modified method
//Matrix A is stored in a way that columns can be accessed consecutively
//modified method prints matrix in correct order
void print_matrix(char* name, double matrix[3][3]) {
	int i, j;
	printf("Matrix %s: \n", name);
	for (j = 0; j < 3; j++) {	// outer loop - not vectorized
		for (i = 0; i < 3; i++) { //not vectorized - vector dependence
			printf(" %f ", matrix[i][j]);	
		}
		printf(" ;\n");
	}
}

void print_vector(char* name, double vec[3]) {
	int i;
	printf("vector %s: \n", name);
	for (i=0;i<3;i++)	//not vectorized due to vector dependence
	{
		printf(" %f \n",vec[i]);
	}
}

