#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

const int N = 2048;
const int B = 32;

int main()
{
    int n = N;
    int i, j, k, i1, j1, k1;
    double *a, *b, *c, *cc;
    a = (double *) malloc(N*N*sizeof(double));
    b = (double *) malloc(N*N*sizeof(double));
    c = (double *) malloc(N*N*sizeof(double));
    cc = (double *) malloc(N*N*sizeof(double));
    if(!a || !b || !c || !cc)
    {
	printf("allocate failed!\n");
	exit(-1);
    }
    for(i = 0; i < N*N; i++)
    {
	a[i] = rand()%1000000/1000.0;
	b[i] = rand()%1000000/1000.0;
	c[i] = 0.0;
	cc[i] = 0.0;
    }

    time_t start, end;
    start = clock();
    for(i = 0; i < n; i++)
	for(j = 0; j < n; j++)
	    for(k = 0; k < n; k++)
	    {
		cc[i*n+j] += a[i*n+k] * b[k*n+j];
	    }
    end = clock();
    printf("ijk version: %.6lf\n", (double)(end-start)/CLOCKS_PER_SEC);

    start = clock();
    for(i = 0; i < N; i += B)
	for(j = 0; j < N; j += B)
	    for(k = 0; k < N; k += B)
		for(i1 = i; i1 < i+B; i1+=2)
		    for(j1 = j; j1 < j+B; j1+=2)
			for(k1 = k; k1 < k+B; k1+=2)
			{
			    register int t=i1*n+j1;
			    register int tt=t+n;
			    register int ta = i1*n+k1;
			    register int tta = ta+n;
			    register int tb = k1*n+j1;
			    register int ttb = tb+n;
			    register double c00= c[t];
			    register double c01= c[t+1];
			    register double c10= c[tt];
			    register double c11= c[tt+1];
			    register double b00= b[tb];
			    register double b01= b[tb+1];
			    register double b10= b[ttb];
			    register double b11= b[ttb+1];
			    register double a00= a[ta];
			    register double a01= a[ta+1];
			    register double a10= a[tta];
			    register double a11= a[tta+1];

			    c00 += a00*b00;
			    c01 += a00*b01;
			    c10 += a10*b00;
			    c11 += a10*b01;
			    c00 += a01*b10;
			    c01 += a01*b11;
			    c10 += a11*b10;
			    c11 += a11*b11;

			    c[t] =c00;
			    c[t+1] = c01;
			    c[tt] = c10;
			    c[tt+1] = c11;
			}
    end = clock();
    printf("Blocking version: %.6lf\n", (double)(end-start)/CLOCKS_PER_SEC);

    /* check error */
    double error = 0;
    for(i = 0; i < N*N; i++)
	error += (c[i] - cc[i]) * (c[i] - cc[i]);

    printf("error = %e\n", sqrt(error));
}
