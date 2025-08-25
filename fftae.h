////////////////////////////////////////////////////////////////////////////////
// 06-11-20 kxsong
// fftae.h
// Ver. 1.0
//note: Fft analyse engine.
////////////////////////////////////////////////////////////////////////////////


#ifndef _fftae_h_
#define _fftae_h_

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fftw3.h"

#define PI 3.14159265358979323846264

//All calculated data.
struct fftadata {
	double *phidata;	/* Input data.*/
	double *miudata;	/* miu data.*/
	double *udata;	/* miu data.*/
	double *vdata;	/* miu data.*/
	double *fftdata;	/* FFT result.*/
	double *corrdata;	/* Corrrelation result.*/
	double *scattI;		/* Scatter function.*/
	double *correD;		/* Correlation length.*/
    double *qindex;		/* Index q for S(q).*/
	double sq;			/* sigma S.*/
	double sqk;			/* sigma Sk.*/
	double sqkk;		/* sigma Skk.*/
	double averphi;		/* average phi.*/
	int dataFinalS;		/* Length of qindex.*/
	double qstar;		/* The q to the highest Sq. */
	double fluxhy;		/* The Flux of hydr. */
	double fluxdi;		/* The Flux of diff. */

	int nx;
	int ny;
	int nn;
};
typedef struct fftadata fftDType;

//define fft variable.
struct fftparam {
    fftw_complex *p_fftwIn, *p_fftwOut;
    fftw_plan d_fftwPlanForward, d_fftwPlanBackward;
    double d_fftwCoefForward, d_fftwCoefBackward;
};
typedef struct fftparam fftPType;

//Allocate the fftDType.
int initfftadata(int dnx, int dny, fftDType *fftData);

//Free the allocated space of fftDType.
int distrfftadata(fftDType *fftData);

//FFt analyse engine.
int fftaengine(fftDType *fftData, fftPType *fftPar);

//Init FFT Engine.
int initfft(fftDType *fftData, fftPType *fftPar);

//Free FFT Engine.
int distrfft(fftPType *fftPar);

//p_fft * Coef.
int pfftMcoef(int nn, fftw_complex *fftwc, double Coef);

//Calculate Powder pattern.
int powderPattern(fftDType *fftData);

//Calculate correlation function.
int correData(fftDType *fftData, fftPType *fftPar);
int ScorreData(fftDType *fftData, fftPType *fftPar);

//Calculate Powder pattern.
int calculateSK(fftDType *fftData);

//measurement of the nematic correlation function.
int nemcorrel(fftDType *fftData, fftPType *fftPar, double *nembx, double *nemby);
int Snemcorrel(fftDType *fftData, fftPType *fftPar, double *nembx, double *nemby);

#endif 	/* _fftae_h_ */
