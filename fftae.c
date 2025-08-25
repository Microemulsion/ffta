////////////////////////////////////////////////////////////////////////////////
// 06-11-20 kxsong
// fftae.c
// Ver. 1.0
//note: Fft analyse engine.
////////////////////////////////////////////////////////////////////////////////


#include "fftae.h"

int initfftadata(int dnx, int dny, fftDType *fftData)
{
	int dnn;

	if( (dnx < 0)||(dnx >10000))
	{
        printf("Wrong nx  Bye, bye!\n");
        exit(1);
	}
	if(  (dny < 0)||(dny >10000))
	{
        printf("Wrong ny  Bye, bye!\n");
        exit(1);
	}
	fftData->nx = dnx;
	fftData->ny = dny;
	fftData->nn = (dnx * dny);

	dnn = fftData->nn;

	fftData->phidata = ( double *)malloc( sizeof(double)*dnn);
	//printf(" fftData->phidata: %f \n", fftData->phidata[2]);
	fftData->miudata = ( double *)malloc( sizeof(double)*dnn);
	fftData->udata = ( double *)malloc( sizeof(double)*dnn);
	fftData->vdata = ( double *)malloc( sizeof(double)*dnn);
	fftData->fftdata = ( double *)malloc( sizeof(double)*dnn);
	fftData->corrdata = ( double *)malloc( sizeof(double)*dnn);
	fftData->correD = ( double *)malloc( sizeof(double)*dnx);
	fftData->scattI = ( double *)malloc( sizeof(double)*dnx);
	fftData->qindex = ( double *)malloc( sizeof(double)*dnx);

    return 0;
}

int distrfftadata(fftDType *fftData)
{
	if ( fftData->phidata !=NULL) free(   fftData->phidata);
	if ( fftData->miudata !=NULL) free(   fftData->miudata);
	if ( fftData->udata !=NULL) free(   fftData->udata);
	if ( fftData->vdata !=NULL) free(   fftData->vdata);
	if ( fftData->fftdata !=NULL) free(   fftData->fftdata);
	if ( fftData->corrdata !=NULL) free(   fftData->corrdata);
	if ( fftData->correD !=NULL)  free(   fftData->correD);
	if ( fftData->scattI !=NULL)  free(   fftData->scattI);
	if ( fftData->qindex !=NULL)  free(   fftData->qindex);

    return 0;
}

int initfft(fftDType *fftData, fftPType *fftPar)
{
	int dnn, dnx;
	dnn = fftData->nn;
	dnx = fftData->nx;

	fftPar->p_fftwIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dnn);
    fftPar->p_fftwOut = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dnn);
//printf(" fftPar->p_fftwIn: %f \n", fftPar->p_fftwIn[2]);

	fftPar->d_fftwCoefForward  = 1.;
	fftPar->d_fftwCoefBackward = 1. / dnn;

	fftPar->d_fftwPlanForward 
			= fftw_plan_dft_2d(dnx, dnx
			                   , fftPar->p_fftwIn, fftPar->p_fftwOut, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftPar->d_fftwPlanBackward 
			= fftw_plan_dft_2d(dnx, dnx
			                   , fftPar->p_fftwIn, fftPar->p_fftwOut, FFTW_FORWARD , FFTW_ESTIMATE);

	return 0;
}

int distrfft(fftPType *fftPar)
{
	fftw_destroy_plan( fftPar->d_fftwPlanForward );
	fftw_destroy_plan( fftPar->d_fftwPlanBackward );
	fftw_free(fftPar->p_fftwIn);
	fftw_free(fftPar->p_fftwOut);

	return 0;
}

int pfftMcoef(int nn, fftw_complex *fftwc, double Coef)
{
	int i;

	for ( i=0; i<nn; i++)
	{
        fftwc[i][0] = fftwc[i][0] * Coef;
		fftwc[i][1] = fftwc[i][1] * Coef;		
	}

	return 0;
}

int fftaengine(fftDType *fftData, fftPType *fftPar)
{
	double tempData = 0.0;
	int i;
	int dnn;
	dnn = fftData->nn;

	for ( i = 0; i< dnn; i++)
	{
        tempData += fftData->phidata[i];
	}

	tempData /= dnn;
	fftData->averphi = tempData;

	for ( i = 0; i< dnn; i++)
	{
        fftData->phidata[i] -= tempData;
		fftPar->p_fftwIn[i][0] = fftData->phidata[i];
		fftPar->p_fftwIn[i][1] = 0.0;
		fftPar->p_fftwOut[i][0] = 0.0;
		fftPar->p_fftwOut[i][1] = 0.0;
	}

	fftw_execute(fftPar->d_fftwPlanForward);
	pfftMcoef(dnn, fftPar->p_fftwOut, fftPar->d_fftwCoefForward);


	for ( i = 0; i< dnn; i++)
	{
		tempData = fftPar->p_fftwOut[i][0];
		fftData->fftdata[i]  = tempData * tempData;
		tempData = fftPar->p_fftwOut[i][1];
		fftData->fftdata[i]  += tempData * tempData;

		//12-01-06 kxsong add a factor.
		fftData->fftdata[i] /= dnn;
	}

	powderPattern(fftData);


	correData(fftData, fftPar);


	calculateSK(fftData);


	return 0;
}

int powderPattern(fftDType *fftData)
{
	int dnn;

	int n, yx, nQ;
	int i, j;
	double maxQ;
	double qx, qy, q;
	int nx = fftData->nx;
	int ny = fftData->nx;
	int xMin = 0;
	int yMin = 0;
	double dx = 1.0;
	double dy = 1.0;
	double dd = 1.0;

	dnn = fftData->nn;


	qy = yMin + dy * ( ny/2 );
	qx = xMin + dx * ( nx/2 );

	maxQ = sqrt( qx * qx + qy * qy );
	
	nQ = (int) ( maxQ / dd ) + 1;

	fftData->dataFinalS = nQ;

	for( i = 0; i < nx; i++ ) {
		fftData->qindex[i] = 0.0;
		fftData->scattI[i] = 0.0;
	}

	for( j = 0; j < ny; j++ ) {

		if ( j < ny/2) {
			qy = yMin + dy * j;
		}
		else {
			qy = yMin + dy * (j - ny);
		}

		for( i = 0; i < nx; i++ ) {
			
			if ( i < nx/2) {
				qx = xMin + dx * i;
			}
			else {
				qx = xMin + dx * (i - nx);
			}

			yx = j * nx + i;


			q = sqrt( qx * qx + qy * qy );
			n = (int) ( ( q + dd / 2. ) / dd );

			fftData->scattI[n] += fftData->fftdata[yx];
			fftData->qindex[n] += 1.0;

		}
	}

	for( i = 0; i < nQ; i++ ) {
		fftData->scattI[i] /= fftData->qindex[i];
	}


	dd *= (2*PI/nx);
	for( i = 0; i < nQ; i++ ) {
		fftData->qindex[i] = (dd * i);
	}

	return 0;
}

int correData(fftDType *fftData, fftPType *fftPar)
{
	int dnn;

	int n, yx, nQ;
	int i, j;
	double qx, qy, q;
	int nx = fftData->nx;
	int ny = fftData->ny;
	int xMin = 0;
	int yMin = 0;
	double dx = 1.0;
	double dy = 1.0;
	double tempData1 = 0.0;
	double tempData2 = 0.0;
	double dd = 1.0;

	dnn = fftData->nn;

	for ( i = 0; i< dnn; i++)
	{
		fftPar->p_fftwIn[i][0] = fftData->fftdata[i];
		fftPar->p_fftwIn[i][1] = 0.0;
		fftPar->p_fftwOut[i][0] = 0.0;
		fftPar->p_fftwOut[i][1] = 0.0;
	}

	fftw_execute(fftPar->d_fftwPlanBackward);
	pfftMcoef(dnn, fftPar->p_fftwOut, fftPar->d_fftwCoefBackward);


	for ( i = 0; i< dnn; i++)
	{
		tempData1 = fftPar->p_fftwOut[i][0];
		fftData->corrdata[i] = tempData1;
	}

	nQ = (nx/2);;

	for( i = 0; i < nx; i++ ) {
		fftData->qindex[i] = 0.0;
		fftData->correD[i] = 0.0;
	}

	for( j = 0; j < ny; j++ ) {

		if ( j < ny/2) {
			qy = yMin + dy * j;
		}
		else {
			qy = yMin + dy * (j - ny);
		}

		for( i = 0; i < nx; i++ ) {
			
			if ( i < nx/2) {
				qx = xMin + dx * i;
			}
			else {
				qx = xMin + dx * (i - nx);
			}

			yx = j * nx + i;

			q = sqrt( qx * qx + qy * qy );
			n = (int) ( ( q + dd / 2. ) / dd );

			fftData->correD[n] += fftData->corrdata[yx];
			fftData->qindex[n] += 1.0;

		}
	}

	for( i = 0; i < nQ; i++ ) {
		fftData->correD[i] /= fftData->qindex[i];
	}

	dd *= (2*PI)/nx;
	for( i = 0; i < nx; i++ ) {
		fftData->qindex[i] = (dd * i);
	}

	return 0;
}

int calculateSK(fftDType *fftData)
{
	int i;
	double S = 0;
	double Sk = 0;
	double Skk = 0;
	double Sfirst, Ssecond;

	double qsta0 = 0.0;
	double sq0 = 0.0;

	for( i = 0; i < (fftData->dataFinalS); i++ ) {
		Sfirst	= fftData->qindex[i];
		Ssecond = fftData->scattI[i];
		S += Ssecond ;
		Sk += Sfirst * Ssecond ;
		Skk += Sfirst * Sfirst * Ssecond;
		
		if ( Ssecond > sq0 )
		{
			sq0 = Ssecond;
			qsta0 = Sfirst;
		}
	}

	fftData->sq = S;
	fftData->sqk = Sk;
	fftData->sqkk = Skk;

	fftData->qstar = qsta0;

	return 0;
}

//measurement of the nematic correlation function.
int nemcorrel(fftDType *fftData, fftPType *fftPar, double *nembx, double *nemby)
{
	double tempData = 0.0;
	int i, j;
	int dnn;

	double tempData1 = 0.0;
	double tempData2 = 0.0;

	int n, yx, nQ;
	double dd = 1.0;
	double qx, qy, q;
	int nx = fftData->nx;
	int ny = fftData->ny;
	int xMin = 0;
	int yMin = 0;
	double dx = 1.0;
	double dy = 1.0;


	//dnn = fftData->nn;

	dnn = fftData->nn;

	//for ( i = 0; i< dnn; i++)
	//{
	//	tempData += nembx[i];
	//}

	//tempData /= dnn;
	//fftData->averphi = tempData;

	for ( i = 0; i< dnn; i++)
	{
        //fftData->phidata[i] -= tempData;
		fftPar->p_fftwIn[i][0] = nembx[i];
		fftPar->p_fftwIn[i][1] = 0.0;
		fftPar->p_fftwOut[i][0] = 0.0;
		fftPar->p_fftwOut[i][1] = 0.0;
	}

	fftw_execute(fftPar->d_fftwPlanForward);
	pfftMcoef(dnn, fftPar->p_fftwOut, fftPar->d_fftwCoefForward);


	for ( i = 0; i< dnn; i++)
	{
		tempData = fftPar->p_fftwOut[i][0];
		fftData->fftdata[i]  = tempData * tempData;
		tempData = fftPar->p_fftwOut[i][1];
		fftData->fftdata[i]  += tempData * tempData;

		//12-01-06 kxsong add a factor.
		fftData->fftdata[i] /= dnn;
	}

	for ( i = 0; i< dnn; i++)
	{
		fftPar->p_fftwIn[i][0] = fftData->fftdata[i];
		fftPar->p_fftwIn[i][1] = 0.0;
		fftPar->p_fftwOut[i][0] = 0.0;
		fftPar->p_fftwOut[i][1] = 0.0;
	}

	fftw_execute(fftPar->d_fftwPlanBackward);
	pfftMcoef(dnn, fftPar->p_fftwOut, fftPar->d_fftwCoefBackward);


	for ( i = 0; i< dnn; i++)
	{
		tempData1 = fftPar->p_fftwOut[i][0];
		fftData->corrdata[i] = tempData1;
	}

	//calculate another part. And add to the first part.
	//for ( i = 0; i< dnn; i++)
	//{
	//	tempData += nembx[i];
	//}

	//tempData /= dnn;
	//fftData->averphi = tempData;

	for ( i = 0; i< dnn; i++)
	{
        //fftData->phidata[i] -= tempData;
		fftPar->p_fftwIn[i][0] = nemby[i];
		fftPar->p_fftwIn[i][1] = 0.0;
		fftPar->p_fftwOut[i][0] = 0.0;
		fftPar->p_fftwOut[i][1] = 0.0;
	}

	fftw_execute(fftPar->d_fftwPlanForward);
	pfftMcoef(dnn, fftPar->p_fftwOut, fftPar->d_fftwCoefForward);


	for ( i = 0; i< dnn; i++)
	{
		tempData = fftPar->p_fftwOut[i][0];
		fftData->fftdata[i]  = tempData * tempData;
		tempData = fftPar->p_fftwOut[i][1];
		fftData->fftdata[i]  += tempData * tempData;

		//12-01-06 kxsong add a factor.
		fftData->fftdata[i] /= dnn;
	}

	for ( i = 0; i< dnn; i++)
	{
		fftPar->p_fftwIn[i][0] = fftData->fftdata[i];
		fftPar->p_fftwIn[i][1] = 0.0;
		fftPar->p_fftwOut[i][0] = 0.0;
		fftPar->p_fftwOut[i][1] = 0.0;
	}

	fftw_execute(fftPar->d_fftwPlanBackward);
	pfftMcoef(dnn, fftPar->p_fftwOut, fftPar->d_fftwCoefBackward);


	for ( i = 0; i< dnn; i++)
	{
		tempData1 = fftPar->p_fftwOut[i][0];
		fftData->corrdata[i] += tempData1;
	}

	nQ = (nx/2);;

	for( i = 0; i < nx; i++ ) {
		fftData->qindex[i] = 0.0;
		fftData->correD[i] = 0.0;
	}

	for( j = 0; j < ny; j++ ) {

		if ( j < ny/2) {
			qy = yMin + dy * j;
		}
		else {
			qy = yMin + dy * (j - ny);
		}

		for( i = 0; i < nx; i++ ) {
			
			if ( i < nx/2) {
				qx = xMin + dx * i;
			}
			else {
				qx = xMin + dx * (i - nx);
			}

			yx = j * nx + i;

			q = sqrt( qx * qx + qy * qy );
			n = (int) ( ( q + dd / 2. ) / dd );

			if ( n < nQ ) {
                fftData->correD[n] += fftData->corrdata[yx];
                fftData->qindex[n] += 1.0;
			}
		}
	}

	for( i = 0; i < nQ; i++ ) {
		fftData->correD[i] /= fftData->qindex[i];
	}

	dd *= (2*PI)/nx;
	for( i = 0; i < nx; i++ ) {
		fftData->qindex[i] = (dd * i);
	}

	return 0;
}

int ScorreData(fftDType *fftData, fftPType *fftPar)
{
	int yx, nQ;
	int i, j, k, l, ki, lj;
	double maxQ;
	double dd = 1.0;
	int nx = fftData->nx;
	int ny = fftData->ny;
	int dnn = fftData->nn;
	int d1nn = nx-1;
	int xMin = 0;
	int yMin = 0;
	double dx = 1.0;
	double dy = 1.0;
	double tempData1 = 0.0;
	double tempData2 = 0.0;
	int NRange=64;
	int NR2=NRange*NRange;

	for ( i = 0; i< dnn; i++)
	{
		fftData->corrdata[i] = 0.0;
	}

	for ( i = 0; i< nx; i++){
											
		for ( j = 0; j< nx; j++){


			for ( k = 0; k< nx; k++){
				for ( l = 0; l< nx; l++){
					ki=(k-i);
					lj=(l-j);
					if ((ki*ki+lj*lj)<NR2) {
                        ki=abs(ki);
                        lj=abs(lj);
						yx = i*nx+j;
						tempData1 = fftData->phidata[yx];
						yx = k*nx+l;
						tempData2 = fftData->phidata[yx];
						yx = ki*nx+lj;
						fftData->corrdata[yx] += tempData1*tempData2;
						yx = (d1nn-ki)*nx+(d1nn-lj);
						fftData->corrdata[yx] += 1.0;
					}
				}
			}
		}
	}

	for ( i = 0; i< NRange; i++){											
		for ( j = 0; j< NRange; j++){

					if ((i*i+j*j)<NR2) {
						yx = (d1nn-i)*nx+(d1nn-j);
						tempData1 = fftData->corrdata[yx];

						yx = i*nx+j;
						fftData->corrdata[yx] /= tempData1;
					}
		}
	}

	for ( i = 0; i< NRange; i++){											
		for ( j = 0; j< NRange; j++){

					if ((i*i+j*j)<NR2) {
						yx = (d1nn-i)*nx+(d1nn-j);
						tempData1 = fftData->corrdata[yx];
						fftData->corrdata[yx] = 0.0;

						yx = i*nx+j;
						fftData->corrdata[yx] /= tempData1;
					}
		}
	}

	for( i = 0; i < nx; i++ ) {
        fftData->qindex[i] = 0.0;
        fftData->correD[i] = 0.0;
	}

	for ( i = 0; i< NRange; i++){											
		for ( j = 0; j< NRange; j++){
			maxQ = sqrt( i * i + j * j );
			nQ = (int) ( ( maxQ + dd / 2. ) / dd );
			yx = i*nx+j;
            if ( nQ < NRange ) {
				fftData->correD[nQ] += fftData->corrdata[yx];
                fftData->qindex[nQ] += 1.0;

			}
		}
	}

	for( i = 0; i < nQ; i++ ) {
		fftData->correD[i] /= fftData->qindex[i];
	}

	dd = (2*PI)/nx;
	for( i = 0; i < nx; i++ ) {
		fftData->qindex[i] = (dd * i);
	}

	return 0;
}

int Snemcorrel(fftDType *fftData, fftPType *fftPar, double *nembx, double *nemby)
{
	int yx, nQ;
	int i, j, k, l, ki, lj;
	double maxQ;
	double dd = 1.0;
	int nx = fftData->nx;
	int ny = fftData->ny;
	int dnn = fftData->nn;
	int d1nn = nx-1;
	int xMin = 0;
	int yMin = 0;
	double dx = 1.0;
	double dy = 1.0;
	double tempData1 = 0.0;
	double tempData2 = 0.0;
	double tempData3 = 0.0;
	double tempData4 = 0.0;
	int NRange=64;
	int NR2=NRange*NRange;

	for ( i = 0; i< dnn; i++)
	{
		fftData->corrdata[i] = 0.0;
	}

	for ( i = 0; i< nx; i++){
											
		for ( j = 0; j< nx; j++){


			for ( k = 0; k< nx; k++){
				for ( l = 0; l< nx; l++){
					ki=(k-i);
					lj=(l-j);
					if ((ki*ki+lj*lj)<NR2) {
                        ki=abs(ki);
                        lj=abs(lj);
						yx = i*nx+j;
						tempData1 = nembx[yx];
						tempData3 = nemby[yx];
						yx = k*nx+l;
						tempData2 = nembx[yx];
						tempData4 = nemby[yx];
						yx = ki*nx+lj;
						fftData->corrdata[yx] += tempData1*tempData2;
						fftData->fftdata[yx] += tempData3*tempData4;
						yx = (d1nn-ki)*nx+(d1nn-lj);
						fftData->corrdata[yx] += 1.0;
						fftData->fftdata[yx] += 1.0;
					}
				}
			}
		}
	}

	for ( i = 0; i< NRange; i++){											
		for ( j = 0; j< NRange; j++){

					if ((i*i+j*j)<NR2) {
						yx = (d1nn-i)*nx+(d1nn-j);
						tempData1 = fftData->corrdata[yx];
						tempData3 = fftData->fftdata[yx];

						yx = i*nx+j;
						fftData->corrdata[yx] /= tempData1;
						fftData->fftdata[yx] /= tempData3;
					}
		}
	}

	for ( i = 0; i< NRange; i++){											
		for ( j = 0; j< NRange; j++){

					if ((i*i+j*j)<NR2) {
						yx = (d1nn-i)*nx+(d1nn-j);
						tempData1 = fftData->corrdata[yx];
						tempData3 = fftData->fftdata[yx];
						fftData->corrdata[yx] = 0.0;
						fftData->fftdata[yx] = 0.0;

						yx = i*nx+j;
						fftData->corrdata[yx] /= tempData1;
						fftData->fftdata[yx] /= tempData1;
						fftData->corrdata[yx] += fftData->fftdata[yx];
					}
		}
	}

	for( i = 0; i < nx; i++ ) {
        fftData->qindex[i] = 0.0;
        fftData->correD[i] = 0.0;
	}

	for ( i = 0; i< NRange; i++){											
		for ( j = 0; j< NRange; j++){            
			maxQ = sqrt( i * i + j * j );
			nQ = (int) ( ( maxQ + dd / 2. ) / dd );
			yx = i*nx+j;
            if ( nQ < NRange ) {
				fftData->correD[nQ] += fftData->corrdata[yx];
                fftData->qindex[nQ] += 1.0;

			}
		}
	}

	for( i = 0; i < nQ; i++ ) {
		fftData->correD[i] /= fftData->qindex[i];
	}

	dd *= (2*PI)/nx;
	for( i = 0; i < nx; i++ ) {
		fftData->qindex[i] = (dd * i);
	}

	return 0;
}
