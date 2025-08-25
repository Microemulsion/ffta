////////////////////////////////////////////////////////////////////////////////
// 06-11-20 kxsong
// ffta.cpp
// Ver. 1.0
// note: main function.
// Copyright (c) 2006 Kxsong
// Copyright (c) 2006 Changchun China
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
////////////////////////////////////////////////////////////////////////////////


#include "ffta.h"


int main (
	int argc,
	char *argv[]
){
	char	filename[1024];
	char	filekind[128];
	int		filenum = 1;
	FILE	*in;
	FILE	*propty;
	FILE	*scatte;
	int i, j;
	

	fftDType fftD, *fftData;
	fftPType fftP, *fftPar;

	fftData = &(fftD);
	fftPar = &(fftP);

	initfftadata(fftData);
	initfft(fftPar);


	printf("/****************************************/\n");
	printf("/*                                      */\n");
	printf("/*                                      */\n");
	printf("/*  FFT analyse Tools Version 1.0       */\n");
	printf("/*  Song Kaixu                          */\n");
	printf("/*  Changchun  China                    */\n");
	printf("/*  Copyright 2006                      */\n");
	printf("/*                                      */\n");
	printf("/*                                      */\n");
	printf("/****************************************/\n\n\n");

	if( !( propty = fopen( "property.txt", "w")))
	{
        printf("ERROR: fopen(\"property.txt\",\"rb\") = NULL.  Bye, bye!\n");
        exit(1);
	}

	if( !( scatte = fopen( "scaterf.txt", "w")))
	{
        printf("ERROR: fopen(\"scaterf.txt\",\"rb\") = NULL.  Bye, bye!\n");
        exit(1);
	}

	fprintf(propty, " filename	AVG	Qtop		phi*S/SK	phi*S/SKK	SKK/SK\n");
	

	if( argc == 2)
        {
            printf("argv = \"%s\"\n", argv[1]);
            }
     else if( argc == 1)
         {
             printf("PLease use \"ffta filename\" format!");
			 exit(1);
             }
         else
             {
                 printf(" Only one parameter \"filename\" is needed!");
                 exit(1);
                 }

	// Open the first input file.
	sprintf(filekind, ".%s", "iph");
	sprintf(filename,"%s%03d%s" , argv[1], filenum, filekind) ;
	if( !( in = fopen( filename, "rb")))
        {
            sprintf(filekind, ".%s", "iqd");
            sprintf(filename,"%s%03d%s" , argv[1], filenum, filekind) ;
            if( !( in = fopen( filename, "rb")))
			{
                printf("ERROR: fopen(\"%s\",\"rb\") = NULL.  Bye, bye!\n", filename);
                exit(1);
			}
         }
	fclose( in);


	for ( filenum = 1; filenum < 1000; filenum++)
	{
		sprintf(filename,"%s%03d%s" , argv[1], filenum, filekind) ;
        if( ( in = fopen( filename, "r")))
		{
			printf("Reading file: \"%s\"\n", filename);

			readata(fftData, in);

			fftaengine(fftData, fftPar);

			fprintf(propty, "%d	%f	%f	%f	%f	%f\n", filenum, fftData->averphi, fftData->qstar, (PI * fftData->sq / fftData->sqk), (PI * fftData->sq / fftData->sqkk), (fftData->sqkk / fftData->sqk));

			fprintf(scatte, "Filename %d\n", filenum);
			fprintf(scatte, "q	powderF\n");
			for ( i = 0; i < fftData->dataFinalS ; i++)
			{
				fprintf(scatte, "%f	%f\n", fftData->qindex[i], fftData->scattI[i]);
			}

			fprintf(scatte, "\n");
			fprintf(scatte, "index	correL\n");
			for ( i = 0; i < fftData->dataFinalS ; i++)
			{
				fprintf(scatte, "%d	%f\n", i, fftData->correD[i]);
			}

			output(fftData, argv[1], filenum);

            fclose( in);
		} else {
			//printf("Note: fopen(\"%s\",\"r\") = NULL.\n", filename);
		}
	}


    distrfftadata(fftData);
	distrfft(fftPar);

	fclose(propty);
	fclose(scatte);

	return 0;
}


int readata(fftDType *fftData, FILE *in)
{
	int i, j, yx, NN, Nnumbt, Ni;
	double tempphi ;
	char	numbt[102] = "";
	char	chart[102] = "";
	char	oneline[10240] = "";
	char   *stopstring;
	char   tempC;


	//Ignore the first two strings.
	fscanf(in, " %s ", chart);
	fscanf(in, " %s ", chart);


	for (j = 0; j < DATADIM; j++) {

		fscanf(in, " %s ", oneline);
		NN = strlen(oneline);
		Nnumbt = 0;
		Ni = 0;

		for (i = 0; i < NN; i++) {

			tempC = oneline[i];
			//printf("%c", tempC);
			numbt[Nnumbt] = tempC;
			Nnumbt += 1;

			if (tempC == 44)
			{
				//printf ("\n");
				Nnumbt = 0;
				tempphi = strtod(numbt, &stopstring);
				//printf(" %s \n", numbt);
				yx = j*DATADIM + Ni;
				fftData->phidata[yx] = tempphi;
				//printf ("tt  %d  %f |||",yx, fftData->phidata[yx]);
				Ni += 1;
			}
		}
		fread(numbt,1,1,in);
	}
	return 0;
}

int output(fftDType *fftData, char	*argv1, int filenum)
{
	int x,y,yx;
	int qx,qy;
	int istatus;
	int max_x = DATADIM;
	int max_y = DATADIM;
	char	filename[1024];

 	FILE *fp;

	
	double phi0;
	double fft0;
	double fft1;
	double core1;
	double core0;


	sprintf(filename,"%s%03d.dat" , argv1, filenum) ;
	fp = fopen(filename, "wb");
	if(!fp) return; /* was return(0) */



	/* Write the first line of plot3d format. */
	fprintf (fp,"TITLE = \"FFT ANALYSIS RESULT\"\n") ;
	/* Write the second line of plot3d format. */
	fprintf (fp,"VARIABLES = \"X\", \"Y\", \"PHI\", \"FFT\", \"FFTO\", \"CORREL\", \"CORREL0\", \n");

	/* Print the first line of the plot3d file. */
	istatus = fprintf(fp, "ZONE	I=%d,	J=%d,	F=POINT\n", max_x, max_y);

	for(y=0;y<max_y;y++)
	{
		for(x=0;x<max_x;x++)
		{
			yx=x+y*max_x;

			phi0 = fftData->phidata[yx];

			fft1 = fftData->fftdata[yx];

			core0 = fftData->corrdata[yx];


			qx = (x - DATADIM/2);
			qy = (y - DATADIM/2);
			if ( y < DATADIM/2)
				qy = (y + DATADIM/2);
			if ( x < DATADIM/2)
				qx = (x + DATADIM/2);
			yx = qx + qy * max_y;
			fft0 = fftData->fftdata[yx];

			core1 = fftData->corrdata[yx];

			istatus = fprintf(fp, "%d	%d	%f  %f	%f  %f  %f\n", x, y, phi0 , fft0, fft1, core1, core0);

			}
		}


	fclose(fp);

	return 0;

}