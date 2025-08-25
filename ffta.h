////////////////////////////////////////////////////////////////////////////////
// 06-11-20 kxsong
// ffta.h
// Ver. 1.0
////////////////////////////////////////////////////////////////////////////////


#ifndef _ffta_h_
#define _ffta_h_


#include "fftae.h"


// function prototype

int readata(fftDType *fftData, FILE *in);

int output(fftDType *fftData, char	*argv1, int filenum);

int main (
        int argc,
        char *argv[]
);



#endif 	/* _ffta_h_ */
