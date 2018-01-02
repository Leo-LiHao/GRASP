/* compile: mex -g -I/LIGO/GRASP/include -I/LIGO/frame/v3.60/include inspfilt.c -L/LIGO/GRASP/lib 	*/
/*	-lgrasp -L/usr/local/lib -lrecipes_c -L/LIGO/frame/v3.60/lib -lFrame				*/
/*													*/
/* example: load data;[snr_max,timestart]=inspfilt(m1(1:10),m2(1:10),data,srate,flo,fri);		*/
/* Steve Drasco												*/
/* summer 1998												*/

#include "grasp.h"
#include "mex.h"

#define HSCALE 1.e21       /* A convenient scaling factor; results independent of it (used to be 1.e21) */
#define SAFETY 1000        /* Padding safety factor to avoid wraparound errors */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        float	*m1, *m2, *data, srate, *snr_max, flo, *fri;
	double	*timestart;
        int     i, ntmplt, npoint, frinum;
	void	multifilter(float*,float*,float*,int,float,int,float,float*,int,float*,double*);

        if (nrhs != 6)
                mexErrMsgTxt("Need 6 inputs (m1, m2, data, srate, flo, fri)");
        if (nlhs != 2)
                mexErrMsgTxt("Need 2 outputs (snr_max, timestart)");

	/* scalars */
	ntmplt = (int) mxGetM(prhs[0]) * mxGetN(prhs[0]);
	npoint = (int) mxGetM(prhs[2]) * mxGetN(prhs[2]);
	srate = (float) mxGetScalar(prhs[3]);
	flo =(float) mxGetScalar(prhs[4]);
	frinum = (int) mxGetM(prhs[5]) * mxGetN(prhs[5]);

	/* arrays */
	m1 = (float *) mxMalloc(ntmplt*sizeof(float));
	m2 = (float *) mxMalloc(ntmplt*sizeof(float));
	snr_max = (float *) mxMalloc(ntmplt*sizeof(float));
	fri = (float *) mxMalloc(frinum*sizeof(float));
	timestart = (double *) mxMalloc(ntmplt*sizeof(double));
	data = (float *) mxMalloc(npoint*sizeof(float));
	for (i = 0; i < ntmplt; i++) {
		m1[i] = (float) mxGetPr(prhs[0])[i];
		m2[i] = (float) mxGetPr(prhs[1])[i];
	}
	for (i = 0; i < npoint; i++) {
		data[i] = (float) mxGetPr(prhs[2])[i];
	}
	for (i = 0; i < frinum; i++) {
                fri[i] = (float) mxGetPr(prhs[5])[i];
        }

	/* call filter function */
        multifilter(m1, m2, data, ntmplt, srate, npoint, flo, fri, frinum, snr_max, timestart);

	/* create and fill output matrices */
        plhs[0] = mxCreateDoubleMatrix(ntmplt, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(ntmplt, 1, mxREAL);
        for (i = 0; i < ntmplt; i++) {
                mxGetPr(plhs[0])[i] = snr_max[i];
		mxGetPr(plhs[1])[i] = timestart[i];
        } 
	
	/* tell me you're done */
	mexPrintf("\n\nDone filtering %d data points through %d filters\n\n",npoint,ntmplt);
}

/* GRASP: Copyright 1997,1998  Bruce Allen */
	void multifilter(float *m1,float *m2,float *data,int ntmplt,float srate,int npoint,float flo,float *fri,int frinum,float *snr_max,double *timestart) 
	{
	void realft(float*,unsigned long,int);
	int i,j,code,maxi,chirplen,impulseoff,chirppoints,indices[8],reset=0;
	short *datas;
	float distance,var,*mean_pow_spec, *response;
	float *htilde,*chirp0,*chirp90,*ch0tilde,*ch90tilde, *output0, *output90;
	float n0,n90,inverse_distance_scale,decaytime,*twice_inv_noise,tc;
	float lin0,lin90,varsplit,stats[8],gammq(float,float);
	double decay,norm=0.0,prob,timeoff;
	float *inject_chirp0,*inject_chirp90,invMpc_inject,inject_tc,*inject_output0;
	int   inject_chirppoints;
	
	/* stores data as short integers */
	datas = (short *) mxMalloc(npoint*sizeof(short));
        for (i=0;i<npoint;i++) {
		datas[i]=data[i];
	}

	/* Memory allocation */
	chirp0=(float *)mxMalloc(sizeof(float)*npoint);
	chirp90=(float *)mxMalloc(sizeof(float)*npoint);
	ch0tilde=(float *)mxMalloc(sizeof(float)*npoint);
	ch90tilde=(float *)mxMalloc(sizeof(float)*npoint);
	response=(float *)mxMalloc(sizeof(float)*(npoint+2));
	htilde=(float *)mxMalloc(sizeof(float)*npoint);
	mean_pow_spec=(float *)mxMalloc(sizeof(float)*(npoint/2+1));
	twice_inv_noise=(float *)mxMalloc(sizeof(float)*(npoint/2+1));
	output0=(float *)mxMalloc(sizeof(float)*(npoint+2));
	output90=(float *)mxMalloc(sizeof(float)*(npoint+2));

	/* set constants */
        inverse_distance_scale=2.0*HSCALE*(TSOLAR*C_LIGHT/MPC);

	decaytime=15.0*npoint/srate;
        decay=exp(-1.0*npoint/(srate*decaytime));

        /* get the response function, and put in scaling factor */
        GRnormalize(fri,frinum,npoint,srate,response);
        for (i=0;i<npoint+2;i++) response[i]*=HSCALE/ARMLENGTH_1994;

	/* ############################### inject chirp ################################# */
        /* If you want to use this you must alter the mass values in the call to 	  */
	/* make_filters. Otherwise the program will hang.

	inject_chirp0=(float *)mxMalloc(sizeof(float)*npoint);
	inject_chirp90=(float *)mxMalloc(sizeof(float)*npoint);
	inject_output0=(float *)mxMalloc(sizeof(float)*npoint);

	make_filters(m1[277],m2[286],inject_chirp0,inject_chirp90,flo,npoint,srate,&inject_chirppoints,&inject_tc,0,4);
	for (i=0;i<inject_chirppoints;i++) {
		inject_chirp0[i]*=inverse_distance_scale;
		inject_chirp90[i]*=inverse_distance_scale;
	}
	invMpc_inject=100.0;
        time_inject_chirp(1.0,0.0,12345,invMpc_inject,inject_chirp0,inject_chirp90,data,response,inject_output0,npoint);

        										  */
        /* ############################################################################## */

	 /* find the FFT of data*/
        realft(data-1,npoint,1);

	/* normalized delta-L/L tilde */
        product(htilde,data,response,npoint/2);

	/* get the inverse of the auto-regressive-mean power-spectrum */
        avg_inv_spec(flo,srate,npoint,decay,&norm,htilde,mean_pow_spec,twice_inv_noise);

	/* main loop to filter data */
	for(i = 0; i < ntmplt; i++) {
		/* manufacture two chirps (dimensionless strain at 1 Mpc distance) */
		make_filters(m1[i],m2[i],chirp0,chirp90,flo,npoint,srate,&chirppoints,&tc,0,4);

		/* normalization */
		for (j=0;j<chirppoints;j++){
			ch0tilde[j]=chirp0[j]*=inverse_distance_scale;
			ch90tilde[j]=chirp90[j]*=inverse_distance_scale;
		}

		/* zero out the unused elements of the tilde arrays */
		for (j=chirppoints;j<npoint;j++)
			ch0tilde[j]=ch90tilde[j]=0.0;

		/* FFT the chirps */
		realft(ch0tilde-1,npoint,1);
		realft(ch90tilde-1,npoint,1);

		/* set length of template including a safety margin */
		chirplen=chirppoints+SAFETY;
		if (chirplen>npoint) {
			snr_max[i] = timestart[i] = 0;
			mexPrintf("Length of Chirp greater than number of data points. Skipping filter number %d and returnning 0 value",i);
			goto endofloop;
		}

		/* orthogonalize the chirps: we never modify ch0tilde, only ch90tilde */
		orthonormalize(ch0tilde,ch90tilde,twice_inv_noise,npoint,&n0,&n90);

		/* distance scale Mpc for SNR=1 */
		distance=sqrt(1.0/(n0*n0)+1.0/(n90*n90));
		

		/* find the moment at which SNR is a maximum */
		find_chirp(htilde,ch0tilde,ch90tilde,twice_inv_noise,n0,n90,output0,output90,
			npoint,chirplen,&maxi,(snr_max+i),&lin0,&lin90,&var);

		/* identify when an impulse would have caused observed filter output */
		impulseoff=(maxi+chirppoints)%npoint;
		timeoff=impulseoff/srate;
		timestart[i]=maxi/srate;

		/* if SNR greater than 5, then print details, else just short message */
		if (snr_max[i] <5.0) {
			mexPrintf("max snr: %.2f offset: %d variance: %.5f\n", snr_max[i],maxi,var);
		} else {
			/* See if the nominal chirp can pass a frequency-space veto test */
			varsplit=splitup_freq2(lin0*n0/sqrt(2.0),lin90*n90/sqrt(2.0),ch0tilde,ch90tilde,2.0/(n0*n0),
				twice_inv_noise,npoint,maxi,8,indices,stats,output0,htilde);
			prob=gammq(4.0,4.0*varsplit);
			mexPrintf("\nMax SNR: %.2f (offset %d) variance %f\n",snr_max[i],maxi,var);
			mexPrintf("   If impulsive event, offset %d or time %.2f\n",impulseoff,timeoff);
			mexPrintf("   If inspiral, template start offset %d (time %.2f) coalescence time %.2f\n",
				maxi,timestart[i],timestart[i]+tc);
			mexPrintf("   Normalization: S/N=1 at %.2f kpc\n",1000.0*distance);
			mexPrintf("   Linear combination of max SNR: %.4f x phase_0 + %.4f x phase_pi/2\n",lin0,lin90);
			if (prob<0.01)
				mexPrintf("   Less than 1%% probability that this is a chirp (p=%f).\n",prob);
			else
				mexPrintf("   POSSIBLE CHIRP!  with > 1%% probability (p=%f).\n",prob);
		
			/* See if the time-domain statistics are unusual or appears Gaussian */
			if (is_gaussian(datas,npoint,-2048,2047,1))
				mexPrintf("   Distribution does not appear to have outliers...\n\n");
			else
				mexPrintf("   Distribution has outliers! Reject\n\n");
		}
		endofloop:;
	}

}
