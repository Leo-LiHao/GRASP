/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
static char *rcsid="$Id: matched.c,v 1.12 1998/10/26 01:37:45 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

void freq_inject_chirp(float c0,float c90,int offset,float invMpc,float* ch0tilde,float* ch1tilde,float* htilde,int n) {

	int i,re,im;
	float sreal,simag,preal,pimag,angle;
        /* test that memory has been allocated */

        if (htilde==NULL) {
            GR_start_error("freq_inject_chirp()",rcsid,__FILE__,__LINE__);
            GR_report_error("Second to last argument can not be NULL pointer.\n");
            GR_end_error();
            abort();
        }
        
	/* if the source is infinitely far away, return */
	if (invMpc==0.0)
		return;

	/* get correct linear combination & normalization */
	c0*=invMpc;
	c90*=invMpc;

	/* phase angle needed to transform to correct "time" offset */
	angle=(-2.0*M_PI*offset)/n;

	/* now loop through, correcting phase and summing into htilde */
	for (i=0;i<n/2;i++) {
		/* subscript of real and imaginary parts */
		im=(re=i+i)+1;

		/* real and imaginary parts of linear combination of chirps */
		sreal=c0*ch0tilde[re]+c90*ch1tilde[re];
		simag=c0*ch0tilde[im]+c90*ch1tilde[im];

		/* real and imaginary parts of time-offset induced phase */
		preal=cos(i*angle);
		pimag=sin(i*angle);

		/* add in the phase-shifted chirp */
		htilde[re]+=sreal*preal+simag*pimag;
		htilde[im]+=simag*preal-sreal*pimag;
	}
	return;
}

void time_inject_chirp(float c0,float c90,int offset,float invMpc,float *chirp0,float *chirp1,
                 float* data,float* response,float* work,int n) {

	void realft(float*,unsigned long,int);
	int i;

	 /* if the source is infinitely far away, return */
	if (invMpc==0.0)
		return;

	/* get correct linear combination & normalization (factor of 2/n for FFT/IFFT pair) */
	invMpc*=2.0/n;
	c0*=invMpc;
	c90*=invMpc;

	/* set the first offset entries to zero */
	clear(work,offset,1);

	/* find chirp in time domain, with offset  */
	for (i=0;i<n-offset;i++)
		work[i+offset]=c0*chirp0[i]+c90*chirp1[i];

	/* FFT chirp into frequency domain */
	realft(work-1,n,1);

	/* now divide by response function (skipping DC, Nyquist) */
	ratio(work+2,work+2,response+2,n/2-1);

	/* deal with DC, Nyquist */
	work[0]=work[1]=0.0;

	/* and return to the time domain */
	realft(work-1,n,-1);

	/* and add into data, without quantization errors */
	for (i=offset;i<n;i++) data[i]+=work[i];

	return;
}

float splitup_freq(float c0,float c1,float *chirp0, float *chirp1,float norm,
		float* twice_inv_noise,int npoint,int offset,int nbands,int* indices,float* stats,
			float* working,float* htilde) {

	float varsplit,sigsum=0.0;
	int i,j;

	/* determine the different frequency bands */
	splitup(working,chirp0,twice_inv_noise,npoint,norm,nbands,indices);

	/* set total variance to zero */
	varsplit=0.0;
	j=0;

	/* step through all the frequency bands */
	for (i=0;i<nbands;i++) {

		/* clear the working array */
		clear(working,npoint,1);

		/* construct filter with correct support */
		while (j<=indices[i]) {
			working[j]=c0*chirp0[j]+c1*chirp1[j];
			j++;
		}
		correlate(working,htilde,working,twice_inv_noise,npoint);
		sigsum+=stats[i]=working[offset];
	}
	for (i=0;i<nbands;i++)
		varsplit+=(sigsum/nbands-stats[i])*(sigsum/nbands-stats[i]);

	return varsplit;
}


/* This is a more efficent version of the previous routine */
float splitup_freq2(float c0,float c1,float *chirp0, float *chirp1,float norm,
		float* twice_inv_noise,int npoint,int offset,int nbands,int* indices,float* stats,
			float* working,float* htilde) {

	float angle,sreal,simag,preal,pimag,sigsum,varsplit,signal;
	int i,j,im,re;
	/* determine the different frequency bands */
	splitup(working,chirp0,twice_inv_noise,npoint,norm,nbands,indices);

	/* now loop through, correcting phase and summing into htilde */
	angle=(-2.0*M_PI*offset)/npoint;
	for (i=0;i<npoint/2;i++) {
		/* subscript of real and imaginary parts */
		im=(re=i+i)+1;

		/* real and imaginary parts of linear combination of chirps */
		sreal=c0*chirp0[re]+c1*chirp1[re];
		simag=c0*chirp0[im]+c1*chirp1[im];

		/* real and imaginary parts of time-offset induced phase */
		preal=cos(i*angle);
		pimag=sin(i*angle);

		/* put in the phase-shifted chirp */
		working[re]=sreal*preal+simag*pimag;
		working[im]=simag*preal-sreal*pimag;
	}

	/* now calculate the variance in the signal contributions */
	sigsum=varsplit=0.0;
	j=0;
	for (i=0;i<nbands;i++) {
		signal=0.0;
		while (j<=indices[i]) {
			signal+=htilde[j]*working[j]*twice_inv_noise[j/2];
			j++;
		}
		stats[i]=signal;
		sigsum+=signal;
	}
	for (i=0;i<nbands;i++)
		varsplit+=(sigsum/nbands-stats[i])*(sigsum/nbands-stats[i]);

	return varsplit;
}


/* This is the splitup_freq routine when the phase of the signal is unknown */ 
float splitup_freq3(float c0,float c1,float *chirp0, float *chirp1,
		    float norm, float* twice_inv_noise,int npoint,
		    int offset,int nbands,int* indices,float* stats, 
		    float* working,float* htilde) {
    
    float angle,preal,pimag,sigsumreal,sigsumimag;
    float varsplit,signalreal,signalimag,oldnorm;
    int i,j,im,re;
    /* determine the different frequency bands */
    splitup(working,chirp0,twice_inv_noise,npoint,norm,nbands,indices);

    /* match old normalization if desired.  If not, set c0=c1=sqrt(2) */
    oldnorm=0.5*sqrt(c0*c0+c1*c1);

    /* initialize variables */
    sigsumreal=sigsumimag=varsplit=0.0;

    /* now loop through, correcting phase and summing into htilde */
    angle=(-2.0*M_PI*offset)/npoint;
    /* Compute the real part */
    for (i=0;i<npoint/2;i++) {
	/* subscript of real and imaginary parts */
	im=(re=i+i)+1;

	/* real and imaginary parts of time-offset induced phase */
	preal=oldnorm*cos(i*angle);
	pimag=oldnorm*sin(i*angle);

	/* compute the phase-shifted chirp. */
	/* Complex conjugate of T0^* exp(-2 pi i f t0) */
	/* See e:twophasesig in the GRASP manual for comparison */
	working[re]=chirp0[re]*preal+chirp0[im]*pimag;
	working[im]=chirp0[im]*preal-chirp0[re]*pimag;
    }
    j=0;
    for (i=0;i<nbands;i++) {
	signalreal=0.0;
	while (j<indices[i]) {
	    /* See GRASP equation e:twophasesig to understand the factor of 2*/
	    signalreal+=(htilde[j]*working[j])*2*twice_inv_noise[j/2];
	    j++;
	}
	/* remember stats[0..2p-1] has 2p elements */
	re=i+i;
	stats[re]=signalreal;
	/* the entire signal */
	sigsumreal+=signalreal;
    }

    /* Compute the imaginary part */
    for (i=0;i<npoint/2;i++) {
	/* subscript of real and imaginary parts */
	im=(re=i+i)+1;

	/* real and imaginary parts of time-offset induced phase */
	preal=oldnorm*cos(i*angle);
	pimag=oldnorm*sin(i*angle);

	/* compute the phase-shifted chirp. */
	/* Complex conjugate of [T90^*] exp(-2 pi i f t0) */
	/* See e:twophasesig in the GRASP manual for comparison */
	working[re]=chirp1[re]*preal+chirp1[im]*pimag;
	working[im]=chirp1[im]*preal-chirp1[re]*pimag;
    }
    j=0;
    for (i=0;i<nbands;i++) {
	signalimag=0.0;
	while (j<indices[i]) {
	    /* See GRASP equation e:twophasesig to understand the factor of 2*/
	    signalimag+=(htilde[j]*working[j])*2*twice_inv_noise[j/2];
	    j++;
	}
	/* remember stats[0..2p-1] has 2p elements */
	im=(i+i+1);
	stats[im]=signalimag;
	/* the entire signal */
	sigsumimag+=signalimag;
    }

    varsplit=0.0;
    for (i=0;i<nbands;i++){
	im=(re=i+i)+1;
	signalreal=(sigsumreal/nbands-stats[re]);
	signalimag=(sigsumimag/nbands-stats[im]);
	varsplit+=(signalreal*signalreal+signalimag*signalimag);
    }
    return varsplit;
}


/* This is the splitup_freq routine when the phase of the signal is unknown */
/* More efficient version by Bruce Allen, 10/21/98 */
float splitup_freq4(float c0,float c1,float *chirp0, float *chirp1,
		    float norm, float* twice_inv_noise,int npoint,
		    int offset,int nbands,int* indices,float* stats, 
		    float* working,float* htilde) {
    
    float preal,pimag,sigsumreal,sigsumimag;
    float varsplit,signalreal,signalimag,oldnorm,*working2,temp;
    int i,j,im,re;
    double cosd,sind,cosi,sini,cosin,angle;
    /* determine the different frequency bands */
    splitup(working,chirp0,twice_inv_noise,npoint,norm,nbands,indices);

    /* match old normalization if desired.  If not, set c0=c1=sqrt(2) */
    oldnorm=0.5*sqrt(c0*c0+c1*c1);

    /* initialize variables */
    sigsumreal=sigsumimag=varsplit=0.0;

    /* now loop through, correcting phase and summing into htilde */
    angle=(-2.0*M_PI*offset)/npoint;

    /* iterate these to find sin and cosine of angle.  Will give float precision */
    cosd=cos(angle);
    sind=sin(angle);
    cosi=1.0;
    sini=0.0;

    /* define start of second work array */
    working2=working+npoint;

    /* Compute the real part */
    for (i=0;i<npoint/2;i++) {
	/* subscript of real and imaginary parts */
	im=(re=i+i)+1;

	/* compute the phase-shifted chirp. */
	/* Complex conjugate of T0^* exp(-2 pi i f t0) */
	/* See e:twophasesig in the GRASP manual for comparison */
	working[re]=chirp0[re]*cosi+chirp0[im]*sini;
	working[im]=chirp0[im]*cosi-chirp0[re]*sini;

	/* compute the phase-shifted chirp. */
	/* Complex conjugate of [T90^*] exp(-2 pi i f t0) */
	/* See e:twophasesig in the GRASP manual for comparison */
        working2[re]=chirp1[re]*cosi+chirp1[im]*sini;
	working2[im]=chirp1[im]*cosi-chirp1[re]*sini;

	/* real and imaginary parts of time-offset induced phase */
        /* cosine and sine of (n+1)delta in terms of cosine and sine of n delta */
        cosin=cosi*cosd-sini*sind;
        sini=cosi*sind+sini*cosd;
        cosi=cosin;
    }

    /* compute signals in each band */
    j=0;
    varsplit=0.0;

    for (i=0;i<nbands;i++) {
	signalreal=0.0;
	signalimag=0.0;
	while (j<indices[i]) {
            temp=htilde[j]*twice_inv_noise[j/2];
	    signalreal+=temp*working[j];
	    signalimag+=temp*working2[j];
	    j++;
	}

	/* See GRASP equation e:twophasesig to understand the factor of 2*/
	signalreal*=2.0*oldnorm;
	signalimag*=2.0*oldnorm;

	/* remember stats[0..2p-1] has 2p elements */
	im=(re=i+i)+1;
	stats[re]=signalreal;
	stats[im]=signalimag;

	/* the total signal sum */
	sigsumreal+=signalreal;
	sigsumimag+=signalimag;
    }

    varsplit=0.0;
    for (i=0;i<nbands;i++){
	im=(re=i+i)+1;
	signalreal=(sigsumreal/nbands-stats[re]);
	signalimag=(sigsumimag/nbands-stats[im]);
	varsplit+=(signalreal*signalreal+signalimag*signalimag);
    }

    return varsplit;
}


/* This is the splitup_freq routine when the phase of the signal is unknown */
/* By Jolien Creighton <jolien@tapir.caltech.edu> 22 October 1998 */
/* and Bruce Allen <ballen@dirac.phys.uwm.edu> 21 October 1998 */
float splitup_freq5(float c0, float c1, float *chirp0, float *chirp1,
		    float norm, float *twice_inv_noise, int npoint,
		    int offset, int nbands, int *indices, float *stats, 
		    float *working, float *htilde)
{
  int j, i = 0, ir = 0, ii = 1;
  float oldnorm, ssre = 0, ssim = 0, varsplit = 0;
  double angle, cosd, sind, cosin, cosi = 1, sini = 0;

  /* determine the different frequency bands */
  splitup(working,chirp0,twice_inv_noise,npoint,norm,nbands,indices);

  /* match old normalization if desired; if not, set c0 = c1 = sqrt(2) */
  oldnorm = 0.5*sqrt(c0*c0 + c1*c1);

  /* phase corresponding to start time of the signal */
  angle = (-2*M_PI*offset)/npoint;
  cosd = cos(angle);
  sind = sin(angle);

  for (j = 0; j < nbands; ++j) { /* compute signals in each band */
    float sre = 0, sim = 0;
    int jr = j + j;
    int ji = 1 + jr;
    int imax = indices[j];
    while (ir < imax && ii < imax) { /* sum over frequencies in band */
      float tin = twice_inv_noise[i];
      float hre = htilde[ir]*tin;
      float him = htilde[ii]*tin;
      
      /* compute phase-shifted chirp (cf. e:twophasesig in GRASP manual) */
      /* complex conjugate of T0^* exp(-2 pi i f t0) */
      float c0re = chirp0[ir]*cosi + chirp0[ii]*sini;
      float c0im = chirp0[ii]*cosi - chirp0[ir]*sini;
      /* complex conjugate of T90^* exp(-2 pi i f t0) */
      float c1re = chirp1[ir]*cosi + chirp1[ii]*sini;
      float c1im = chirp1[ii]*cosi - chirp1[ir]*sini;
      sre += hre*c0re + him*c0im;
      sim += hre*c1re + him*c1im;
      
      /* increment indices and apply recurrence relation to trig functions */
      i++;
      ii = 1 + (ir = i + i);
      cosin= cosi*cosd - sini*sind;
      sini = cosi*sind + sini*cosd;
      cosi = cosin;
    }

    /* record the total signal in the frequency band */
    sre *= 2*oldnorm; /* factor of two: cf. e:twophasesig in GRASP manual */
    sim *= 2*oldnorm;
    stats[jr] = sre;
    stats[ji] = sim;
    
    /* the total signal sum */
    ssre += sre;
    ssim += sim;
  }

  /* mean signal in each frequency band */
  ssre /= nbands;
  ssim /= nbands;

  /* sum the variances of the signal in each frequency band */
  for (j = 0; j < nbands; ++j) {
    int jr = j + j;
    int ji = 1 + jr;
    float sre = ssre - stats[jr];
    float sim = ssim - stats[ji];
    varsplit += sre*sre + sim*sim;
  }

  return varsplit;
}

void avg_inv_spec(float flo,float srate,int npoint,double decay,double *norm,
                  float *htilde, float* mean_pow_spec,float* twice_inv_noise) {
		int icut,im,re,i;
		double lastnorm,delta;
		float spec;

		/* lower cutoff frequency */
		icut=npoint*flo/srate;
		if (icut<1) icut=1;

		/* remember the last normalization constant */
		lastnorm=(*norm);

		/* set the new value */
		(*norm)=1.0+decay*lastnorm;
		delta=decay*lastnorm/(*norm);

		/* average the (one-sided) power spectra */
		for (i=icut;i<npoint/2;i++) {
			im=(re=i+i)+1;
			spec=2.0*(htilde[re]*htilde[re]+htilde[im]*htilde[im]);
			mean_pow_spec[i]*=delta;
			mean_pow_spec[i]+=spec/(*norm);
		}

		/* inverse of the spectral noise f<fcut */
		for (i=0;i<icut;i++)
			twice_inv_noise[i]=0.0;

		/* inverse of the spectral noise f<fcut */
		for (i=icut;i<npoint/2;i++)
			twice_inv_noise[i]=2.0/mean_pow_spec[i];

		/* inverse of the spectral noise at Nyquist freq */
		twice_inv_noise[npoint/2]=0.0;

		return;
}


void orthonormalize(float* ch0tilde,float* ch1tilde,float* twice_inv_noise,int npoint,float* n0,float* n1) {

		float c00,c10,c11,real0,imag0,real1,imag1,shinv,eps;
		int i,im,re;

		/* Note <noise^2> = 1/2 (Q,Q) = 1/2 c00 = 1/2 c11, i.e. c11=(Q,Q) */
		/* inner product of ch1tilde with ch0tilde and ch1tilde with ch0tilde */
		c11=c10=c00=0.0;
		for (i=0;i<npoint/2;i++) {
			im=(re=i+i)+1;
			real0=ch0tilde[re];
			imag0=ch0tilde[im];
			real1=ch1tilde[re];
			imag1=ch1tilde[im];
			shinv=twice_inv_noise[i];
			c00+=shinv*(real0*real0+imag0*imag0);
			c10+=shinv*(real0*real1+imag0*imag1);
			c11+=shinv*(real1*real1+imag1*imag1);
		}
		/* This loop equivalent to (but faster than!) 
		correlate(output0,ch0tilde,ch0tilde,twice_inv_noise,npoint);
		c00=output0[0];
		correlate(output0,ch1tilde,ch0tilde,twice_inv_noise,npoint);
		c10=output0[0];
		correlate(output0,ch1tilde,ch1tilde,twice_inv_noise,npoint);
		c11=output0[0];
		*/

		/* compute new orthogonal mode */
		eps=c10/c00;
		/* Note: this may result in rounding errors after many of iterations... */
		for (i=0;i<npoint;i++)
			ch1tilde[i]-=eps*ch0tilde[i];

		/* and its inner product with itself */
		c11=c11-2.0*eps*c10+eps*eps*c00;

		/* for <noise^2> = 1, would multiply ch0tilde, ch1tilde by n0, n1 */
		*n0=1.0/sqrt(0.5*c00);
		*n1=1.0/sqrt(0.5*c11);

		return;
}


void find_chirp(float* htilde,float* ch0tilde,float* ch1tilde,float* twice_inv_noise,float n0,float n1,
                float* output0,float* output1,int npoint,int chirplen,int* offset,float* snr_max,float* lin0,float* lin1,float *var) {

		double mean2;
		float max,out0,out1,snr2;
		int maxi=0,i;

		/* compute two signal values */
 		correlate(output0,htilde,ch0tilde,twice_inv_noise,npoint);
		correlate(output1,htilde,ch1tilde,twice_inv_noise,npoint);

		/* look for maximum S/N over two different phase signals */
		mean2=0.0;
		max=0.0;
 		for (i=0;i<npoint-chirplen;i++) {
			out0=n0*output0[i];
			out1=n1*output1[i];

			/* compute snr at that offset */
 			snr2=out0*out0+out1*out1;
			mean2+=snr2;

			/* if larger than previous value, save it, and the offset */
			if (snr2>max) {
				max=snr2;
				maxi=i;
			}
		}
		max*=0.5;
		*var=0.5*mean2/(npoint-chirplen);

		/* maximum signal-to-noise, offset achieved */
		*snr_max=sqrt(max);
		*offset=maxi;

		/* coefficients of best-fit linear combination */
		*lin0=n0*output0[maxi]/((*snr_max)*sqrt(2.0));
		*lin1=n1*output1[maxi]/((*snr_max)*sqrt(2.0));

		return;
}

/* here a and b are FT's of real functions, hence arrays of
length n, c is length n, s is length n/2+1 (all frequency values).
This routine returns:

\int c = a(f) b^*(f) exp(-2 pi i f t) s_real(f)

*/

void correlate(float *c,float *a,float *b,float *s,int n){
	void realft(float*,unsigned long,int);
	float *savec;
	int m;

	/* save pointer to c */
	savec=c;

	/* set m to half of n */
	m=n/2;

	/* put c = a x b* */
	productc(c,a,b,m);

	/* DC term */
	(*c++)*=(*s);

	/* Nyquist term */
	(*c++)*=*(s+m);

	/* all other frequencies */
	s++;

	while (--m>0) {
		(*c++)*=(*s);
		(*c++)*=(*s++);
	}

	realft(savec-1,n,-1);
	return;
}
		
void splitup(float *c,float *a,float *s,int n,float total,int segments,int *indices){
	float *savec;
	double sum,increment;
	int m,which;

	/* save pointer to c */
	savec=c;

	/* set m to half of n */
	m=n/2;

	/* put c = a x a* */
	productc(c,a,a,m);

	/* DC term */
	(*c++)*=(*s);

	/* Nyquist term */
	(*c++)*=*(s+m);

	/* all other frequencies */
	s++;

	while (--m>0) {
		(*c++)*=(*s);
		(*c++)*=(*s++);
	}

	sum=0.0;
	total/=segments;
	increment=total;
	which=0;
	for (m=0;m<n;m+=2) {
		sum+=savec[m];
		if (sum>total) {
			indices[which++]=m;
			total+=increment;
			if (which==segments)
				break;
		}
	}
	indices[segments-1]=n-1;
	return;
}


