/* GRASP: Copyright 1997,1998  Bruce Allen */
static char *rcsid="$Id: pN_chirp.c,v 1.18 1999/04/19 17:41:43 ballen Exp $\
\n$Name: RELEASE_1_9_8 $";

#include "grasp.h"
/* #include "grasp.h" */
/* routine to compute two chirp filters */

#ifdef INLINE_CUBEROOT
float cuberoot(float x);
#endif


/* This code written by Benjamin Owen, 1997 owen@tapir.caltech.edu */

#define MAX_ORDER 5

float sp_phase(float f, float t_c, float m1, float m2, int order);

void sp_filters(float m1, float m2, float *ch1, float *ch2, float fstart,
    int n, float srate, float f_c, float t_c, int order) {
  float M, eta;
  int i;

  /* Check post-Newtonian order */
  if (order < 0 || order > MAX_ORDER) {
    GR_start_error("sp_filters()", rcsid, __FILE__, __LINE__);
    GR_report_error("Can't do post-Newtonian order > %d/2\n", MAX_ORDER);
    GR_end_error();
    abort();
  }

  /* Compute chirp parameters */
  M = m1+m2;
  eta = m1*m2/M/M;

  /* Sanity check on frequencies */
  if (fstart<=0. || fstart>=srate/2.) {
    GR_start_error("sp_filters()", rcsid, __FILE__, __LINE__);
    GR_report_error("Must have 0 < fstart < srate/2\n");
    GR_end_error();
    abort();
  }
  if (f_c<fstart) {
    GR_start_error("sp_filters()", rcsid, __FILE__, __LINE__);
    GR_report_error("f_c < fstart for %.1f+%.1f solar masses, ", m1, m2);
    GR_report_error("chirps are arrays of zeroes\n");
    GR_end_error();
    clear(ch1, n, 1);
    clear(ch2, n, 1);
    return;
  }

  /* Now generate the chirp */
  for (i=0; i<n; i+=2) {
    float f;
    f = i*srate/n/2.;
    if (f < fstart || f > f_c)
      ch1[i] = ch1[i+1] = ch2[i] = ch2[i+1] = 0.;
    else {
      float Psi;
      Psi = sp_phase(f, t_c, m1, m2, order)-M_PI/4.;
      f = pow(f,-7./6.)*pow(TSOLAR,-1./6.)*pow(M/M_PI/M_PI,1./3.)*sqrt(5.*M*eta/96.);
      ch1[i] = f*cos(Psi);
      ch1[i+1] = f*sin(Psi);
      ch2[i] = ch1[i+1];
      ch2[i+1] = -ch1[i];
    }
  }
}

/* This phase includes only frequency-dependent terms */
float sp_phase(float f, float t_c, float m1, float m2, int order) {
  float x, eta, Psi;

  x = pow(M_PI*(m1+m2)*TSOLAR*f, 1./3.);
  eta = m1*m2/pow(m1+m2,2.);
  Psi = pow(x,-5.);
  if (order >= 2)
    Psi += (4.9140212+6.1111111*eta)*pow(x,-3.);
  if (order >= 3)
    Psi += -16.*M_PI*pow(x,-2.);
  if (order >= 4)
    Psi += (30.103153+53.859127*eta+42.847222*eta*eta)/x;
  if (order >= 5)
    Psi += (153.35317+5.*eta)*M_PI*log(x);
  Psi *= 0.0234375/eta;
  Psi += 2.*M_PI*f*t_c;
  return Psi;
}


/* End of code written by Benjamin Owen, 1997 owen@tapir.caltech.edu */


/*  
     In order to form the gravitational wave strain for
     a restricted PN template.

     h_+ = (TSOLAR*CLIGHT/L) * [-(1+cos(incl)^2] *ptrptrCos[i]
     h_x = (TSOLAR*CLIGHT/L) * [-2 cos(incl) ] *ptrptrCos[i]

     where L is the luminostity distance to the source
     and incl  is the inclination of the source relative to the line 
     of sight.

                 PHASE AND FREQUENCY MODULE

   This routine computes the phase and frequency evolution
   of for and inspiral.

*/

#define N_MEMORY_CHUNK 4096  /*chunks of memory (in floats) allocated */
static float loc_freq_coef[MAX_PHS_TERMS],loc_phaseterms[MAX_PHS_TERMS];
static float loc_lhs,loc_eta,loc_m_tot,loc_phase_coal;
static float cfreq0,cfreq2,cfreq3,cfreq4,cfreq5,loc_theta_scl;
static float cphse0,cphse2,cphse3,cphse4,cphse5;
static int L_52PN;
float freq_of_theta(float theta),init_theta_of_freq(float freq);
float phase_of_theta(float theta);
int error_code(int code,int id1,float fd1, float fd2);
int loc_trmntn_code,loc_err_cd_sprs;
void compute_phase_freq_coefs(float m1,float m2,float spin1,float spin2,
                                     float *loc_freq_coef,int N);

int phase_frequency(float m1,float m2,float spin1,float spin2,int n_phaseterms,
               float *phaseterms, float Initial_Freq,float Max_Freq_Rqst,
               float *Max_Freq_Actual,float Sample_Time,float **ptrptrphase,
               float **ptrptrfrequency,int *N_alloc,int *N_filld,
               float *clscnc_time,int err_cd_sprs)

{  /* begin phase_frequency */
float  thetai,initial_time,initial_phase;
float  phase,dumph,freq,theta,time,theta_check,freq_check;
float  freq_prev_step,fdum1=0.0,fdum2=0.0;
int i,L_alloc,L_continue,trmntn_code,idum1=0;
float cfreqscale,timescale;
float p4,p2,p1,pm8;

     /* Sufficent memory for the storing the phase/freq information? */
     if (n_phaseterms>MAX_PHS_TERMS) {
         GR_start_error("phase_frequency()",rcsid,__FILE__,__LINE__);
         GR_report_error("number of phaseterms %d ",n_phaseterms);
         GR_report_error("is greater than that allowed (%d)\n",MAX_PHS_TERMS);
         GR_end_error();
         abort();
     }

        /* Initialization */
	loc_err_cd_sprs=err_cd_sprs;
	trmntn_code=3201;    /* if not reset, there is a logic problem */
	*N_filld=0;

        /* Compute mass parameters (MASSES in units of solar masses) */
        loc_m_tot=m1+m2;
	loc_eta=m1*m2/loc_m_tot/loc_m_tot;
        timescale=loc_eta/loc_m_tot/TSOLAR/5.;
        cfreqscale=1.0/TSOLAR/8./loc_m_tot/2./M_PI;

        /* Set local logic for which terms included in phase/freq  evolution */
        for (i=0;i<n_phaseterms;i++) loc_phaseterms[i]=phaseterms[i];
        /* L_52PN, whether to bother with 2.5PN calculations */
        L_52PN=0;
        if(loc_phaseterms[5]!=0.) L_52PN=1;
       
        /* ... and compute the phase and frequency coeficients. */
        compute_phase_freq_coefs(m1,m2,spin1,spin2,loc_freq_coef,n_phaseterms);

        /* Determine the start-time for the loop */
        thetai=init_theta_of_freq(Initial_Freq);
        if(loc_trmntn_code==3204) return 3204;
        loc_theta_scl=thetai;  /* scale for logarithm in post2.5 corrections */
        loc_lhs=0.;
        initial_time=thetai*(5.*loc_m_tot*TSOLAR/loc_eta);
        *clscnc_time=initial_time;
        loc_phase_coal=0.;
        loc_phase_coal= -phase_of_theta(thetai); /* sets phase=0 at start time */
        initial_phase=phase_of_theta(thetai); /* convoluted, but flexibility often used in dvlpmnt */

        /* determine if chirp starts with increasing frequency (sanity check): */
        theta_check=(initial_time-Sample_Time)*loc_eta/loc_m_tot/TSOLAR/5.;
        freq_check=freq_of_theta(theta_check);
        if((freq_check-Initial_Freq)<0.){ 
            trmntn_code=3203;
            error_code(trmntn_code,idum1,Initial_Freq,loc_m_tot);
            return trmntn_code;
        }
              

/* 
   Much of the  logic in this routine depends on whether we are to allocate
   the memory here (L_alloc=1) or whether the memory has been allocated
   elsewhere (L_alloc=0).  Determine which mode we are opperating in.
*/
     L_alloc=0;
     if(*ptrptrphase==NULL || *ptrptrfrequency==NULL ) L_alloc=1;

     /* Allocate first chunk of memory, if so requested */
     if(L_alloc==1){
        *N_alloc=N_MEMORY_CHUNK;
        *ptrptrphase=(float *)malloc(sizeof(float)*(*N_alloc));
        *ptrptrfrequency=(float *)malloc(sizeof(float)*(*N_alloc));
        /* if unable to allocate this first chunk of memory print error */
        if((*ptrptrphase==NULL) ||  (*ptrptrfrequency==NULL) ){
           trmntn_code=2102;
           error_code(trmntn_code,*N_alloc,fdum1,fdum2);
           return trmntn_code;
        }
     }

     /* Initialize arrays. This takes care of zero padding. */
     for (i=0;i<*N_alloc;i++) (*ptrptrphase)[i]=0.0;
     for (i=0;i<*N_alloc;i++) (*ptrptrfrequency)[i]=0.0;

     /* store the initial phase and frequency and initialize the loop */
     (*ptrptrfrequency)[0]=Initial_Freq;
     (*ptrptrphase)[0]=initial_phase;
     time=initial_time;
     freq_prev_step=freq_of_theta(thetai);
     L_continue=1;



/* Finally!  The Phase-and-Frequency Loop:  
  This is the crankshaft of the engine that powers coalescing-binary
  gravitation-wave signal analysis.  If it breaks you are stranded.  */


for (i=1;L_continue!=0;i++) {
 
    /* If we are out of allocated memory ... */
    if( i>(*N_alloc-1) ) {
        /* and we are not allocating more in this routine: TERMINATE CHIRP */
        if( L_alloc==0 ) {
            trmntn_code=2001;
            error_code(trmntn_code,i,freq_prev_step,Max_Freq_Rqst);
            *N_filld=i;
            *Max_Freq_Actual=freq_prev_step;
            return trmntn_code;
        }
        /* or if we are allocating here: then allocate more ...  */
        if( L_alloc==1 ) {
            *ptrptrphase=(float *)realloc(*ptrptrphase,
                sizeof(float)*(*N_alloc+N_MEMORY_CHUNK));
            *ptrptrfrequency=(float *)realloc(*ptrptrfrequency,
                sizeof(float)*(*N_alloc+N_MEMORY_CHUNK));
            /* unless there is no more to allocate: so TERMINATE CHIRP */
            if((*ptrptrphase==NULL) ||  (*ptrptrfrequency==NULL) ){
                 trmntn_code=2103;
                 error_code(trmntn_code,i,freq_prev_step,Max_Freq_Rqst);
                 *N_filld=i;
                 *Max_Freq_Actual=freq_prev_step;
                 return trmntn_code;
            }
            *N_alloc+=N_MEMORY_CHUNK;
        }
    }

    time=initial_time-((float)i)*Sample_Time;
    theta=time*timescale;

    /* If we have not reached the coalescence time ... */
    if(theta>0.0){
          /* compute the phase and frequency ... */
          /* THE LINES BELOW ARE THE EXACT (BUT OPTIMIZED) EQUIVALENT OF 
             phase=phase_of_theta(theta);  */

          p4=sqrt(theta);
          p2=sqrt(p4);
          p1=sqrt(p2);
          dumph=p1*(p4*cphse0+p2*cphse2+p1*cphse3+cphse4);
             /* if requested the post 2.5 correction */
          if(L_52PN==1) dumph += cphse5*log(theta/loc_theta_scl);
          phase= -dumph/loc_eta+loc_phase_coal;

	  /* THE LINES BELOW ARE THE EXACT (BUT OPTIMIZED) EQUIVALENT OF  
             freq=freq_of_theta(theta); */
          pm8=1.0/theta;
          freq=cfreqscale*pm8*(p1*(p4*cfreq0+p2*cfreq2)+p2*cfreq3+p1*cfreq4);
             /* if requested the post 2.5 correction */
          if(L_52PN==1) freq += cfreqscale*cfreq5/theta;
    }

    /* terminate loop if we have reached coalescence time. */
    else {
          trmntn_code=1202;
          *N_filld=i;
          *Max_Freq_Actual=freq_prev_step;
          error_code(trmntn_code,i,freq_prev_step,Max_Freq_Rqst);
          return trmntn_code;
    }          

    /* LOOK FOR A REASON TO QUIT COMPUTING THE PHASE AND FREQUENCY    */

    /* terminate loop if the frequency is not monotonic increasing ... */
    if(freq<freq_prev_step){
          /* if not, terminate the chirp and send a warning.  */
          *N_filld=i;
          *Max_Freq_Actual=freq_prev_step;
          trmntn_code=1201;
          error_code(trmntn_code,i,freq_prev_step,Max_Freq_Rqst);
          return trmntn_code;
    }

    /* terminate loop if the frequency is increasing to rapidly ... */
    if(freq>1.5*freq_prev_step){
          *N_filld=i;
          *Max_Freq_Actual=freq_prev_step;
          trmntn_code=1203;
          error_code(trmntn_code,i,freq_prev_step,Max_Freq_Rqst);
          return trmntn_code;
    }

    /* terminate loop when PN expansion clearly no good: ie when r<2m_tot... */
    /* at 1-pn order theta^(-1/4)=4(m/r) and p2 is set to theta^(1/4)  */
    if(p2<.5){
          *N_filld=i;
          *Max_Freq_Actual=freq_prev_step;
          trmntn_code=1204;
          error_code(trmntn_code,i,freq_prev_step,Max_Freq_Rqst);
          return trmntn_code;
    }

    /* Normal Termination  if  the maximum requested frequency is reached. */
    if(freq>=Max_Freq_Rqst) {
          /* this is the only place the the termination code gets set to 0 */
          *N_filld=i;
          *Max_Freq_Actual=freq_prev_step;
          trmntn_code=0;
          return trmntn_code;
    }

    /* Couldn't find a reason to discontinue loop, so store info and repeat.*/
    (*ptrptrfrequency)[i]=freq;
    (*ptrptrphase)[i]=phase;
    freq_prev_step=freq;

} /* End Phase and Frequency Loop */


return trmntn_code;
} /* end phase_frequency */

/*  routines for finding time as function of freq. */
float init_theta_of_freq(float freq)
{       float rtbis(float (*func)(float),float x1,float x2,float xacc);
        float  th,oldth=0.0,g,final,acc,theta,gold;
        int id1=0;
        loc_lhs=freq;
        /* First, a loop to find a ROUGH BRACKET of the  root.
        Start with Netonian value of theta which gives the frequency/8 */
        th=pow(1./TSOLAR/loc_m_tot/(freq/8.)/16./M_PI,8./3.);
        g=freq_of_theta(th);
        gold=g;
        if (g==0.0) return th;  /* lucky! guess the root right off */
        /* robust, but rough estimator of the root. Start way to the right
        and come in looking for it. See fig21 in Grasp manual. 
        Rough bracketing comes from the right or the left, BUT NOT BOTH,
        thus the strange logical structure of this section.
        There should be no logic fault even if the root happens to
        be local min, max or inflection point. */
        if (g<0.0) while (g<0.0) {
                         oldth=th;
                         th/=1.05;
                         g=freq_of_theta(th);
                         if(gold>g){
                           loc_trmntn_code=3204;
                           error_code(loc_trmntn_code,id1,loc_m_tot,freq);
                         }
                         else gold=g;
                   }
        /* initial guess was so conservative that root shouldn't be to
        the right, but we will try it anyway */
        else       while (g>0.0) {
                         oldth=th;
                         th*=1.05;
                         g=freq_of_theta(th);
                         if(gold<g){
                           loc_trmntn_code=3204;
                           error_code(loc_trmntn_code,id1,loc_m_tot,freq);
                         }
                         else gold=g;
                     }
        /*End ROUGH BRACKET of root. Root now bracketed between (oldth,th).
        Call the root finder to finish it off. */
        acc=1.e-12*(oldth+th);
        final=rtbis(&freq_of_theta,oldth,th,acc);
        theta=final;
return theta;}

/* ################################################################### */
/*   Routine freq_theta requires modification to incorpotate higher    */
/*               post-Newtonian corrections (beyond 2.5PN).            */
/* ################################################################### */
/* This is a direct implementaion of Equation 8 in BIWW */
float freq_of_theta(float theta)
{
float freq,dumfr,pm1,pm2,pm3,pm4,pm5,pm6,pm7;
int trmntn_code;
	if(theta<=0.){
		trmntn_code=3002;
		error_code(trmntn_code,0,theta,0.0);
		exit(1);
	}
	pm4=sqrt(1./theta);
	pm2=sqrt(pm4);
	pm1=sqrt(pm2);
	pm3=pm1*pm2;
	pm5=pm1*pm4;
	pm6=pm2*pm4;
	pm7=pm6*pm1;
	dumfr=pm3*(loc_freq_coef[0])*(loc_phaseterms[0])
	     +pm5*(loc_freq_coef[2])*(loc_phaseterms[2])
	     +pm6*(loc_freq_coef[3])*(loc_phaseterms[3])
	     +pm7*(loc_freq_coef[4])*(loc_phaseterms[4]);
             if(L_52PN==1) dumfr += 
                   (loc_freq_coef[5])*(loc_phaseterms[5])/theta;
	     /* additional correction to phase would be added here */
	freq= (dumfr/TSOLAR/8./loc_m_tot/2./M_PI);
        /* note: this routine is used both to compute the freq as a 
        fcn of time, but also it is called by the root finder to find the
        time (ie theta) when the freq =, say, 60Hz. Therefore it returns 
        the freq minus the static variable loc_lhs. If you want the frequency,
        set loc_lhs =0. */
return freq-loc_lhs;
}

/* ################################################################### */
/*   Routine phase_of_theta requires modification to incorpotate       */
/*           higher post-Newtonian corrections (beyond 2.5PN).         */
/* ################################################################### */
/* This is a direct implementaion of Equation 7 in BIWW */
float phase_of_theta(float theta)
{
float phase,dumph,p1,p2,p3,p4,p5;
int trmntn_code;
	if(theta<=0.){
	       trmntn_code=3003;
	       error_code(trmntn_code,0,theta,0.0);
	       exit(1);
	}

	p4=sqrt(theta);
	p2=sqrt(p4);
	p1=sqrt(p2);
	p3=p1*p2;
	p5=p1*p4;
	dumph=p5*(loc_freq_coef[0])*(loc_phaseterms[0])
	     +p3*(loc_freq_coef[2])*(loc_phaseterms[2])*(5./3.)
	     +p2*(loc_freq_coef[3])*(loc_phaseterms[3])*(5./2.)
	     +p1*(loc_freq_coef[4])*(loc_phaseterms[4])*(5.);
             /* if requested the post 2.5 correction */
              if(L_52PN==1) dumph += 
                (loc_freq_coef[5])*(loc_phaseterms[5])*log(theta/loc_theta_scl);
	     /* additional correction to phase would be added here */

	phase= -dumph/loc_eta +loc_phase_coal;
        /* note: this routine is used both to compute the phase as a 
        fcn of time, but also could be used by a root finder to find the
        time (ie theta) when the phase  =, say, 2pi. Therefore it returns 
        the phase minus the static variable loc_phase_coal.  During the
        actual chirp calculation the loc_phase_col is adjusted to so
        that phase(t=0)=0. */
return phase;
}


/* ################################################################### */
/*  Routine phase_of_freq requires modification to incorpotate higher  */
/*               post-Newtonian corrections (beyond 2.5PN).            */
/* ################################################################### */
void compute_phase_freq_coefs(float m1,float m2,float spin1,float spin2,
                              float *loc_freq_coef,int N)
{
float eta,dmom,chi_s,chi_a;
int i;
      eta=m1*m2/(m1+m2)/(m1+m2);
      dmom=(m1-m2)/(m1+m2);         /* (m1-m2)/mtot, GRASP convention m1<=m2 */
      chi_s=.5*(spin1+spin2);       /* symmetric spin quantity */
      chi_a=.5*(spin1-spin2);       /* antisymmetric spin quantity */

      /* initialize the coeficients to zero  */
      for (i=0;i<N;i++) loc_freq_coef[i]=0;

/* 
     Enter the coeficients in the frequecny formula:  This is version
     is intended to be a faithful representation of the formulae in
     BIWW equation 8.  The indexing goes as the O[1/c^n] correction.
     The spin corrections can be found in Will and Wiseman PRD 54,pp
     4813, equation F22a.
*/

   loc_freq_coef[0]=1.0;                           /* Leading order             */
   loc_freq_coef[1]=0.0;                           /* No O(1/c) correction      */
   loc_freq_coef[2]=743.0/2688.0+11.0*eta/32.0;    /* 1-PN correction O(1/c^2)  */
   loc_freq_coef[3]= -3.0*M_PI/10.0                 /* 1.5-PN Tail Term O(1/c^3) */
   +(113./160.)*(chi_s+dmom*chi_a)-19.*eta*chi_s/40.; /* plus spin-orbit terms  */
   loc_freq_coef[4]=1855099.0/14450688.0
       +56975.0*eta/258048.0+371.0*eta*eta/2048.0  /* 2-PN correction O(1/c^4)  */
       -(237./512.)*eta*(chi_s*chi_s-chi_a*chi_a); /* plus 2PN spin-spin term   */
   loc_freq_coef[5]= -(7729./21504.+3.*eta/256.)*M_PI;  /* 2.5PN order O(1/c^5)  */
/*
      + ... 
      Put additional Post-Newtonian and spin corections here.
*/


       /* The following coefficients are redundant, 
       but used for more efficiant in-line code      */
       cphse0=(loc_freq_coef[0])*(loc_phaseterms[0]);
       cphse2=(loc_freq_coef[2])*(loc_phaseterms[2])*(5./3.);
       cphse3=(loc_freq_coef[3])*(loc_phaseterms[3])*(5./2.);
       cphse4=(loc_freq_coef[4])*(loc_phaseterms[4])*(5.);
       cphse5=(loc_freq_coef[5])*(loc_phaseterms[5])*(5./8.);

       cfreq0=(loc_freq_coef[0])*(loc_phaseterms[0]);
       cfreq2=(loc_freq_coef[2])*(loc_phaseterms[2]);
       cfreq3=(loc_freq_coef[3])*(loc_phaseterms[3]);
       cfreq4=(loc_freq_coef[4])*(loc_phaseterms[4]);
       cfreq5=(loc_freq_coef[5])*(loc_phaseterms[5]);
 
return;
}


/*              Error Code Printing Routine                */

int error_code(int code,int idum1,float fdum1, float fdum2)
{

/*  
     Index for termination codes.

     Normal termination code:                            0
     Termination with a warning (w/o memory allocation): 1001-1099
     Termination with a warning (w   memory allocation): 1101-1199
     Termination with a warning (indep of memory alloc): 1201-1199
     Serious Problems (w/o memory allocation):           2001-2099
     Serious Problems (w   memory allocation):           2101-2199
     Serious Problems (indep of memory allocation):      2201-2299
     Pure bone-head FATAL terminations:                  3001-...

*/

  if(code==1201){
    if(code>loc_err_cd_sprs){
    GR_start_error("phase_frequency()",rcsid,__FILE__,__LINE__);
    GR_report_error("Frequency evolution no longer monotonic.\n");
    GR_report_error("Phase evolution terminated at frequency and step: %f\t%i\n",
            fdum1,idum1);
    GR_report_error("Terminating chirp. Termination code set to:\t%i\n",code);
    GR_report_error("Returning to calling routine.\n");
    GR_end_error();
    }
    return 0;
   }

  if(code==1202){
    if(code>loc_err_cd_sprs){
    GR_start_error("phase_frequency()",rcsid,__FILE__,__LINE__);
    GR_report_error("Time parameter is less than zero.\n");
    GR_report_error("Indicates that next time step is beyond coalescence.\n");
    GR_report_error("Terminating chirp. Termination code set to:\t%i\n",code);
    GR_report_error("Returning to calling routine.\n");
    GR_end_error();
    }
    return 0;
  }

  if(code==1203){
    if(code>loc_err_cd_sprs){
    GR_start_error("phase_frequency()",rcsid,__FILE__,__LINE__);
    GR_report_error("Frequency rising too rapidly.\n");
    GR_report_error("Indicates post Newtonian calculation not valid.\n");
    GR_report_error("Terminating chirp. Termination code set to:\t%i\n",code);
    GR_report_error("Returning to calling routine.\n");
    GR_end_error();
    }
    return 0;
  }

  if(code==1204){
    if(code>loc_err_cd_sprs){
    GR_start_error("phase_frequency()",rcsid,__FILE__,__LINE__);
    GR_report_error("Post Newtonian Equations no longer\n");
    GR_report_error("valid. [(Gm/rc^2)>1]\n");
    GR_report_error("Terminating chirp. Termination code set to:\t%i\n",code);
    GR_report_error("Returning to calling routine.\n");
    GR_end_error();
    }
    return 0;
  }

  if(code==2001){
    if(code>loc_err_cd_sprs){
    GR_start_error("phase_frequency()",rcsid,__FILE__,__LINE__);
    GR_report_error("Allocated memory is filled up before\n");
    GR_report_error("reaching the maximum frequency reqested for this chirp.\n");
    GR_report_error("Freq Reached: %f, Number of points: %i\n",fdum1,idum1);
    GR_report_error("Terminating chirp. Termination code set to:\t%i\n",code);
    GR_report_error("Returning to calling routine.\n");
    GR_end_error();
    }
    return 0;
  }

  if(code==2102){
    if(code>loc_err_cd_sprs){
    GR_start_error("phase_frequency()",rcsid,__FILE__,__LINE__);
    GR_report_error("Insufficiant memory to store requested\n");
    GR_report_error("first memory chunk. Floats required for chunk:\t%i\n",2*idum1);
    GR_report_error("Terminating chirp. Termination code set to:\t%i\n",code);
    GR_report_error("Returning to calling routine.\n");
    GR_end_error();
    }
    return 0;
  }

  if(code==2103){
    if(code>loc_err_cd_sprs){
    GR_start_error("phase_frequency()",rcsid,__FILE__,__LINE__);
    GR_report_error("Insufficiant memory to store requested\n");
    GR_report_error("chirp. Floats required so far:\t%i\n",2*idum1);
    GR_report_error("Reached frequency:\t%f\n",fdum1);
    GR_report_error("Terminating chirp. Termination code set to:\t%i\n",code);
    GR_report_error("Returning to calling routine.\n");
    GR_end_error();
    }
    return 0;
  }

  if(code==3201){
    if(code>loc_err_cd_sprs){
    GR_start_error("phase_frequency()",rcsid,__FILE__,__LINE__);
    GR_report_error("Error code not set in phase_frequency().\n");
    GR_report_error("This indicates logical flaw.\n");
    GR_report_error("Termination code set to:\t%i\n",code);
    GR_report_error("Exiting to system.\n");
    GR_end_error();
    }
    return 0;
  }

  if(code==3002){
    if(code>loc_err_cd_sprs){
    GR_start_error("freq_of_theta()",rcsid,__FILE__,__LINE__);
    GR_report_error("Invalid time parameter (negative theta)\n");
    GR_report_error("has been sent. theta:\t%f\n",fdum1);
    GR_report_error("Death is nigh. Termination code set to:\t%i\n",code);
    GR_report_error("Exiting to system.\n");
    GR_end_error();
    }
    return 0;
  }

  if(code==3003){
    if(code>loc_err_cd_sprs){
    GR_start_error("phase_of_theta()",rcsid,__FILE__,__LINE__);
    GR_report_error("Invalid time parameter (negative theta)\n");
    GR_report_error("has been sent. theta:\t%f\n",fdum1);
    GR_report_error("Death is nigh. Termination code set to:\t%i\n",code);
    GR_report_error("Exiting to system.\n");
    GR_end_error();
    }
    return 0;
  }

  if(code==3203){
    if(code>loc_err_cd_sprs){
    GR_start_error("phase_frequency()",rcsid,__FILE__,__LINE__);
    GR_report_error("This chirp does not even start with a \n");
    GR_report_error("monotonically increasing frequency. Indicates too large initial\n");
    GR_report_error("for masses or problems with root finder for initial time.\n");
    GR_report_error("Initial Frequency:\t%f\t Total Mass\t%f\n",fdum1,fdum2);
    GR_report_error("Termination code set to:\t%i\n",code);
    GR_report_error("Returning to calling routine.\n");
    GR_end_error();
    }
    return 0;
  }

  if(code==3204){
    if(code>loc_err_cd_sprs){
    GR_start_error("phase_frequency()",rcsid,__FILE__,__LINE__);
    GR_report_error("No meaningful start time for this \n");
    GR_report_error("chirp for this mass and initial frequency.:\t%i\n",2*idum1);
    GR_report_error("Initial Frequency:\t%f\t Total Mass\t%f\n",fdum2,fdum1);
    GR_report_error("Termination code set to:\t%i\n",code);
    GR_report_error("Returning to calling routine.\n");
    GR_end_error();
    }
    return 0;
  }

GR_start_error("error_code()",rcsid,__FILE__,__LINE__);
GR_report_error("Error code has not been implemented:\t%i.\n",code);
GR_end_error();
exit(1);
return 0;
}


void make_filters(float m1,float m2,float *ch1,float* ch2,float
	fstart,int steps_alloc,float srate,int *filled,float *t_coal,int err_cd_sprs,int order) {

	float spin1,spin2,phaseterms[MAX_PHS_TERMS],Initial_Freq,Max_Freq_Rqst,Max_Freq_Actual,Sample_Time;
	int steps_filld,n_phaseterms,filter_ok,i;

	if (order<0 || order >5) {
		GR_start_error("make_filters()",rcsid,__FILE__,__LINE__);
		GR_report_error("The order argument of make_filters() is the order in (v/c) past the quadrupole approximation.\n");
		GR_report_error("If the order is set to zero, this gives the quadrupole approximation, and\n"); 
		GR_report_error("if set to four, this gives the second-post-Newtonian approximation.\n"); 
		GR_report_error("The order is currently set to %d; the allowed range is 0 to 5.\n",order);
		GR_end_error();
		abort();
	}

	/* Set physical parameters of the orbital system:  */
	spin1=0.;
	spin2=0.;

	/* Set ORBITAL frequency range of the chirp and sample time:  */
	Initial_Freq=fstart/2.0;			/* in cycles/second */
	Max_Freq_Rqst=srate/4.0;			/* in cycles/second */
	Sample_Time=1.0/srate;			  /* in seconds		 */

	/* Establish what post-Newtonian [O(1/c^n)] terms you wish to include
		in the phase and frequency evolution: */
	n_phaseterms=MAX_PHS_TERMS;
	for (i=0;i<n_phaseterms;i++) phaseterms[i]=0.0;
	for (i=0;i<=order;i++) phaseterms[i]=1.0;
	phaseterms[1]=0.;		 /* There is no O(1/c) correction */

	/* Use chirp_filters() to compute the two filters: */
	 filter_ok=chirp_filters(m1,m2,spin1,spin2,n_phaseterms,phaseterms,
		Initial_Freq,Max_Freq_Rqst,&Max_Freq_Actual,Sample_Time,
		&ch1,&ch2,&steps_alloc,&steps_filld,t_coal,err_cd_sprs);

	/* And print out the results: */
	/*
	printf("\nsteps_filld %i\t steps_alloc %i\n\n",steps_filld,steps_alloc);
	printf("Termination code: %i\n\n",filter_ok);
	*/
 	*filled=steps_filld;

	/* zero pad the remainder of the array */
	for (i=steps_filld;i<steps_alloc;i++) ch1[i]=0.0;
	for (i=steps_filld;i<steps_alloc;i++) ch2[i]=0.0;
	return;
}


#ifdef INLINE_TRIGS
/* Chebyshev polynomials for expansion of sin/cos on -pi < x < pi */
static const float cs0=1.000000e+00, cs2=-4.999999e-01, cs4=4.166650e-02,
cs6=-1.388787e-03, cs8=2.477176e-05, cs10=-2.709995e-07,
cs12=1.732713e-09;

static const float cs1=1.000000e+00, cs3=-1.666669e-01, cs5=8.333445e-03,
cs7=-1.984473e-04, cs9=2.760671e-06, cs11=-2.530742e-08,
cs13=1.544217e-10;
#endif


int chirp_filters(float m1,float m2,float spin1,float spin2,int n_phaseterms,
              float *phaseterms, float Initial_Freq,float Max_Freq_Rqst,
              float *ptrMax_Freq_Actual,float Sample_Time,
              float **ptrptrCos,float **ptrptrSin, int *N_alloc,
              int *N_filld,float *clscnc_time,int err_cd_sprs)
{

	float *ptrphase,*ptrfrequency,m_red,m_tot,normalize,phs,phs_const,factor,y,si,co;
	int i,chirp_ok;
#ifdef INLINE_CUBEROOT
	static unsigned int psave=123456789;
	unsigned int p;
	static float scale1,scale2;
	extern const float gcheb0,gcheb1,gcheb2,gcheb3,gcheb4,gcheb5;
	extern const float scaling_cuberoot[];
	float x,z;
#endif
#ifdef INLINE_TRIGS
	float phs2;
	double offset=0.0;
#endif
        ptrphase=*ptrptrCos;
        ptrfrequency=*ptrptrSin;

        /* compute phase and frequency */
        /* Use phase_frequency() to compute phase and frequency evolution: */
        chirp_ok=phase_frequency(m1,m2,spin1,spin2,n_phaseterms,phaseterms,
           Initial_Freq,Max_Freq_Rqst,ptrMax_Freq_Actual,Sample_Time,
           &ptrphase,&ptrfrequency,N_alloc,N_filld,
           clscnc_time,err_cd_sprs);

        /* compute the filters */
        m_tot=m1+m2;
        m_red=m1*m2/m_tot;
        normalize=2.*m_red*pow((TSOLAR*m_tot*2.*M_PI),2./3.);
        phs_const=4.*TSOLAR*m_tot;
     
        /* notice this loop destroys phase and frequency information */
        for (i=0;i<*N_filld;i++){
#ifdef INLINE_CUBEROOT
		/* HARDWIRED INLINE CUBE ROOT (see cuberoot()) */
		/* y=cuberoot(ptrfrequency[i]); */
		x=ptrfrequency[i];
		p=(*(int *)&x) >> 23;
		if (p!=psave) {
			psave=p;
			scale1=scaling_cuberoot[2*p];
			scale2=scaling_cuberoot[2*p+1];
		}
		x*=scale1;
		z=gcheb0+x*(gcheb1+x*(gcheb2+x*(gcheb3+x*(gcheb4+x*gcheb5))));
		y=x*z;
		y*=(4.0f-z*y*y)*scale2;
		/* END OF HARDWIRED INLINE CUBE ROOT (see cuberoot()) */
#else
		y=pow(ptrfrequency[i],1.0/3.0);
#endif
		factor=y*y*normalize;
		phs=ptrphase[i];
		/* log correction for phase commented out */
		/*phs-=phs_const*(ptrfrequency[i])*log(ptrfrequency[i]/Initial_Freq);*/
#ifdef INLINE_TRIGS
		/*  THE FOLLOWING LINES ARE APPROXIMATIONS OF SIN() and COS()  */
		phs=2.0*phs+offset;
		if (phs>=M_PI) {
			offset-=2.0*M_PI;
			phs-=2.0*M_PI;
		}
		phs2=phs*phs;

		si=phs*(cs1+phs2*(cs3+phs2*(cs5+phs2*(cs7+phs2*(cs9+phs2*(cs11+phs2*cs13))))));
		co=cs0+phs2*(cs2+phs2*(cs4+phs2*(cs6+phs2*(cs8+phs2*(cs10+phs2*cs12)))));
#else
		si=sin(2.0*phs);
		co=cos(2.0*phs);
#endif
		ptrphase[i]=factor*co;
		ptrfrequency[i]=factor*si;
        }
	*ptrptrCos=ptrphase;
	*ptrptrSin=ptrfrequency;
return chirp_ok;
}


#ifdef INLINE_CUBEROOT
float cuberoot(float x) {
	static unsigned int psave=4096;
	static float scale1,scale2;
	unsigned int p;
	float y,z;
	extern const float gcheb0,gcheb1,gcheb2,gcheb3,gcheb4,gcheb5;
	extern const float scaling_cuberoot[];

	/* the exponent of the float rep of x */
	p=(*(int *)&x) >> 23;

	/* if exponent hasn't changed, neither have scale factors */
	if (p!=psave) {
		psave=p;
		scale1=scaling_cuberoot[2*p];
		scale2=scaling_cuberoot[2*p+1];
	}
	
	/* rescale x to the range from 1.0->2.0 inclusive */
	x*=scale1;

	/* fifth-order fit to z=x^(-2/3) with uniform-error Chebyshevs */
	z=gcheb0+x*(gcheb1+x*(gcheb2+x*(gcheb3+x*(gcheb4+x*gcheb5))));

	/* y ~ x^(-2/3) x ~ x^(1/3) */
	y=x*z;

	/* Now do one pass of Newton-Raphson on f(z)=z^-3-x^2
	   so dz=-f/f' = -(z^-3 - x^2)/(-3 z^-4)
           so z'= z+dz = z(4 - z^3 x^2)/3 = z(4- y^2 z)/3
           This pass doubles # of decimal places of precision.
        */
	y*=(4.0f-z*y*y)*scale2;

	return y;
}

/* Coefficients of Chebyshev polynomial approx for x^(-2/3) on interval
   1<x<2.  Error about part in 10^4.
*/
const float 
	gcheb0=2.8567583422f,
	gcheb1=-4.0127830607f,
	gcheb2=3.4804896648f,
	gcheb3=-1.7355624528f,
	gcheb4=0.4620750592f,
	gcheb5=-0.05099705856f;

/* Table to reduce domain of input for Chebychev polynomial approx.
   First column is (1/2)^k.  Second column 1/3 2^(-k/3) */
const float scaling_cuberoot[2*256] = {
1.701412e+38,	6.015553e-14,
8.507059e+37,	7.579123e-14,
4.253530e+37,	9.549096e-14,
2.126765e+37,	1.203111e-13,
1.063382e+37,	1.515825e-13,
5.316912e+36,	1.909819e-13,
2.658456e+36,	2.406221e-13,
1.329228e+36,	3.031649e-13,
6.646140e+35,	3.819638e-13,
3.323070e+35,	4.812443e-13,
1.661535e+35,	6.063298e-13,
8.307675e+34,	7.639277e-13,
4.153837e+34,	9.624886e-13,
2.076919e+34,	1.212660e-12,
1.038459e+34,	1.527855e-12,
5.192297e+33,	1.924977e-12,
2.596148e+33,	2.425319e-12,
1.298074e+33,	3.055711e-12,
6.490371e+32,	3.849954e-12,
3.245186e+32,	4.850639e-12,
1.622593e+32,	6.111421e-12,
8.112964e+31,	7.699908e-12,
4.056482e+31,	9.701277e-12,
2.028241e+31,	1.222284e-11,
1.014120e+31,	1.539982e-11,
5.070602e+30,	1.940255e-11,
2.535301e+30,	2.444569e-11,
1.267651e+30,	3.079963e-11,
6.338253e+29,	3.880511e-11,
3.169127e+29,	4.889137e-11,
1.584563e+29,	6.159927e-11,
7.922816e+28,	7.761022e-11,
3.961408e+28,	9.778274e-11,
1.980704e+28,	1.231985e-10,
9.903520e+27,	1.552204e-10,
4.951760e+27,	1.955655e-10,
2.475880e+27,	2.463971e-10,
1.237940e+27,	3.104409e-10,
6.189700e+26,	3.911310e-10,
3.094850e+26,	4.927941e-10,
1.547425e+26,	6.208817e-10,
7.737125e+25,	7.822619e-10,
3.868563e+25,	9.855883e-10,
1.934281e+25,	1.241763e-09,
9.671407e+24,	1.564524e-09,
4.835703e+24,	1.971177e-09,
2.417852e+24,	2.483527e-09,
1.208926e+24,	3.129048e-09,
6.044629e+23,	3.942353e-09,
3.022315e+23,	4.967054e-09,
1.511157e+23,	6.258095e-09,
7.555786e+22,	7.884706e-09,
3.777893e+22,	9.934108e-09,
1.888947e+22,	1.251619e-08,
9.444733e+21,	1.576941e-08,
4.722366e+21,	1.986822e-08,
2.361183e+21,	2.503238e-08,
1.180592e+21,	3.153882e-08,
5.902958e+20,	3.973643e-08,
2.951479e+20,	5.006476e-08,
1.475740e+20,	6.307765e-08,
7.378698e+19,	7.947286e-08,
3.689349e+19,	1.001295e-07,
1.844674e+19,	1.261553e-07,
9.223372e+18,	1.589457e-07,
4.611686e+18,	2.002591e-07,
2.305843e+18,	2.523106e-07,
1.152922e+18,	3.178914e-07,
5.764608e+17,	4.005181e-07,
2.882304e+17,	5.046212e-07,
1.441152e+17,	6.357829e-07,
7.205759e+16,	8.010362e-07,
3.602880e+16,	1.009242e-06,
1.801440e+16,	1.271566e-06,
9.007199e+15,	1.602072e-06,
4.503600e+15,	2.018485e-06,
2.251800e+15,	2.543132e-06,
1.125900e+15,	3.204145e-06,
5.629500e+14,	4.036970e-06,
2.814750e+14,	5.086263e-06,
1.407375e+14,	6.408290e-06,
7.036874e+13,	8.073939e-06,
3.518437e+13,	1.017253e-05,
1.759219e+13,	1.281658e-05,
8.796093e+12,	1.614788e-05,
4.398047e+12,	2.034505e-05,
2.199023e+12,	2.563316e-05,
1.099512e+12,	3.229576e-05,
5.497558e+11,	4.069011e-05,
2.748779e+11,	5.126632e-05,
1.374390e+11,	6.459151e-05,
6.871948e+10,	8.138021e-05,
3.435974e+10,	1.025326e-04,
1.717987e+10,	1.291830e-04,
8.589935e+09,	1.627604e-04,
4.294967e+09,	2.050653e-04,
2.147484e+09,	2.583661e-04,
1.073742e+09,	3.255208e-04,
5.368709e+08,	4.101305e-04,
2.684355e+08,	5.167321e-04,
1.342177e+08,	6.510417e-04,
6.710886e+07,	8.202611e-04,
3.355443e+07,	1.033464e-03,
1.677722e+07,	1.302083e-03,
8.388608e+06,	1.640522e-03,
4.194304e+06,	2.066928e-03,
2.097152e+06,	2.604167e-03,
1.048576e+06,	3.281044e-03,
5.242880e+05,	4.133857e-03,
2.621440e+05,	5.208333e-03,
1.310720e+05,	6.562089e-03,
6.553600e+04,	8.267714e-03,
3.276800e+04,	1.041667e-02,
1.638400e+04,	1.312418e-02,
8.192000e+03,	1.653543e-02,
4.096000e+03,	2.083333e-02,
2.048000e+03,	2.624835e-02,
1.024000e+03,	3.307085e-02,
5.120000e+02,	4.166667e-02,
2.560000e+02,	5.249671e-02,
1.280000e+02,	6.614171e-02,
6.400000e+01,	8.333334e-02,
3.200000e+01,	1.049934e-01,
1.600000e+01,	1.322834e-01,
8.000000e+00,	1.666667e-01,
4.000000e+00,	2.099868e-01,
2.000000e+00,	2.645668e-01,
1.000000e+00,	3.333333e-01,
5.000000e-01,	4.199737e-01,
2.500000e-01,	5.291337e-01,
1.250000e-01,	6.666667e-01,
6.250000e-02,	8.399473e-01,
3.125000e-02,	1.058267e+00,
1.562500e-02,	1.333333e+00,
7.812500e-03,	1.679895e+00,
3.906250e-03,	2.116535e+00,
1.953125e-03,	2.666667e+00,
9.765625e-04,	3.359789e+00,
4.882812e-04,	4.233069e+00,
2.441406e-04,	5.333333e+00,
1.220703e-04,	6.719579e+00,
6.103516e-05,	8.466139e+00,
3.051758e-05,	1.066667e+01,
1.525879e-05,	1.343916e+01,
7.629395e-06,	1.693228e+01,
3.814697e-06,	2.133333e+01,
1.907349e-06,	2.687831e+01,
9.536743e-07,	3.386456e+01,
4.768372e-07,	4.266667e+01,
2.384186e-07,	5.375663e+01,
1.192093e-07,	6.772911e+01,
5.960464e-08,	8.533334e+01,
2.980232e-08,	1.075133e+02,
1.490116e-08,	1.354582e+02,
7.450581e-09,	1.706667e+02,
3.725290e-09,	2.150265e+02,
1.862645e-09,	2.709164e+02,
9.313226e-10,	3.413333e+02,
4.656613e-10,	4.300530e+02,
2.328306e-10,	5.418329e+02,
1.164153e-10,	6.826667e+02,
5.820766e-11,	8.601061e+02,
2.910383e-11,	1.083666e+03,
1.455192e-11,	1.365333e+03,
7.275958e-12,	1.720212e+03,
3.637979e-12,	2.167332e+03,
1.818989e-12,	2.730667e+03,
9.094947e-13,	3.440424e+03,
4.547474e-13,	4.334663e+03,
2.273737e-13,	5.461333e+03,
1.136868e-13,	6.880849e+03,
5.684342e-14,	8.669326e+03,
2.842171e-14,	1.092267e+04,
1.421085e-14,	1.376170e+04,
7.105427e-15,	1.733865e+04,
3.552714e-15,	2.184533e+04,
1.776357e-15,	2.752339e+04,
8.881784e-16,	3.467730e+04,
4.440892e-16,	4.369067e+04,
2.220446e-16,	5.504679e+04,
1.110223e-16,	6.935461e+04,
5.551115e-17,	8.738134e+04,
2.775558e-17,	1.100936e+05,
1.387779e-17,	1.387092e+05,
6.938894e-18,	1.747627e+05,
3.469447e-18,	2.201872e+05,
1.734723e-18,	2.774184e+05,
8.673617e-19,	3.495253e+05,
4.336809e-19,	4.403743e+05,
2.168404e-19,	5.548369e+05,
1.084202e-19,	6.990507e+05,
5.421011e-20,	8.807486e+05,
2.710505e-20,	1.109674e+06,
1.355253e-20,	1.398101e+06,
6.776264e-21,	1.761497e+06,
3.388132e-21,	2.219348e+06,
1.694066e-21,	2.796203e+06,
8.470329e-22,	3.522994e+06,
4.235165e-22,	4.438695e+06,
2.117582e-22,	5.592406e+06,
1.058791e-22,	7.045989e+06,
5.293956e-23,	8.877390e+06,
2.646978e-23,	1.118481e+07,
1.323489e-23,	1.409198e+07,
6.617445e-24,	1.775478e+07,
3.308722e-24,	2.236962e+07,
1.654361e-24,	2.818396e+07,
8.271806e-25,	3.550956e+07,
4.135903e-25,	4.473924e+07,
2.067952e-25,	5.636791e+07,
1.033976e-25,	7.101912e+07,
5.169879e-26,	8.947849e+07,
2.584939e-26,	1.127358e+08,
1.292470e-26,	1.420382e+08,
6.462349e-27,	1.789570e+08,
3.231174e-27,	2.254716e+08,
1.615587e-27,	2.840765e+08,
8.077936e-28,	3.579140e+08,
4.038968e-28,	4.509433e+08,
2.019484e-28,	5.681530e+08,
1.009742e-28,	7.158279e+08,
5.048710e-29,	9.018866e+08,
2.524355e-29,	1.136306e+09,
1.262177e-29,	1.431656e+09,
6.310887e-30,	1.803773e+09,
3.155444e-30,	2.272612e+09,
1.577722e-30,	2.863312e+09,
7.888609e-31,	3.607546e+09,
3.944305e-31,	4.545224e+09,
1.972152e-31,	5.726623e+09,
9.860761e-32,	7.215093e+09,
4.930381e-32,	9.090447e+09,
2.465190e-32,	1.145325e+10,
1.232595e-32,	1.443019e+10,
6.162976e-33,	1.818089e+10,
3.081488e-33,	2.290649e+10,
1.540744e-33,	2.886037e+10,
7.703720e-34,	3.636179e+10,
3.851860e-34,	4.581299e+10,
1.925930e-34,	5.772074e+10,
9.629650e-35,	7.272358e+10,
4.814825e-35,	9.162597e+10,
2.407412e-35,	1.154415e+11,
1.203706e-35,	1.454472e+11,
6.018531e-36,	1.832519e+11,
3.009266e-36,	2.308830e+11,
1.504633e-36,	2.908943e+11,
7.523164e-37,	3.665039e+11,
3.761582e-37,	4.617659e+11,
1.880791e-37,	5.817886e+11,
9.403955e-38,	7.330078e+11,
4.701977e-38,	9.235319e+11,
2.350989e-38,	1.163577e+12,
1.175494e-38,	1.466016e+12,
5.877472e-39,	1.847064e+12,
2.938736e-39,	2.327155e+12};
#endif

