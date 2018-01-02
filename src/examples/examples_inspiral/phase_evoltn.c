/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"
int main() {
   float m1,m2,spin1,spin2,phaseterms[MAX_PHS_TERMS],clscnc_time,*ptrphase,*ptrfrequency;
   float time,Initial_Freq,Max_Freq_Rqst,Max_Freq_Actual,Sample_Time,time_in_band;
   int steps_alloc,steps_filld,i,n_phaseterms,err_cd_sprs,chirp_ok;

   /* Set masses and spins of the orbital system:  */
   m1=m2=1.4;
   spin1=spin2=0.;

   /* Set ORBITAL frequency range of the chirp and sample time:  */
   Initial_Freq=60.;                  /* in cycles/second */
   Max_Freq_Rqst=2000.;               /* in cycles/second */
   Sample_Time=1./9868.4208984375;    /* in seconds       */

/* Use this block to compare with Will & Wiseman, PRD 54, 4813 (1996) Figure 10, page 4846. 
   spin1=0.1;           spin2=0.5;
   m1=1.4;              m2=12.0;
   Initial_Freq=75.0;   Max_Freq_Rqst=180.0;       */

   /* post-Newtonian [O(1/c^n)] terms you wish to include (or supress) 
      in the phase and frequency evolution: */
   n_phaseterms=5;             /* the number of entries in phaseterms */
   phaseterms[0] =1.;          /* The Newtonian piece            */
   phaseterms[1] =0.;          /* There is no O(1/c) correction  */
   phaseterms[2] =1.;          /* The post-Newtonian correction  */
   phaseterms[3] =1.;          /* The 3/2 PN correction          */
   phaseterms[4] =1.;          /* The 2 PN correction            */

   /* Set memory-allocation and error-code supression logic: */
   ptrphase=ptrfrequency=NULL;
   err_cd_sprs=0;

   /* Use phase_frequency() to compute phase and frequency evolution: */
   chirp_ok=phase_frequency(m1,m2,spin1,spin2,n_phaseterms,phaseterms,
      Initial_Freq,Max_Freq_Rqst,&Max_Freq_Actual,Sample_Time,&ptrphase,
      &ptrfrequency,&steps_alloc,&steps_filld,&clscnc_time,err_cd_sprs);

   /* ... and print out the results: */
   time_in_band=(float)(steps_filld-1)*Sample_Time;
   fprintf(stderr,"\nm1=%f  m2=%f  Initial_Freq=%f\n", m1,m2,Initial_Freq);
   fprintf(stderr,"steps_filld=%i  steps_alloc=%i  Max_Freq_Actual=%f\n",
            steps_filld,steps_alloc,Max_Freq_Actual);
   fprintf(stderr,"time_in_band=%f clscnc_time=%f\n",time_in_band,clscnc_time);
   fprintf(stderr,"Termnination code: %i\n\n",chirp_ok);

   for (i=0;i<steps_filld;i++){
      time=i*Sample_Time;
      printf("%i\t%f\t%f\t%f\n",i,time,ptrphase[i],ptrfrequency[i]);
   }
return 0;
}
