/* tmwave.c
   This example calculates the two waveforms using the testmass formulas. 
   It saves these in the file waveorm.dat without doing much else. The purpose
   of this example is to demonstrate how to calculate a waveform using GRASP.
     
   Author: S. Droz (droz at physics.uoguelph.ca)
*/


#include "grasp.h"

#ifdef __MACOS__
/* If we run this on a Macintosh compatible machine use
   the console package SIOUX. This is ignored on any other 
   platform. */
#include <sioux.h>

/* Prototypes: */
int PlayAudio(float *f, double Rate, int n);
#endif 

int main()
{
   int NoPo =0;               /* Read in as many data points as possible  */    
   float *x = NULL;           /* We let GRASP take care of all the memory */
   float *Phase = NULL;       /* allocation.                              */
   float *hplus = NULL;
   float *hcross = NULL;
   int NoOfWavePoints = 10000;    /* Calculate as many points as needed. */
   float *f = NULL;
   float fend,fstart,dt;
   int   error, i;
   FILE  *fp;
   float m1 = 1.4;      /* Mass of the first body in units of the solar mass */ 
   float m2 = 1.4;      /* Mass of the second body in units of the solar mass */ 
   float theta = 1.2;   /* The inclination angle in radians */      
   float phi = 0.0;     /* The azimuthal angle in radians */
   
   /* Which modes should we include? (1 include, 0 omit) */         
   /*		m = -5 -4 -3 -2 -1  1  2  3  4  5	   l     */
   int modes[28] = {         1, 1, 1, 1,            /* l = 2 */
                          1, 1, 1, 1, 1, 1,         /* l = 3 */
                       0, 0, 0, 0, 0, 0, 0, 0,      /* l = 4 */
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; /* l = 5 */
   
#ifdef __MACOS__ /* Mac stuff, don't worry */
   SIOUXSettings.asktosaveonclose = 0;
#endif
   
   /* First we have to read in the data files. This will only work if you've 
      set the environment variable GRASP_PARAMETERS. We just read in the default
      files, so we can give NULL as filenames. x[0..NoPo-1] will contain all the
      v - values.  This routine sets up memory for the A_{lm}(v)'s and P(v)
      internally. It will also calculate V(t) and save it.                 */   
   printf(" Reading data ...\n");
   error = ReadData(NULL, NULL, &x, &NoPo);
   if ( error ) return error;
   
   /* We now have to calculate the phase function 
      Phi(f_0,v). This function already knows
      how many points to calculate; the same number 
      as the number of datapoints. Since we want to plot
      the wave function around ~ 100 Hz we set fo = 400.0 */
   printf(" Calculating the phase ...\n");
   error = calculate_testmass_phase(400.0, (m1+m2) ,&Phase);
   if ( error ) printf("Error calculating the phase\n");
   
   
   /* Uncomment the following code if you want to save Phi(f_0,v) */
   /*
   fp = fopen("Phase.dat","w");
   for (i = 0; i < NoPo; i++) 
     fprintf(fp,"%f %f %20.18f\n",x[i], pow(x[i],3.0)/(2.0*(m1+m2)*TSOLAR*Pi), Phase[i]);
   fclose(fp);
   */
   
   
   /* We're now ready to calculate the chirp itself. */
   printf(" Calculating the chirp ...\n");

   dt =Get_Duration(60.0, 785.0,m1,m2)/(1.0*NoOfWavePoints-1.0); /* Set the timestep in seconds */
    
   testmass_chirp(m1, m2, theta, phi , Phase, 60.0, 785.0, &fstart, &fend, 
               dt,  &hplus, &hcross, &f,  &NoOfWavePoints, 3,  modes);
               
   printf(" Calculated %d  data points\nin the frequency intervall [%f, %f].\n",
           NoOfWavePoints,fstart, fend);
   printf(" The chirp lasted %f seconds.\n",dt*NoOfWavePoints);	  
   printf(" Writing data to disk. This might take a few seconds.\n");            
   fp = fopen("waveform.dat","w");
   for (i = 0; i < NoOfWavePoints; i++) 
     fprintf(fp,"%f %f %f %f %f\n", i*dt, f[i] , hplus[i],
     hcross[i],pow(f[i]*2.0*(m1+m2)*TSOLAR*M_PI,1.0/3.0)  );
   fclose(fp);

#ifdef __MACOS__
   printf("Playing wave .... \n");  /* Play the wave */
   error = PlayAudio(hplus, 1.0/dt, Noofcp-1);
#endif
   
   Clean_Up_Memory(Phase); /* Clean up all the memory which was used internally. */
   free(hplus);            /* Get rid of the waveforms */
   free(hcross);
   printf("Goodbye...\n"); /* That's all folks. */
   return error;
}
