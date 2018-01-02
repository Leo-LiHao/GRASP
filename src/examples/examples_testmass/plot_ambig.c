/* plot_ambig.c
   Calculate a series of values of the ambiguity function,
   using 2 pN waveforms as templates and a mode calculated from
   black hole perturbation theory as signal.
     
   Author: S. Droz (droz at physics.uoguelph.ca)
*/

#include "grasp.h"

/* Prototypes: */
void realft(float *, unsigned long, int); 
float norm(float* T,float* twice_inv_noise,int npoint);

float norm(float* That,float* twice_inv_noise,int npoint)
/* Calculate $\int df / S(|f|) T(f) * T^*(f)$  or, in the notaion
   of the manual (T/Sh , T/Sh) = <T,T>. */
{
    int i,im,re;
	float real, imag,c=0;
	
	/* This loop is equivalent to (but faster than!) 
		correlate(output0,That,That,twice_inv_noise,npoint);
		c=output0[0]; */
     		
	for (i=1;i<npoint/2;i++)  /* Negclect the DC and fc values */
	{
		im=(re=i+i)+1;
		real=That[re];
		imag=That[im];
		c+=twice_inv_noise[i]*(real*real+imag*imag);
	}
	return sqrt(c);  /* Note that the 2 from 2/S compensates for the fact that
	                    we only sum over positive frequencies. */
}

	 
int main()
{
   int NoPo =0;               /* Read in as many data points as possible  */    
   float *x = NULL;           /* We let GRASP take care of all the memory */
   float *Phase = NULL;       /* allocation.                              */
   float *hplus = NULL;
   float *hcross = NULL;
   float hNorm;
   int Npoints = 32768;               /* 2^20 points */
   int NoOfWavePoints = 10000;        /* The number of points we want saved  */
   int NoOfPointsGen;
   float *f = NULL;
   float fend,dt;
   int   error, i,j;
   FILE  *fp;
   float m1 = 4.5;            /* Mass of the first body in solar masses */ 
   float m2 = 4.5;            /* Mass of the second body in solar masses */ 
   float eta,Mc,e,m,xx,yy;
   float fstart = 70.0;        /* Starting ORBITAL frequency. */   
   float theta = 1.2;          /* Pick an angle */
   float phi = 0.0;              
   float *pN0, *pN90,n0,n90,c0,c90,*output0,*output90;
   float SNR,var;
   int offset;
   float *twice_inv_noise;
   double *temp;
   float t_coal=0;
   float MaxAmb = 0.0; 
   float Mx=0.0,My=0.0;  
   
   /* Which modes should we include? (1 include, 0 omit) */         
   /*          m = -5 -4 -3 -2 -1  1  2  3  4  5 	   l     */
   int modes[28] = {         1, 1, 1, 1,			/* l = 2 */
                          1, 1, 1, 1, 1, 1,         /* l = 3 */
                       0, 0, 0, 0, 0, 0, 0, 0,      /* l = 4 */
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; /* l = 5 */
   
   /* First we have to read in the data file. This will only work if you've
     set the environment variable GRASP_PARAMETERS. We just read in the default
     files, so we can give NULL as filenames. x[0..NoPo-1] will contain all the
     v - values. This routine sets up memory for the A_{lm}(v)'s and 
     P(v) internally. It will also calculate V(t) and save it.           */   
   printf(" Reading data ...\n");
   error = ReadData(NULL, NULL, &x, &NoPo);
   if ( error ) return error;
   
   /* We now have to calculate the phase function 
      Phi(f_0,v). This function already knows
      how many points to calculate; the same number 
      as we've read datapoints. Since we are interested in freqencies 
      of a couple of  100 Hz we set f_0 = 200.0 */
   printf(" Calculating the phase ...\n");
   error = calculate_testmass_phase(200.0, (m1+m2) ,&Phase);
   if ( error ) printf("Error calculating the phase\n");

   /* We're now ready to calculate the chirp itself. */
   printf(" Calculating the chirp ...\n");
   printf(" MaxF = %f -> T =  %e\n",Get_Fmax(m1,m2),Get_Duration(fstart, Get_Fmax(m1,m2),m1,m2));
   dt = Get_Duration(fstart, Get_Fmax(m1,m2),m1,m2)/(NoOfWavePoints-1); /* Set the timestep in seconds */
   printf(" dt = %e\n", dt);
   NoOfWavePoints = Npoints;
   
   testmass_chirp(m1, m2, theta, phi , Phase, fstart ,Get_Fmax(m1,m2)-10, &fstart, &fend, 
               dt,  &hplus, &hcross, &f,  &NoOfWavePoints, 3,  modes);      
   Clean_Up_Memory(Phase); /* Clean up all the memory which was used internally. */
   free(hcross); /* We don't need htimes. */
   printf(" Calculated %d  data points\n in the frequency intervall [%f, %f].\n",
           NoOfWavePoints,fstart, fend);
   printf(" The cirp lasted %f seconds.\n",dt*NoOfWavePoints);	 
   /* Zero out the remaining points */
   clear(hplus +  NoOfWavePoints, Npoints-NoOfWavePoints,1);
   
   /* Get the spectral desnity */
   
   twice_inv_noise = (float *)malloc((Npoints/2+1)*sizeof(float));
   temp = (double *)malloc((Npoints/2+1)*sizeof(double));
   if ( ! ( temp && twice_inv_noise )) return -1; /* Not enough memory */

   noise_power("noise_40smooth.dat", Npoints/2, 1.0/(Npoints*dt), temp);
   for (i=0; i< Npoints/2 ; i++) twice_inv_noise[i] = (float)(2.0e-31/temp[i]);
   free(temp);
   
   /* Allocate memory for the templates, etc.  */
   pN0  = (float *)malloc(Npoints*sizeof(float));
   pN90  = (float *)malloc(Npoints*sizeof(float));
   output0  = (float *)malloc(Npoints*sizeof(float));
   output90 = (float *)malloc(Npoints*sizeof(float));
   if ( ! ( pN0 && pN90 && output0 && output90)) return -1; /* Not enough memory */
     
   
   realft(hplus-1,Npoints,1);   /* FFT the signal */
   hNorm = norm(hplus, twice_inv_noise, Npoints); /* Get the signal's norm */

 
   /* Now get ready to loop over the mass range. We use the chirp mass and
      mass ratio eta as parameters. */
   Mc = pow(m1*m1*m1*m2*m2*m2/(m1+m2),1.0/5.0);
   eta = m1*m2/pow(m1+m2,2); 
   printf(" Chirpmass Mc = %e Msun, eta = %e\n",Mc,eta);
   fp = fopen("scan.dat","w");
   for (i = 0; i<=50; i++)
   {
     for (j = 0; j<=50;j++)
     {
        xx = (0.25 + j*(1.0 - 0.25)/50.0);    /* Deviation from the `true' value */
        yy = (1.00 + i*(1.3 - 1.0)/50.0);
        e = eta*xx;
        m = Mc*yy;
        m1 = 0.5*m*pow(e,-3.0/5.0)*(1-sqrt(1-4.0*e));
        m2 = 0.5*m*pow(e,-3.0/5.0)*(1+sqrt(1-4.0*e)); 
        /* Use make_filters to make the templates, then FFT and orthonormalize. */
		make_filters(m1, m2, pN0, pN90, 2.0*fstart, Npoints, 1.0/dt ,&NoOfPointsGen, 
          &t_coal, 2000, 4);  
        realft(pN0-1,Npoints,1);
        realft(pN90-1,Npoints,1);
     	orthonormalize(pN0,pN90, twice_inv_noise, Npoints,&n0,&n90); 
        
		find_chirp(hplus,pN0,pN90,twice_inv_noise,n0,n90,output0, output90 , Npoints, NoOfWavePoints, 
              &offset, &SNR, &c0,&c90,&var);
        if ( SNR > MaxAmb) 
        {
          MaxAmb = SNR;
          Mx = xx;
          My = yy;
        }
        fprintf(fp, "%e %e %e\n",xx,yy, SNR/hNorm); /* Save {\cal A} */		
      }
      fprintf(fp,"\n");fflush(fp);  
      printf(".");fflush(stdout);
    }  
   fclose(fp);
   MaxAmb /= hNorm;
      
   /* Clean up the remaining memory and exit */
   free(hplus);            /* Get rid of the waveforms */
   free(pN0);
   free(pN90);
   free(output0);
   free(output90);
   free(twice_inv_noise);
   printf("\n The maximum of A_i in the scaned intervall was %4.1f%% and occured at eta*=%4.3f, Mc*=%4.3f \n",100*MaxAmb,Mx,My); 
   printf("\nGoodbye...\n"); /* That's it folks. */
   return error;
}
