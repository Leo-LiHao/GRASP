#include "grasp.h"
#include "tfmain.h"

extern struct_signalparameters snpar;
extern struct_tfparam tfparam;
extern dl_options dlopt;

void picturerawprint(float **pic)
{
        /* output the TF map as as a linear array to the file "picture"
         the matrix is scanned from top to bottom and from left to right.*/  
    char *rawfile="picture";
    FILE *fp;
    int i,j;
    
    fp = fopen(rawfile,"w");
    for(i=0;i<PDIM;i++)
        for(j=0;j<PDIM;j++)
            fprintf(fp,"%f\n",*(pic[i] + j));
    fclose(fp);
}

void getshf(int n,float *shf,float deltaf)
/* initialize the power spectrum shf[]. shf[0] contains the PSD at 0 frequency
   and shf[i] contains the power spectrum at frequency i*deltaf */
{
    double mini;
    
    if((tfparam.noisetype==NOISE_WHITE))
    {
    int i;   
    for (i=0;i<=n;i++)
        shf[i] = 10000.0;
    }
    if(tfparam.noisetype==NOISE_LIGO_INI)
    {
        double *shfd;
        int i;
        
        shfd = (double *)malloc((n+1) * sizeof(double));
        noise_power("noise_ligo_init.dat",n+1,deltaf,shfd);
        for (i=0;i<=n;i++)shfd[i] *= tfparam.hscale;
        for (i=0;i<=n;i++)shfd[i] *= tfparam.hscale;
        mini=shfd[0];
        for(i=0;i<n;i++) if(shfd[i]<mini) mini = shfd[i];
        for(i=0;i<n;i++) if(shfd[i]>(100000.0*mini)) shfd[i]=(100000.0*mini);
        for (i=0;i<=n;i++) shf[i]=(float)shfd[i]*1600000.0;
        free (shfd);
    }
}

void get_coalescence(float *ch, int np,  float fa , float fs ,int *filled)
{
    FILE *fp;
    int i=0;
    
    if((fp = fopen("MergeSig.dat","r"))==NULL){
        printf("the file containing the signal is not present\n");
        exit(-1);
    }
    while((fscanf(fp,"%f",&ch[i])==1)&&(i<NDIM))
        i++;
    *filled=i;
    fclose(fp);
    for (i=*filled;i<NDIM;i++)
        ch[i] = 0.0;
    return;
}

void fill_data_with_signal(int n,float *a,float *b,float tempfloat)
{
    int i,j,noofsub;

    for(i=0;i<DATADIM;i++) *(b+i) = 0.0;
        /* number of subsegments in each segment of data */
    noofsub = (DATADIM-(POSTSAFETY+PRESAFETY))/NDIM;
    for(i=1;i<=noofsub;i++){
        for(j=0;j<NDIM;j++){
            *(b + NDIM*(i-1) + PRESAFETY + j) = tempfloat  *  (*(a + PRESAFETY + j));
        }
    }
}

void gettfparameters()
{
    FILE *fpin;
    char dummy[1024];
    
        /* get the parameters from the file tfmain.in */
    fpin = fopen("tfmain.in","r");
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%d\n",&(tfparam.run_number));
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%f\n",&(tfparam.f_lower));
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%d\n",&(tfparam.start_segment));
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%d\n",&(tfparam.transformtype));
    fscanf(fpin,"%s\n",dummy);    
    fscanf(fpin,"%d\n",&(tfparam.windowwidth));
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%d\n",&(tfparam.offset_step_size));
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%d\n",&(snpar.signaltype));
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%d\n",&(snpar.signaloffset));
    fscanf(fpin,"%s\n",dummy);    
    fscanf(fpin,"%f\n",&(snpar.m1));
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%f\n",&(snpar.m2));
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%f\n",&(snpar.pind));
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%f\n",&(snpar.ampind));
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%f\n",&(snpar.timfrac));
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%f\n",&(snpar.snr));
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%d\n",&(tfparam.num_of_segments));
    tfparam.DIM = DATADIM;
    tfparam.ND = NDIM;
    tfparam.PD = PDIM;
    tfparam.PRE = PRESAFETY;
    tfparam.POST = POSTSAFETY;
    tfparam.TD = (NDIM/PDIM);
    tfparam.FD = (NDIM/(4*PDIM));
    fscanf(fpin,"%s\n",dummy);    
    fscanf(fpin,"%lf\n",&(dlopt.sigma));
    fscanf(fpin,"%s\n",dummy);
    fscanf(fpin,"%lf\n",&(dlopt.high));
    fclose(fpin);
    dlopt.low = dlopt.high/3.0;
    snpar.signaloffset += PRESAFETY;
    tfparam.hscale = 1.e21;
    snpar.cons_noise_pow = 1.e8;
    tfparam.noisetype = NOISE_TYPE;
    tfparam.srate = SRATE;
    tfparam.rescale_factor = RESCALE_FACTOR;
    snpar.freqfrac = 0.1;
}
 
void noise_gau_col_fr(long *idum,unsigned long nr, float *randarr, float *shf_root)
{
    int i;
    float mygasdev(long *);
    
    for(i=0;i<nr;i++)
        randarr[i] = mygasdev(idum);
    for(i=1;i<nr/2;i++){
        randarr[2*i]   *= shf_root[i];
        randarr[2*i+1] *= shf_root[i];
    }
    randarr[0] *= shf_root[0];
    randarr[1] *= shf_root[nr/2];
}

float mygasdev(long *idum)
{
        float ran2(long *idum);
        static int iset=0;
        static float gset;
        float fac,rsq,v1,v2;

        if  (iset == 0) {
                do {
                    v1=2.0*ran2(idum)-1.0;
                    v2=2.0*ran2(idum)-1.0;
                    rsq=v1*v1+v2*v2;
                } while (rsq >= 1.0 || rsq == 0.0);
                fac=sqrt(-2.0*log(rsq)/rsq);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
            iset=0;
            return gset;
        }
}

void timstat(int seg, FILE *fp, float *h)
{
    int i,j,noofsub,k,kmax;
    float max,tmp;
    
        /* the number of subsegments in each segment of data */
    noofsub = (DATADIM-(POSTSAFETY+PRESAFETY))/NDIM;
    for (i=0;i<noofsub;i++){
        j = i*NDIM + PRESAFETY;
        max = 0.0;
        kmax = 0;
        for(k=0;k<NDIM;k++){
            tmp = *(h + j + k) * *(h + j + k);
            if(tmp>max) {
                max = tmp;
                kmax = k;
            }
            
        }
        max = sqrt(max);
        fprintf(fp,"%d %d %f\n", seg*noofsub + i, kmax, max);
    }
    return;
}

void whiten_filter(float *htilde,long npoint, float *inv_noise)
{
        /* A filter to whiten the noise */
    int i;

    htilde[0] *= sqrt(inv_noise[0]);
    htilde[1] *= sqrt(inv_noise[npoint/2]);
    for (i=1;i<npoint/2;i++){
            htilde[2*i] *= sqrt(inv_noise[i]);
            htilde[2*i+1] *= sqrt(inv_noise[i]);
    }
}

void over_whiten_filter(float *htilde, long npoint, float *inv_noise)
{
            /* A filter to whiten the noise */
    int i;

    htilde[0] *= inv_noise[0];
    htilde[1] *= inv_noise[npoint/2];
    for (i=1;i<npoint/2;i++){
            htilde[2*i] *= inv_noise[i];
            htilde[2*i+1] *= inv_noise[i];
    }
    
}
