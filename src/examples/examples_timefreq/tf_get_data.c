#include "grasp.h"
#include "tfmain.h"

float *twice_inv_noise,*ch0,*ch1,*buff,*shf,*shf_root,deltaf;
extern float *htilde,srate;
extern int npoint;
extern long longn;
extern struct_tfparam tfparam;
extern struct_signalparameters snpar;
long idummy=-89473884;      /* initialized to a random seed*/

int get_time_series_data()
{
    return gethtilde();
}

int gethtilde()
{
    int i,order=1,err_cd_sprs=4000,filled;
    static int first=1;
    float t_coal,tempfloat;
    FILE *randfp;
    
    if(first){
        first = 0;
            /* read in the random number seed */
        if((randfp = fopen("randomseeds","r"))==NULL){
            printf("the randomseeds file is not present\n");
            printf("Please create a file called randomseeds which contains\n");
            printf("a column of negative random numbers seeds \n");
            exit(-1);
        }
        for(i=0;i<tfparam.run_number;i++)
            if(fscanf(randfp,"%ld\n",&idummy)<0){
                printf("the randomseeds file does not contain enough random numbers\n");
                exit(-1);
            }
        fclose(randfp);
            /* allocate memory for the inverse power spectrum, signal arrays and buffers*/
        twice_inv_noise = (float *) malloc(sizeof(float)*(npoint/2 + 1));
        ch0 = (float *) malloc(sizeof(float)*npoint);
        ch1 = (float *) malloc(sizeof(float)*npoint);
        buff = (float *) malloc(sizeof(float)*npoint);
        shf = (float *)malloc (sizeof(float)*(npoint/2+1));
        shf_root = (float *)malloc (sizeof(float)*(npoint/2+1));
            /* initialize the signal arrays to zero */
        for(i=0;i<npoint;i++) ch0[i]=ch1[i]=0.0;
            /* switch between the signal types */
        switch(snpar.signaltype){
            case INSERT_QUAS_PER:
                gen_quasiperiodic_signal(ch1,NDIM,tfparam.f_lower,srate,snpar.pind, snpar.ampind,snpar.timfrac,
                                         snpar.freqfrac,&filled);
                break;
            case INSERT_INSPIRAL:
                make_filters(snpar.m1,snpar.m2,ch0,ch1,tfparam.f_lower,npoint,
                             srate,&filled,&t_coal, err_cd_sprs, order);
                break;
            case INSERT_COALESCENCE:
                get_coalescence(ch1,npoint,tfparam.f_lower,srate,&filled);
                break;
        }
            /* offset the signal to the right by tfparam.signaloffset points */
        for(i=filled-1;i>=0;i--) ch1[i+snpar.signaloffset] = ch1[i];
        for(i=0;i<snpar.signaloffset;i++) ch1[i] = 0.0;
            /* copy ch1 to ch0 */
        for(i=0;i<npoint;i++)
            ch0[i] = ch1[i];
            /* take the Fourier transform of the signal*/
        realft(ch0-1,longn,1);
            /* get the power spectrum for noise */
        deltaf = srate/npoint;
        getshf(npoint/2+1,shf,deltaf);
        for(i=0;i<=npoint/2;i++) twice_inv_noise[i] = 1./shf[i];
        for(i=0;i<=npoint/2;i++) shf_root[i] = sqrt(shf[i]);
            /* normalize the signal to a particular SNR*/
        correlate(buff,ch0,ch0,twice_inv_noise,npoint);
        tempfloat = snpar.snr/sqrt(buff[0]);
        fill_data_with_signal(npoint,ch1,ch0,tempfloat);
        realft(ch0-1,longn,1);
    }
    noise_gau_col_fr(&idummy, npoint, htilde, shf_root);
        /* add the signal to the noise */
    if(snpar.addsignal)
        for(i=0;i<npoint;i++) htilde[i] += ch0[i];
    over_whiten_filter(htilde,npoint,twice_inv_noise);
        /* zero out the higher frequency to avoid aliasing */
    if((tfparam.transformtype==WIGNERTF)||(tfparam.transformtype==CHOIWILLIAMS))
        for(i=npoint/2;i<npoint;i++) htilde[i] = 0.0;
    realft(htilde - 1, longn, -1);
    return DATADIM;
}
