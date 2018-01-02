#include "grasp.h"
void realft(float *, unsigned long, int);

/* this file contains code for computing the various time frequency distributions */

void time_freq_map(float *htilde, struct_tfparam *tfs, int ind, float **tfgeneral,float **pic)
{
    int i,j;
    
    for(i=0;i<(*tfs).ND;i+=(*tfs).TD){
      for(j=0; j < (*tfs).TD/(*tfs).offset_step_size; j++){
	/* branch depending on the type of the distribution */
	switch((*tfs).transformtype){
	case WIGNERTF:
	  wigner((*tfs).ND, htilde + ind*(*tfs).ND + (*tfs).PRE, tfgeneral[j], \
		 i+(*tfs).offset_step_size*j);
	  break;
	case WFFTWTF:
	  windowfft((*tfs).ND, htilde + ind*(*tfs).ND + (*tfs).PRE, tfgeneral[j],\
		    i+(*tfs).offset_step_size*j,tfs);
	  break;
	case CHOIWILLIAMS:
	  choiwill((*tfs).ND, htilde + ind*(*tfs).ND + (*tfs).PRE,tfgeneral[j],\
		   i+(*tfs).offset_step_size*j,tfs);
	  break;  
	case WIGNERTF_NP:
	  wigner_np((*tfs).ND, htilde + ind*(*tfs).ND + (*tfs).PRE, tfgeneral[j], \
		 i+(*tfs).offset_step_size*j);
	  break;	  
	}  
      }
      /* average the time frequency region to a raster of dimension (*tfs).PD*/
      averageblock(tfgeneral,pic[i/(*tfs).TD],tfs);
    }   
}

void wigner(unsigned long n, float *inp, float *out, int wig_offset)
/* computes the wigner tranform of the array inp
   assuming it is the time representation
   Only alternate elements of array inp are used and so the array
   "out" contains only n/2 elements.
   Zeropadding is done.
   */
{
    int i,j,k;
        /* out(t,tau) = inp(t-tau/2) * inp(t+tau/2)*/
    for(i=0;i<n;i+=2){
        j = wig_offset - i/2;
        k = wig_offset + i/2;
        if((j<0)||(k>n))
            out[i/2] =0.0;
        else
            out[i/2] = inp[j]*inp[k];
    }
    realft(out-1,n/2,1);
        /* zero out the negative elements */
    for(i=0;i<=n/2;i+=1)
        if(out[i]<0) out[i]=0.0;
}

void wigner_np(unsigned long n, float *inp, float *out, int wig_offset)
/* computes the wigner tranform of the array inp
   assuming it is the time representation
   Only alternate elements of array inp are used and so the array
   "out" contains only n/2 elements.
   Zero padding is NOT done
   */
{
    int i,j,k;
        /* out(t,tau) = inp(t-tau/2) * inp(t+tau/2)*/
    for(i=0;i<n;i+=2){
        j = wig_offset - i/2;
        k = wig_offset + i/2;
	out[i/2] = inp[j]*inp[k];
    }
    realft(out-1,n/2,1);
        /* zero out the negative elements */
    for(i=0;i<=n/2;i+=1)
        if(out[i]<0) out[i]=0.0;
}



void windowfft(long n, float *wave, float *waveft, int offset,struct_tfparam *tfs)
{
/*  Compute the windowed fourier transform with a Welch Window of width
    "winwidth"
 */
    float window;
    int j,winwidth;
    
    winwidth = (*tfs).windowwidth;
        /* window data */
    for(j=0;j<offset;j++)
        waveft[j]=0;
    for(j=offset;((j<offset+winwidth) && (j<n));j++){   
        window=1.0-4.0*((j-offset)-winwidth/2.0)*((j-offset)-winwidth/2.0)
            /(float)(winwidth*winwidth);
        waveft[j]=wave[j]*window;
    }
    for(j=offset+winwidth;j<n;j++)
        waveft[j]=0;
        /* FFT the data  */
    realft(waveft-1,n,1);
}



void choiwill(unsigned long n, float *inp, float *out, int offset, struct_tfparam *tfs)
{
/* Computes the Choiwilliams of the array "inp".
   only alternate elements of array inp are used and so the array
   "out" contains only n/2 elements.
   */
    
    int j,k,l,m,width;
    float temp,tmp1,tmp2,tmp3,tmp4;
    double sum=0;

    width = (*tfs).windowwidth;
    tmp2 = -1.0*PI*PI*width/(*tfs).srate;
    tmp3 = sqrt(PI*width/(*tfs).srate);
    tmp4 = tmp3/(*tfs).srate;
    for(j=0;j<n;j+=2){
        sum=0.0;
        l=j/2;
        m=n-l;
        tmp1 = 1.0/(j+1);
        for(k=l;k<m;k++){
            temp =  (k-offset)*tmp1;
            temp *= temp*tmp2;
            if(temp>-60.0)  /* other wise the exp() term is too small */
                sum += inp[k-l]*inp[k+l]*exp(temp);
        }
        out[l] = sum*tmp4*tmp1;
    }
    realft(out-1,n/2,1);
       /* zero out the negative elements */
    for(j=2;j<=n/2;j++)
        if(out[j]<0) out[j]=0.0;
}


void averageblock(float **tfb, float *raster,struct_tfparam *tfs)
{
    switch((*tfs).transformtype){
        case WIGNERTF:
        {
            int i,j,k,l;
            double sum;
            for(i=0;i<(*tfs).PD;i++){
                sum=0.0;
                for(j=0;j<(*tfs).TD/(*tfs).offset_step_size;j++){
                    for(k=0;k<(*tfs).FD;k++){
                        l = 2*(i*(*tfs).FD + k);
                        sum += *(tfb[j] + l);
                    }
                }
                raster[i] = (sum*(*tfs).offset_step_size)/((*tfs).FD*(*tfs).TD);
            }
        }
        break;

        case WIGNERTF_NP:
        {
            int i,j,k,l;
            double sum;
            for(i=0;i<(*tfs).PD;i++){
                sum=0.0;
                for(j=0;j<(*tfs).TD/(*tfs).offset_step_size;j++){
                    for(k=0;k<(*tfs).FD;k++){
                        l = 2*(i*(*tfs).FD + k);
                        sum += *(tfb[j] + l);
                    }
                }
                raster[i] = (sum*(*tfs).offset_step_size)/((*tfs).FD*(*tfs).TD);
            }
        }
        break;
                
        case WFFTWTF:
        {
            int i,j,k,l,m;
            double sum,val;
            
            for(i=0;i<(*tfs).PD;i++){
                sum=0.0;
                for(j=0;j<(*tfs).TD/(*tfs).offset_step_size;j++){
                    for(k=0;k<(*tfs).FD;k++){
                        l = 2*(i*(*tfs).FD + k);
                        m=l+1;
                        val=((*(tfb[j]+l))*(*(tfb[j]+l)) + (*(tfb[j]+m))*(*(tfb[j]+m)));
                        sum += val;
                    }
               }
                raster[i] = (sum*((*tfs).offset_step_size))/((*tfs).FD*(*tfs).TD);
            }
        }
        break;
        case CHOIWILLIAMS:
        {
            int i,j,k,l;
            double sum;
            for(i=0;i<(*tfs).PD;i++){
                sum=0.0;
                for(j=0;j<(*tfs).TD/(*tfs).offset_step_size;j++){
                    for(k=0;k<(*tfs).FD;k++){
                        l = 2*(i*(*tfs).FD + k);
                        sum += *(tfb[j] + l);
                    }
                }
                raster[i] = (sum*(*tfs).offset_step_size)/((*tfs).FD*(*tfs).TD);
            }
        }
        break;
                
    }    
}

void normalize_picture(float **pic,int PDIM)
{
    float max,min;
    int i,j;
    
        /* uniformly scale the values in the tf picture to values between 0 and 1 */
    max=min=*(pic[0]);
    for(i=0;i<PDIM;i++)
        for(j=0;j<PDIM;j++){
            if(*(pic[i] + j)>max) max = *(pic[i] + j);
            if(*(pic[i] + j)<min) min = *(pic[i] + j);
        }
    if(max==min){
        for(i=0;i<PDIM;i++)
            for(j=0;j<PDIM;j++)
                *(pic[i] + j) = 0.5;
        return;
    }        
    for(i=0;i<PDIM;i++)
        for(j=0;j<PDIM;j++)
            *(pic[i] + j) = (*(pic[i] + j) - min)/(max-min);   
}

void rescale(float **pic,int PDIM, float rescaleval)
{
        /* to rescale the picture to fit within the range of a character variable */
    int i,j;
    
    for(i=0;i<PDIM;i++){
        for(j=0;j<PDIM;j++){
            *(pic[i]+j) /= rescaleval;
        }
    }
    for(i=0;i<PDIM;i++)
        for(j=0;j<PDIM;j++)
            if(*(pic[i] + j)>1.0)  *(pic[i] + j) = 1.0;   
}


float compute_scalefactor(float **pic,float rescale_factor, int pdim)
{
    float max=0.0;
    int i,j;

    for (i=0;i<pdim;i++)
        for(j=0;j<pdim;j++)
            if(max<*(pic[i] + j)) max = *(pic[i] + j);
    return (rescale_factor*max);
}


void gen_quasiperiodic_signal(float *ch0,int np,float fa, float fs, float pind, float ampind,
                              float timfrac, float freqfrac, int *filled)
{
    float fup,ph,freqc,cons,tim,amp;
    int i;
    
    fup = fs * freqfrac;
        /* freq = const * t^(pind) */
        /* amplitude goes as f^(ampind) */
        /* where the const is decided by np, fa, fup fsamp*/
    cons = (np/fs)*timfrac;
    cons=1./cons;
    cons = (fup-fa)*pow(cons,pind);
    *filled = timfrac*np - 1;
    for(i=0;i<np;i++){
        if(i<*filled){
            tim=i/fs;
            if(tim<=0.0) freqc=fa;
            else
                freqc = fa + cons * pow(tim,pind);
            amp = pow(freqc,ampind);
            if(tim<=0.0) ph = 0.0;
            else
                ph = 2.*PI*(fa*tim + cons*pow(tim,pind+1.0)/(pind + 1.));
            ch0[i] = amp * sin(ph);
        }
        else {
            ch0[i] = 0.0;
        }
    }
}


#define RGBABS1 .3
#define RGBABS2 .6
#define RGBABS3 .73
#define RGBABS4 .93
#define RGBABS5  1.0



void huetorgb(float magn, float *colr, float *colg, float *colb)
{
    
    static float l10=1./RGBABS1,l21=1./(RGBABS2-RGBABS1),l32=1./(RGBABS3-RGBABS2);
    static float l43=1./(RGBABS2-RGBABS1),l54=1./(RGBABS5-RGBABS4);
    
     
    if(magn<RGBABS1){
        *colr = 1;
        *colg = magn*l10;
        *colb = 0.0;
        return;
    }
    if(magn<RGBABS2){
        *colr = 1 - (magn-RGBABS1)*l21;
        *colg = 1.0;
        *colb = 0.0;
        return;
    }
    if(magn<RGBABS3){
        *colr = 0;
        
        *colg = 1.;
        *colb=(magn-RGBABS2)*l32;
        return;
    }
    if(magn<RGBABS4){
        *colr = 0;
        *colg = 1.0 - (magn-RGBABS3)*l43;
        *colb = 1.0;
        return;
    }    
    {
        *colr=(magn-RGBABS4)*l54;
        *colg=0.0;
        *colb = 1.0;
        return;
    }
}



void pgmprint(float **pic,char *pgmfile, int pdim)
{
        /* output the TF map as a PGM file to file picture.pgm */
        /* assumes the picture contains values between 0 and 1 */
    FILE *fp;
    int i,j;
    unsigned char p1,p2;
    
    fp = fopen (pgmfile,"w");
    fprintf(fp,"%s\n","P2");
    fprintf(fp,"%s\n","#Creator program: GRASP:timefrequency graph");
    fprintf(fp,"%d %d \n",pdim,pdim);
    fprintf(fp,"%d \n",255);
    for(i=pdim-1;i>=0;i--){
        for(j=0;j<pdim;j++){
            p1 =*(pic[j]+i)*255;
                /*p2 = 255 - p1;*/
            p2=p1; /* color inversion */
            fprintf(fp,"%3d ",p2);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

void ppmprint(float **pic,char *ppmfile, int pdim)
{
     /* output the TF map as a PPM file to file picture.ppm  */
    FILE *fp;
    int i,j;
    float colr,colg,colb;
    unsigned char icolr,icolg,icolb;
    
    fp = fopen ( ppmfile,"w");
    fprintf(fp,"%s\n","P3");
    fprintf(fp,"%s\n","#Creator program: GRASP: timefrequency graph");
    fprintf(fp,"%d %d \n",pdim,pdim);
    fprintf(fp,"%d \n",255);
    for(i=pdim-1;i>=0;i--){
        for(j=0;j<pdim;j++){
            huetorgb(*(pic[j]+i),&colr,&colg,&colb);
            icolr=255*colr;
            icolg=255*colg;
            icolb=255*colb;
            fprintf(fp,"%3d %3d %3d\n",icolr,icolg,icolb);
        }
    }
    fclose(fp);
}
















