/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main() {
  short *datas;
  int i,num_points,num_win,num_freq,padded_length,max_lines,num_removed,remain;
  float nwdt,*data,*mtap_spec_init,*mtap_spec_final,freq,tstart,srate,*initial_data,amp,phi;
  struct removed_lines *line_list;
  FILE *fpifo,*fplock,*fpout1,*fpout2;

  /* open the IFO output file and lock file */
  fpifo =grasp_open("GRASP_DATAPATH","channel.0","r");
  fplock=grasp_open("GRASP_DATAPATH","channel.10","r");

  /* data length, padded length, num frequencies including DC, Nyquist */
  num_points=8192;
  padded_length=65536;
  num_freq=1+padded_length/2;

  /* number of taper windows to use, and time-freq bandwidth */
  num_win=5;
  nwdt=3.0;

  /* maximum number of lines to remove */
  max_lines=100;
	
  /* allocate arrays */
  datas=(short *)malloc(sizeof(short)*num_points);
  data=(float *)malloc(sizeof(float)*num_points);
  mtap_spec_init=(float *)malloc(sizeof(float)*num_freq);
  mtap_spec_final=(float *)malloc(sizeof(float)*num_freq);
  line_list=(struct removed_lines *)malloc(sizeof(struct removed_lines)*max_lines);
  initial_data=(float *)malloc(sizeof(float)*num_points);

  /* get a section of data... */
  get_data(fpifo,fplock,&tstart,num_points,datas,&remain,&srate,0);

  /* copy short data to float data,and save initial data set  */
  for (i=0;i<num_points;i++) initial_data[i]=data[i]=datas[i];

  /* remove the spectral lines from the data set */
  remove_spectral_lines(data,num_points,padded_length,nwdt,num_win,
     max_lines,500,&num_removed,line_list,mtap_spec_init,mtap_spec_final,1,0,num_freq);
    
  /* print out a list of lines removed */
    printf("Total number of lines removed: %d\n",num_removed);
    for (i=0;i<num_removed;i++) {
    freq=0.5*line_list[i].index   *srate/num_freq;
    amp=2.0*sqrt(line_list[i].re*line_list[i].re+line_list[i].im*line_list[i].im);
    phi=180*atan2(line_list[i].im,line_list[i].re)/M_PI;
    printf("Removed line frequency %.3f Hz amplitude %.2f phase %.2f (F-test %.1f)\n",
		freq,amp,phi,line_list[i].fvalue);
  }

  /* now output a file containing the initial and final data... */
  fpout1=fopen("ifo_clean_data.out","w");
  fprintf(fpout1,"# Three columns are:\n# Time (sec)  Initial data  Final Data\n");
  for (i=0;i<num_points;i++)
    fprintf(fpout1,"%f\t%f\t%f\n",i/srate,initial_data[i],data[i]);
  fclose(fpout1);

  /* ... and the initial and final spectra, for graphing by xmgr */
  fpout2=fopen("ifo_clean_spec.out","w");
  fprintf(fpout2,"# Three columns are:\n# Freq (Hz)  Initial spectrum  Final spectrum\n");
  for (i=0;i<num_freq;i++)
    fprintf(fpout2,"%f\t%f\t%f\n",0.5*i*srate/num_freq,mtap_spec_init[i],mtap_spec_final[i]);
  fclose(fpout2);

  return 0;
}
