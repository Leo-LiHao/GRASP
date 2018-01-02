/* GRASP: Copyright 1997,1998  Bruce Allen */
#include "grasp.h"

int main(int argc, char *argv[]) {   
   FILE *fp;
   struct fgetinput fgetinput;
   struct fgetoutput fgetoutput;
   float  *data,delta_f;
   int *npoints;
   int i,j,check,min_points;
   char fname[256],detector[256],fft_dir[256],cmd[256];
     
   if ( argc != 2) { 
      printf("Usage: corr_init configuration-file\n");
      exit(1);
   } 

   fp = fopen(argv[1],"r"); 
   if ( fp == NULL ) {
    fprintf(stderr,"Problems opening %s\n",argv[1]);
    exit(1);    
   }   
   fprintf(stderr,"Reading %s\n",argv[1]); 

   while (1) {
     fgets(detector,sizeof(detector),fp);
     if (detector[0] != '#') break;
   }
   detector[strlen(detector)-1]='\0';
   check=fscanf(fp,"%d",&fgetinput.nchan);
    
   /* storage for channel names, data locations, points returned, ratios */
   fgetinput.chnames=(char **)malloc(fgetinput.nchan*sizeof(char *));
   for (i=0;i<fgetinput.nchan;i++)   
     fgetinput.chnames[i]=(char *)malloc(256*sizeof(char));
   fgetinput.locations=(short **)malloc(fgetinput.nchan*sizeof(short *));
   fgetoutput.npoint=(int *)malloc(fgetinput.nchan*sizeof(int));
   fgetoutput.ratios=(int *)malloc(fgetinput.nchan*sizeof(int));
   fgetinput.datatype=(char *)malloc(fgetinput.nchan*sizeof(char));
   npoints=(int *)malloc(fgetinput.nchan*sizeof(int));
 
   for (i=0;i<fgetinput.nchan;i++) {
     check=fscanf(fp,"%s %c %d",fgetinput.chnames[i],&fgetinput.datatype[i],&npoints[i]);
     /* the next fgetinput.nchan lines of the configuration file should contain 3 columns */
     /*                 - if not print an error message                                   */
     if (check != 3) {
       fprintf(stderr,"Problems reading data from %s\n",argv[1]); 
       exit(1);   
     }
   }

   fclose(fp);

   /* number of points to get */
    fgetinput.npoint=npoints[0];

   /* allocate storage space for data */
   for(i=0;i<fgetinput.nchan;i++) {
     switch (fgetinput.datatype[i]) {
     case 'S': /* short data */
     case 'u': /* unsigned short data */
       fgetinput.locations[i]=(short *)malloc(npoints[i]*sizeof(short)); 
       break;
     case 'I':  /* integer data */
     case 'i':  /* unsigned integer data */
       fgetinput.locations[i]=(short *)malloc(npoints[i]*sizeof(int));
       break;
     case 'L': /* long data */
     case 'l': /* unsigned long data */
       fgetinput.locations[i]=(short *)malloc(npoints[i]*sizeof(long));
       break;
     case 'F': /* float data */
       fgetinput.locations[i]=(short *)malloc(npoints[i]*sizeof(float));
       break;
     case 'D': /* double data */
       fgetinput.locations[i]=(short *)malloc(npoints[i]*sizeof(double));
       break;
     case 'C': /* character data */
     case 'f': /* complex float data */ 
     case 'd': /* complex double data */
     case 's': /* character string data */
     case 'c': /* unsigned character data */
       fprintf(stderr,"Data type %c cannot be plotted\n",fgetinput.datatype[0]);
       break;
     default:
       fprintf(stderr,"Unknown data type %c\n",fgetinput.datatype[0]);
       exit(1);
     }      
   }  

   /* don't have inlock channel if not 40m */
   if (detector != "C1") { 
     fgetinput.inlock=0;
   }

   /* but we don't need calibration information */
   fgetinput.calibrate=0;

   /* don't seek, we need the sample values! */
   fgetinput.seek=0;

   /* source of files */
   fgetinput.files=framefiles;

   check=fget_ch(&fgetoutput,&fgetinput); 
   if ( check == 0) {
     fprintf(stderr,"not enough data!!!!!!!!");
     exit(1);
   }
   fprintf(stderr,"Signal (%s) sample rate is %f\n",fgetinput.chnames[0],fgetoutput.srate);

   min_points=npoints[0];
   for(i=1;i<fgetinput.nchan;i++) {
     if (npoints[i]<min_points) min_points=npoints[i];
   }

   delta_f=fgetoutput.srate/npoints[0];

   for (i=0;i<256;i++) {
     if (argv[1][i] == '\0' || argv[1][i] == '.') {
       fft_dir[i] = '\0';
       break;
     }
     fft_dir[i]=argv[1][i];
   }

   strcat(fft_dir,"_fft");
   sprintf(cmd,"mkdir %s 2> /dev/null",fft_dir);
   system(cmd);

   for (i=0;i<fgetinput.nchan;i++) {
     /* check frames are consistent with configuration file */
     if (npoints[i] != fgetoutput.srate/fgetoutput.ratios[i]) {
       fprintf(stderr,"Sample rates in %s is %f\n",
	       fgetinput.chnames[i],fgetoutput.srate/fgetoutput.ratios[i]);
     }

     /* assign memory for storing data in floats */
     data=(float *)malloc(sizeof(float)*npoints[i]);
 
     switch (fgetinput.datatype[i]) {
     case 'S': /* short data */
     case 'u': /* unsigned short data */
      for(j=0;j<npoints[i];j++)  data[j]=fgetinput.locations[i][j];  
      break;
     case 'I':  /* integer data */
     case 'i':  /* unsigned integer data */
        for(j=0;j<npoints[i];j++)  data[j]=((int *)(fgetinput.locations[i]))[j];
       break;
     case 'L': /* long data */
     case 'l': /* unsigned long data */
       for(j=0;j<npoints[i];j++)  data[j]=((long *)(fgetinput.locations[i]))[j];
       break;
     case 'F': /* float data */
       for(j=0;j<npoints[i];j++)  data[j]=((float *)(fgetinput.locations[i]))[j];
       break;
     case 'D': /* double data */
       for(j=0;j<npoints[i];j++)  data[j]=((double *)(fgetinput.locations[i]))[j];
       break;
     case 'C': /* character data */
     case 'f': /* complex float data */ 
     case 'd': /* complex double data */
     case 's': /* character string data */
     case 'c': /* unsigned character data */
       fprintf(stderr,"Data type %c cannot be converted to float\n",fgetinput.datatype[0]);
       break;
     default:
       fprintf(stderr,"Unknown data type %c\n",fgetinput.datatype[0]);
       exit(1);
     }          
 
     /* We can now free fgetinput.locations[i] */
     free(fgetinput.locations[i]); 
     /* Take the Fourier transforms */
     realft(data-1,npoints[i],1);
     data[1]=0;
 
     /* Print the Fourier transform to file */
     sprintf(fname,"%s/%s-%s_fft.b",fft_dir,detector,fgetinput.chnames[i]);
     fp = fopen(fname,"wb"); 
     if ( fp == NULL) {
       printf("Cannot open %s\n",fname);
       exit(1);
     } 
     fprintf(stderr,"writing %s\n",fname);
     fwrite(&delta_f,sizeof(float),1,fp);
     fwrite(data,sizeof(float),min_points,fp);	    
     fclose(fp);

     
     free(data);
   }
     
   /* free memory */
   free(fgetinput.locations);
   free(fgetinput.chnames);
   free(fgetinput.locations);
   free(fgetoutput.npoint);
   free(fgetoutput.ratios);
   free(npoints);

   return 0;
}
