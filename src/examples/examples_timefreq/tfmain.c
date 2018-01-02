      /* To investigate the Time Frequency  distribution of chirps */
      /* using MPI and GL calls */
      /* prototype program */

#include "mpi.h"
#include "grasp.h"
#include "tfmain.h"

/**************************************************
  GLOBAL VARIABLES TO BE USED ACROSS FILES
  ************************************************/

struct_tfparam tfparam; /* the parameters of the timefrequency program defined in file timefreq.h */
dl_options dlopt;       /* the parameters for the linerecognition algorithms*/
struct_signalparameters snpar; /* the parameters of the signal and noise; defined in tfmain.h */
MPI_Comm comm = 91;
long longn=DATADIM;
float *htilde,srate=SRATE;
int reply=1,count=0,counter;
int npoint=DATADIM,numprocs,myid,new_lock=0;

int main(int argc,char **argv)
{

    char processor_name[256],command[256];
    int namelen,mypid;
    
        /* initialize the MPI processes*/
    MPI_Init ( &argc, &argv);
    MPI_Comm_size ( comm, &numprocs );
    MPI_Comm_rank ( comm, &myid );
    MPI_Get_processor_name(processor_name,&namelen);

        /* Renice the processes */
    mypid = getpid();
    sprintf(command,"renice 10  %d &\n",mypid);
    system(command);
         /* For debugging purposes */
#if(DEBUG==1)
    if(myid==0){
        sprintf(command,"xxgdb tf %d &\n",mypid);
        system(command);
    }
    sleep(20);
#endif
        /*allocate space various common arrays */
    htilde = (float *) malloc(sizeof(float)*npoint);
       /* get the parameters from the file tfmain.in */
    gettfparameters();
        /* branch off to either the master or the slaves */
    if(myid==0)
        master();
    else
        slave();
        /* exit */
    MPI_Finalize();
    return 0;
}

void master()
{
    int slaves,recv=0,from,mdata=0,i;
    MPI_Status status;
    char tmp_str[256];
    float *scale;
    FILE *fp,*fp1,*fp2;
    
        /* create the output directory and copy the input file there */
    sprintf(tmp_str,"mkdir run%02d\n",tfparam.run_number);
    system(tmp_str);
    sprintf(tmp_str,"cp tfmain.in run%02d\n",tfparam.run_number);
    system(tmp_str);
    sprintf(tmp_str,"cp tfmain.h run%02d\n",tfparam.run_number);
    system(tmp_str);    
    sprintf(tmp_str,"run%02d/timstat",tfparam.run_number);
    fp = fopen(tmp_str,"w");
    sprintf(tmp_str,"run%02d/rescale",tfparam.run_number);
    fp1 = fopen(tmp_str,"w");
    if(snpar.signaltype==INSERT_COALESCENCE){
        sprintf(tmp_str,"cp MergeSig.dat run%02d\n",tfparam.run_number);
        system(tmp_str);
    }
    sprintf(tmp_str,"run%02d/segments",tfparam.run_number);
    fp2 = fopen(tmp_str,"w");
    sprintf(tmp_str,"printenv > run%02d/environment\n",tfparam.run_number);
    system(tmp_str);
    sprintf(tmp_str,"echo number of processes = %d >> run%02d/prog",numprocs,tfparam.run_number);
    system(tmp_str);
        /* set the rescale parameter */
    snpar.addsignal = 0;
    scale = (float *) malloc(sizeof(float)*numprocs);
    for (slaves=1;slaves<numprocs;slaves++){
        mdata = get_time_series_data();
        MPI_Send(&count, 1, MPI_INT, slaves , 1001, comm);
        MPI_Send(htilde,npoint,MPI_FLOAT,slaves,1002,comm);
    }
    for(i=0;i<tfparam.start_segment;i++){
        mdata = get_time_series_data();
        printf("Skipping segment %d, mdata = %d\n",i,mdata);
        fflush(stdout); 
    }
    for(slaves=1;slaves<numprocs;slaves++){
        MPI_Recv(scale+slaves, 1,MPI_FLOAT, slaves , 1004, comm, &status);
        fprintf(fp1,"%d %f\n",slaves,scale[slaves]);
    }
    fclose(fp1);
        /* Compute average of the rescale values */
    tfparam.maxpixelval = 0.0;
    for(slaves=1;slaves<numprocs;slaves++) tfparam.maxpixelval += scale[slaves];
    tfparam.maxpixelval /= (numprocs-1);
    printf(" The average rescale value : %f\n",tfparam.maxpixelval);
    fflush(stdout);
    free(scale);
        /* send the rescale value to all the slaves */
    MPI_Bcast(&tfparam.maxpixelval, 1, MPI_FLOAT, 0, comm);
    for(slaves=1;slaves<numprocs;slaves++){
        MPI_Recv(&reply, 1,MPI_INT, slaves , 1003, comm, &status);
        MPI_Recv(&counter , 1,MPI_INT, slaves , 1003, comm, &status);
    }

    count= tfparam.start_segment;
    recv = tfparam.start_segment;
    snpar.addsignal = 1;
        /* all set to start the actual simulation loop over
           slaves and send the datasegments  to the slaves */
    sprintf(tmp_str,"date >> run%02d/prog\n",tfparam.run_number);
    system(tmp_str);
    for (slaves=1;slaves<numprocs;slaves++){
            /* get the data */
        mdata = get_time_series_data();
#if(DEBUG1)
        graph(htilde,npoint,1);
        sleep(5);
#endif
        timstat(count, fp, htilde);
        MPI_Send(&count, 1, MPI_INT, slaves , 1001, comm);
        MPI_Send(htilde,npoint,MPI_FLOAT,slaves,1002,comm);
        printf("master   : sent segment %d to slave %d mdata = %d\n",count,slaves,mdata);
        fflush(stdout);
        count++;
    }
        /* wait for the slaves to send the done message and send fresh datasegments */
    while(recv<count){
        MPI_Recv(&reply , 1, MPI_INT, MPI_ANY_SOURCE, 1003, comm, &status);
        from = status.MPI_SOURCE;
        MPI_Recv(&counter , 1,MPI_INT, from , 1003, comm, &status);
        printf("master   : Received reply for segment %d back from slave %d\n",counter,from);
        fflush(stdout);
        recv++;
        if(mdata) mdata = get_time_series_data();
        printf("master   :count=%d,recv=%d,mdata=%d\n",count,recv,mdata);
        if((count<tfparam.num_of_segments)&&(mdata)) {
            fprintf(fp2,"%d %d %d \n",count, from,new_lock);
            MPI_Send(&count, 1, MPI_INT, from , 1001, comm);
            MPI_Send(htilde, npoint, MPI_FLOAT, from , 1002, comm);
            printf("master   : sent segment %d to slave %d\n",count,from);
            fflush(stdout);
            timstat(count, fp,htilde);
            count++;
        }
            /* send slave the termination message */
        else{
            counter=-1;
            printf("master   :Sending termination signal to slave  %d \n",from);
            MPI_Send(&counter, 1, MPI_INT, from , 1001, comm);
        }
    }
    fclose(fp);
    fclose(fp2);
    sprintf(tmp_str,"date >> run%02d/prog\n",tfparam.run_number);
    system(tmp_str);
    return;
}


void slave()
{
    MPI_Status status;
    float *tfgeneral[TDIM],*pic[PDIM],*out_img[PDIM],scalefactor=0.0;
    int i,rows=PDIM,cols=PDIM,noofsub,ind;
    char filen[256];
    static int first = 1;
    
        /* the number of subsegments in each segment of data */
    noofsub = (DATADIM-(POSTSAFETY+PRESAFETY))/NDIM;
        /* allocate memory of the timefrequency block and for the picture matrices*/
    for(i=0;i<TDIM/tfparam.offset_step_size;i++)
        tfgeneral[i] = (float *)malloc(sizeof(float)*NDIM);
    for(i=0;i<PDIM;i++){
        pic[i] = (float *)malloc(sizeof(float)*PDIM);
        out_img[i] = (float *)malloc(sizeof(float)*PDIM);
    }
        /* loop for receiving the data */
    while(1){
        MPI_Recv(&count,1,MPI_INT,0,1001,comm,&status);
            /* Terminate if count is less that 0*/
        if(count<0) return;
            /* Recv the htilde array from the master*/
        MPI_Recv(htilde, npoint, MPI_FLOAT, 0, 1002, comm, &status);
        printf("slave %d: received segment %d from master\n",myid,count);
        fflush(stdout);
            /* compute the time-frequency maps for the data segment */
        for(ind=0;ind<noofsub;ind++){
            /* compute the time-frequency maps for each subsegment */
            time_freq_map(htilde, &tfparam, ind,tfgeneral,pic);            
            if(first){
                scalefactor += compute_scalefactor(pic,tfparam.rescale_factor,PDIM);
#if(DEBUG1)                
                printf("slave %d: nofsub = %d  scalefactor = %f\n",myid,ind,compute_scalefactor(pic,tfparam.rescale_factor,PDIM));
                fflush(stdout);
#endif                
                if(ind==(noofsub-1)){
                    first = 0;
                    scalefactor /= noofsub;
                    MPI_Send(&scalefactor,1,MPI_FLOAT, 0, 1004, comm);
                    MPI_Bcast(&tfparam.maxpixelval, 1, MPI_FLOAT, 0, comm);
                }
            }
            else{
                rescale(pic,PDIM,tfparam.maxpixelval);
                sprintf(filen,"./run%02d/out_%d.%02d",tfparam.run_number,count,ind);
                get_line_lens(dlopt.sigma,dlopt.high,dlopt.low,rows,cols,pic,filen);
#if(DEBUG2)
                if((ind==0)&&(count==0))
                    ppmprint(pic,"picture.ppm",PDIM);
#ifdef HAVE_GL
                    /* if you want to display the TF map on the screen */
                if(myid==1)
                    plottf(pic,PDIM);
#endif
                printf("slave %d: nofsub = %d\n",myid,ind);
                fflush(stdout);
#endif
            }    
        }
            /* inform master that you have finished the current segment */
        MPI_Send(&reply, 1, MPI_INT, 0, 1003, comm);
            /* send the counter back to the master as a check */
        MPI_Send(&count, 1, MPI_INT, 0, 1003, comm);
    }
}
