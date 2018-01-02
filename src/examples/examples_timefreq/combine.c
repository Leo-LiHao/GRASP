/* PROGRAM TO COMBINE THE OUTPUT FILES FROM THE TF PROGRAM*/

#include <stdio.h>
#include <unistd.h>


int main(int argc, char **argv)
{
    char fil[256];
    int i,j,k,no_cur,len;
    float linestr;
    FILE *fp,*fpout;
    int NOFSEG,NOFSUBSEG;

        /* PROGRAM USAGE */
    if(argc<4){
        printf("Usage: %s directory_name no_of_segments no_of_subsegments\n",argv[0]);
        exit(1);
    }
        /* GET THE NUMBER OF SEGMENTS AND SUBSEGMENTS */
    NOFSEG = atoi(argv[2]);
    NOFSUBSEG = atoi(argv[3]);
        /* OPEN OUTPUT FILE */
    sprintf(fil,"%s/allcurv.dat",argv[1]);
    fpout = fopen(fil,"w");
    for(i=0;i<NOFSEG;i++){
        for(j=0;j<NOFSUBSEG;j++){
            fprintf(fpout,"%d ",i*NOFSUBSEG + j);
                /* OPEN EACH FILE PRODUCED BY TF PROGRAM */
            sprintf(fil,"%s/out_%d.%02d",argv[1],i,j);
            if((fp = fopen(fil,"r")) != NULL){
                fscanf(fp,"%d\n",&no_cur);
                fprintf(fpout,"%d ",no_cur);
                for(k=0;k<no_cur;k++){
                    fscanf(fp,"%d %f",&len,&linestr);
                    fprintf(fpout,"%d %f ",len,linestr);
                    
                }
                fclose(fp);
                sprintf(fil,"/bin/rm %s/out_%d.%02d",argv[1],i,j);
                system(fil);
                fprintf(fpout,"\n");
            }
            else{
                fprintf(fpout,"%d \n",0);
            }
            
        }
    }
    fclose(fpout);
    exit(0);
}
