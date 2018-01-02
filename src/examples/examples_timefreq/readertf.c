/* PROGRAM TO READ THE ALLCURV.DAT FILE AND FIND THE LONGEST CURVE IN EACH SEGMENT */

#include <stdio.h>
#include <unistd.h>

int main(int argc, char **argv)
{
    char fil[256];
    int i,k,no_cur,len,mlen;
    float linstr,mlinstr,avlinstr;
    FILE *fp,*fpout;

    if(argc<2){
        printf("Usage: %s directory_name \n",argv[0]);
        exit(1);
    }
    no_cur = 0;
    sprintf(fil,"%s/allcurv.dat",argv[1]);
    fp = fopen(fil,"r");
    fpout = fopen("curve.dat","w");
    while(fscanf(fp,"%d %d\n",&i,&no_cur)==2){
        mlen=0;
        mlinstr = 0.0;
        for(k=0;k<no_cur;k++){
            fscanf(fp,"%d %f\n",&len,&linstr);
            if(len>mlen) {
                mlen=len;
                mlinstr =linstr;
            }
        }
        
        if(mlen!=0)
            avlinstr = mlinstr/mlen;
        else
            avlinstr = 0.0;
        fprintf(fpout,"%d %d %d %f %f\n",i,no_cur,mlen,mlinstr,avlinstr);
    }
    fclose(fpout);
    fclose(fp);
    exit(0);
}
