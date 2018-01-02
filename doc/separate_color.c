#include <stdio.h>
#define LEN 1024

int pairs[2*LEN];
int lastpage=-1;
int retval;


int makepairs(int *pairs,int *pages) {
int i=0,j=0;

	while (i<lastpage+1) {
		if (pages[i]==0 && i<lastpage+1) {
			i++;
		}
		else {
			pairs[2*j]=i;
			while (pages[i]==1 && i<lastpage+1)
				i++;
			pairs[2*j+1]=i-1;
			j++;
		}
	}
	return j;
}

main() {
FILE *fp;
int page=0;
int i,start,stop,colorlist[LEN],bwlist[LEN];
int colorpages=0,num;

	for (i=0;i<LEN;i++) {
		colorlist[i]=0;
		bwlist[i]=1;
	}
	fp=fopen("manual.idx","r");
	while (1) {
		retval=fscanf(fp,"\\indexentry{colorpage|hyperpage}{%d}\n",&page);
		if (retval==EOF) break;
		if (retval==1) {
			colorlist[page]=1;
			bwlist[page]=0;
			colorpages++;
		}
		if (retval==0)
			retval=fscanf(fp,"alastpage|hyperpage}{%d}\n",&lastpage);
		if (retval==EOF) break;
	}

	if (lastpage==-1) {
		fprintf(stderr,"No final page found in index!\n");
		exit(1);
	}
	if (lastpage>=LEN) {
		fprintf(stderr,"MACRO LEN = %d too small!  Recompile!\n",LEN);
		exit(1);
	}

	printf("#!\n\n# First make a postscript file of everything.\n");
	printf("dvips -z -o junk.ps manual\n");
	printf("# Filter to put %%%%BeginDocument lines at start of new lines\n");
	printf("sed -e 's/@setspecial%%%%BeginDocument/@setspecial\\\\\n%%%%BeginDocument/' junk.ps >! manual.ps\n");
	printf("rm -f junk.ps\n\n");


	printf("#\n# Now a postscript file of the color pages only.\n");
	printf("dvips -pp ");
	num=makepairs(pairs,colorlist);
	for (i=0;i<num;i++) {
		if (pairs[2*i]==pairs[2*i+1])
			printf("%d",pairs[2*i]);
		else
			printf("%d-%d",pairs[2*i],pairs[2*i+1]);
		if (i!=num-1)
			printf(",");
	}
	printf(" -o manual_color.ps manual\n\n");


	printf("#\n# Now a postscript file of the black and white pages only.\n");
	printf("dvips -pp ");
	num=makepairs(pairs,bwlist);
	for (i=0;i<num;i++) {
		if (pairs[2*i]==pairs[2*i+1])
			printf("%d",pairs[2*i]);
		else
			printf("%d-%d",pairs[2*i],pairs[2*i+1]);
		if (i!=num-1)
			printf(",");
	}
	printf(" -o manual_bw.ps manual\n\n");


	return 0;
}
		

