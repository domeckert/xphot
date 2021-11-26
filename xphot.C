#include <stdlib.h>                                                        
#include <stdio.h>                                                         
#include <cstring>
#include "math.h"
#include "fitsio.h"
#include "wcslib/wcs.h"
#include "sbprofs.h"

//const double BIG=1e50;
//const double pi=3.14159265358979312;

int line_num(char *filename){
	FILE *f=fopen(filename, "r");;
	char c;                       
	int lines = 0;                
	if(f == NULL) return 0;       
	while((c = fgetc(f)) != EOF){ 
		if(c == '\n') lines++;
	}                             
	fclose(f);                    
	return lines;                 
}                                     

void initwcs(struct wcsprm *m_wcs, double* crpix,
             double* crval, double* cdelt,char **ctype)
{
    // Call wcsini() in WCSLIB.
    int naxis = 2;
    m_wcs->flag = -1;
    wcsini(1, naxis, m_wcs);
    
    for( int i=0; i<naxis; i++){
        strcpy(m_wcs->ctype[i], ctype[i] ); //axis type
        m_wcs->crval[i] = crval[i]; // reference value
        m_wcs->crpix[i] = crpix[i]; // pixel coordinate
        m_wcs->cdelt[i] = cdelt[i]; // scale factor
    }
    int status=0;
    if ((status = wcsset(m_wcs))) {
        fprintf(stderr, "wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
    }
}

void region(int nsrc,double *xpos,double *ypos,double *radexc,double *exposure,long *axes){
    //printf("    %d regions will be ignored\n",nsrc);
    // Modified 01-24-2017 ; improved performance
    for (int ns=0;ns<nsrc;ns++){
        int boxsize=(int)round(radexc[ns]+0.5);
        int cx=(int)round(xpos[ns]);
        int cy=(int)round(ypos[ns]);
        for (int i=cx-boxsize; i<cx+boxsize+1; i++) {
            for (int j=cy-boxsize; j<cy+boxsize+1; j++) {
                if (i>=0 && j>=0 && i<axes[0] && j<axes[1]) {
                    double posx=(i-xpos[ns]);
                    double posy=(j-ypos[ns]);
                    double dist=sqrt(posx*posx+posy*posy);
                    if (dist<radexc[ns]){
                        exposure[j*axes[0]+i]=0.0;
                    }
                }
            }
        }
    }
}

int readimg(int band,char *fimg,char *expnam,double ra,double dec,double rad,int nsrc,double *rasrc,double *decsrc,double *radsrc,double &sb,double &esb,double &counts,double &expo,double &area,double &sb_bg,double &esb_bg,double &counts_bg,double &expo_bg,double &area_bg){
    // Reading files
    printf("Reading images for band %d...\n",band);
    int status=0;
    long axn[2];
    printf("Reading image file %s...\n",fimg);
    fitsfile *x;
    fits_open_file(&x,fimg,READONLY,&status);
    if (status!=0) {
        printf("    Error %d\n",status);
        return status;
    }
    int bitpix,naxis;
    fits_get_img_param(x,2,&bitpix,&naxis,axn,&status);
    long start[2]={1,1};
    int anynul;
    double *imgn=new double[axn[0]*axn[1]];
    fits_read_pix(x,TDOUBLE,start,axn[0]*axn[1],NULL,imgn,&anynul,&status);
    if (status!=0) {
        printf("    Error %d\n",status);
        return status;
    }
    double cdelt1,cdelt2,crval1,crval2,crpix1,crpix2;
    char **ctype=new char*[2];
    ctype[0]=new char[20];
    ctype[1]=new char[20];
    fits_read_key(x,TDOUBLE,(char *)"CDELT2",&cdelt2,NULL,&status);
    fits_read_key(x,TDOUBLE,(char *)"CDELT1",&cdelt1,NULL,&status);
    fits_read_key(x,TDOUBLE,(char *)"CRVAL1",&crval1,NULL,&status);
    fits_read_key(x,TDOUBLE,(char *)"CRVAL2",&crval2,NULL,&status);
    fits_read_key(x,TDOUBLE,(char *)"CRPIX1",&crpix1,NULL,&status);
    fits_read_key(x,TDOUBLE,(char *)"CRPIX2",&crpix2,NULL,&status);
    fits_read_key(x,TSTRING,(char *)"CTYPE1",ctype[0],NULL,&status);
    fits_read_key(x,TSTRING,(char *)"CTYPE2",ctype[1],NULL,&status);
    fits_close_file(x,&status);
    double crpix[2]={crpix1,crpix2};
    double cdelt[2]={cdelt1,cdelt2};
    double crval[2]={crval1,crval2};
    struct wcsprm *wcs_n=new struct wcsprm;
    initwcs(wcs_n,crpix,crval,cdelt,ctype);
    fitsfile *y;
    printf("Reading exposure map %s...\n",expnam);
    fits_open_file(&y,expnam,READONLY,&status);
    double *expn=new double[axn[0]*axn[1]];
    fits_read_pix(y,TDOUBLE,start,axn[0]*axn[1],NULL,expn,&anynul,&status);
    if (status!=0) {
        printf("    Error %d\n",status);
        return status;
    }
    fits_close_file(y,&status);
    printf("Images successfully read for band %d\n",band);
    //End reading data
    
    //Computing aperture photometry
    double world[2],imgcrd[2],pixcrd[2];
    world[0]=ra;
    world[1]=dec;
    double inra,indec;
    int stat;
    status = wcss2p(wcs_n, 1, 2, world, &inra,&indec,imgcrd, pixcrd,&stat);
    //printf("Pixel coordinates for the source: %g , %g\n",pixcrd[0],pixcrd[1]);
    
    //Extract source spectrum
    calc_xphot(imgn,expn,rad,axn,pixcrd[0],pixcrd[1],cdelt2,sb,esb,counts,area,expo);
    printf("Number of counts and count rate for source: %g counts ; CR = %g +/- %g cts/s\n",counts,sb,esb);
    
    //Mask sources
    double *xpos=new double[nsrc];
    double *ypos=new double[nsrc];
    double *radexc=new double[nsrc];
    for (int i=0; i<nsrc; i++) {
        world[0]=rasrc[i];
        world[1]=decsrc[i];
        status = wcss2p(wcs_n, 1, 2, world, &inra,&indec,imgcrd, pixcrd,&stat);
        xpos[i]=pixcrd[0];
        ypos[i]=pixcrd[1];
        radexc[i]=radsrc[i]/cdelt2/3600.; //From arcsec to pixel
		  //printf("i, xpos, ypos, radexc: %d %g %g %g\n",i,xpos[i],ypos[i],radexc[i]);
    }
    region(nsrc,xpos,ypos,radexc,expn,axn);
    
    //Extract background spectrum
    calc_xphot_bkg(imgn,expn,rad,axn,pixcrd[0],pixcrd[1],cdelt2,sb_bg,esb_bg,counts_bg,area_bg,expo_bg);
    printf("Number of counts and count rate for background: %g counts ; CR = %g +/- %g cts/s\n",counts_bg,sb_bg,esb_bg);
    
    printf("Band %d successfully processed\n",band);
    delete [] imgn;
    delete [] expn;
    delete [] xpos;
    delete [] ypos;
    delete [] radexc;
    return status;
}

int write_table(char *outfile,int nbin,double *counts,double *sb,double *esb,double *area,double *expo,char *backfile){
    int status=0;
    char outf[200];
    sprintf(outf,"!%s",outfile);
    fitsfile *x;
    fits_create_file(&x,outf,&status);
    if (status!=0) {
        printf("Error %d\n",status);
        return status;
    }
    int *channel=new int[nbin];
    for (int i=0; i<nbin; i++) {
        channel[i]=i;
    }
    /*char *ttype[2]={(char *)"CHANNEL",(char *)"COUNTS"};
    char *tform[2]={(char *)"1I",(char *)"1E"};
    char *tunit[2]={(char *)"channel",(char *)"count"};
    int tfields=2;*/
    char *ttype[4]={(char *)"CHANNEL",(char *)"COUNTS",(char *)"RATE",(char *)"STAT_ERR"};
    char *tform[4]={(char *)"1I",(char *)"1E",(char *)"1E",(char *)"1E"};
    char *tunit[4]={(char *)"channel",(char *)"count",(char *)"count/s",(char *)"count/s"};
    int tfields=4;
    char extname[]="SPECTRUM";
    fits_create_tbl(x,BINARY_TBL,nbin,tfields,(char **)ttype,(char **)tform,(char **)tunit,extname,&status);
    fits_write_col(x,TINT,1,1,1,nbin,channel,&status);
    fits_write_col(x,TDOUBLE,2,1,1,nbin,counts,&status);
    fits_write_col(x,TDOUBLE,3,1,1,nbin,sb,&status);
    fits_write_col(x,TDOUBLE,4,1,1,nbin,esb,&status);
    if (status!=0) {
        printf("    Error %d\n",status);
    }
    fits_write_key(x,TSTRING,(char *)"CREATOR",(void *)"xphot",NULL,&status);
    fits_write_key(x,TSTRING,(char *)"EXTNAME",(char *)"SPECTRUM",(char *)"Spectrum",&status);
    fits_write_key(x,TSTRING,(char *)"TELESCOP",(char *)"XMM",(char *)"Telescope",&status);
    fits_write_key(x,TSTRING,(char *)"INSTRUME",(char *)"EMOS1",(char *)"Instrument",&status);
    fits_write_key(x,TSTRING,(char *)"HDUCLASS",(char *)"OGIP",(char *)"Format conforms to OGIP/GSFC conventions",&status);
    fits_write_key(x,TSTRING,(char *)"HDUCLAS1",(char *)"SPECTRUM",(char *)"File contains a spectrum",&status);
    fits_write_key(x,TSTRING,(char *)"HDUCLAS2",(char *)"TOTAL",(char *)"File contains gross countst",&status);
    //fits_write_key(x,TSTRING,(char *)"HDUCLAS3",(char *)"COUNTS",(char *)"Spectrum is stored as counts",&status);
    int tl=1;
    fits_write_key(x,TLOGICAL,(char *)"POISSERR",&tl,(char *)"Poisson errors appropriate",&status);
    fits_write_key(x,TSTRING,(char *)"HDUVERS1",(char *)"1.1.0   ",(char *)"Version of format",&status);
    double syserr=0.0;
    fits_write_key(x,TDOUBLE,(char *)"SYS_ERR",&syserr,(char *)"Global systematic error",&status);
    fits_write_key(x,TSTRING,(char *)"CHANTYPE",(char *)"PI",(char *)"Type of channel data",&status);
    fits_write_key(x,TINT,(char *)"DETCHANS",&nbin,(char *)"Total number of detector channels available",&status);
    int quality=0;
    fits_write_key(x,TINT,(char *)"QUALITY",&quality,(char *)"All channels have good quality",&status);
    fits_write_key(x,TINT,(char *)"GROUPING",&quality,(char *)"No data grouping done",&status);
    if (nbin>1) {
        fits_write_key(x,TDOUBLE,(char *)"EXPOSURE",&expo[1],(char *)"Weighted exposure time for band 0.7-1.2",&status);
        fits_write_key(x,TDOUBLE,(char *)"BACKSCAL",&area[1],(char *)"Scaling factor for background",&status);
    }
    else{
        fits_write_key(x,TDOUBLE,(char *)"EXPOSURE",&expo[0],(char *)"Weighted exposure time for band 0.7-1.2",&status);
        fits_write_key(x,TDOUBLE,(char *)"BACKSCAL",&area[0],(char *)"Scaling factor for background",&status);
    }
    double scal=1.0;
    fits_write_key(x,TDOUBLE,(char *)"AREASCAL",&scal,(char *)"Nominal scaling factor for data",&status);
    fits_write_key(x,TDOUBLE,(char *)"CORRSCAL",&scal,(char *)"Nominal scaling factor for correction file",&status);
    fits_write_key(x,TSTRING,(char *)"BACKFILE",backfile,(char *)"Background FITS file",&status);
    fits_write_key(x,TSTRING,(char *)"RESPFILE",(char *)"m1_rebin.rmf",(char *)"redistribution matrix",&status);
    fits_write_key(x,TSTRING,(char *)"ANCRFILE",(char *)"m1_unvig.arf",(char *)"ancillary response",&status);
    char expband[20];
    char expcomm[100];
    for (int i=0;i<nbin;i++){
        sprintf(expband,"EXPO%d",i+1);
        sprintf(expcomm,"Weighted exposure time for band %d",i+1);
        fits_write_key(x,TDOUBLE,expband,&expo[i],expcomm,&status);
    }
	if (status!=0) {
        printf("    Error %d\n",status);
    }
    fits_close_file(x,&status);
    return status;
}

int get_num_sources(char *srclist){
    fitsfile *fsrc;
	 int status=0;
    char name[200];
    sprintf(name,"%s[1]",srclist);
    fits_open_file(&fsrc,name,READONLY,&status);
    if (status!=0){
        printf("Error %d\n",status);
        return status;
    }
    long nsrctemp;
    fits_get_num_rows(fsrc,&nsrctemp,&status);
    if (status!=0) {
        printf("Error %d\n",status);
        return status;
    }
	 int nsrc=nsrctemp;
	 return nsrc;
}

int read_srclist(char *srclist,int nsrc,double *raexc,double *decexc,double *radecerr,int status){
    fitsfile *fsrc;
    char name[200];
    sprintf(name,"%s[1]",srclist);
    fits_open_file(&fsrc,name,READONLY,&status);
    if (status!=0){
        printf("Error %d\n",status);
        return status;
    }
    int colnum;
    int anynul;
    fits_get_colnum(fsrc,CASEINSEN,(char *)"RA",&colnum,&status);
    if (status!=0) {
        printf("Error %d\n",status);
        return status;
    }
    fits_read_col(fsrc,TDOUBLE,colnum,1,1,nsrc,NULL,raexc,&anynul,&status);
    if (status!=0) {
        printf("Error %d\n",status);
        return status;
    }
    fits_get_colnum(fsrc,CASEINSEN,(char *)"DEC",&colnum,&status);
    if (status!=0) {
        printf("Error %d\n",status);
        return status;
    }
    fits_read_col(fsrc,TDOUBLE,colnum,1,1,nsrc,NULL,decexc,&anynul,&status);
    if (status!=0) {
        printf("Error %d\n",status);
        return status;
    }
    fits_get_colnum(fsrc,CASEINSEN,(char *)"RADEC_ERR",&colnum,&status);
    if (status!=0) {
        printf("Error %d\n",status);
        return status;
    }
    fits_read_col(fsrc,TDOUBLE,colnum,1,1,nsrc,NULL,radecerr,&anynul,&status);
    if (status!=0) {
        printf("Error %d\n",status);
        return status;
    }
    fits_close_file(fsrc,&status);
    if (status!=0) {
        printf("Error %d\n",status);
        return status;
    }
    else {
        printf("Source list successfully read\n");
    }
    return status;
}

void select_sources(double ra,double dec,double rad,int nexc,double *raexc,double *decexc,double *radecerr,int &nsrc,double *rasrc,double *decsrc,double *radsrc){
    // Routine to select only sources close to the source (speeds up masking process)
    nsrc=0;
    for (int i=0; i<nexc; i++) {
        double dist=sqrt((ra-raexc[i])*(ra-raexc[i])+(dec-decexc[i])*(dec-decexc[i]));
        if (dist*3600.<3.*rad) {
            rasrc[nsrc]=raexc[i];
            decsrc[nsrc]=decexc[i];
            radsrc[nsrc]=radecerr[i];
            nsrc++;
        }
    }
}

int main(int argc, char **argv){
	do {
		int status=0;
		if (argc!=9){
			printf("Usage: \n");
			printf("xphot imglist explist ra dec radius srcfile outspec outbkg\n");
            printf("imglist,explist: ASCII files containing paths to count images and exposure maps (one per line)\n");
            printf("ra, dec: source coordinates\n");
            printf("radius: size of extraction region (in arcsec); an annulus with the same size around the source is used for background estimation\n");
            printf("srcfile: file containing sources to exclude when calculating background\n");
            printf("outspec: output source spectrum file\n");
            printf("outbkg: output background spectrum file\n");
			break;
		}
		*argv++;
		char *flist=*argv;
		*argv++;
		char *explist=*argv;
		*argv++;
		double ra=atof(*argv);
		*argv++;
        double dec=atof(*argv);
        *argv++;
        double rad=atof(*argv);
        *argv++;
        char *srcfile=new char[200];
        sprintf(srcfile,"%s",*argv);
        *argv++;
		char *outfile=new char[200];
		sprintf(outfile,"%s",*argv);
        *argv++;
        char *outbkg=new char[200];
        sprintf(outbkg,"%s",*argv);

        int nband=line_num(flist);
        printf("%d files found in input list\n",nband);
		FILE *ff=fopen(flist,"r");
		FILE *fexp=fopen(explist,"r");
		char **allfiles=new char*[nband];
		char **allexp=new char*[nband];
		for (int i=0;i<nband;i++){
			allfiles[i]=new char[200];
			fscanf(ff,"%s\n",allfiles[i]);
			printf("File %d: %s\n",i+1,allfiles[i]);
			allexp[i]=new char[200];
			fscanf(fexp,"%s\n",allexp[i]);
			printf("Exposure map %d: %s\n",i+1,allexp[i]);
		}
		fclose(ff);
		fclose(fexp);
		int nexc=get_num_sources(srcfile);
        double *raexc=new double[nexc];
        double *decexc=new double[nexc];
        double *radecerr=new double[nexc];
        status=read_srclist(srcfile,nexc,raexc,decexc,radecerr,status);
        if (status!=0) {
            printf("Exiting with status %d\n",status);
            break;
        }
        double *rasrc=new double[nexc];
        double *decsrc=new double[nexc];
        double *radsrc=new double[nexc];
        int nsrc;
        select_sources(ra,dec,rad,nexc,raexc,decexc,radecerr,nsrc,rasrc,decsrc,radsrc);
        
        //Reading images and exposure maps
        double *sb=new double[nband];
        double *esb=new double[nband];
        double *counts=new double[nband];
        double *expo=new double[nband];
        double *area=new double[nband];
        
        double *sb_bg=new double[nband];
        double *esb_bg=new double[nband];
        double *counts_bg=new double[nband];
        double *expo_bg=new double[nband];
        double *area_bg=new double[nband];
        
		for (int i=0;i<nband;i++){
			status=readimg(i+1,allfiles[i],allexp[i],ra,dec,rad,nsrc,rasrc,decsrc,radsrc,sb[i],esb[i],counts[i],expo[i],area[i],sb_bg[i],esb_bg[i],counts_bg[i],expo_bg[i],area_bg[i]);
            if (status!=0) {
                printf("Exiting with status %d\n",status);
            }
		}
        
        char dummybkg[20];
        sprintf(dummybkg,"NONE");
        write_table(outfile,nband,counts,sb,esb,area,expo,outbkg);
        write_table(outbkg,nband,counts_bg,sb_bg,esb_bg,area_bg,expo_bg,dummybkg);
        if (status==0) {
            printf("Results written in file %s (source) ; %s (background)\n",outfile,outbkg);
        }
        else{
            printf("Files could not be written");
        }
        delete [] sb;
        delete [] esb;
        delete [] counts;
        delete [] expo;
        delete [] area;
        delete [] sb_bg;
        delete [] esb_bg;
        delete [] counts_bg;
        delete [] expo_bg;
        delete [] area_bg;
        delete [] raexc;
        delete [] decexc;
        delete [] radecerr;
        delete [] rasrc;
        delete [] decsrc;
        delete [] radsrc;
	}while (0);
}
