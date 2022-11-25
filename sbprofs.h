/*
 *  sbprofs.h
 *  
 *
 *  Created by Dominique Eckert on 23.09.10.
 *  Copyright 2010 INAF/IASF-Milano. All rights reserved.
 *
 */

double max(double *vec, int nel){
	double tm=0.0;
	for (int i=0;i<nel;i++){
		if (vec[i]>tm) {
			tm=vec[i];		
		}
	}
	return tm;
}

void calc_xphot(double *img,double *exposure,double srcrad,long *axes,double centroid_ra,double centroid_dec,double pixsize,double &sb,double &esb,double &counts,double &area,double &expo){
    sb=0.0;
    esb=0.0;
    counts=0.0;
    area=0.0;
    expo=0.0;
    
	//double maxexp=max(exposure, axes[0]*axes[1]);

    int npix=0;
    //Size of the box
    int nbox=(int)floor(srcrad/(pixsize*60.*60.)+1.);
    int cx=(int)round(centroid_ra-1.);
    int cy=(int)round(centroid_dec-1.);
    
	for (int i=cx-nbox;i<cx+nbox+1;i++){
		for (int j=cy-nbox;j<cy+nbox+1;j++){
            double posx=(i-centroid_ra)*pixsize*3600;//arcsec
            double posy=(j-centroid_dec)*pixsize*3600;
            double dist=sqrt(posx*posx+posy*posy);
            if ((i>=0)&&(j>=0)&&(i<axes[0])&&(j<axes[1])){
                if ((dist<srcrad)&&(exposure[j*axes[0]+i]>0.0)){
                    sb+=img[j*axes[0]+i]/exposure[j*axes[0]+i];
                    esb+=img[j*axes[0]+i]/exposure[j*axes[0]+i]/exposure[j*axes[0]+i];
                    expo+=exposure[j*axes[0]+i];
                    counts+=img[j*axes[0]+i];
                    npix++;
					//printf("i, j, dist, counts, expo, npix: %d %d %g %g %g %d\n", i, j, dist, img[j*axes[0]+i], exposure[j*axes[0]+i], npix);
                }
            }
		}
	}
    esb=sqrt(esb);
    expo/=npix;
    area=npix*pixsize*3600.*pixsize*3600./0.05/0.05;//backscal parameter
}

void calc_xphot_bkg(double *img,double *exposure,double srcrad,long *axes,double centroid_ra,double centroid_dec,double pixsize,double &sb,double &esb,double &counts,double &area,double &expo){
    sb=0.0;
    esb=0.0;
    counts=0.0;
    area=0.0;
    expo=0.0;
    
    int npix=0;
    //Size of the box
    int nbox=(int)floor(4.*srcrad/(pixsize*3600.)+1.);
    int cx=(int)round(centroid_ra-1.);
    int cy=(int)round(centroid_dec-1.);
    
    for (int i=cx-nbox;i<cx+nbox+1;i++){
        for (int j=cy-nbox;j<cy+nbox+1;j++){
            double posx=(i-centroid_ra)*pixsize*3600;//arcsec
            double posy=(j-centroid_dec)*pixsize*3600;
            double dist=sqrt(posx*posx+posy*posy);
            if ((i>=0)&&(j>=0)&&(i<axes[0])&&(j<axes[1])){
                if ((dist>=2.*srcrad)&&(dist<3.*srcrad)&&(exposure[j*axes[0]+i]>0.0)){ //Select pixels between 2 and 3 times srcrad
					sb+=img[j*axes[0]+i]/exposure[j*axes[0]+i];
                    esb+=img[j*axes[0]+i]/exposure[j*axes[0]+i]/exposure[j*axes[0]+i];
                    expo+=exposure[j*axes[0]+i];
                    counts+=img[j*axes[0]+i];
                    npix++;
					//printf("npix, counts: %d , %g\n", npix, counts);       
                }
            }
        }
    }
    esb=sqrt(esb);
    expo/=npix;
    area=npix*pixsize*3600.*pixsize*3600./0.05/0.05;//backscal parameter
}

int save_img(char *temp,double *modimg,long *axes,double crpix1,double crval1,double crpix2,double crval2,double cdelt1,double pixsize){
	int status=0;
	fitsfile *x;
	char nfn[200];
	sprintf(nfn,"!%s",temp);
	fits_create_file(&x,nfn,&status);
	if (status!=0){
		printf("    Error %d\n",status);
		return status;
	}
	fits_create_img(x,-32,2,axes,&status);
	long start[2]={1,1};
	fits_write_pix(x,TDOUBLE,start,axes[0]*axes[1],modimg,&status);
	if (status!=0){
		printf("    Error %d\n",status);
		return status;
	}
	fits_write_key(x,TSTRING,(char *)"CREATOR",(void *)"proffit 1.5",NULL,&status);
	fits_write_key(x,TSTRING,(char *)"CTYPE1",(void *)"RA---TAN",(char *)"LONGPROJ where LONG can be RA, GLON, ELON and PROJ can be CAR, TAN or AIT",&status);
	fits_write_key(x,TDOUBLE,(char *)"CRPIX1",&crpix1,(char *)"Pixel at reference point",&status);
	fits_write_key(x,TDOUBLE,(char *)"CRVAL1",&crval1,(char *)"LONG at the reference value",&status);
	fits_write_key(x,TSTRING,(char *)"CUNIT1",(void *)"deg",(char *)"Physical units of axis 1",&status);
	fits_write_key(x,TSTRING,(char *)"CTYPE2",(void *)"DEC--TAN",(char *)"LAT-PROJ where LAT can be DEC, GLAT, ELAT and PROJ can be CAR, TAN or AIT",&status);
	fits_write_key(x,TDOUBLE,(char *)"CRPIX2",&crpix2,(char *)"Pixel at reference point",&status);
	fits_write_key(x,TDOUBLE,(char *)"CRVAL2",&crval2,(char *)"LAT at the reference value",&status);
	fits_write_key(x,TSTRING,(char *)"CUNIT2",(void *)"deg",(char *)"Physical units of axis 2",&status);
	fits_write_key(x,TDOUBLE,(char *)"CDELT1",&cdelt1,(char *)"Element (1,1) of coordinate transf. matrix (default 1)",&status);
	fits_write_key(x,TDOUBLE,(char *)"CDELT2",&pixsize,(char *)"Element (2,2) of coordinate transf. matrix (default 1)",&status);
	fits_write_key(x,TSTRING,(char *)"RADECSYS",(void *)"FK5",(char *)"Stellar reference frame",&status);
	fits_write_key(x,TSTRING,(char *)"EQUINOX",(void *)"2000.0",(char *)"Coordinate system equinox",&status);
	fits_close_file(x,&status);
	if (status!=0){
		printf("    Error %d\n",status);
		return status;
	}
	else {
		printf("    Image succesfully written\n");
	}
	return status;
}	

