/*
 *  sbprofs.h
 *  
 *
 *  Created by Dominique Eckert on 23.09.10.
 *  Copyright 2010 INAF/IASF-Milano. All rights reserved.
 *
 */

void calc_xphot(double *img,double *exposure,double srcrad,long *axes,double centroid_ra,double centroid_dec,double pixsize,double &sb,double &esb,double &counts,double &area,double &expo){
    sb=0.0;
    esb=0.0;
    counts=0.0;
    area=0.0;
    expo=0.0;
    
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
                }
            }
        }
    }
    esb=sqrt(esb);
    expo/=npix;
    area=npix*pixsize*3600.*pixsize*3600./0.05/0.05;//backscal parameter
}
