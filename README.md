# xphot

A C++ tool which reads a set of X-ray images and exposure maps in narrow bands and extracts XSPEC-compatible source and background spectra using aperture photometry

Usage:
    xphot imglist explist ra dec radius srcfile outspec outbkg
    
- imglist,explist: ASCII files containing paths to count images and exposure maps (one per line)
- ra, dec: source coordinates
- radius: size of extraction region (in arcsec); an annulus with the same size around the source is used for background estimation
- srcfile: file containing sources to exclude when calculating background
- outspec: output source spectrum file
- outbkg: output background spectrum file

