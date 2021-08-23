# function to create array containing all wavelength values 
# of the grid where a spectrum has been sampled to
# from a fits header object

from astropy.io import fits
import numpy as np

def wavel(header,naxis=3,cdel_key='CD3_3'):
   """
   xax = wavel(header,naxis=3,cdel_key='CD3_3')
   header ... a pyfits header object containing naxis, crval, crpix & crdel values
              needed for creation of wavelength grid array
   naxis ... the wavelength axis (used to determine CRVALn, NAXISn & CRPIXn values)
   cdel_key ... the key word that contains the wavelength increment array
   """
   nax = naxis

   naxis = header['NAXIS'+str(nax)]
   crval = header['CRVAL'+str(nax)]
   crpix = header['CRPIX'+str(nax)]
   crdel = header[cdel_key]

   xax = crval + (np.arange(naxis) - (crpix - 1))*crdel

   return xax


def wavel_to_index(xax,lambda_in):
   """
   index_out = wavel_to_index(xax, lambda_in)
   """
   return np.abs(xax - lambda_in).argmin()
   

def wavel_old(header,dim=3):
   """
   xax = wavel_old(header,dim=3)
   header ... a pyfits header object containing naxis, crval, crpix & crdel values
              needed for creation of wavelength grid array
   dim ... dimension which correpsonds to wavelength (default: 3)
   """
   assert type(dim) is int
   assert type(header) is fits.header.Header

   naxis = header['NAXIS'+str(dim)]
   crval = header['CRVAL'+str(dim)]
   crpix = header['CRPIX'+str(dim)]
   crdel = header['CDELT'+str(dim)]

   xax = crval + (np.arange(naxis) - (crpix - 1))*crdel

   return xax
