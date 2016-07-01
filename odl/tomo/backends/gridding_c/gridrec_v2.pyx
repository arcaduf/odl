import cython

import numpy as np
cimport numpy as np


cdef extern void gridrec_v2_backproj( float *S , int npix , int nang , float *angles ,
                                   float *param , float *I )

cdef extern void gridrec_v2_forwproj( float *S , int npix , int nang , float *angles ,
                                   float *param , float *I )


                    

@cython.boundscheck( False )
@cython.wraparound( False )
def backproj( np.ndarray[ float , ndim=2 , mode="c" ] sinogram not None ,
              np.ndarray[ float , ndim=1 , mode="c" ] angles not None ,
              np.ndarray[ float , ndim=1 , mode="c" ] param = None ):

    cdef int nang, npix

    nang , npix = sinogram.shape[0], sinogram.shape[1]
    myFloat = sinogram.dtype

    image = np.zeros( ( npix , npix ) , dtype=myFloat, order='C' )
    
    cdef float [:,::1] cimage = image

    if param is None:
        param = np.array( [ 0.0 , 0.0 ] ).astype( np.float32 )

    gridrec_v2_backproj( &sinogram[0,0] , npix , nang , &angles[0] , &param[0] ,
                         &cimage[0,0] )

    return image



@cython.boundscheck( False )
@cython.wraparound( False )
def forwproj( np.ndarray[ float , ndim=2 , mode="c" ] image not None ,
              np.ndarray[ float , ndim=1 , mode="c" ] angles not None ,
               ):

    cdef int nang, npix

    npix = image.shape[0]
    nang = len( angles )
    myFloat = image.dtype

    #pdim = 2 * int( 2**( int( np.ceil( np.log10( npix )/np.log10(2) ) ) ) )
    #print pdim
    sino = np.zeros( ( nang , npix ) , dtype=myFloat, order='C' )
    param = np.zeros( 3 , dtype=myFloat , order='C' )
    
    cdef float [:,::1] csino = sino
    cdef float [:] cparam = param

    gridrec_v2_forwproj( &csino[0,0] , npix , nang , &angles[0] , &cparam[0] , &image[0,0] )

    return sino

