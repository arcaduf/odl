###########################################################
###########################################################
####                                                   ####
####              TEST ADJOINT OPERATOR                ####
####                                                   ####
###########################################################
###########################################################



####  PYTHON MODULES
from __future__ import division , print_function
import sys
import numpy as np
import gridrec_v2 as grid
import myPrint as pp
import myImageIO as io
import myImageDisplay as dis


####  MY VARIABLE FORMAT
myFloat = np.float32




###########################################################
###########################################################
####                                                   ####
####                         MAIN                      ####
####                                                   ####
###########################################################
###########################################################  

def main():
    print('\n')

    ##  Choose dimension
    dim = int( sys.argv[1] )

    ##  Set parameters
    param = np.array([0.0, 0], dtype=myFloat )

    ##  Random array of values
    angles = np.random.random( dim ).astype(myFloat)

    ##  Random matrix representing the input sinogram
    sino1 = np.random.random(( dim , dim ) ).astype(myFloat)

    ##  Random matrix representing the input image
    img2 = np.random.random(( dim , dim ) ).astype(myFloat)

    
    ##  Output of GRIDREC( sino1 )
    img3 = grid.backproj( sino1 , angles , param )

    ##  Output of ADJOINT_GRIDREC( img2 )
    sino3 = grid.forwproj( img2 , angles )

    ##  Computing  < b , T( a ) >
    innerProduct = np.sum( np.conjugate( img2 ) * img3 )
    print("\n\nCython version I scalar products:")
    print("< b , Tt( a ) >  = ", innerProduct )
                
    ##  Computing  < Tt( b ) , a >
    innerProduct = np.sum( np.conjugate( sino3 ) * sino1 )
    print("< T( b ) , a >  = ", innerProduct )
    print('\n')   
   


if __name__ == '__main__':
    main()
