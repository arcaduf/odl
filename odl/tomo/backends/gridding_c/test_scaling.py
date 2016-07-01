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
myfloat = np.float32




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
    n1 = 256
    n2 = 40

    ##  Set parameters
    param = np.array([0.0, 0], dtype=myfloat )

    ##  Random array of values
    angles = np.linspace( 0 , np.pi , 2000 ).astype( myfloat )

    ##  Create square
    phantom = np.zeros( ( n1 , n1 ) ).astype( myfloat )
    n1h      = np.int( n1 * 0.5 )
    phantom[n1h-n2:n1h+n2,n1h-n2:n1h+n2] = 1.0

    ##  Test
    sum1 = np.sum( phantom )
    
    sino = grid.forwproj( phantom , angles )
    sum2 = np.mean( np.sum( sino , axis=1 ) )
    
    dis.plot( phantom , 'Phantom' )    
    dis.plot( grid.forwproj( phantom , angles ) , 'Forward projection' )

    print( 'sum1: ' , sum1 )
    print( 'sum2: ' , sum2 )
    print( 'ratio: ' , sum2 / sum1 )

    print('\n')   
   


if __name__ == '__main__':
    main()
