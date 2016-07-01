####  PYTHON MODULES
import numpy as np
import sys
#sys.path.append( 'gridding_c/' )
from .gridding_c import gridrec_v2 as gr


####  MY FORMAT VARIABLES & CONSTANTS
myfloat = np.float32
myint = np.int
pi = np.pi
eps = 1e-8 




##########################################################
##########################################################
####                                                  ####
####                CALL OF THE PROJECTORS            ####
####                                                  ####
##########################################################
##########################################################   
    
##  Forward projector
def radon_forward( x , geometry ): 
    #scale = x.space.cell_sides[0]
    return gr.forwproj( x.asarray().astype( myfloat ) , geometry.angles.astype( myfloat ) )
    

    
##  Backprojector
def back_projector( x , geometry ):
    return gr.backproj( x.asarray().astype( myfloat ) , geometry.angles.astype( myfloat ) )
