#include <stdio.h>
#include <stdlib.h>


//**********************************************
// Define one of the filter names:

#define  MANYFILTERS
#define Haar_2
#include "orth1d.c"
#undef Haar_2
#define Daub_4
#include "orth1d.c" 
#undef Daub_4 
#define Daub_6  
#include "orth1d.c"
#undef Daub_6  
#define Daub_8 
#include "orth1d.c"
#undef Daub_8 
#define Daub_10 
#include "orth1d.c"
#undef Daub_10 
#define Daub_12
#include "orth1d.c"
#undef Daub_12
#define Daub_14 
#include "orth1d.c"
#undef Daub_14 
#define Daub_16 
#include "orth1d.c"
#undef Daub_16 
#define Daub_18 
#include "orth1d.c"
#undef Daub_18 
#define Daub_20 
#include "orth1d.c"
#undef Daub_20 
#define Coif_6 
#include "orth1d.c"
#undef Coif_6 
#define Coif_8
#include "orth1d.c"
#undef Coif_8
#define Coif_12 
#include "orth1d.c"
#undef Coif_12 
#define Coif_18 
#include "orth1d.c"
#undef Coif_18 
#define Coif_24 
#include "orth1d.c"
#undef Coif_24 
#define Coif_30   
#include "orth1d.c"
#undef Coif_30   

/**********************************************

#define Daub_12

****************************************************
the corresposding procedures will get the name:

haar_2(), daub_4() ,.. and so on
           
        with arguments (invector,step in_difference,highpass,lowpass,step difference in output,length,normptr)
            
             and for the inverse proceedure:

invhaar_2(),invdaub_4(),... and so on

       with arguments (highpass,lowpass,outvector,step differnce in output,length,normptr)
 ***********************************************************************************/   


