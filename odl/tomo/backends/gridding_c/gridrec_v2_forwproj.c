/**********************************************************
 **********************************************************
 ***                                                    ***
 ***             FORWARD REGRIDDING PROJECTOR           ***
 ***                                                    ***
 ***        Written by F. Arcadu on the 19/03/2014      ***
 ***                                                    ***
 **********************************************************
 **********************************************************/




/**********************************************************
 **********************************************************
 ***                                                    ***
 ***                       HEADERS                      ***
 ***                                                    ***
 **********************************************************
 **********************************************************/

#include <math.h>

#ifndef _PWSF_LIB
#define _PSWF_LIB
#include "pswf.h"
#endif

#include <fftw3.h>




/**********************************************************
 **********************************************************
 ***                                                    ***
 ***                        MACROS                      ***
 ***                                                    ***
 **********************************************************
 **********************************************************/

#define Cnvlvnt(X) ( wtbl[ (int)(X+0.5) ] )
#define myAbs(X) ((X)<0 ? -(X) : (X))
#define C 7.0
#define PI 3.141592653589793




/**********************************************************
 **********************************************************
 ***                                                    ***
 ***           FORWARD REGRIDDING PROJECTOR             *** 
 ***                                                    ***
 **********************************************************
 **********************************************************/

void gridrec_v2_forwproj( float *S , int npix , int nang , float *angles , float *param , float *I , char *fftwfn  ) {
    
  /*
   *   Define variables
   */

  int pdim;             //  Number of operative pixels, either equal to "npix" or to "npix" + zero-padding
  int padleft;          //  Pixel position of the last padding zero on the left of the actual projection
  int padright;         //  Pixel position of the first padding zero on the right of the actual projection  
  int pdim_h;           //  Half number of "pdim" pixels
  int npix1, npix2;     //  Pixels at the beginning and at the end of cproj to be saved in the final sinogram
   
  unsigned long j, k, w, n, index;
  unsigned long iul, ivl, iuh, ivh, iv, iu;
  unsigned long ltbl;   // Number of PSWF elements
  
  long ustart, vstart, ufin, vfin;
  
  float Cdata1R, Cdata1I, Cdata2R, Cdata2I, CtmpR, CtmpI;
  float U, V;
  float lconv;        //  Size of the convolution window
  float lconv_h;       //  Half size of the convolution window
  float rtmp;
  float scaling;      //  Convert frequencies to grid units
  float tblspcg;
  float convolv;
  float corrn_u, corrn;
  float ctr;
  float *filter_no;  
  fftwf_complex *cproj;    
  float *SINE, *COSE;
  float *wtbl;
  float *dwtbl;
  float *winv;
  float *work;
  fftwf_complex *H;
  float x;
  float tmp;
  float norm;



  FILE *fp = fopen(fftwfn,"r");
  if(fp){
    fftwf_import_wisdom_from_file(fp); // Load wisdom file for faster FFTW
    fclose(fp);
  }



  /*
   *   Calculate number of operative padded pixels
   */

  pdim = (int) pow( 2 , (int)( ceil( log10( npix )/log10(2) ) ) ) * 2;
  padleft = (int) ( ( pdim - npix ) * 0.5 );
  padright = (int) ( padleft + npix );
  pdim_h = 0.5 * pdim;
  npix1 = (int)( npix * 0.5 ); 
  npix2 = 2 * (int)( pdim - npix1 );


  /*
   *   Allocate memory for complex array storing each
   *   projection and its FFT at time
   */

  cproj = (fftwf_complex *)fftwf_malloc( pdim*sizeof(fftwf_complex));

  for(n=0;n<pdim;n++){
      cproj[n][0]=0.0;
      cproj[n][1]=0.0;
  }    

  /*
   *   Enable scaling
   */

  //rtmp, scaling = npix / (float)pdim ;
  rtmp, scaling = 1.0;
  
  
  
  /*
   *   Initialize look-up-table for sin and cos of the projection angles
   */

  SINE = (float *)calloc( nang , sizeof(float) );
  COSE = (float *)calloc( nang , sizeof(float) );  
  lutTrig(  nang , angles , SINE , COSE );
  
  

  /*
   *   Initialize look-up-table for the PSWF interpolation kernel
   *   and the correction matrix
   */
  
  lconv = (float)( 2*C*1.0 / PI ); 
  lconv_h = lconv * 0.5; 
  ltbl = 2048;
  tblspcg = 2 * ltbl / lconv;  
  wtbl = (float *)calloc( ltbl + 1 , sizeof(float) );
  dwtbl = (float *)calloc( ltbl + 1 , sizeof(float) );  
  winv = (float *)calloc( pdim+1 , sizeof(float) );    
  work = (float *)calloc( (int)lconv + 1 , sizeof(float) );   
  lutPswf( ltbl , pdim_h , wtbl , dwtbl , winv ); 


    
  /*
   *   Allocate memory for cartesian Fourier grid
   */

  H = (fftwf_complex*)fftwf_malloc(pdim*pdim*sizeof(fftwf_complex));
  for(n=0;n<pdim*pdim;n++){
      H[n][0]=0;
      H[n][1]=0;
  }


 
  /*
   *   Correct for center of rotation
   */    

  ctr = npix * 0.5;
  if( pdim != npix )
    ctr += (float) padleft;

  int pdim_d = 2 * pdim;

  filter_no = (float *)calloc( pdim_d , sizeof(float) );
  tmp = (float)( 2 * PI * ctr / pdim );
  norm = (float)( PI / pdim / nang );
  norm = 1.0;

  float tmp1;

  for( j=0,k=0 ; j<pdim ; j+=2,k++ ){
    x = k * tmp;
    float fValue = 1.0;
    filter_no[j] = fValue*norm * cos(x);
    filter_no[j+1] = fValue*norm * sin(x);
  }



  /*
   *   Multiplication for a correction matrix
   */

  j = 0;
  ustart = pdim - pdim_h;
  ufin = pdim;
  
  while (j < pdim) {
    for( iu = ustart ; iu < ufin ; j++ , iu++ ) {
      corrn_u = winv[j];
      k = 0;
      vstart = pdim - pdim_h ;
      vfin = pdim;
      
      while( k < pdim )	{
	    for( iv = vstart ; iv < vfin ; k++ , iv++ ) {
	      corrn = corrn_u * winv[k];

          if( padleft <= j && j < padright && padleft <= k && k < padright ){
	        index = ( npix - 1 - (k-padleft) ) * npix + (j-padleft);
	        H[iu * pdim + iv][0] = corrn * I[ index ];
          }
	    }

        if (k < pdim) {
	      vstart = 0; 
	      vfin = pdim_h;
	    }
      }
    }

    if (j < pdim) {
      ustart = 0;
      ufin = pdim_h;
    }
  }



  /*
   *   Perform 2D FFT of the cartesian Fourier Grid
   */

   fftwf_plan p = fftwf_plan_dft_2d(pdim, pdim,
                                    H, H,
                                    FFTW_BACKWARD, FFTW_ESTIMATE);
   fftwf_execute(p);
 


  /*
   *   Interpolation of the cartesian Fourier grid with PSWF
   */
  
  fftwf_plan p1 = fftwf_plan_dft_1d(pdim, cproj, cproj, FFTW_FORWARD, FFTW_ESTIMATE);

  tmp = 1.0/(float)pdim;
  //float tmp1;

  for( n=0 ; n<nang ; n++ ){ 
    
    /*
     *   Loop on half number of harmonics, because hermitianity
     *   is exploited
     */

    for( j=0, w=0 ; j < pdim_h ; j++, w++ ) {  
      
      Cdata1R = 0.0;
      Cdata1I = 0.0;
      Cdata2R = 0.0;
      Cdata2I = 0.0; 
      
      
      /*
       *   Get cartesian neighbouring points for each polar point
       */

      U = ( rtmp = scaling * w ) * COSE[n] + pdim_h; 
      V = rtmp * SINE[n] + pdim_h;           
      
      iul = (long)ceil( U - lconv_h ); iuh = (long)floor( U + lconv_h );
      ivl = (long)ceil( V - lconv_h ); ivh = (long)floor( V + lconv_h );
      
      if ( iul<0 ) iul = 0; if ( iuh >= pdim ) iuh = pdim-1;   
      if ( ivl<0 ) ivl = 0; if ( ivh >= pdim ) ivh = pdim-1;
      

      for (iv = ivl, k=0; iv <= ivh; iv++, k++)
	    work[k] = Cnvlvnt( myAbs( V - iv ) * tblspcg );

      
      /*
       *   Calculate the contribution of all neighbouring cartesian points
       *   to each polar Fourier sample
       */

      for( iu=iul ; iu<=iuh ; iu++ ){
	    rtmp = Cnvlvnt( myAbs( U - iu ) * tblspcg );
	    
        for( iv=ivl , k=0 ; iv<= ivh ; iv++,k++ ){
	      convolv = rtmp * work[k];
              
          if (iu!=0 && iv!=0 && w!=0){ 
		    Cdata1R += H[iu * pdim + iv][0] * convolv;
            Cdata1I += H[iu * pdim + iv][1] * convolv;
            Cdata2R += H[(pdim-iu) * pdim + pdim - iv][0] * convolv;
            Cdata2I += H[(pdim-iu) * pdim + pdim - iv][1] * convolv;

          } 
           
          else {
            Cdata1R += H[iu * pdim + iv][0] * convolv;
            Cdata1I += H[iu * pdim + iv][1] * convolv;
          }
        }
	  }

      
      CtmpR = filter_no[2*j];
      CtmpI = -filter_no[2*j+1];

      if( j!=0 ){
        cproj[pdim-j][0] = CtmpR * Cdata1R  -  CtmpI * Cdata1I;
        cproj[pdim-j][1] = CtmpR * Cdata1I  +  CtmpI * Cdata1R;  

        CtmpI = -CtmpI;
        
        cproj[j][0] = CtmpR * Cdata2R  -  CtmpI * Cdata2I;
        cproj[j][1] = CtmpR * Cdata2I  +  CtmpI * Cdata2R; 
      }

      else {
        cproj[j][0] = CtmpR * Cdata1R  -  CtmpI * Cdata1I;
        cproj[j][1] = CtmpR * Cdata1I  +  CtmpI * Cdata1R;
      }
    } // End loop on (half) transform data
    
    
    /*
     *   Perform 1D IFFT of cproj
     */
    
    fftwf_execute(p1); 

    
    /*
     *   Assign sinogram elements
     */
    //  Without adjoint of the zero-padding
    for( k=0 , j=padleft ; k<npix ; j++,k++ )
        S[ n*npix + k ] =  tmp * cproj[j][0];
    
    cproj[pdim_h][0] = 0.0;
    cproj[pdim_h][1] = 0.0;
    
  } // End for loop on angles
 
   fftwf_destroy_plan(p1); 
   fftwf_destroy_plan(p); 



  /*
   *  Free memory
   */

  free( SINE );
  free( COSE );
  free( wtbl );
  free( dwtbl );
  free( winv );
  free( work );
  fftwf_free( cproj );
  fftwf_free( H );
  free( filter_no );

  
  return;
}  
