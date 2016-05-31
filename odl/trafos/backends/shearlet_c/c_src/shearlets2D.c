#define WITH_MAIN
#define Daub_6
#define NO_INTERPOLATION

#define MAXNSECTOR 200


//#define REAL   double
#define REAL   float


#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "filter_select.h"


#define operator_skip_H(wavelet) wavelet##_skip_H
#define operator_skip_L(wavelet) wavelet##_skip_L
#define operator(wavelet) wavelet

#define operator_par_skip_H(wavelet) wavelet##_par_skip_H
#define operator_par_skip_L(wavelet) wavelet##_par_skip_L
#define operator_par(wavelet) wavelet##_par

#define MULTILOWPASS_xLINE(vptr,xlen,lowpassoffset){	\
outvptr=vptr;\
L=tempvector;\
clen=xlen;\
 for(level=0;level<Level1;level++){\
   if(level==Level1-1)L=lowpassoffset;\
clen_2 =(clen+1)/2;\
 operator_skip_H(daub_6)(outvptr,1,H,L,1,clen,normptr);	\
     outvptr=L;\
     clen=clen_2;}\
     lowpassoffset +=clen;\
  }


#define MULTILOWPASS_parallell_yLINEs(vptr,len,lowpassoffset,parallell){\
outvptr=vptr;\
L=tempvector;\
clen=len;\
 for(level=0;level<Level1;level++){\
  if(level==Level1-1)L=lowpassoffset;\
clen_2 =(clen+1)/2;\
H=L+clen_2;\
 operator_par_skip_H(daub_6)(outvptr,1,H,L,1,clen,normptr,parallell);	\
     outvptr=L;\
     clen=clen_2;\
    }\
    lowpassoffset +=clen*parallell;\
}


REAL **a_memoryvector;
REAL **b_memoryvector;
char *dofilter;
REAL *tempvector;

int filter6(REAL * invector,
	    REAL *filter,
	    REAL **yfilter_alloc,
	    char direction,
	    int xlen,
	    int ylen,
	    int zlen,   // to be set to 1  for 2-dim case
	    REAL *outvector);

int  operator(daub_6)(REAL* invector,
		 int in_diff,
		 REAL *H,
		 REAL *L,
		 int out_diff,
		 int length,
		 REAL *normptr);

 int  operator_par(daub_6)(REAL* invector,
		 int in_diff,
		 REAL *H,
		 REAL *L,
		 int out_diff,
		 int length,
		 REAL *normptr,
		 int parallell);

int  daub_6_skip_H(REAL* invector,
		 int in_diff,
		 REAL *H_none,
		 REAL *L,
		 int out_diff,
		 int length,
		 REAL *normptr);

int  operator_par_skip_H(daub_6)(REAL* invector,
		 int in_diff,
		 REAL *H_none,
		 REAL *L,
		 int out_diff,
		 int length,
		 REAL *normptr,
		 int parallell);

int  operator_skip_L(daub_6)(REAL* invector,
		 int in_diff,
		 REAL *H,
		 REAL *L_none,
		 int out_diff,
		 int length,
		 REAL *normptr);

 int  operator_par_skip_L(daub_6)(REAL* invector,
		 int in_diff,
		 REAL *H,
		 REAL *L_none,
		 int out_diff,
		 int length,
		 REAL *normptr,
		 int parallell);

        

       			 int setup_interpolatevector(REAL** interpolatevptr,int scalelevels);


int shearlets2D(REAL *Image,
                int Lx,
                int Ly,
                int wfilterlength,
                REAL *wfilter,
                int Level0,
                int Level1,
                int nsector,
                REAL *coefficients){
    REAL  *A0=NULL;
    REAL  *B0=NULL;
    REAL  *A=NULL;
    REAL  *B=NULL;
    REAL  *C=NULL;
    REAL *Outvector=NULL;
    REAL *outvector=NULL;
    REAL *filter=NULL;
    REAL *fptr;
    REAL   *ptra, *ptrb,*ptraend;
    REAL *invptr,*invptr2;
    int  *Sh=NULL; 
    int  *iSh=NULL,*iiSh=NULL;
    int *ishp, *iishp;
    int  *newSh=NULL;
    int  *lastSh=NULL;
    int *SH=NULL,*iSH=NULL,*iiSH=NULL;    
    REAL  *Astart; 
    REAL  *Bstart; 
    REAL  *Aend;
    REAL  *Bend;
    REAL *Ap,*Ap0;
    REAL *Bp;
    REAL  a,b;
    int direction;
    unsigned int size;
    int astart, bstart,aend,bend;
    int n,m, last_i;
    int Lmax; 
    int w;
    int alen;
    int current_i;
    REAL *L,*H;  
    REAL *L0;
    REAL *Ltemp;
    REAL *L_none,*H_none;
    REAL * LowLevelCoeff;
    REAL *LowLeveloffset;
    REAL *Invectoroffset;
    REAL *Coeffoffset;
    char scaletype;
    int i,j,k;
    int i0=0,j0=0;
    int l;
    int J[MAXNSECTOR];
    REAL *Interpolatevector=NULL; 
    REAL *normptr;
    REAL Norm=1.0;
    int clen,clen_2;
    REAL *outvptr;
    REAL *yfilter_alloc[6];
    REAL  *fp_3,*fp_2,*fp_1,*fp0,*fp1,*fp2;
    int nsector_2;
    extern char *dofilter;
    char *Dofilter=NULL; 
    int level;
    int cylength0;
    int cxlength0;
    int cxymaxlength0; 
    int xlength = 0;
    int hxlength;
    int lxlength;
    int ylength = 0;
    int hylength;
    int lylength;

    extern REAL **a_memoryvector;
    extern REAL **b_memoryvector;
    extern    REAL *tempvector;

    REAL *tempvector0;
    int u;
    int fl_u,fract_u;
    int n_interpolate=0;
    int  n_noninterpolate=0;
    int Length1,Length2;
    int two20 = 1<<20;
    int two14 = 1<<14;

    nsector_2=(nsector+1)/2;
    for(l=0;l<nsector_2;l++)J[l]=l+nsector_2-1;
    for(l=nsector_2;l<nsector;l++)J[l]=-l+nsector-1;

    tempvector0=NULL;

    for(m=0;m<6;m++)yfilter_alloc[m]=NULL;

    size=Lx * Ly;
    Lmax=(Lx>Ly?Lx:Ly);
    cxlength0=Lx;cylength0=Ly;  

    for(level=0;level<Level1;level++){
      cxlength0=(cxlength0+1)/2;
      cylength0=(cylength0+1)/2;
    }

    cxymaxlength0=cylength0*Lx;

    cxymaxlength0=(cxlength0*Ly>cxymaxlength0?cxlength0*Ly:cxymaxlength0);

 
    /**************MEMORY ALLOCATION**************/
    a_memoryvector=(REAL**)malloc((size_t)((FILTERLENGTH+1)/2)*sizeof(REAL*));
    if(a_memoryvector==NULL)
        return -1;

    b_memoryvector=(REAL**)malloc((size_t)((FILTERLENGTH+1)/2)*sizeof(REAL*));
    if(b_memoryvector==NULL)
        return -1;
  
    for(m=0;m<(FILTERLENGTH+1)/2;m++){
        a_memoryvector[m] = (REAL*)malloc((size_t)Lmax*sizeof(REAL));
        if(a_memoryvector[m]==NULL)
            return -1;
    }

    for(m=0;m<(FILTERLENGTH+1)/2;m++){
        b_memoryvector[m]=(REAL*)malloc((size_t)Lmax*sizeof(REAL));
        if(b_memoryvector[m]==NULL)
            return -1;
    }

    tempvector0=(REAL*)malloc((size_t)20*Lmax*sizeof(REAL));
    if(tempvector0==NULL)
        return -1;

    Dofilter=(char*)malloc((size_t)nsector*Lmax*sizeof(char)); 

    for(m=0;m<6;m++){
        yfilter_alloc[m]=( REAL*)malloc((size_t)(Lmax*sizeof( REAL)));
        if(yfilter_alloc[m]==NULL)
            return -1;
    }
    tempvector=tempvector0;

    A0=( REAL*)malloc((size_t)(size*sizeof( REAL)));
    if(A0==NULL)
        return -1;
    A=A0;

    B0=( REAL*)malloc((unsigned int)(size*sizeof( REAL)));
    if(B0==NULL)
        return -1;
    B=B0;
    C=B;

    SH=(int*)malloc((size_t)Lmax*nsector*sizeof(int));
    iSH=(int*)malloc((size_t)Lmax*nsector*sizeof(int));
    iiSH=(int*)malloc(((size_t)(2*Lmax+1)*nsector*sizeof(int)));
    lastSh=(int*)malloc((size_t)Lmax*sizeof(int));
    newSh=(int*)malloc((size_t)Lmax*sizeof(int));

    #ifndef NO_INTERPOLATION
    Outvector = ( REAL*)malloc((size_t)(size*sizeof( REAL)));
    if(Outvector == NULL)
        return -1;
    #endif

    LowLevelCoeff = ( REAL*)malloc((size_t)nsector*cxymaxlength0 *sizeof( REAL));
    if(LowLevelCoeff == NULL)
        return -1;

    normptr=&Norm;
    LowLeveloffset= LowLevelCoeff;

    #ifndef NO_INTERPOLATION
    Interpolatevector = (REAL*)malloc((size_t)(1000*sizeof( REAL)));
    if(LowLevelCoeff==NULL)
        return -1;

    filter=( REAL*)malloc((size_t)(6*Lmax*sizeof( REAL)));
    if(filter==NULL)
        return -1;

    setup_interpolatevector(&Interpolatevector,6);
    #endif

    i0 = (int)Lx/2;  // (i0,j0) corresponds to origo
    j0 = nsector_2-1;

    for(w=0;w<2;w++){
        if(w==0){
            direction=12;   
            Length1=Lx;Length2=Ly;
        }
        if(w==1){
            direction=21;
            Length1=Ly;Length2=Lx;
        }

        for(i=0;i<Lx;i++)
            lastSh[i]=0.0;

        B=Image;
        j=32;
        for(i=0;i<Lx;i++)
            lastSh[i]=0.0;
        B=Image; 

        for(l=0;l<nsector;l++){
            j=J[l];  //reorder the cone-directions

            Sh=SH+Length2*j;
            iSh=iSH+Length2*j;
            iiSh=iiSH+Length2*j;
            dofilter=Dofilter+j*Length2;

            #ifndef NO_INTERPOLATION
            fptr=filter+j*Length2;
            #endif

            if(j==nsector_2-1){
                for(i=0;i<Length2;i++) {
                    newSh[i]=0;
                    fract_u=0;
                    dofilter[i]=0;
                }
            }
            else{
                if(j==nsector_2 -2){
                    for(i=0;i<Length2;i++)
                        lastSh[i]=0;
                }

                fptr=filter;

                for(i=0;i<Length2;i++){ 
                    u=(i-i0)*(j-j0)+two20;
                    fl_u= ((u&(0x7fffffffffffffe0))>>6) - two14;
                    fract_u = u&(0x3f);
                    newSh[i] = fl_u; 
                    if(fract_u)
                        dofilter[i]=1;
                    else
                        dofilter[i]=0; 

                    #ifndef NO_INTERPOLATION
                    *fptr++ = Interpolatevector[fract_u];
                    *fptr++ = Interpolatevector[fract_u+64];
                    *fptr++ = Interpolatevector[fract_u+128];
                    *fptr++ = Interpolatevector[fract_u+192];
                    *fptr++ = Interpolatevector[fract_u+256];
                    *fptr++ = Interpolatevector[fract_u+320];
                    #endif
                }

                for(i=0;i<Lx;i++)
                    Sh[i] = newSh[i]-lastSh[i];
                for(i=0;i<Lx;i++)
                    lastSh[i]=newSh[i];
                i=1;
                k=1;
                iiSh[0]=0;
                last_i=0;
                while(i<Lx){
                    while((i<Lx)&&(Sh[i]==Sh[i-1]))
                        i++;
                    if(i<Lx){
                        iSh[k]=i;
                        iiSh[k] = i-last_i;
                        last_i=i;
                        iSh[k+1] = (Sh[i]>Sh[i-1]?1:-1);
                        iiSh[k+1] = iiSh[k-1]+iiSh[k] + (Sh[i]>Sh[i-1]?1:-1) * Lx;
                        k +=2;
                        i++;
                    }
                }
            }
        }

        for(l=0;l<nsector;l++){
            j=J[l];  //reorder the cone-directions

            if(j==nsector_2-1){
                /*****NO SHEARING *************/
                if(direction==12){ 
                    outvptr=Image; 
                    LowLeveloffset=LowLevelCoeff+ j*cylength0*Lx; 
                    Coeffoffset=coefficients;	
                    if(Level0>Level1){
                        invptr=B;
                        L=LowLeveloffset;
                        clen=Ly;				
                        for(level=0;level<Level0;level++){		      
                            clen_2 =(clen+1)/2;
                            if(level==Level1)
                                L=Coeffoffset;	
                            operator_par_skip_H(daub_6)(invptr,1,H_none,L,1,clen,normptr,Lx);
                            invptr=L;						
                            clen=clen_2;
                        }
                    }

                    if(Level0<=Level1){
                        invptr=B;
                        L=Coeffoffset;
                        clen=Ly;				
                        for(level=0;level<Level1;level++){		      
                            clen_2 =(clen+1)/2;
                            if(level==(Level0-1))
                                L=LowLeveloffset;	
                            operator_par_skip_H(daub_6)(invptr,1,H_none,L,1,clen,normptr,Lx);
                            invptr=L;						
                            clen=clen_2;
                        }
                    }
         
                    if(Level0==Level1)
                        Invectoroffset = LowLeveloffset;
                    else
                        Invectoroffset = Coeffoffset;

                    for(k=0;k<cylength0;k++){
                        invptr=Invectoroffset;			
                        L=Coeffoffset;
                        clen=Lx;			
                        for(level=0;level<Level0;level++){
                            clen_2 =(clen+1)/2;
                            operator_skip_H(daub_6)(invptr,1,H_none,L,1,Lx,normptr);
                            invptr=L;
                            clen=clen_2;
                        }
                        LowLeveloffset +=clen;
                        Invectoroffset += Lx;	      
                    }
                }
            
                if(direction==21){
                    Invectoroffset=B;	      
                    LowLeveloffset=LowLevelCoeff+(65+j)*cxlength0*Ly; 
                    for(k=0;k<Ly;k++){
                        invptr=Invectoroffset;
                        L=LowLeveloffset;
                        clen=Lx;
                        for(level=0;level<Level1;level++){
                            clen_2 =(clen+1)/2;
                            operator_skip_H(daub_6)(invptr,1,H_none,L,1,clen,normptr);
                            invptr=L;
                            clen=clen_2;
                        }
                        LowLeveloffset +=clen;
                        Invectoroffset += Lx;	      
                    }
                }
                /********END NO SHEARING ********/
            }
            else{
                if(j==nsector_2-2)
                    B=Image;
                //INNER LOOPS STARTS HERE 
                if(direction==21){
                    outvector = LowLevelCoeff + j * cylength0 * Lx; 
                    for(k=0;k<Ly;k++){
                        astart = k * Lx;
                        bstart = k * Lx;
                        current_i = k;
                        ishp = iSh;       
                        for(i=0;i<Lx;i++){
                            if(i == *ishp++){    
                                if(*ishp++ > 0){
                                    astart += Lx;
                                    current_i++;
                                    if(current_i == Ly){
                                        astart -=size;
                                        current_i=0;
                                    }
                                }
                                else{
                                    astart -= Lx;
                                    current_i--;

                                    if(current_i<0){
                                        astart +=size;
                                        current_i=Lx - 1;
                                    }
                                }
                            }    
                            A[astart+i] = B[bstart+i];
                        }
                    }
                }

                if(direction==12){
                    outvector=Outvector;
                    outvptr=outvector;
                    for(k=0;k<Ly;k++){
                        if(Sh[k]>=0)
                            astart = k * Lx + Sh[k];
                        else
                            astart = (k+1) * Lx + Sh[k];
                        bstart = k * Lx;
                        aend = (k+1) * Lx;
                        bend = (k+1) * Lx;
                        m = bstart;
                        for(n=astart;n<aend;n++)
                            A[n] = B[m++];
                        n = bstart;
                        while(m<bend)
                            A[n++] = B[m++];
             
                        #ifndef NO_INTERPOLATION
                        filter6(B + k * Lx, filter + 6 * k, yfilter_alloc, direction, Lx, 1, 1, outvptr);
                        outvptr += Lx;
                        #endif
                    }

                }

                B = A;
                A = C;
                C = B;
                tempvector = A;
                #ifndef NO_INTERPOLATION
                if(direction==21)
                    filter6(B,filter,yfilter_alloc,direction,Lx,Ly,1,outvector);
                if(direction==12){
                    LowLeveloffset=LowLevelCoeff+j*cylength0*Lx; 
                    MULTILOWPASS_parallell_yLINEs(outvptr,Ly,LowLeveloffset,Lx)
                }
                #else
                if(direction==12){
                    LowLeveloffset=LowLevelCoeff+j*cylength0*Lx; 
                    MULTILOWPASS_parallell_yLINEs(B,Ly,LowLeveloffset,Lx)
                }
                if (direction==21){
                    LowLeveloffset=LowLevelCoeff+j*cxlength0*Ly; 
                    for(k=0;k<Ly;k++){
                        MULTILOWPASS_xLINE(B,Ly,LowLeveloffset)
                    }
                }
                #endif

                /**********************/
            }
        }

        if(direction==12){
            xlength=Lx;
            ylength=cylength0;
        }

        if(direction==21){
            xlength=cxlength0;
            ylength=Ly;
        }

        hxlength=xlength/2;
        lxlength=(xlength+1)/2;
        hylength=ylength/2;
        lylength=(ylength+1)/2; 

        Ltemp=tempvector; 
        Invectoroffset=LowLevelCoeff;
        LowLeveloffset=LowLevelCoeff;
        Coeffoffset=coefficients;


        for(level=0;level<Level0;level++){
            if(level==0)
                scaletype='A';
            if(level==1)
                scaletype='B';
            if(level==2)
                scaletype='A';
            if(level==3)
                scaletype='B';

            H=Coeffoffset;
            L=LowLeveloffset;
            invptr=Invectoroffset;  

            if(scaletype=='A'){
                if(direction==12){
                    for(j=0;j<nsector;j++){  
                        if((j&1)==0){
                            for(k=0;k<ylength;k++){
                                operator(daub_6)(invptr,1,H,L,1,xlength,normptr);
                                invptr +=xlength;
                                H +=hxlength;
                                L +=lxlength;
                            }
                        }
                        if((j&1)==1){
                            for(k=0;k<ylength;k++){
                                operator_skip_L(daub_6)( invptr,1,H,L,1,xlength,normptr);
                                invptr +=xlength;
                                H +=hxlength;
                            }
                        }	 	 
                    }
                    xlength=lxlength; 
                    hxlength=xlength/2;
                    lxlength=(xlength+1)/2;
                }

                if(direction==21){
                    for(j=0;j<nsector;j++){
                        if((j&1)==0){
                            operator_par(daub_6)( invptr,1,H,L,1,ylength,normptr,xlength);
                            invptr +=xlength*ylength;
                            H +=hylength*xlength;
                            L +=lylength*xlength;             
                        }
                        if((j&1)==1){
                            operator_par_skip_L(daub_6)( invptr,1,H,L_none,1,ylength,normptr,xlength);
                            invptr +=xlength*ylength;
                            H +=hylength*xlength;
                        }
                    }
                    ylength=lylength; 
                    hylength=ylength/2;
                    lylength=(ylength+1)/2; 
                }
                nsector=(nsector+1)/2;
            }

            if(scaletype=='B'){
                if(direction==12){
                    for(j=0;j<nsector;j++){
                        L0=Ltemp;    
                        for(k=0;k<ylength;k++){
                            operator(daub_6)( invptr,1,H,L0,1,xlength,normptr);
                            invptr +=xlength;
                            H +=hxlength;
                            L0 +=lxlength;             
                        }
                        invptr2=Ltemp;
                        operator_par_skip_H(daub_6)(invptr2,1,H_none,L,1,ylength,normptr,lxlength);
                        L += lxlength * lylength;             
                    }     
                    xlength = lxlength;
                    hxlength = xlength/2;
                    lxlength = (xlength+1)/2;
                    ylength = lylength; 
                    hylength = ylength/2;
                    lylength = (ylength+1)/2; 
                }
            if(direction==21){
                for(j=0;j<nsector;j++){
                    operator_par(daub_6)( invptr,1,H,Ltemp,1,ylength,normptr,xlength);
                    invptr +=xlength*ylength;
                    H +=xlength*hylength;
                    invptr2=Ltemp;   
                    for(k=0;k<lylength;k++){
                        operator_skip_H(daub_6)( invptr2,1,H_none,L,1,xlength,normptr);
                        invptr2 +=xlength;
                        L +=lxlength;             
                    }    
                }
                xlength = lxlength;
                hxlength = xlength/2;
                lxlength = (xlength+1)/2;
                ylength = lylength; 
                hylength = ylength/2;
                lylength = (ylength+1)/2; 
            }
        }
        Coeffoffset=H;
        Invectoroffset=LowLeveloffset; 
        LowLeveloffset=L;  
        }       


    }
    //INNER LOOPS ENDS  HERE 


    /* Free stuff */
    if(A0!=NULL)
        free(A0);
    if(B0!=NULL)
        free(B0);

    if(filter!=NULL)
        free(filter);

    #ifndef NO_INTERPOLATION
    if(Outvector!=NULL)
        free(Outvector);

    for(m=0;m<6;m++){
        if(yfilter_alloc[m]!=NULL)
            free(yfilter_alloc[m]);
        }
    #endif

    if(a_memoryvector!=NULL){
        for(m=0;m<(FILTERLENGTH+1)/2;m++){
            if(a_memoryvector[m]!=NULL)
                free(a_memoryvector[m]);
        }
        free(a_memoryvector);
    }

    if(Interpolatevector!=NULL)
        free(Interpolatevector);

    if(b_memoryvector!=NULL){
        for(m=0;m<(FILTERLENGTH+1)/2;m++){
            if(b_memoryvector[m]!=NULL)
                free(b_memoryvector[m]);
        }
    }

    if(tempvector0!=NULL)
        free(tempvector0);

    if(Sh!=NULL)
        free(Sh); 
    if(iSh!=NULL)
        free(iSh); 
    if(iiSh!=NULL)
        free(iiSh); 
    if(lastSh!=NULL)
        free(lastSh); 
    if(newSh!=NULL)
        free(newSh); 

    return 0;
}







#ifdef WITH_MAIN

 int main(){
    REAL *Image;  
    REAL *coefficients;
    int xlen = (1<<12);
    int ylen = (1<<12);
    int wfilterlength=0;
    REAL *wfilter=NULL;     
   
 int nsector=65;
    int shearinglevels=6;
    int lowpasslevels=5;
    int i;


Image=( REAL*)malloc((size_t)(xlen*sizeof( REAL)));
if(Image==NULL){
   printf("can not allocate Image  with size %dx%d\n",xlen,ylen);
return -1;
 }


coefficients=( REAL*)malloc((unsigned int)(4*xlen*ylen*sizeof( REAL)));
if(coefficients==NULL){
   printf("can not allocate MyMatrix with size %dx%d\n",xlen,ylen);
return -1;
 }

  for(i=0;i<xlen*ylen;i++){
    Image[i]=(REAL)(random()&(0xffff));
}


 shearlets2D(Image,
	     xlen,
	     ylen,
	     wfilterlength,
	     wfilter,
	     lowpasslevels,
	     shearinglevels,
	     nsector,
	     coefficients);
     
  

  if(Image!=NULL)free(Image);
  if(coefficients!=NULL)free(coefficients);


			  return 0;

}

#endif















