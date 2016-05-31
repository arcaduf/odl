#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAXROTLEVELS 16

/**********************************************
 Define one of the filter names:

Haar_2  Daub_4  Daub_6  Daub_8 Daub_10 Daub_12
Daub_14 Daub_16 Daub_18 Daub_20 Coif_6 Coif_8
Coif_12 Coif_18 Coif_24 Coif_30   

**********************************************

#define Daub_12

****************************************************
the corresposding procedures will get the name:

haar_2(), daub_4() ,.. and so on
           
        with arguments (invector,step in_difference,highpass,lowpass,step difference in output,length,normptr)
            
             and for the inverse proceedure:

invhaar_2(),invdaub_4(),... and so on

       with arguments (highpass,lowpass,outvector,step differnce in output,length,normptr)
 ***********************************************************************************/   

#include "undefmacros.h"
#include "filter_select.h"
#include "orth1d.h" 

void  FILTER_PROCEEDURE(invector,indiff,H,L,diff,length,normptr)

REAL *invector;
int indiff;
REAL *L;
REAL *H;
int diff;
int length;
REAL *normptr;
//REAL *filterrotations;
{
  REAL  norm;
  int rotlevels;
int filterlength;
  int halflength;
  int halfrotlevels;
  REAL *invectorend;
 REAL *invectornear1_end;
 REAL *invectornear_end;
 REAL *invectornear_start;
  REAL *Lend;
  REAL *Lnear_end;
  REAL * Hend;
  REAL * Hnear_end;
  REAL x[MAXROTLEVELS];
  REAL ix[MAXROTLEVELS];
  REAL a[MAXROTLEVELS];
  REAL b[MAXROTLEVELS];
  REAL tt,t[4*MAXROTLEVELS];
  int j,k;
  REAL a_0=0,b_0=0;
  DECLARE_LOOPVARIABLES(a)
  DECLARE_LOOPVARIABLES(b) 
  REAL *ptr;
  REAL *ptr0;
  REAL *ptr1;
REAL Htemp[9];
REAL Ltemp[9];
 REAL *Ltempend;
 REAL *Htempend;
 REAL s0,d0;


  REAL *ptr0temp;
    REAL *ptr1temp;
  int ifrotlevelsodd;
  int indiff2;  


  indiff2 =2*indiff;
  halflength=length/2;
  rotlevels= FILTERLENGTH/2;
  halfrotlevels=(rotlevels)>>1;
 filterlength =FILTERLENGTH;
 ifrotlevelsodd =rotlevels%2;
  invectorend=invector+length*indiff;  
  invectornear1_end=invector+(length-1)*indiff;  
  invectornear_end=invector+(length-filterlength)*indiff;  
  invectornear_start=invector+filterlength*indiff;  
Lend=L+halflength*diff;
  Lnear_end=L+(halflength-halfrotlevels)*diff;
  Hend=H+halflength*diff;
  Hnear_end=H+(halflength-halfrotlevels)*diff;
if(rotlevels>MAXROTLEVELS){
      printf("MAXROTLEVELS is exceeded\n");
    return  ;
} 



 norm=*normptr;
  INIT_x
    if(norm==1.0){
  for(j=0;j<rotlevels;j++){
      norm *=(1+x[j]*x[j]);
  }
  if(norm!=1.0)norm=1.0/sqrt((double)norm);
  *normptr=norm;      
 }



if(4*rotlevels<length){ 

   ptr=invector;
#ifndef SKIP_L 
 ptr0temp=Ltemp;
#endif
#ifndef SKIP_H
  ptr1temp=Htemp;
#endif

  while(ptr<invectornear_start){
if((ptr!=invector)||ifrotlevelsodd){
    a_0=(*ptr);   
ptr +=indiff;     
}       
  
b_0=(*(ptr));
ptr +=indiff;          
MULTADDLOOP(b,a,-,x,b) 
#ifndef SKIP_L
*ptr0temp++ =norm*LASTMULTADD(a,-,x,b) ;
#endif
#ifndef SKIP_H
*ptr1temp++=norm*LASTMULTADD(b,+,x,a) 
#endif
MULTADD_reverse_LOOP(a,b,+,x,a)     
  }


  ptr0=L+halfrotlevels;
  ptr1=H+halfrotlevels; 

  while(ptr<invectornear1_end){     /*the main loops starts*/
    a_0=(*(ptr));
      ptr +=indiff;        
      b_0=(*(ptr));
      ptr +=indiff;        
       MULTADDLOOP(b,a,-,x,b) 
#ifndef SKIP_L
	 *ptr0=norm*LASTMULTADD(a,-,x,b)
#endif 
	 ptr0+=diff;
#ifndef SKIP_H
      *ptr1=norm*LASTMULTADD(b,+,x,a)
       ptr1+=diff;
#endif
      MULTADD_reverse_LOOP(a,b,+,x,a) 
        }
 

                             /*the main loops ends */
  b_0=0;
#ifndef SKIP_L
  Ltempend=ptr0temp;  
  ptr0temp=Ltemp;
#endif
#ifndef SKIP_H  
  Htempend=ptr1temp;
  ptr1temp=Htemp;
  while(ptr1<Hend)
#else
  while(ptr0<Lend)
#endif
    {
      if(ptr<invectorend){
	a_0=*(ptr);
	ptr +=indiff ;       
      }else {
	a_0=0;
      }
      MULTADDLOOP(b,a,-,x,b) 
#ifndef SKIP_L
	*ptr0 =*(ptr0temp++) +norm*LASTMULTADD(a,-,x,b)
	ptr0+=diff;
#endif
#ifndef SKIP_H
        *ptr1=(*ptr1temp++)+norm*LASTMULTADD(b,+,x,a)
	ptr1+=diff;
#endif
      MULTADD_reverse_LOOP(a,b,+,x,a) 
  a_0=0; 
	}

	
  
    ptr0=L;
    ptr1=H; 
    
#ifndef SKIP_L
    while(ptr0temp <Ltempend){
#else
   while(ptr1temp <Htempend){
	
#endif                               
	MULTADDLOOP(b,a,-,x,b) 
#ifndef SKIP_L
	  *ptr0 =*(ptr0temp++) +norm*LASTMULTADD(a,-,x,b)  
	  ptr0+=diff;
#endif
	
#ifndef SKIP_H
	*ptr1=(*ptr1temp++)+norm*LASTMULTADD(b,+,x,a)
	  ptr1+=diff;
#endif
	MULTADD_reverse_LOOP(a,b,+,x,a) 
	  }
 }  

if(4*rotlevels>=length){
  for(k=0;k<length;k++){
    t[k]=invector[k*indiff] ;
}
  for(j=0;j<rotlevels;j++){
    if((rotlevels-j)&1){
      for(k=0;k<length;k+=2){
	tt=x[j]*t[k];
	t[k]-=x[j]*t[k+1];
	t[k+1]+=tt;
      }
    }
    if(!((rotlevels-j)&1)){
      for(k=1;k<length-1;k+=2){
	tt=x[j]*t[k];
	t[k]-=x[j]*t[k+1];
	t[k+1]+=tt;
      }
      tt=x[j]*t[length-1];
      t[length-1]-=x[j]*t[0];
      t[0]+=tt;
      
    }
  }
  ptr0=L;
  ptr1=H;
  for(k=0;k<length;k+=2){
#ifndef SKIP_L  
  *ptr0=norm*t[k];
   ptr0+=diff;
#endif
  #ifndef SKIP_H
    *ptr1=norm*t[k+1];
     ptr1+=diff;
#endif
  }
 }
}

    

void  INVERSE_FILTER_PROCEEDURE(H,L,outvector,diff,length,normptr)

REAL *L;
REAL *H;
REAL *outvector;
int diff;
int length;
REAL *normptr;
/*int filterlength;
REAL *filterrotations;*/
{
  REAL  norm;
  int rotlevels;
  int halflength;
  int halfrotlevels;
  REAL * outvectorend;
  REAL *Lend;
  REAL *Hend;
  REAL x[MAXROTLEVELS];
  REAL ix[MAXROTLEVELS];
  REAL a[MAXROTLEVELS];
  REAL b[MAXROTLEVELS];
  REAL a_0 = 0, b_0 = 0;
  DECLARE_LOOPVARIABLES(a)
  DECLARE_LOOPVARIABLES(b) 
  REAL tt,t[4*MAXROTLEVELS];
  int j,k;
  REAL *ptr;
  REAL *ptr0;
  REAL *ptr1;
  int ifrotlevelsodd;


  halflength=length/2;
  rotlevels= FILTERLENGTH/2;
  halfrotlevels=rotlevels>>1;
  ifrotlevelsodd =rotlevels%2;
  outvectorend=outvector+length*diff;  
  Lend=L+halflength;
  Hend=H+halflength;
  rotlevels= FILTERLENGTH/2;
  if(rotlevels>MAXROTLEVELS){
  printf("MAXROTLEVELS is exceeded\n");
    return  ;
} 

 norm=*normptr;
  INIT_x
    if(norm==1.0){
  for(j=0;j<rotlevels;j++){
  /*  x[j]=filterrotations[j];*/
    norm *=(1+x[j]*x[j]);
  }
  if(norm!=1.0)norm=1.0/sqrt((REAL)norm);
  *normptr=norm;   
 }

if(4*rotlevels<length){

  ptr0=L;
  ptr1=H;
  ptr=outvectorend-(rotlevels-1)*diff;
  while((ptr<outvectorend)&&(ptr>outvector+diff)){
       a_0=(*(ptr0++))*norm;      
  b_0=(*(ptr1++))*norm;         
       MULTADDLOOP(b,a,+,ix,b) 
      *ptr=LASTMULTADD(a,+,ix,b) ;
      ptr+=diff;
    if(ptr==outvectorend) ptr=outvector;
     *ptr=LASTMULTADD(b,-,ix,a) ;
      ptr+=diff;
      MULTADD_reverse_LOOP(a,b,-,ix,a) 
  }
  if(ptr==outvectorend) ptr=outvector;

  while(ptr0<Lend){
  a_0=(*(ptr0++))*norm;       
  b_0=(*(ptr1++))*norm;       
      MULTADDLOOP(b,a,+,ix,b) 
      *ptr=LASTMULTADD(a,+,ix,b) ;
      ptr+=diff;
      *ptr=LASTMULTADD(b,-,ix,a) ;
      ptr+=diff;
      MULTADD_reverse_LOOP(a,b,-,ix,a) 
     }
  a_0=0;            
  b_0=0;            
  while((ptr<outvectorend)&&(ptr>outvector+diff)){
     MULTADDLOOP(b,a,+,ix,b) 
      *ptr+=LASTMULTADD(a,+,ix,b) ;
      ptr+=diff;
  if(ptr==outvectorend) ptr=outvector;
      *ptr+=LASTMULTADD(b,-,ix,a) ;
      ptr+=diff;
      MULTADD_reverse_LOOP(a,b,-,ix,a) 
   }
  if(ptr==outvectorend) ptr=outvector;
  ptr0=Lend-halfrotlevels;                             
while(ptr0<Lend){                               
     MULTADDLOOP(b,a,+,ix,b) 
      *ptr+=LASTMULTADD(a,+,ix,b) ;
      ptr+=diff;
       *ptr+=LASTMULTADD(b,-,ix,a) ;
      ptr+=diff;
      MULTADD_reverse_LOOP(a,b,-,ix,a) 
     
  ptr0++;
   }
 }  
if(4*rotlevels>=length){

  ptr0=L;
  ptr1=H;
  for(k=0;k<length;k+=2){
    t[k]= *(ptr0++) * norm;
    t[k+1]=*(ptr1++) * norm;
}
   for(j=0;j<rotlevels;j++){
    if(!(j&1)){
      for(k=0;k<length;k+=2){
	tt=ix[j]*t[k];
	t[k]+=ix[j]*t[k+1];
	t[k+1]-=tt;
      }
    }
    if(j&1){
      for(k=1;k<length-1;k+=2){
	tt=ix[j]*t[k];
	t[k]+=ix[j]*t[k+1];
	t[k+1]-=tt;
      }
      tt=ix[j]*t[length-1];
      t[length-1]+=ix[j]*t[0];
      t[0]-=tt;
      
    }
  }
  for(k=0;k<length;k++){
    outvector[k*diff]=t[k];
}
}
}









































