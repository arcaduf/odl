#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXROTLEVELS 16


/*This file do the process parallell in the fastest direction */

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



void  FILTER_PROCEEDURE_parallell(invector,indiff,H,L,diff,length,normptr,parallell)





REAL *invector;
int indiff;
REAL *L;
REAL *H;
int diff;
int length;
REAL *normptr;
int parallell;
/*
REAL *filterrotations;*/
{
  REAL  norm;
  int rotlevels;
  int halflength;
  int halfrotlevels;
int filterlength;  
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
  REAL tt;
  REAL *t=NULL;
  int j,k;
  REAL *a_0=NULL,*b_0=NULL;
  DECLARE_LOOPVARIABLES_par(a)
  DECLARE_LOOPVARIABLES_par(b) 
  REAL *ptr;
  REAL *ptr0;
  REAL *ptr1;
  REAL *ptr0temp;
  REAL *ptr1temp;
  int ifrotlevelsodd;
  int indiffx2;  
  int v,m;
  REAL *Ltemp;
  REAL *Htemp;
  REAL *Ltempend;
  REAL *Htempend;

  extern   REAL **a_memoryvector;
  extern   REAL **b_memoryvector;
  extern REAL *tempvector;

  indiffx2 =2*indiff;
  halflength=length/2;
  rotlevels= FILTERLENGTH/2;
  halfrotlevels=(rotlevels)>>1;
  filterlength=FILTERLENGTH;
  ifrotlevelsodd =rotlevels%2;
  invectorend=invector+length*indiff*parallell;  
  invectornear1_end=invector+(length-1)*indiff*parallell;  
  invectornear_end=invector+(length-2*halfrotlevels)*indiff*parallell;    
  invectornear_start=invector+ filterlength*indiff*parallell;
  Lend=L+halflength*diff*parallell;
  Lnear_end=L+(halflength-halfrotlevels)*diff*parallell;
  Hend=H+halflength*diff*parallell;
  Hnear_end=H+(halflength-halfrotlevels)*diff*parallell;




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


  Ltemp=tempvector;
#ifndef SKIP_H
  Htemp=tempvector+ 10*parallell;
#endif
if(4*rotlevels<length){ 

 
   ptr=invector;
#ifndef SKIP_L
  ptr0temp= Ltemp;
#endif
#ifndef SKIP_H
  ptr1temp=Htemp;
#endif

ALLOCATE_LOOPVARIABLES_par(a,a_memoryvector);
ALLOCATE_LOOPVARIABLES_par(b,b_memoryvector);

NULLSET_LOOPVARIABLES_par(parallell,a_memoryvector)			\
NULLSET_LOOPVARIABLES_par(parallell,b_memoryvector)			\
a_0=a_memoryvector[0];
b_0=b_memoryvector[0];

  while(ptr<invectornear_start){
if((ptr!=invector)||ifrotlevelsodd){
   a_0=ptr;   
}       
ptr +=indiff*parallell;          
 b_0=ptr;
ptr +=indiff*parallell;          
MULTADDLOOP_par(b,a,-,x,b,parallell) 
#ifndef SKIP_L
LASTMULTADD_par(ptr0temp,a,-,x,b,parallell) ;
ptr0temp +=parallell;
#endif
#ifndef SKIP_H
 LASTMULTADD_par(ptr1temp,b,+,x,a,parallell) ;
ptr1temp +=parallell;
#endif
MULTADD_reverse_LOOP_par(a,b,+,x,a,parallell)     
  }

  ptr0=L+halfrotlevels*parallell;
  ptr1=H+halfrotlevels*parallell;; 

  while(ptr<invectornear1_end){     /*the main loops starts*/
    a_0=ptr;
     ptr +=indiff*parallell;        
     b_0=ptr;    
     ptr +=indiff*parallell;        
     MULTADDLOOP_par(b,a,-,x,b,parallell) 
#ifndef SKIP_L
       LASTMULTADD_par(ptr0,a,-,x,b,parallell) ;
     ptr0+=diff*parallell;
#endif
#ifndef SKIP_H      
     LASTMULTADD_par(ptr1,b,+,x,a,parallell) ;
     ptr1+=diff*parallell;
#endif
      MULTADD_reverse_LOOP_par(a,b,+,x,a,parallell) 
        }                             /*the main loops ends */
 

memset(b_0,0,parallell*sizeof(REAL)); 

#ifndef SKIP_L
  Ltempend=ptr0temp;
  ptr0temp= Ltemp;
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
     a_0=ptr;
    ptr +=indiff*parallell ;       
   }else {
  memset(a_0,0,parallell*sizeof(REAL)); 
   } 
   MULTADDLOOP_par(b,a,-,x,b,parallell) 
#ifndef SKIP_L
     LASTMULTADD_ACUM_par(ptr0,ptr0temp,a,-,x,b,parallell) ;
   ptr0+=diff*parallell;
   ptr0temp+=parallell;
#endif

#ifndef SKIP_H
   LASTMULTADD_ACUM_par(ptr1,ptr1temp,b,+,x,a,parallell) ;
   ptr1+=diff*parallell;
   ptr1temp+=parallell;
#endif  
  MULTADD_reverse_LOOP_par(a,b,+,x,a,parallell) 
  }
 memset(a_0,0,parallell*sizeof(REAL)); 

ptr0=L;
ptr1=H; 
#ifndef SKIP_L
 while(ptr0temp<Ltempend)
#else
   while(ptr1temp<Ltempend)
#endif
{                               
    MULTADDLOOP_par(b,a,-,x,b,parallell) 
#ifndef SKIP_L
    LASTMULTADD_ACUM_par(ptr0,ptr0temp,a,-,x,b,parallell) ;
    ptr0+=diff*parallell;
    ptr0temp +=parallell;
#endif
    
#ifndef SKIP_H 
    LASTMULTADD_ACUM_par(ptr1,ptr1temp,b,+,x,a,parallell) ;
    ptr1+=diff*parallell;
    ptr1temp +=parallell;

#endif
    MULTADD_reverse_LOOP_par(a,b,+,x,a,parallell) 
      }
 }


if(4*rotlevels>=length){
  t=tempvector;
  memset(t,0,length*parallell*sizeof(REAL));
  for(k=0;k<length;k++){
    for(v=0;v<parallell;v++)
      t[k*parallell+v]=invector[k*indiff*parallell+v] ;
  }
  for(j=0;j<rotlevels;j++){
    if((rotlevels-j)&1){
      for(k=0;k<length;k+=2){
	for(v=0;v<parallell;v++){
	  tt=x[j]*t[k*parallell+v];
	  t[k*parallell+v]-=x[j]*t[(k+1)*parallell+v];
	  t[(k+1)*parallell +v]+=tt;
	} }
    }
    if(!((rotlevels-j)&1)){
      for(k=1;k<length-1;k+=2){
	for(v=0;v<parallell;v++){
	  tt=x[j]*t[k*parallell+v];
	  t[k]-=x[j]*t[(k+1)*parallell+v];
	  t[(k+1)*parallell+v]+=tt;
	}}
      
      tt=x[j]*t[length-1];
      t[length-1]-=x[j]*t[0];
      t[0]+=tt;
      
    }
  }
  
	 
	 ptr0=L;
	 ptr1=H;
	 for(k=0;k<length;k+=2){
	   for(v=0;v<parallell;v++){
#ifndef SKIP_L  
    ptr0[v]=norm*t[k*parallell+v];
#endif
#ifndef SKIP_H  
	     ptr1[v]=norm*t[(k+1)*parallell+v];
#endif
	   }
#ifndef SKIP_L  
	     ptr1+=diff*parallell;
#endif
#ifndef SKIP_H  	
     ptr0+=diff*parallell;
	#endif

	 }	 
 }
}










































