//#define REAL double
#define REAL float


#if(FILTERLENGTH==2)
#define ix_0 x_0	  

#endif
#if(FILTERLENGTH==4)
#define ix_1 x_0	  
#define ix_0 x_1	  

#endif
#if(FILTERLENGTH==6)
#define ix_2  x_0	  
#define ix_1  x_1	  
#define ix_0  x_2	  

#endif
#if(FILTERLENGTH==8)
#define ix_3 x_0	  
#define ix_2 x_1	  
#define ix_1 x_2	  
#define ix_0 x_3	  

#endif
#if(FILTERLENGTH==10)

#define ix_4 x_0	  
#define ix_3 x_1	  
#define ix_2 x_2	  
#define ix_1 x_3	  
#define ix_0 x_4	  

#endif
#if(FILTERLENGTH==12)

#define ix_5 x_0	  
#define ix_4 x_1	  
#define ix_3 x_2  
#define ix_2 x_3	  
#define ix_1 x_4	  
#define ix_0 x_5  

#endif
#if(FILTERLENGTH==14)

#define ix_6  x_0	  
#define ix_5  x_1	  
#define ix_4  x_2	  
#define ix_3  x_3	  
#define ix_2  x_4	  
#define ix_1  x_5	  
#define ix_0  x_6	  

#endif
#if(FILTERLENGTH==16)

#define ix_7  x_0	  
#define ix_6  x_1	  
#define ix_5  x_2	  
#define ix_4  x_3	  
#define ix_3  x_4	  
#define ix_2  x_5	  
#define ix_1  x_6	  
#define ix_0  x_7	  

#endif
#if(FILTERLENGTH==18)

#define ix_8  x_0	  
#define ix_7  x_1	  
#define ix_6  x_2	  
#define ix_5  x_3	  
#define ix_4  x_4	  
#define ix_3  x_5	  
#define ix_2  x_6	  
#define ix_1  x_7	  
#define ix_0  x_8	  

#endif
#if(FILTERLENGTH==20)

#define ix_9  x_0	  
#define ix_8  x_1	  
#define ix_7  x_2	  
#define ix_6  x_3	  
#define ix_5  x_4	  
#define ix_4  x_5	  
#define ix_3  x_6	  
#define ix_2  x_7	  
#define ix_1  x_8	  
#define ix_0  x_9	  

#endif
#if(FILTERLENGTH==22)


#define ix_10  x_0                   	  
#define ix_9   x_1	  
#define ix_8   x_2	  
#define ix_7   x_3	  
#define ix_6   x_4	  
#define ix_5   x_5	  
#define ix_4   x_6	  
#define ix_3   x_7	  
#define ix_2   x_8	  
#define ix_1   x_9	  
#define ix_0   x_10	  

#endif
#if(FILTERLENGTH==24)

#define ix_11     x_0 
#define ix_10     x_1              	  
#define ix_9	  x_2 
#define ix_8	  x_3
#define ix_7	  x_4
#define ix_6	  x_5
#define ix_5	  x_6
#define ix_4	  x_7
#define ix_3	  x_8 
#define ix_2	  x_9
#define ix_1	  x_10
#define ix_0	  x_11

#endif
#if(FILTERLENGTH==26)

#define ix_12   x_0
#define ix_11   x_1
#define ix_10   x_2                	  
#define ix_9	x_3
#define ix_8	x_4  
#define ix_7	x_5  
#define ix_6	x_6  
#define ix_5	x_7  
#define ix_4	x_8  
#define ix_3	x_9  
#define ix_2	x_10  
#define ix_1	x_11  
#define ix_0	x_12 

#endif
#if(FILTERLENGTH==28)

#define ix_13  x_0
#define ix_12  x_1
#define ix_11  x_2 
#define ix_10  x_3                 	  
#define ix_9   x_4  
#define ix_8   x_5  
#define ix_7   x_6	  
#define ix_6   x_7	  
#define ix_5   x_8	  
#define ix_4   x_9	  
#define ix_3   x_10	  
#define ix_2   x_11	  
#define ix_1   x_12	  
#define ix_0   x_13	  

#endif

#if(FILTERLENGTH==30)

#define ix_14 x_0
#define ix_13 x_1
#define ix_12 x_2
#define ix_11 x_3
#define ix_10 x_4
#define ix_9  x_5
#define ix_8  x_6
#define ix_7  x_7
#define ix_6  x_8
#define ix_5  x_9
#define ix_4  x_10
#define ix_3  x_11
#define ix_2  x_12
#define ix_1  x_13
#define ix_0  x_14

#endif
      
#if(FILTERLENGTH ==2)
#define DECLARE_LOOPVARIABLES(a) 
#define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0 ;

#define LASTMULTADD  LASTMULTADD_0
#define MULTADDLOOP MULTADD_0
#define MULTADD_reverse_LOOP MULTADD_rev_0
#endif


#if(FILTERLENGTH ==4)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0;

#define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[1]=x_1 ;  ix[1]=ix_1 ;

#define LASTMULTADD LASTMULTADD_1
#define MULTADDLOOP MULTADD_1
#define MULTADD_reverse_LOOP MULTADD_rev_1
#endif


#if(FILTERLENGTH ==6)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0,a##_2=0;

#define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[1]=x_1 ;  ix[1]=ix_1 ;\
 x[2]=x_2 ;  ix[2]=ix_2;

#define LASTMULTADD LASTMULTADD_2
#define MULTADDLOOP MULTADD_2
#define MULTADD_reverse_LOOP MULTADD_rev_2
 #endif


#if(FILTERLENGTH ==8)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0,a##_2=0,a##_3=0;

 #define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[1]=x_1 ;  ix[1]=ix_1 ;\
 x[2]=x_2 ;  ix[2]=ix_2;\
 x[3]=x_3 ;  ix[3 ]=ix_3 ;

#define LASTMULTADD LASTMULTADD_3
#define MULTADDLOOP MULTADD_3
#define MULTADD_reverse_LOOP MULTADD_rev_3
#endif

#if(FILTERLENGTH ==10)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0,a##_2=0,a##_3=0,a##_4=0;
 
 #define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0;\
 x[1]=x_1 ;  ix[1]=ix_1 ;\
 x[2]=x_2 ;  ix[2]=ix_2;\
 x[3]=x_3 ;  ix[3 ]=ix_3 ;\
 x[4]=x_4 ;  ix[4 ]=ix_4 ;

#define LASTMULTADD LASTMULTADD_4
#define MULTADDLOOP MULTADD_4
#define MULTADD_reverse_LOOP MULTADD_rev_4
#endif

#if(FILTERLENGTH ==12)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0,a##_2=0,a##_3=0,a##_4=0,a##_5=0;

 #define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[1]=x_1 ;  ix[1]=ix_1 ;\
 x[2]=x_2 ;  ix[2]=ix_2;\
 x[3]=x_3 ;  ix[3 ]=ix_3 ;\
 x[4]=x_4 ;  ix[4 ]=ix_4 ;\
 x[5]=x_5 ;  ix[5 ]=ix_5 ;

#define LASTMULTADD LASTMULTADD_5
#define MULTADDLOOP MULTADD_5
#define MULTADD_reverse_LOOP MULTADD_rev_5
#endif

#if(FILTERLENGTH ==14)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0,a##_2=0,a##_3=0,a##_4=0,a##_5=0,a##_6=0;

 #define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[1]=x_1 ;  ix[1]=ix_1 ;\
 x[2]=x_2 ;  ix[2]=ix_2;\
 x[3]=x_3 ;  ix[3 ]=ix_3 ;\
 x[4]=x_4 ;  ix[4 ]=ix_4 ;\
 x[5]=x_5 ;  ix[5 ]=ix_5 ;\
 x[6]=x_6 ;  ix[6 ]=ix_6 ;\

#define LASTMULTADD LASTMULTADD_6
#define MULTADDLOOP MULTADD_6
#define MULTADD_reverse_LOOP MULTADD_rev_6
#endif

#if(FILTERLENGTH ==16)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0,a##_2=0,a##_3=0,a##_4=0,a##_5=0,a##_6=0,a##_7=0;

 #define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[1]=x_1 ;  ix[1]=ix_1 ;\
 x[2]=x_2 ;  ix[2]=ix_2;\
 x[3]=x_3 ;  ix[3 ]=ix_3 ;\
 x[4]=x_4 ;  ix[4 ]=ix_4 ;\
 x[5]=x_5 ;  ix[5 ]=ix_5 ;\
 x[6]=x_6 ;  ix[6 ]=ix_6 ;\
 x[7]=x_7 ;  ix[7 ]=ix_7 ;

#define LASTMULTADD LASTMULTADD_7
#define MULTADDLOOP MULTADD_7
#define MULTADD_reverse_LOOP MULTADD_rev_7
#endif

#if(FILTERLENGTH ==18)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0,a##_2=0,a##_3=0,a##_4=0,a##_5=0,a##_6=0,a##_7=0,a##_8=0;

#define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[1]=x_1 ;  ix[1]=ix_1 ;\
 x[2]=x_2 ;  ix[2]=ix_2;\
 x[3]=x_3 ;  ix[3 ]=ix_3 ;\
 x[4]=x_4 ;  ix[4 ]=ix_4 ;\
 x[5]=x_5 ;  ix[5 ]=ix_5 ;\
 x[6]=x_6 ;  ix[6 ]=ix_6 ;\
 x[7]=x_7 ;  ix[7 ]=ix_7 ;\
 x[8]=x_8 ;  ix[8 ]=ix_8 ;

#define LASTMULTADD LASTMULTADD_8
#define MULTADDLOOP MULTADD_8
#define MULTADD_reverse_LOOP MULTADD_rev_8
#endif


#if(FILTERLENGTH ==20)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0,a##_2=0,a##_3=0,a##_4=0,a##_5=0,a##_6=0,a##_7=0,a##_8=0,a##_9=0;

#define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[1]=x_1 ;  ix[1]=ix_1 ;\
 x[2]=x_2 ;  ix[2]=ix_2;\
 x[3]=x_3 ;  ix[3 ]=ix_3 ;\
 x[4]=x_4 ;  ix[4 ]=ix_4 ;\
 x[5]=x_5 ;  ix[5 ]=ix_5 ;\
 x[6]=x_6 ;  ix[6 ]=ix_6 ;\
 x[7]=x_7 ;  ix[7 ]=ix_7 ;\
 x[8]=x_8 ;  ix[8 ]=ix_8 ;\
 x[9]=x_9 ;  ix[9 ]=ix_9 ;\


#define LASTMULTADD LASTMULTADD_9
#define MULTADDLOOP MULTADD_9
#define MULTADD_reverse_LOOP MULTADD_rev_9
#endif

#if(FILTERLENGTH ==22)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0,a##_2=0,a##_3=0,a##_4=0,a##_5=0,a##_6=0,a##_7=0,a##_8=0,a##_9=0,\
a##_10=0;

#define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[1]=x_1 ;  ix[1]=ix_1 ;\
 x[2]=x_2 ;  ix[2]=ix_2;\
 x[3]=x_3 ;  ix[3 ]=ix_3 ;\
 x[4]=x_4 ;  ix[4 ]=ix_4 ;\
 x[5]=x_5 ;  ix[5 ]=ix_5 ;\
 x[6]=x_6 ;  ix[6 ]=ix_6 ;\
 x[7]=x_7 ;  ix[7 ]=ix_7 ;\
 x[8]=x_8 ;  ix[8 ]=ix_8 ;\
 x[9]=x_9 ;  ix[9 ]=ix_9 ;\
 x[10]=x_10 ;  ix[10 ]=ix_10 ;

#define LASTMULTADD LASTMULTADD_10
#define MULTADDLOOP MULTADD_10
#define MULTADD_reverse_LOOP MULTADD_rev_10
#endif

#if(FILTERLENGTH ==24)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0,a##_2=0,a##_3=0,a##_4=0,a##_5=0,a##_6=0,a##_7=0,a##_8=0,a##_9=0,\
a##_10=0,a##_11=0;

#define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[1]=x_1 ;  ix[1]=ix_1 ;\
 x[2]=x_2 ;  ix[2]=ix_2;\
 x[3]=x_3 ;  ix[3 ]=ix_3 ;\
 x[4]=x_4 ;  ix[4 ]=ix_4 ;\
 x[5]=x_5 ;  ix[5 ]=ix_5 ;\
 x[6]=x_6 ;  ix[6 ]=ix_6 ;\
 x[7]=x_7 ;  ix[7 ]=ix_7 ;\
 x[8]=x_8 ;  ix[8 ]=ix_8 ;\
 x[9]=x_9 ;  ix[9 ]=ix_9 ;\
 x[10]=x_10 ;  ix[10 ]=ix_10 ;\
 x[11]=x_11 ;  ix[11 ]=ix_11 ;

#define LASTMULTADD LASTMULTADD_11
#define MULTADDLOOP MULTADD_11
#define MULTADD_reverse_LOOP MULTADD_rev_11
#endif

#if(FILTERLENGTH ==26)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0,a##_2=0,a##_3=0,a##_4=0,a##_5=0,a##_6=0,a##_7=0,a##_8=0,a##_9=0,\
a##_10=0,a##_11=0,a##_12=0;

#define INIT_x \ 
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[1]=x_1 ;  ix[1]=ix_1 ;\
 x[2]=x_2 ;  ix[2]=ix_2;\
 x[3]=x_3 ;  ix[3 ]=ix_3 ;\
 x[4]=x_4 ;  ix[4 ]=ix_4 ;\
 x[5]=x_5 ;  ix[5 ]=ix_5 ;\
 x[6]=x_6 ;  ix[6 ]=ix_6 ;\
 x[7]=x_7 ;  ix[7 ]=ix_7 ;\
 x[8]=x_8 ;  ix[8 ]=ix_8 ;\
 x[9]=x_9 ;  ix[9 ]=ix_9 ;\
 x[10]=x_10 ;  ix[10 ]=ix_10 ;\
 x[11]=x_ ;  ix[11 ]=ix_11 ;\
 x[12]=x_12 ;  ix[12]=ix_12 ;

#define LASTMULTADD LASTMULTADD_12
#define MULTADDLOOP MULTADD_12
#define MULTADD_reverse_LOOP MULTADD_rev_12
#endif

#if(FILTERLENGTH ==28)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0,a##_2=0,a##_3=0,a##_4=0,a##_5=0,a##_6=0,a##_7=0,a##_8=0,a##_9=0,\
a##_10=0,a##_11=0,a##_12=0,a##_13=0;

#define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[1]=x_1 ;  ix[1]=ix_1 ;\
 x[2]=x_2 ;  ix[2]=ix_2;\
 x[3]=x_3 ;  ix[3 ]=ix_3 ;
 x[4]=x_4 ;  ix[4 ]=ix_4 ;\
 x[5]=x_5 ;  ix[5 ]=ix_5 ;\
 x[6]=x_6 ;  ix[6 ]=ix_6 ;\
 x[7]=x_7 ;  ix[7 ]=ix_7 ;\
 x[8]=x_8 ;  ix[8 ]=ix_8 ;\
 x[9]=x_9 ;  ix[9 ]=ix_9 ;\
 x[10]=x_10 ;  ix[10 ]=ix_10 ;\
 x[11]=x_ ;  ix[11 ]=ix_11 ;\
 x[12]=x_12 ;  ix[12]=ix_12 ;\
 x[13]=x_13 ;  ix[13 ]=ix_13 ;

#define LASTMULTADD LASTMULTADD_13
#define MULTADDLOOP MULTADD_13
#define MULTADD_reverse_LOOP MULTADD_rev_13
#endif

#if(FILTERLENGTH ==30)
#define DECLARE_LOOPVARIABLES(a) \
REAL a##_1=0,a##_2=0,a##_3=0,a##_4=0,a##_5=0,a##_6=0,a##_7=0,a##_8=0,a##_9=0,\
a##_10=0,a##_11=0,a##_12=0,a##_13=0,a##_14=0;

#define INIT_x \
 x[0]=x_0 ;  ix[0]=ix_0 ;\
 x[1]=x_1 ;  ix[1]=ix_1 ;\
 x[2]=x_2 ;  ix[2]=ix_2;\
 x[3]=x_3 ;  ix[3 ]=ix_3 ;\
 x[4]=x_4 ;  ix[4 ]=ix_4 ;\
 x[5]=x_5 ;  ix[5 ]=ix_5 ;\
 x[6]=x_6 ;  ix[6 ]=ix_6 ;\
 x[7]=x_7 ;  ix[7 ]=ix_7 ;\
 x[8]=x_8 ;  ix[8 ]=ix_8 ;\
 x[9]=x_9 ;  ix[9 ]=ix_9 ;\
 x[10]=x_10 ;  ix[10 ]=ix_10 ;\
 x[11]=x_11 ;  ix[11 ]=ix_11 ;\
 x[12]=x_12 ;  ix[12]=ix_12 ;\
 x[13]=x_13 ;  ix[13 ]=ix_13 ;\
 x[14]=x_14 ;  ix[14]=ix_14 ;


#define LASTMULTADD LASTMULTADD_14
#define MULTADDLOOP MULTADD_14
#define MULTADD_reverse_LOOP MULTADD_rev_14
#endif


#define MULTADD_0(ab,a,sign,x,b) ;

#define MULTADD_1(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;

#define MULTADD_2(ab,aG,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;\
     ab##_2=a##_1  sign  x##_1*b##_1;

#define MULTADD_3(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_3=a##_2  sign  x##_2*b##_2;

#define MULTADD_4(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_4=a##_3  sign  x##_3*b##_3;

#define MULTADD_5(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_5=a##_4  sign  x##_4*b##_4;

#define MULTADD_6(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_6=a##_5  sign  x##_5*b##_5;

#define MULTADD_7(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_7=a##_6  sign  x##_6*b##_6;

#define MULTADD_8(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_8=a##_7  sign  x##_7*b##_7;

#define MULTADD_9(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_8=a##_7  sign  x##_7*b##_7;\
     ab##_9=a##_8  sign  x##_8*b##_8;

#define MULTADD_10(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_8=a##_7  sign  x##_7*b##_7;\
     ab##_9=a##_8  sign  x##_8*b##_8;\
     ab##_10=a##_9  sign  x##_9*b##_9;

#define MULTADD_11(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_8=a##_7  sign  x##_7*b##_7;\
     ab##_9=a##_8  sign  x##_8*b##_8;\
     ab##_10=a##_9  sign  x##_9*b##_9;\
     ab##_11=a##_10  sign  x##_10*b##_10;

#define MULTADD_12(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_8=a##_7  sign  x##_7*b##_7;\
     ab##_9=a##_8  sign  x##_8*b##_8;\
     ab##_10=a##_9  sign  x##_9*b##_9;\
     ab##_11=a##_10  sign  x##_10*b##_10;\
     ab##_12=a##_11  sign  x##_11*b##_11;

#define MULTADD_13(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_8=a##_7  sign  x##_7*b##_7;\
     ab##_9=a##_8  sign  x##_8*b##_8;\
     ab##_10=a##_9  sign  x##_9*b##_9;\
     ab##_11=a##_10  sign  x##_10*b##_10;\
     ab##_12=a##_11  sign  x##_11*b##_11;\
     ab##_13=a##_12  sign  x##_12*b##_12;


#define MULTADD_14(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_8=a##_7  sign  x##_7*b##_7;\
     ab##_9=a##_8  sign  x##_8*b##_8;\
     ab##_10=a##_9  sign  x##_9*b##_9;\
     ab##_11=a##_10  sign  x##_10*b##_10;\
     ab##_12=a##_11  sign  x##_11*b##_11;\
     ab##_13=a##_12  sign  x##_12*b##_12;\
     ab##_14=a##_13  sign  x##_13*b##_13;

#define LASTMULTADD_0(a,sign,x,b) \
  (a##_0 sign x##_0*b##_0);

#define LASTMULTADD_1(a,sign,x,b) \
  (a##_1 sign x##_1*b##_1);

#define LASTMULTADD_2(a,sign,x,b) \
  (a##_2 sign x##_2*b##_2);

#define LASTMULTADD_3(a,sign,x,b) \
  (a##_3 sign x##_3*b##_3);

#define LASTMULTADD_4(a,sign,x,b)		\
  (a##_4 sign x##_4*b##_4);

#define LASTMULTADD_5(a,sign,x,b) \
  (a##_5 sign x##_5*b##_5);

#define LASTMULTADD_6(a,sign,x,b) \
  (a##_6 sign x##_6*b##_6);

#define LASTMULTADD_7(a,sign,x,b) \
  (a##_7 sign x##_7*b##_7);


#define LASTMULTADD_8(a,sign,x,b) \
  (a##_8 sign x##_8*b##_8);

#define LASTMULTADD_9(a,sign,x,b) \
  (a##_9 sign x##_9*b##_9);

#define LASTMULTADD_10(a,sign,x,b) \
  (a##_10 sign x##_10*b##_10);

#define LASTMULTADD_11(a,sign,x,b) \
  (a##_11 sign x##_11*b##_11);

#define LASTMULTADD_12(a,sign,x,b) \
  (a##_12 sign x##_12*b##_12);

#define LASTMULTADD_13(a,sign,x,b) \
 (a##_13 sign x##_13*b##_13);

#define LASTMULTADD_14(a,sign,x,b) \
  (a##_14 sign x##_14*b##_14);
         
         
#define MULTADD_rev_14(ab,a,sign,x,b) \
     ab##_14=a##_13  sign  x##_13*b##_13;\
     ab##_13=a##_12  sign  x##_12*b##_12;\
     ab##_12=a##_11  sign  x##_11*b##_11;\
     ab##_11=a##_10  sign  x##_10*b##_10;\
     ab##_10=a##_9  sign  x##_9*b##_9;\
     ab##_9=a##_8  sign  x##_8*b##_8;\
     ab##_8=a##_7  sign  x##_7*b##_7;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_1=a##_0  sign  x##_0*b##_0;


#define MULTADD_rev_13(ab,a,sign,x,b) \
     ab##_13=a##_12  sign  x##_12*b##_12;\
     ab##_12=a##_11  sign  x##_11*b##_11;\
     ab##_11=a##_10  sign  x##_10*b##_10;\
     ab##_10=a##_9  sign  x##_9*b##_9;\
     ab##_9=a##_8  sign  x##_8*b##_8;\
     ab##_8=a##_7  sign  x##_7*b##_7;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_1=a##_0  sign  x##_0*b##_0;


#define MULTADD_rev_12(ab,a,sign,x,b) \
     ab##_12=a##_11  sign  x##_11*b##_11;\
     ab##_11=a##_10  sign  x##10*b##_10;\
     ab##_10=a##_9  sign  x##_9*b##_9;\
     ab##_9=a##_8  sign  x##_8*b##_8;\
     ab##_8=a##_7  sign  x##_7*b##_7;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_1=a##_0  sign  x##_0*b##_0;


#define MULTADD_rev_11(ab,a,sign,x,b) \
     ab##_11=a##_10  sign  x##_10*b##_10;\
     ab##_10=a##_9  sign  x##_9*b##_9;\
     ab##_9=a##_8  sign  x##_8*b##_8;\
     ab##_8=a##_7  sign  x##_7*b##_7;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_1=a##_0  sign  x##_0*b##_0;


#define MULTADD_rev_10(ab,a,sign,x,b) \
     ab##_10=a##_9  sign  x##_9*b##_9;\
     ab##_9=a##_8  sign  x##_8*b##_8;\
     ab##_8=a##_7  sign  x##_7*b##_7;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_1=a##_0  sign  x##_0*b##_0;


#define MULTADD_rev_9(ab,a,sign,x,b) \
     ab##_9=a##_8  sign  x##_8*b##_8;\
     ab##_8=a##_7  sign  x##_7*b##_7;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_1=a##_0  sign  x##_0*b##_0;


#define MULTADD_rev_8(ab,a,sign,x,b) \
     ab##_8=a##_7  sign  x##_7*b##_7;\
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_1=a##_0  sign  x##_0*b##_0;


#define MULTADD_rev_7(ab,a,sign,x,b) \
     ab##_7=a##_6  sign  x##_6*b##_6;\
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_1=a##_0  sign  x##_0*b##_0;


#define MULTADD_rev_6(ab,a,sign,x,b) \
     ab##_6=a##_5  sign  x##_5*b##_5;\
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_1=a##_0  sign  x##_0*b##_0;


#define MULTADD_rev_5(ab,a,sign,x,b) \
     ab##_5=a##_4  sign  x##_4*b##_4;\
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_1=a##_0  sign  x##_0*b##_0;


#define MULTADD_rev_4(ab,a,sign,x,b) \
     ab##_4=a##_3  sign  x##_3*b##_3;\
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_1=a##_0  sign  x##_0*b##_0;


#define MULTADD_rev_3(ab,a,sign,x,b) \
     ab##_3=a##_2  sign  x##_2*b##_2;\
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_1=a##_0  sign  x##_0*b##_0;


#define MULTADD_rev_2(ab,a,sign,x,b) \
     ab##_2=a##_1  sign  x##_1*b##_1;\
     ab##_1=a##_0  sign  x##_0*b##_0;


#define MULTADD_rev_1(ab,a,sign,x,b) \
     ab##_1=a##_0  sign  x##_0*b##_0;

#define MULTADD_rev_0(ab,a,sign,x,b) ;




/**************  Macros for parallell filters *************/


#define MULTADD_par_0(ab,a,sign,x,b,parallell);

#define MULTADD_par_1(ab,a,sign,x,b,parallell)\
  for(v=0;v<parallell;v++)ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];

#define MULTADD_par_2(ab,a,sign,x,b,parallell)\
  for(v=0;v<parallell;v++)ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];\
  for(v=0;v<parallell;v++)ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];



#define MULTADD_par_3(ab,a,sign,x,b,parallell)\
  for(v=0;v<parallell;v++)ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];\
  for(v=0;v<parallell;v++)ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];\
  for(v=0;v<parallell;v++)ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];
  

#define MULTADD_par_4(ab,a,sign,x,b,parallell)				\
   for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];\
    for(v=0;v<parallell;v++) ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];\
    for(v=0;v<parallell;v++) ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];\
    for(v=0;v<parallell;v++) ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];
  
#define MULTADD_par_5(ab,a,sign,x,b,parallell) \
   for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];\
    for(v=0;v<parallell;v++) ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];\
    for(v=0;v<parallell;v++) ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];\
    for(v=0;v<parallell;v++) ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];\
    for(v=0;v<parallell;v++)ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];
  
#define MULTADD_par_6(ab,a,sign,x,b,parallell)					\
   for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];\
    for(v=0;v<parallell;v++) ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];\
    for(v=0;v<parallell;v++) ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];\
    for(v=0;v<parallell;v++) ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];\
    for(v=0;v<parallell;v++)ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];\
    for(v=0;v<parallell;v++)ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];

  


#define MULTADD_par_7(ab,a,sign,x,b,parallell) \
   for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];\
    for(v=0;v<parallell;v++) ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];\
    for(v=0;v<parallell;v++) ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];\
    for(v=0;v<parallell;v++) ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];\
    for(v=0;v<parallell;v++)ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];\
    for(v=0;v<parallell;v++)ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];\
   for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];
  
#define MULTADD_par_8(ab,a,sign,x,b,parallell) \
   for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];\
    for(v=0;v<parallell;v++) ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];\
    for(v=0;v<parallell;v++) ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];\
    for(v=0;v<parallell;v++) ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];\
    for(v=0;v<parallell;v++) ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];\
    for(v=0;v<parallell;v++) ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];\
    for(v=0;v<parallell;v++) ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];	\
    for(v=0;v<parallell;v++) ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];
  



#define MULTADD_par_9(ab,a,sign,x,b,parallell) \
   for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];\
    for(v=0;v<parallell;v++) ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];\
    for(v=0;v<parallell;v++) ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];\
    for(v=0;v<parallell;v++) ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];\
    for(v=0;v<parallell;v++)ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];\
    for(v=0;v<parallell;v++)ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];\
   for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];\
    for(v=0;v<parallell;v++) ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];\
    for(v=0;v<parallell;v++) ab##_9[v]=a##_8[v]  sign  x##_8*b##_8[v];

 

#define MULTADD_par_10(ab,a,sign,x,b,parallell) \
   for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];\
    for(v=0;v<parallell;v++) ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];\
    for(v=0;v<parallell;v++) ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];\
    for(v=0;v<parallell;v++) ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];\
    for(v=0;v<parallell;v++)ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];\
    for(v=0;v<parallell;v++)ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];\
   for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];\
    for(v=0;v<parallell;v++) ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];\
    for(v=0;v<parallell;v++) ab##_9[v]=a##_8[v]  sign  x##_8*b##_8[v];\
    for(v=0;v<parallell;v++) ab##_10[v]=a##_9[v]  sign  x##_9*b##_9[v];

 


#define MULTADD_par_11(ab,a,sign,x,b,parallell) \
   for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];\
    for(v=0;v<parallell;v++) ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];\
    for(v=0;v<parallell;v++) ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];\
    for(v=0;v<parallell;v++) ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];\
    for(v=0;v<parallell;v++)ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];\
    for(v=0;v<parallell;v++)ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];\
   for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];\
    for(v=0;v<parallell;v++) ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];\
    for(v=0;v<parallell;v++) ab##_9[v]=a##_8[v]  sign  x##_8*b##_8[v];\
    for(v=0;v<parallell;v++) ab##_10[v]=a##_9[v]  sign  x##_9*b##_9[v];\
    for(v=0;v<parallell;v++)ab##_11[v]=a##_10[v]  sign  x##_10*b##_10[v];
 



#define MULTADD_par_12(ab,a,sign,x,b,parallell) \
   for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];\
    for(v=0;v<parallell;v++) ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];\
    for(v=0;v<parallell;v++) ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];\
    for(v=0;v<parallell;v++) ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];\
    for(v=0;v<parallell;v++)ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];\
    for(v=0;v<parallell;v++)ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];\
   for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];\
    for(v=0;v<parallell;v++) ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];\
    for(v=0;v<parallell;v++) ab##_9[v]=a##_8[v]  sign  x##_8*b##_8[v];\
    for(v=0;v<parallell;v++) ab##_10[v]=a##_9[v]  sign  x##_9*b##_9[v];\
    for(v=0;v<parallell;v++)ab##_11[v]=a##_10[v]  sign  x##_10*b##_10[v];\
    for(v=0;v<parallell;v++)ab##_12[v]=a##_11[v]  sign  x##_11*b##_11[v];


#define MULTADD_par_13(ab,a,sign,x,b,parallell) \
   for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];\
    for(v=0;v<parallell;v++) ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];\
    for(v=0;v<parallell;v++) ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];\
    for(v=0;v<parallell;v++) ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];\
    for(v=0;v<parallell;v++)ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];\
    for(v=0;v<parallell;v++)ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];\
   for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];\
    for(v=0;v<parallell;v++) ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];\
    for(v=0;v<parallell;v++) ab##_9[v]=a##_8[v]  sign  x##_8*b##_8[v];\
    for(v=0;v<parallell;v++) ab##_10[v]=a##_9[v]  sign  x##_9*b##_9[v];\
    for(v=0;v<parallell;v++)ab##_11[v]=a##_10[v]  sign  x##_10*b##_10[v];\
    for(v=0;v<parallell;v++)ab##_12[v]=a##_11[v]  sign  x##_11*b##_11[v];\
    for(v=0;v<parallell;v++)ab##_13[v]=a##_12[v]  sign  x##_12*b##_10[v];


#define MULTADD_par_14(ab,a,sign,x,b,parallell) \
   for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];\
    for(v=0;v<parallell;v++) ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];\
    for(v=0;v<parallell;v++) ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];\
    for(v=0;v<parallell;v++) ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];\
    for(v=0;v<parallell;v++)ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];\
    for(v=0;v<parallell;v++)ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];\
   for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];\
    for(v=0;v<parallell;v++) ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];\
    for(v=0;v<parallell;v++) ab##_9[v]=a##_8[v]  sign  x##_8*b##_8[v];\
    for(v=0;v<parallell;v++) ab##_10[v]=a##_9[v]  sign  x##_9*b##_9[v];\
    for(v=0;v<parallell;v++)ab##_11[v]=a##_10[v]  sign  x##_10*b##_10[v];\
    for(v=0;v<parallell;v++)ab##_12[v]=a##_11[v]  sign  x##_11*b##_11[v];\
    for(v=0;v<parallell;v++)ab##_13[v]=a##_12[v]  sign  x##_12*b##_10[v];\
    for(v=0;v<parallell;v++)ab##_14[v]=a##_13[v]  sign  x##_13*b##_11[v];






#define MULTADD_rev_par_0(ab,a,sign,x,b,parallell) ;	\


#define MULTADD_rev_par_1(ab,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];




#define MULTADD_rev_par_2(ab,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)  ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];


#define MULTADD_rev_par_3(ab,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)  ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];	\
  for(v=0;v<parallell;v++)  ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];


#define MULTADD_rev_par_4(ab,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)  ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];  \
  for(v=0;v<parallell;v++)  ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];	\
  for(v=0;v<parallell;v++)  ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];


#define MULTADD_rev_par_5(ab,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)  ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];	\
  for(v=0;v<parallell;v++)  ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];	\
  for(v=0;v<parallell;v++)  ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];	\
  for(v=0;v<parallell;v++)  ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];

#define MULTADD_rev_par_6(ab,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)  ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];	\
  for(v=0;v<parallell;v++)  ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];	\
  for(v=0;v<parallell;v++)  ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];	\
  for(v=0;v<parallell;v++)  ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];	\
  for(v=0;v<parallell;v++)  ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];



#define MULTADD_rev_par_7(ab,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];	\
  for(v=0;v<parallell;v++)  ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];	\
  for(v=0;v<parallell;v++)  ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];	\
  for(v=0;v<parallell;v++)  ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];	\
  for(v=0;v<parallell;v++)  ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];	\
  for(v=0;v<parallell;v++)  ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];



#define MULTADD_rev_par_8(ab,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)  ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];	\
  for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];	\
  for(v=0;v<parallell;v++)  ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];	\
  for(v=0;v<parallell;v++)  ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];	\
  for(v=0;v<parallell;v++)  ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];	\
  for(v=0;v<parallell;v++)  ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];	\
  for(v=0;v<parallell;v++)  ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];



#define MULTADD_rev_par_9(ab,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)  ab##_9[v]=a##_8[v]  sign  x##_8*b##_8[v];	\
  for(v=0;v<parallell;v++)  ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];	\
  for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];	\
  for(v=0;v<parallell;v++)  ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];	\
  for(v=0;v<parallell;v++)  ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];	\
  for(v=0;v<parallell;v++)  ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];	\
  for(v=0;v<parallell;v++)  ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];	\
  for(v=0;v<parallell;v++)  ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];



#define MULTADD_rev_par_10(ab,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)  ab##_10[v]=a##_9[v]  sign  x##_9*b##_9[v];	\
  for(v=0;v<parallell;v++)  ab##_9[v]=a##_8[v]  sign  x##_8*b##_8[v];	\
  for(v=0;v<parallell;v++)  ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];	\
  for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];	\
  for(v=0;v<parallell;v++)  ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];	\
  for(v=0;v<parallell;v++)  ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];	\
  for(v=0;v<parallell;v++)  ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];	\
  for(v=0;v<parallell;v++)  ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];	\
  for(v=0;v<parallell;v++)  ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];




#define MULTADD_rev_par_11(ab,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)  ab##_11[v]=a##_10[v]  sign  x##_10*b##_10[v];\
  for(v=0;v<parallell;v++)  ab##_10[v]=a##_9[v]  sign  x##_9*b##_9[v];	\
  for(v=0;v<parallell;v++)  ab##_9[v]=a##_8[v]  sign  x##_8*b##_8[v];	\
  for(v=0;v<parallell;v++)  ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];	\
  for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];	\
  for(v=0;v<parallell;v++)  ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];	\
  for(v=0;v<parallell;v++)  ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];	\
  for(v=0;v<parallell;v++)  ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];	\
  for(v=0;v<parallell;v++)  ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];	\
  for(v=0;v<parallell;v++)  ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];




#define MULTADD_rev_par_12(ab,a,sign,x,b,parallell)	\
 for(v=0;v<parallell;v++)  ab##_12[v]=a##_11[v]  sign  x##_11*b##_11[v]; \
  for(v=0;v<parallell;v++)  ab##_11[v]=a##_10[v]  sign  x##_10*b##_10[v];\
  for(v=0;v<parallell;v++)  ab##_10[v]=a##_9[v]  sign  x##_9*b##_9[v];	\
  for(v=0;v<parallell;v++)  ab##_9[v]=a##_8[v]  sign  x##_8*b##_8[v];	\
  for(v=0;v<parallell;v++)  ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];	\
  for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];	\
  for(v=0;v<parallell;v++)  ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];	\
  for(v=0;v<parallell;v++)  ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];	\
  for(v=0;v<parallell;v++)  ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];	\
  for(v=0;v<parallell;v++)  ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];	\
  for(v=0;v<parallell;v++)  ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];



#define MULTADD_rev_par_13(ab,a,sign,x,b,parallell)	\
 for(v=0;v<parallell;v++)  ab##_13[v]=a##_12[v]  sign  x##_12*b##_12[v]; \
 for(v=0;v<parallell;v++)  ab##_12[v]=a##_11[v]  sign  x##_11*b##_11[v]; \
  for(v=0;v<parallell;v++)  ab##_11[v]=a##_10[v]  sign  x##_10*b##_10[v];\
  for(v=0;v<parallell;v++)  ab##_10[v]=a##_9[v]  sign  x##_9*b##_9[v];	\
  for(v=0;v<parallell;v++)  ab##_9[v]=a##_8[v]  sign  x##_8*b##_8[v];	\
  for(v=0;v<parallell;v++)  ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];	\
  for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];	\
  for(v=0;v<parallell;v++)  ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];	\
  for(v=0;v<parallell;v++)  ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];	\
  for(v=0;v<parallell;v++)  ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];	\
  for(v=0;v<parallell;v++)  ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];	\
  for(v=0;v<parallell;v++)  ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];


#define MULTADD_rev_par_14(ab,a,sign,x,b,parallell)	\
 for(v=0;v<parallell;v++)  ab##_14[v]=a##_13[v]  sign  x##_13*b##_13[v];\
 for(v=0;v<parallell;v++)  ab##_13[v]=a##_12[v]  sign  x##_12*b##_12[v]; \
 for(v=0;v<parallell;v++)  ab##_12[v]=a##_11[v]  sign  x##_11*b##_11[v]; \
  for(v=0;v<parallell;v++)  ab##_11[v]=a##_10[v]  sign  x##_10*b##_10[v];\
  for(v=0;v<parallell;v++)  ab##_10[v]=a##_9[v]  sign  x##_9*b##_9[v];	\
  for(v=0;v<parallell;v++)  ab##_9[v]=a##_8[v]  sign  x##_8*b##_8[v];	\
  for(v=0;v<parallell;v++)  ab##_8[v]=a##_7[v]  sign  x##_7*b##_7[v];	\
  for(v=0;v<parallell;v++)  ab##_7[v]=a##_6[v]  sign  x##_6*b##_6[v];	\
  for(v=0;v<parallell;v++)  ab##_6[v]=a##_5[v]  sign  x##_5*b##_5[v];	\
  for(v=0;v<parallell;v++)  ab##_5[v]=a##_4[v]  sign  x##_4*b##_4[v];	\
  for(v=0;v<parallell;v++)  ab##_4[v]=a##_3[v]  sign  x##_3*b##_3[v];	\
  for(v=0;v<parallell;v++)  ab##_3[v]=a##_2[v]  sign  x##_2*b##_2[v];	\
  for(v=0;v<parallell;v++)  ab##_2[v]=a##_1[v]  sign  x##_1*b##_1[v];	\
  for(v=0;v<parallell;v++)  ab##_1[v]=a##_0[v]  sign  x##_0*b##_0[v];




#define LASTMULTADD_par_0(ptr,a,sign,x,b,parallell) \
  for(v=0;v<parallell;v++)ptr[v]=norm*( a##_0[v] sign x##_0*b##_0[v]);

#define LASTMULTADD_ACUM_par_0(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_0[v] sign x##_0*b##_0[v]);


#define LASTMULTADD_par_1(ptr,a,sign,x,b,parallell) \
  for(v=0;v<parallell;v++)ptr[v]=norm*( a##_1[v] sign x##_1*b##_1[v]);

#define LASTMULTADD_ACUM_par_1(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_1[v] sign x##_1*b##_1[v]);


#define LASTMULTADD_par_2(ptr,a,sign,x,b,parallell) \
  for(v=0;v<parallell;v++)ptr[v]=norm*( a##_2[v] sign x##_2*b##_2[v]);

#define LASTMULTADD_ACUM_par_2(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_2[v] sign x##_2*b##_2[v]);


#define LASTMULTADD_par_3(ptr,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)ptr[v]= norm*(a##_3[v] sign x##_3*b##_3[v]);

#define LASTMULTADD_ACUM_par_3(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_3[v] sign x##_3*b##_3[v]);


#define LASTMULTADD_par_4(ptr,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)ptr[v]= norm*(a##_4[v] sign x##_4*b##_4[v]);

#define LASTMULTADD_ACUM_par_4(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_4[v] sign x##_4*b##_4[v]);



#define LASTMULTADD_par_5(ptr,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)ptr[v]= norm*(a##_5[v] sign x##_5*b##_5[v]);

#define LASTMULTADD_ACUM_par_5(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_5[v] sign x##_5*b##_5[v]);



#define LASTMULTADD_par_6(ptr,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)ptr[v]= norm*(a##_6[v] sign x##_6*b##_6[v]);


#define LASTMULTADD_ACUM_par_6(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_6[v] sign x##_6*b##_6[v]);

#define LASTMULTADD_par_7(ptr,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)ptr[v]= norm*(a##_7[v] sign x##_7*b##_7[v]);

#define LASTMULTADD_ACUM_par_7(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_7[v] sign x##_7*b##_7[v]);


#define LASTMULTADD_par_8(ptr,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)ptr[v]= norm*(a##_8[v] sign x##_8*b##_8[v]);

#define LASTMULTADD_ACUM_par_8(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_8[v] sign x##_8*b##_8[v]);


#define LASTMULTADD_par_9(ptr,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)ptr[v]= norm*(a##_9[v] sign x##_9*b##_9[v]);

#define LASTMULTADD_ACUM_par_9(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_9[v] sign x##_9*b##_9[v]);


#define LASTMULTADD_par_10(ptr,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)ptr[v]= norm*(a##_10[v] sign x##_10*b##_10[v]);

#define LASTMULTADD_ACUM_par_10(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_10[v] sign x##_10*b##_10[v]);


#define LASTMULTADD_par_11(ptr,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)ptr[v]= norm*(a##_11[v] sign x##_11*b##_11[v]);

#define LASTMULTADD_ACUM_par_11(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_11[v] sign x##_11*b##_11[v]);

#define LASTMULTADD_par_12(ptr,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)ptr[v]= norm*(a##_12[v] sign x##_12*b##_12[v]);

#define LASTMULTADD_ACUM_par_12(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_12[v] sign x##_12*b##_12[v]);


#define LASTMULTADD_par_13(ptr,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)ptr[v]= norm*(a##_13[v] sign x##_13*b##_13[v]);

#define LASTMULTADD_ACUM_par_13(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_13[v] sign x##_13*b##_13[v]);



#define LASTMULTADD_par_14(ptr,a,sign,x,b,parallell)	\
  for(v=0;v<parallell;v++)ptr[v]= norm*(a##_14[v] sign x##_14*b##_14[v]);

#define LASTMULTADD_ACUM_par_14(ptr,ptrtemp,a,sign,x,b,parallell)		\
  for(v=0;v<parallell;v++)ptr[v]=ptrtemp[v]+norm*( a##_14[v] sign x##_14*b##_14[v]);






#if(FILTERLENGTH ==30)
#define DECLARE_LOOPVARIABLES_par(a)\
    REAL *a##_1,*a##_2,*a##_3,*a##_4,*a##_5,*a##_6,*a##_7,*a##_8,\
      *a##_9,*a##_10,*a##_11,*a##_12,*a##_13=0,*a##_14;

#define ALLOCATE_LOOPVARIABLES_par(a,mem)\
    a##_1=mem[1];a##_2=mem[2];a##_3=mem[3];a##_4=mem[4];a##_5=mem[5];a##_6=mem[6];a##_7=mem[7];a##_8=mem[8];a##_9=mem[9];a##_10=mem[10];a##_11=mem[11];a##_12=mem[12];a##_13=mem[13];a##_14=mem[14];

#define NULLSET_LOOPVARIABLES_par(parallell,mem) \
for(m=1;m<=14;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par LASTMULTADD_par_14
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_14
#define MULTADDLOOP_par MULTADD_par_14
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_14
#endif






#if(FILTERLENGTH ==28)
#define DECLARE_LOOPVARIABLES_par(a)\
    REAL *a##_1,*a##_2,*a##_3,*a##_4,*a##_5,*a##_6,*a##_7,*a##_8,\
      *a##_9,*a##_10,*a##_11,*a##_12,*a##_13;

#define ALLOCATE_LOOPVARIABLES_par(a,mem)\
    a##_1=mem[1];a##_2=mem[2];a##_3=mem[3];a##_4=mem[4];a##_5=mem[5];a##_6=mem[6];a##_7=mem[7];a##_8=mem[8];a##_9=mem[9];a##_10=mem[10];a##_11=mem[11];a##_12=mem[12];a##_13=mem[13];

#define NULLSET_LOOPVARIABLES_par(parallell,mem) \
for(m=1;m<=13;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par LASTMULTADD_par_13
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_13
#define MULTADDLOOP_par MULTADD_par_13
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_13
#endif






#if(FILTERLENGTH ==26)
#define DECLARE_LOOPVARIABLES_par(a)\
    REAL *a##_1,*a##_2,*a##_3,*a##_4,*a##_5,*a##_6,*a##_7,*a##_8,\
      *a##_9,*a##_10,*a##_11,*a##_12;

#define ALLOCATE_LOOPVARIABLES_par(a,mem)\
    a##_1=mem[1];a##_2=mem[2];a##_3=mem[3];a##_4=mem[4];a##_5=mem[5];a##_6=mem[6];a##_7=mem[7];a##_8=mem[8];a##_9=mem[9];a##_10=mem[10];a##_11=mem[11];a##_12=mem[12];

#define NULLSET_LOOPVARIABLES_par(parallell,mem) \
for(m=1;m<=12;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par LASTMULTADD_par_12
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_12
#define MULTADDLOOP_par MULTADD_par_12
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_12
#endif




#if(FILTERLENGTH ==24)
#define DECLARE_LOOPVARIABLES_par(a)\
    REAL *a##_1,*a##_2,*a##_3,*a##_4,*a##_5,*a##_6,*a##_7,*a##_8,*a##_9,*a##_10,*a##_11;

#define ALLOCATE_LOOPVARIABLES_par(a,mem)\
    a##_1=mem[1];a##_2=mem[2];a##_3=mem[3];a##_4=mem[4];a##_5=mem[5];a##_6=mem[6];a##_7=mem[7];a##_8=mem[8],a##_9=mem[9];a##_10=mem[10];a##_11=mem[11];

#define NULLSET_LOOPVARIABLES_par(parallell,mem) \
for(m=1;m<=11;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par LASTMULTADD_par_11
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_11
#define MULTADDLOOP_par  MULTADD_par_11
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_11
#endif

#if(FILTERLENGTH ==22)
#define DECLARE_LOOPVARIABLES_par(a)\
    REAL *a##_1,*a##_2,*a##_3,*a##_4,*a##_5,*a##_6,*a##_7,*a##_8,*a##_9,*a##_10;

#define ALLOCATE_LOOPVARIABLES_par(a,mem)\
    a##_1=mem[1];a##_2=mem[2];a##_3=mem[3];a##_4=mem[4];a##_5=mem[5];a##_6=mem[6];a##_7=mem[7];a##_8=mem[8],a##_9=mem[9];a##_10=mem[10];

#define NULLSET_LOOPVARIABLES_par(parallell,mem) \
for(m=1;m<=10;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par LASTMULTADD_par_10
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_10
#define MULTADDLOOP_par  MULTADD_par_10
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_10
#endif

#if(FILTERLENGTH ==20)
#define DECLARE_LOOPVARIABLES_par(a)\
    REAL *a##_1,*a##_2,*a##_3,*a##_4,*a##_5,*a##_6,*a##_7,*a##_8,*a##_9;

#define ALLOCATE_LOOPVARIABLES_par(a,mem)\
    a##_1=mem[1];a##_2=mem[2];a##_3=mem[3];a##_4=mem[4];a##_5=mem[5];a##_6=mem[6];a##_7=mem[7];a##_8=mem[8];

#define NULLSET_LOOPVARIABLES_par(parallell,mem) \
  for(m=1;m<=9;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par LASTMULTADD_par_9
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_9
#define MULTADDLOOP_par MULTADD_par_9
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_9
#endif

#if(FILTERLENGTH ==18)
#define DECLARE_LOOPVARIABLES_par(a)\
    REAL *a##_1,*a##_2,*a##_3,*a##_4,*a##_5=0,*a##_6,*a##_7,*a##_8;

#define ALLOCATE_LOOPVARIABLES_par(a,mem)\
    a##_1=mem[1];a##_2=mem[2];a##_3=mem[3];a##_4=mem[4];a##_5=mem[5];a##_6=mem[6];a##_7=mem[7];a##_8=mem[8];

#define NULLSET_LOOPVARIABLES_par(parallell,mem) \
for(m=1;m<=8;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par  LASTMULTADD_par_8



#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_8
#define MULTADDLOOP_par MULTADD_par_8
#define MULTADD_reverse_LOOP_par  MULTADD_rev_par_8
#endif


#if(FILTERLENGTH ==16)
#define DECLARE_LOOPVARIABLES_par(a)\
    REAL *a##_1,*a##_2,*a##_3,*a##_4,*a##_5=0,*a##_6,*a##_7;

#define ALLOCATE_LOOPVARIABLES_par(a,mem)\
    a##_1=mem[1];a##_2=mem[2];a##_3=mem[3];a##_4=mem[4];a##_5=mem[5];a##_6=mem[6];a##_7=mem[7];

#define NULLSET_LOOPVARIABLES_par(parallell,mem) \
for(m=1;m<=7;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par LASTMULTADD_par_7
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_7
#define MULTADDLOOP_par MULTADD_par_7
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_7
#endif



#if(FILTERLENGTH ==14)
#define DECLARE_LOOPVARIABLES_par(a)\
    REAL *a##_1,*a##_2,*a##_3,*a##_4,*a##_5=0,*a##_6;

#define ALLOCATE_LOOPVARIABLES_par(a,mem)\
    a##_1=mem[1];a##_2=mem[2];a##_3=mem[3];a##_4=mem[4];a##_5=mem[5];a##_6=mem[6];

#define NULLSET_LOOPVARIABLES_par(parallell,mem) \
for(m=1;m<=6;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par LASTMULTADD_par_6
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_6
#define MULTADDLOOP_par MULTADD_par_6
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_6
#endif




#if(FILTERLENGTH ==12)
#define DECLARE_LOOPVARIABLES_par(a)\
    REAL *a##_1,*a##_2,*a##_3,*a##_4,*a##_5;

#define ALLOCATE_LOOPVARIABLES_par(a,mem)\
    a##_1=mem[1];a##_2=mem[2];a##_3=mem[3];a##_4=mem[4];a##_5=mem[5];

#define NULLSET_LOOPVARIABLES_par(parallell,mem) \
for(m=1;m<=5;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par LASTMULTADD_par_5
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_5
#define MULTADDLOOP_par MULTADD_par_5
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_5
#endif



#if(FILTERLENGTH ==10)
#define DECLARE_LOOPVARIABLES_par(a)\
    REAL *a##_1,*a##_2,*a##_3,*a##_4;

#define ALLOCATE_LOOPVARIABLES_par(a,mem)\
    a##_1=mem[1];a##_2=mem[2];a##_3=mem[3];a##_4=mem[4];

#define NULLSET_LOOPVARIABLES_par(parallell,mem) \
for(m=1;m<=4;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par LASTMULTADD_par_4
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_4
#define MULTADDLOOP_par MULTADD_par_4
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_4
#endif




#if(FILTERLENGTH == 8)
#define DECLARE_LOOPVARIABLES_par(a)\
    REAL *a##_1,*a##_2,*a##_3;

#define ALLOCATE_LOOPVARIABLES_par(a,mem)\
    a##_1=mem[1];a##_2=mem[2];a##_3=mem[3];

#define NULLSET_LOOPVARIABLES_par(parallell,mem)\
for(m=1;m<=3;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par LASTMULTADD_par_3
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_3
#define MULTADDLOOP_par MULTADD_par_3
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_3
#endif



#if(FILTERLENGTH == 6)
#define DECLARE_LOOPVARIABLES_par(a) REAL *a##_1,*a##_2;

#define ALLOCATE_LOOPVARIABLES_par(a,mem) a##_1=mem[1];a##_2=mem[2];

#define NULLSET_LOOPVARIABLES_par(parallell,mem)\
for(m=1;m<=2;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par LASTMULTADD_par_2
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_2
#define MULTADDLOOP_par MULTADD_par_2
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_2
#endif

#if(FILTERLENGTH == 4)
#define DECLARE_LOOPVARIABLES_par(a) REAL *a##_1;

#define ALLOCATE_LOOPVARIABLES_par(a,mem) a##_1=mem[1];;

#define NULLSET_LOOPVARIABLES_par(parallell,mem)\
for(m=1;m<=1;m++)memset(mem[m],0,parallell*sizeof(REAL));


#define LASTMULTADD_par LASTMULTADD_par_1
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_1
#define MULTADDLOOP_par MULTADD_par_1
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_1
#endif



#if(FILTERLENGTH == 2)
#define DECLARE_LOOPVARIABLES_par(a);

#define ALLOCATE_LOOPVARIABLES_par(a,mem);	


#define NULLSET_LOOPVARIABLES_par(parallell,mem);


#define LASTMULTADD_par LASTMULTADD_par_0
#define LASTMULTADD_ACUM_par LASTMULTADD_ACUM_par_0
#define MULTADDLOOP_par MULTADD_par_0
#define MULTADD_reverse_LOOP_par MULTADD_rev_par_0
#endif
