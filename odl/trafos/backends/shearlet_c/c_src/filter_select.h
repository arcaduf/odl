/*    Rotation-values for some DAUBECHIES-filter:
      Calculation done by Eirik Fossgaard.
*/

#ifdef Haar_2

#ifdef SKIP_L
#define FILTER_PROCEEDURE           haar_2_skip_L
#define FILTER_PROCEEDURE_parallell haar_2_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invvhaar_2_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           haar_2_skip_H
#define FILTER_PROCEEDURE_parallell haar_2_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invhaar_2_skip_H
#else
#define FILTER_PROCEEDURE           haar_2
#define FILTER_PROCEEDURE_parallell haar_2_par
#define INVERSE_FILTER_PROCEEDURE   invhaar_2
#endif
#define wavelet  haar_2
#endif

#define FILTERLENGTH 2

#define x_0 -1.0

#endif



#ifdef Daub_4          

#ifdef SKIP_L
#define FILTER_PROCEEDURE           daub_4_skip_L
#define FILTER_PROCEEDURE_parallell daub_4_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invdaub_4_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           daub_4_skip_H
#define FILTER_PROCEEDURE_parallell daub_4_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invdaub_4_skip_H
#else
#define FILTER_PROCEEDURE           daub_4
#define FILTER_PROCEEDURE_parallell daub_4_par
#define INVERSE_FILTER_PROCEEDURE   invdaub_4
#endif
#define wavelet  daub_4
#endif




#define FILTERLENGTH 4

#define x_0   .5773502691
#define x_1  -.2679491923


#endif

#ifdef Daub_6

#ifdef SKIP_L
#define FILTER_PROCEEDURE           daub_6_skip_L
#define FILTER_PROCEEDURE_parallell daub_6_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invdaub_6_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           daub_6_skip_H
#define FILTER_PROCEEDURE_parallell daub_6_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invdaub_6_skip_H
#else
#define FILTER_PROCEEDURE           daub_6
#define FILTER_PROCEEDURE_parallell daub_6_par
#define INVERSE_FILTER_PROCEEDURE   invdaub_6
#endif
#define wavelet  daub_6
#endif






#define FILTERLENGTH  6

#define x_0   .4122865951
#define x_1  1.831178514
#define x_2  -.1058894200

#endif

#ifdef Daub_8


#ifdef SKIP_L
#define FILTER_PROCEEDURE           daub_8_skip_L
#define FILTER_PROCEEDURE_parallell daub_8_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invdaub_8_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           daub_8_skip_H
#define FILTER_PROCEEDURE_parallell daub_8_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invdaub_8_skip_H
#else
#define FILTER_PROCEEDURE           daub_8
#define FILTER_PROCEEDURE_parallell daub_8_par
#define INVERSE_FILTER_PROCEEDURE   invdaub_8
#endif
#define wavelet  daub_8
#endif




#define FILTERLENGTH 8


#define x_0   .3222758836
#define x_1  1.233150027
#define x_2  3.856627874
#define x_3  -.04600009616

#endif

#ifdef Daub_10



#ifdef SKIP_L
#define FILTER_PROCEEDURE           daub_10_skip_L
#define FILTER_PROCEEDURE_parallell daub_10_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invdaub_10_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           daub_10_skip_H
#define FILTER_PROCEEDURE_parallell daub_10_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invdaub_10_skip_H
#else
#define FILTER_PROCEEDURE           daub_10
#define FILTER_PROCEEDURE_parallell daub_10_par
#define INVERSE_FILTER_PROCEEDURE   invdaub_10
#endif
#define wavelet  daub_10
#endif



#define FILTERLENGTH  10

#define x_0   .2651451339
#define x_1   .9398995872
#define x_2  2.353886784
#define x_3  7.508378888
#define x_4  -.02083494630

#endif

#ifdef Daub_12

#ifdef SKIP_L
#define FILTER_PROCEEDURE           daub_12_skip_L
#define FILTER_PROCEEDURE_parallell daub_12_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invdaub_12_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           daub_12_skip_H
#define FILTER_PROCEEDURE_parallell daub_12_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invdaub_12_skip_H
#else
#define FILTER_PROCEEDURE           daub_12
#define FILTER_PROCEEDURE_parallell daub_12_par
#define INVERSE_FILTER_PROCEEDURE   invdaub_12
#endif
#define wavelet  daub_12
#endif



#define FILTERLENGTH 12

#define x_0    .2255061720
#define x_1    .7643296306
#define x_2   1.696013010
#define x_3   4.114979257
#define x_4  14.28573961
#define x_5   -.009658362993

#endif

#ifdef Daub_14

#ifdef SKIP_L
#define FILTER_PROCEEDURE           daub_14_skip_L
#define FILTER_PROCEEDURE_parallell daub_14_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invdaub_14_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           daub_14_skip_H
#define FILTER_PROCEEDURE_parallell daub_14_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invdaub_14_skip_H
#else
#define FILTER_PROCEEDURE           daub_14
#define FILTER_PROCEEDURE_parallell daub_14_par
#define INVERSE_FILTER_PROCEEDURE   invdaub_14
#endif
#define wavelet  daub_14
#endif



#define FILTERLENGTH  14


#define x_0    .1963287126
#define x_1    .6466065217
#define x_2   1.333037518
#define x_3   2.764759661
#define x_4   7.035232916
#define x_5  27.00281769
#define x_6   -.004543409641

#endif

#ifdef Daub_16


#ifdef SKIP_L
#define FILTER_PROCEEDURE           daub_16_skip_L
#define FILTER_PROCEEDURE_parallell daub_16_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invdaub_16_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           daub_16_skip_H
#define FILTER_PROCEEDURE_parallell daub_16_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invdaub_16_skip_H
#else
#define FILTER_PROCEEDURE           daub_16
#define FILTER_PROCEEDURE_parallell daub_16_par
#define INVERSE_FILTER_PROCEEDURE   invdaub_16
#endif
#define wavelet  daub_16
#endif

#define FILTERLENGTH 16

#define x_0    .1739238836
#define x_1    .5617332940
#define x_2   1.103629937
#define x_3   2.074598026
#define x_4   4.380557848
#define x_5  12.05139151
#define x_6  49.52666172
#define x_7   -.002443028170

#endif

#ifdef Daub_18



#ifdef SKIP_L
#define FILTER_PROCEEDURE           daub_18_skip_L
#define FILTER_PROCEEDURE_parallell daub_18_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invdaub_18_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           daub_18_skip_H
#define FILTER_PROCEEDURE_parallell daub_18_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invdaub_18_skip_H
#else
#define FILTER_PROCEEDURE           daub_18
#define FILTER_PROCEEDURE_parallell daub_18_par
#define INVERSE_FILTER_PROCEEDURE   invdaub_18
#endif
#define wavelet  daub_18
#endif


#define FILTERLENGTH 18

#define x_0    .1561629731
#define x_1    .4973943657
#define x_2    .9452416623
#define x_3   1.664172294
#define x_4   3.11401686
#define x_5   6.915226655
#define x_6  20.60043019
#define x_7  96.49772819
#define x_8   -.001033336055

#endif

#ifdef Daub_20

#ifdef SKIP_L
#define FILTER_PROCEEDURE           daub_20_skip_L
#define FILTER_PROCEEDURE_parallell daub_20_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invdaub_20_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           daub_20_skip_H
#define FILTER_PROCEEDURE_parallell daub_20_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invdaub_20_skip_H
#else
#define FILTER_PROCEEDURE           daub_20
#define FILTER_PROCEEDURE_parallell daub_20_par
#define INVERSE_FILTER_PROCEEDURE   invdaub_20
#endif
#define wavelet  daub_20
#endif


#define FILTERLENGTH 20

#define x_0      .1417287200
#define x_1      .4467987788
#define x_2      .8289658876
#define x_3     1.394189716
#define x_4     2.402966640
#define x_5     4.635603726
#define x_6    10.98508401
#define x_7    35.63003753
#define x_8   183.0054911
#define x_9    -.0004973444230

#endif

/*Rotation-values for some COIFLET-filter: */

#ifdef Coif_6

#ifdef SKIP_L
#define FILTER_PROCEEDURE           coif_6_skip_L
#define FILTER_PROCEEDURE_parallell coif_6_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invcoif_6_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           coif_6_skip_H
#define FILTER_PROCEEDURE_parallell coif_6_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invcoif_6_skip_H
#else
#define FILTER_PROCEEDURE           coif_6
#define FILTER_PROCEEDURE_parallell coif_6_par
#define INVERSE_FILTER_PROCEEDURE   invcoif_6
#endif
#define wavelet  coif_6
#endif



#define FILTERLENGTH 6

#define x_0  -.2152504427
#define x_1   .3779644639
#define x_2  -.2152504427

#endif

/*Coif_8 computed by use of Maple:  */


#ifdef Coif_8


#ifdef SKIP_L
#define FILTER_PROCEEDURE           coif_8_skip_L
#define FILTER_PROCEEDURE_parallell coif_8_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invcoif_8_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           coif_8_skip_H
#define FILTER_PROCEEDURE_parallell coif_8_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invcoif_8_skip_H
#else
#define FILTER_PROCEEDURE           coif_8
#define FILTER_PROCEEDURE_parallell coif_8_par
#define INVERSE_FILTER_PROCEEDURE   invcoif_8
#endif
#define wavelet  coif_8
#endif

#define FILTERLENGTH 8

#define x_0  2.215250436
#define x_1   .1504720765
#define x_2  -.7384168123
#define x_3  -.451416230

#endif

#ifdef Coif_12


#ifdef SKIP_L
#define FILTER_PROCEEDURE           coif_12_skip_L
#define FILTER_PROCEEDURE_parallell coif_12_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invcoif_12_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           coif_12_skip_H
#define FILTER_PROCEEDURE_parallell coif_12_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invcoif_12_skip_H
#else
#define FILTER_PROCEEDURE           coif_12
#define FILTER_PROCEEDURE_parallell coif_12_par
#define INVERSE_FILTER_PROCEEDURE   invcoif_12
#endif
#define wavelet  coif_12
#endif


#define FILTERLENGTH 12


#define x_0  -.3952094767
#define x_1  -.5625481503
#define x_2   .1165449040
#define x_3  1.317233974
#define x_4  6.198029576
#define x_5  -.04396989341

#endif

#ifdef Coif_18



#ifdef SKIP_L
#define FILTER_PROCEEDURE           coif_18_skip_L
#define FILTER_PROCEEDURE_parallell coif_18_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invcoif_18_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           coif_18_skip_H
#define FILTER_PROCEEDURE_parallell coif_18_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invcoif_18_skip_H
#else
#define FILTER_PROCEEDURE           coif_18
#define FILTER_PROCEEDURE_parallell coif_18_par
#define INVERSE_FILTER_PROCEEDURE   invcoif_18
#endif
#define wavelet  coif_18
#endif

#define FILTERLENGTH 18

#define x_0    -.4874353702
#define x_1  - 1.119071133
#define x_2    -.2570708497
#define x_3     .1290348165
#define x_4     .4411074710
#define x_5    2.215422179
#define x_6    8.338120664
#define x_7  15.03636438
#define x_8    -.009120773147 

#endif

#ifdef Coif_24

#ifdef SKIP_L
#define FILTER_PROCEEDURE           coif_24_skip_L
#define FILTER_PROCEEDURE_parallell coif_24_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invcoif_24_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           coif_24_skip_H
#define FILTER_PROCEEDURE_parallell coif_24_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invcoif_24_skip_H
#else
#define FILTER_PROCEEDURE           coif_24
#define FILTER_PROCEEDURE_parallell coif_24_par
#define INVERSE_FILTER_PROCEEDURE   invcoif_24
#endif
#define wavelet  coif_24
#endif

#define FILTERLENGTH  24

#define x_0  -.5476023581
#define x_1 -1.457533881
#define x_2  -.7720754411
#define x_3  -.1309276144
#define x_4   .1710353887
#define x_5   .2957793746
#define x_6   .8070747686
#define x_7  3.126528296
#define x_8 11.27596534
#define x_9 12.66598170
#define x_10 53.96686137 
#define x_11  - .002000409650

#endif


 
#ifdef Coif_30

#ifdef SKIP_L
#define FILTER_PROCEEDURE           coif_30_skip_L
#define FILTER_PROCEEDURE_parallell coif_30_par_skip_L
#define INVERSE_FILTER_PROCEEDURE   invcoif_30_skip_L
#else
#ifdef SKIP_H
#define FILTER_PROCEEDURE           coif_30_skip_H
#define FILTER_PROCEEDURE_parallell coif_30_par_skip_H
#define INVERSE_FILTER_PROCEEDURE   invcoif_30_skip_H
#else
#define FILTER_PROCEEDURE           coif_30
#define FILTER_PROCEEDURE_parallell coif_30_par
#define INVERSE_FILTER_PROCEEDURE   invcoif_30
#endif
#define wavelet  coif_30
#endif


#define FILTERLENGTH 30

#define x_0   -.5914303923
#define x_1  -1.718001035
#define x_2  -1.195010469
#define x_3   -.4056552189
#define x_4   -.1316532923
#define x_5    .1205373016
#define x_6   .3671126852
#define x_7    .4678947012
#define x_8   1.165968370
#define x_9   4.100416655
#define x_10  15.61099604
#define x_11 11.59905847
#define x_12 37.56973541
#define x_13 197.1316159
#define x_14  -.0004543371650

#endif


