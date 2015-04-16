/*-------------------------------------------------------------------------
  Scientific Computation Research Center, RPI, Troy NY
  (C) Copyright 1995, RPI-SCOREC
 
  Project  : shapeFuntions
  Author(s): Saikat Dey
  Creation : Oct., 95
  Modifi.  : 
  Function :
             mode shape for a simplex edge of order ip, Maple generated.
-------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

double En(int ip, double r, double s) {
   double f = 0.0;
   double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18;
   double t19,t20,t21,t22,t23,t24,t25,t26,t27,t28;

   /* p=2 */
   if( ip==0 ) {
      f = 1.0;
   /* p=3 */
   } else if( ip==1 ) {
      f = s-r;
   /* p=4 */
   } else if( ip==2 ) {
      t1 = s*s;
      t3 = r*r;
      f = t1-3.0*r*s+t3;
   /* p=5 */
   } else if( ip==3 ) {
      t1 = s*s;
      t4 = r*r;
      f = t1*s-6.0*r*t1+6.0*t4*s-t4*r;
   /* p=6 */
   } else if( ip==4 ) {
      t1 = s*s;
      t2 = t1*t1;
      t5 = r*r;
      t9 = t5*t5;
      f = t2-10.0*r*t1*s+20.0*t5*t1-10.0*t5*r*s+t9;
   /* p=7 */
   } else if( ip==5 ) {
      t1 = s*s;
      t2 = t1*t1;
      t5 = r*r;
      t10 = t5*t5;
      f = t2*s-15.0*t2*r+50.0*t5*t1*s-50.0*t5*r*t1+15.0*t10*s-t10*r;
   /* p=8 */
   } else if( ip==6 ) {
      t1 = s*s;
      t2 = t1*t1;
      t6 = r*r;
      t11 = t6*t6;
      f = t2*t1-21.0*r*t2*s+105.0*t6*t2-175.0*t6*r*t1*s+105.0*t11*t1-21.0*t11
*r*s+t11*t6;
   /* p=9 */
   } else if( ip==7 ) {
      t1 = s*s;
      t2 = t1*s;
      t3 = t1*t1;
      t7 = r*r;
      t10 = t7*r;
      t12 = t7*t7;
      f = t3*t2-28.0*r*t3*t1+196.0*t7*t3*s-490.0*t10*t3+490.0*t12*t2-196.0*
t12*r*t1+28.0*t12*t7*s-t12*t10;
   /* p=10 */
   } else if( ip==8 ) {
      t1 = s*s;
      t2 = t1*t1;
      t3 = t2*t2;
      t4 = t1*s;
      t7 = r*r;
      t10 = t7*r;
      t13 = t7*t7;
      t21 = t13*t13;
      f = t3-36.0*r*t2*t4+336.0*t7*t2*t1-1176.0*t10*t2*s+1764.0*t13*t2-1176.0
*t13*r*t4+336.0*t13*t7*t1-36.0*t13*t10*s+t21;
   /* p=11 */
   } else if( ip==9 ) {
      t1 = s*s;
      t2 = t1*t1;
      t3 = t2*t2;
      t6 = r*r;
      t7 = t1*s;
      t10 = t6*r;
      t13 = t6*t6;
      t22 = t13*t13;
      f = t3*s-45.0*t3*r+540.0*t6*t2*t7-2520.0*t10*t2*t1+5292.0*t13*t2*s
-5292.0*t13*r*t2+2520.0*t13*t6*t7-540.0*t13*t10*t1+45.0*t22*s-t22*r;
   /* p=12 */
   } else if( ip==10 ) {
      t1 = s*s;
      t2 = t1*t1;
      t3 = t2*t2;
      t7 = r*r;
      t9 = t7*r;
      t10 = t1*s;
      t13 = t7*t7;
      t23 = t13*t13;
      f = t3*t1-55.0*r*t3*s+825.0*t7*t3-4950.0*t9*t2*t10+13860.0*t13*t2*t1
-19404.0*t13*r*t2*s+13860.0*t13*t7*t2-4950.0*t13*t9*t10+825.0*t23*t1-55.0*t23*r
*s+t23*t7;
   /* p=13 */
   } else if( ip==11 ) {
      t1 = s*s;
      t2 = t1*s;
      t3 = t1*t1;
      t4 = t3*t3;
      t8 = r*r;
      t11 = t8*r;
      t13 = t8*t8;
      t24 = t13*t13;
      f = t2*t4-66.0*r*t4*t1+1210.0*t8*t4*s-9075.0*t11*t4+32670.0*t13*t3*t2
-60984.0*t13*r*t3*t1+60984.0*t13*t8*t3*s-32670.0*t13*t11*t3+9075.0*t24*t2
-1210.0*t24*r*t1+66.0*t24*t8*s-t24*t11;
   /* p=14 */
   } else if( ip==12 ) {
      t1 = s*s;
      t2 = t1*t1;
      t3 = t2*t2;
      t5 = t1*s;
      t8 = r*r;
      t11 = t8*r;
      t14 = t8*t8;
      t25 = t14*t14;
      f = t3*t2-78.0*r*t3*t5+1716.0*t8*t3*t1-15730.0*t11*t3*s+70785.0*t14*t3
-169884.0*t14*r*t2*t5+226512.0*t14*t8*t2*t1-169884.0*t14*t11*t2*s+70785.0*t25*
t2-15730.0*t25*r*t5+1716.0*t25*t8*t1-78.0*t25*t11*s+t25*t14;
   /* p=15 */
   } else if( ip==13 ) {
      t1 = s*s;
      t2 = t1*t1;
      t3 = t2*s;
      t4 = t2*t2;
      t8 = r*r;
      t9 = t1*s;
      t12 = t8*r;
      t15 = t8*t8;
      t18 = t15*r;
      t26 = t15*t15;
      f = t4*t3-91.0*r*t2*t4+2366.0*t8*t4*t9-26026.0*t12*t4*t1+143143.0*t15*
t4*s-429429.0*t18*t4+736164.0*t15*t8*t2*t9-736164.0*t15*t12*t2*t1+429429.0*t26*
t3-143143.0*t26*r*t2+26026.0*t26*t8*t9-2366.0*t26*t12*t1+91.0*t26*t15*s-t26*t18
;
    }
    return f ;
}

#ifdef __cplusplus
}
#endif
