#ifndef __CPSCALC_H
#define __CPSCALC_H


#ifdef __cplusplus
 extern "C" {
#endif 

int  CpsCalc(char  *expr,char  *fileName,long  lineNr);
int  MultiCpsCalc(char  *expr, int n, char  *fileName,long  lineNr);

#define calc(expr)	CpsCalc(expr, __FILE__, __LINE__)
#define mcalc(expr, n)	MultiCpsCalc(expr, n, __FILE__, __LINE__)

#ifdef __cplusplus
 }
#endif 

#endif /*__CPSCALC_H */

