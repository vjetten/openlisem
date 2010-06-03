int  CpsCalc(char  *expr,char  *fileName,long  lineNr);

#define calc(expr)		CpsCalc(expr, __FILE__, __LINE__)
