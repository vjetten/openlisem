/*---------------------------------------------------------------------------
project: openLISEM
name: lisModel.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website SVN: http://sourceforge.net/projects/lisem

---------------------------------------------------------------------------*/
/**************************************************************************/
/*  lookup.c                                                              */ 
/*   Computes Theta from head, K from head or Diff Moist Cap from head    */
/*   currently from linear interpolation in tables, but Van Genuchten     */
/*   possible                                                             */
/*                                                                        */
/**************************************************************************/

//#include "swatremisc.h"
#include "swatre_p.h"
#include "swatresoillut.h"
#include "swatrelookup.h"

//-----------------------------------------------------------------------------------
/* theta from head */
double TheNode(
		double head,           /* current head value of this node */
		const  HORIZON *hor)   /* parameters of horizon this node belongs to */
{
	head = max( head, -1e-10);
	if (head >= -1.0E-2)
		return LUT_Highest(hor->lut, THETA_COL);
	/* SWITCH BETWEEN  VAN_GENUG AND LUT HERE */
	return LUT_LinIntPol(hor->lut, THETA_COL, head, H_COL);
}
//-----------------------------------------------------------------------------------
/*  conductivity from head */
double HcoNode(
		double head,
		const HORIZON *hor,
		double calib,
		double SEC)
{
	if (head >= -1.0E-2)
		return (LUT_Highest(hor->lut, K_COL)*calib/SEC);

	/* SWITCH BETWEEN  VAN_GENUG AND LUT HERE */
	return (LUT_LinIntPol(hor->lut, K_COL, head, H_COL)/SEC);
}
//-----------------------------------------------------------------------------------
/* Differential Moisture Capacity from head */
double DmcNode(
		double head,           /* current head value of this node           */
		const  HORIZON *hor)   /* parameters of horizon this node belongs to */
{
	int i;         /* index in LUT where dmch[i] <= head <= dmch[i+1] */
	const LUT *l;  /* lut of this horizon */

	/* dit gaat niet goed als profiel van verzadigd naar onversazigd swithched:
	if (head >= 0) return 0; */

	if (head >= -1.0E-2)
		return LUT_LinIntPol(hor->lut, DMCC_COL, head, DMCH_COL);

	/* SWITCH BETWEEN  VAN_GENUG AND LUT HERE */

	l = hor->lut;
	i = LUT_Index_LE(l, head, DMCH_COL);
	i = min(LUT_nrRows(l)-2, i);
	i = max(i, 0);

	return LUT_ValueAt(l, DMCC_COL, i) +
			(head - LUT_ValueAt(l, DMCH_COL, i)) *
			(LUT_ValueAt(l,DMCC_COL, i+1)-LUT_ValueAt(l,DMCC_COL, i))/
			(LUT_ValueAt(l,DMCH_COL, i+1)-LUT_ValueAt(l,DMCH_COL, i));
}
//-----------------------------------------------------------------------------------
/*deal with evap sinkterm
            //Sink[i] = SinkNode(h[i], Horizon(p, i), SinkTerm);
double SinkNode(
	double head,
	const HORIZON *hor,
  double SinkTerm)
{
  double sink = 0;
  
  if (head >= -0.1)
     return (sink);  //no Sinkterm when soil nearly saturated

 RootFraction = 
  return (1/(1+(head/-500)**3) * SinkTerm * RootFraction);

	return LUT_LinIntPol(hor->lut, K_COL, head, H_COL);
}
*/
//-----------------------------------------------------------------------------------
