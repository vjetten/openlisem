
/*!
  \file lookup.cpp
  \brief SWATRE: computes Theta from head, K from head or Diff Moist Cap from head.

  functions:\n
- double HNode(double theta, const  HORIZON *hor)\n
- double TheNode(double head, const  HORIZON *hor)\n
- double HcoNode(double head, const HORIZON *hor, double calib, double SEC)\n
- double DmcNode(double head, const  HORIZON *hor) \n

*/

//#include <algorithm>
#include "model.h"

//-----------------------------------------------------------------------------------
double TWorld::FindNode(double head, const  HORIZON *hor, int column)
{
    LUT *l = hor->lut;

    if (head >= 0) {
        return l->hydro[column].last();
    }

    auto it = std::lower_bound(l->hydro[H_COL].begin(), l->hydro[H_COL].end(), head);

    if (it == l->hydro[H_COL].begin()) {
        return(l->hydro[column][0]);
    } else if (it == l->hydro[H_COL].end()) {
        return(l->hydro[column].last());
    } else {
        int lowerIndex = std::distance(l->hydro[H_COL].begin(), it - 1); // Index of lower bound
        int upperIndex = std::distance(l->hydro[H_COL].begin(), it);     // Index of upper bound
        double lV = *(it - 1);
        double uV = *it;

        if (uV == lV) {
            return l->hydro[H_COL][lowerIndex]; // or some default value
        }

        double lTh = l->hydro[column][lowerIndex];
        double uTh = l->hydro[column][upperIndex];

        return lTh + (head - lV) * (uTh - lTh) / (uV-lV);
    }
}
//-----------------------------------------------------------------------------------
/// head from theta
double TWorld::HNode(
        double theta,           // current theta value of this node
        const  HORIZON *hor)    // parameters of horizon this node belongs to
{
    LUT *l = hor->lut;

    auto it = std::lower_bound(l->hydro[THETA_COL].begin(), l->hydro[THETA_COL].end(), theta);
    if (it == l->hydro[THETA_COL].begin()) {
        return(l->hydro[H_COL][0]);
    } else if (it == l->hydro[THETA_COL].end()) {
        return(l->hydro[H_COL].last());
    } else {
        double lH = *(it - 1);
        double uH = *it;
        double f = (theta-lH)/(uH-lH);

        int lowerIndex = std::distance(l->hydro[THETA_COL].begin(), it - 1); // Index of lower bound
        int upperIndex = std::distance(l->hydro[THETA_COL].begin(), it);     // Index of upper bound

        double lTh = l->hydro[H_COL][lowerIndex];
        double uTh = l->hydro[H_COL][upperIndex];
        return (lTh + f*(uTh-lTh));
    }
}
//-----------------------------------------------------------------------------------
/// theta from head
//OBSOLETE
double TWorld::TheNode(
        double head,           // current head value of this node
        const  HORIZON *hor)   // parameters of horizon this node belongs to
{
    LUT *l = hor->lut;
    if (head >= 0) {
        return l->hydro[THETA_COL].last();
    }
    auto it = std::lower_bound(l->hydro[H_COL].begin(), l->hydro[H_COL].end(), head);
    if (it == l->hydro[H_COL].begin()) {
        return(l->hydro[THETA_COL][0]);
    } else if (it == l->hydro[H_COL].end()) {
        return(l->hydro[THETA_COL].last());
    } else {
        double lH = *(it - 1);
        double uH = *it;
        double f = (head-lH)/(uH-lH);

        int lowerIndex = std::distance(l->hydro[H_COL].begin(), it - 1); // Index of lower bound
        int upperIndex = std::distance(l->hydro[H_COL].begin(), it);     // Index of upper bound

        double lTh = l->hydro[THETA_COL][lowerIndex];
        double uTh = l->hydro[THETA_COL][upperIndex];
        return (lTh + f*(uTh-lTh));
    }
}
//-----------------------------------------------------------------------------------
/// hydraulic conductivity from head
/// //OBSOLETE
double TWorld::HcoNode(double head,const HORIZON *hor)
{
    LUT *l = hor->lut;

    if (head >= 0)
        return l->hydro[K_COL].last();

    auto it = std::lower_bound(l->hydro[H_COL].begin(), l->hydro[H_COL].end(), head);
    if (it == l->hydro[H_COL].begin()) {
        return(l->hydro[K_COL][0]);
    } else if (it == l->hydro[H_COL].end()) {
        return(l->hydro[K_COL].last());
    } else {
        double lH = *(it - 1);
        double uH = *it;
        double f = (head-lH)/(uH-lH);
        int lowerIndex = std::distance(l->hydro[H_COL].begin(), it - 1); // Index of lower bound
        int upperIndex = std::distance(l->hydro[H_COL].begin(), it);     // Index of upper bound

        double lK = l->hydro[K_COL][lowerIndex];
        double uK = l->hydro[K_COL][upperIndex];
        return (lK+f*(uK-lK));
    }
}
//-----------------------------------------------------------------------------------
//OBSOLETE
/// Differential Moisture Capacity from head
double TWorld::DmcNode(
        double head,           // current head value of this node
        const  HORIZON *hor, bool on_dmch)   // parameters of horizon this node belongs to
{
    LUT *l = hor->lut;

    if (head >= 0)
        return l->hydro[DMCC_COL].last();

    // DMCC depends on DMCH
    if (on_dmch) {
        auto it = std::lower_bound(l->hydro[DMCH_COL].begin(), l->hydro[DMCH_COL].end(), head);
        if (it == l->hydro[DMCH_COL].begin()) {
            return(l->hydro[DMCC_COL][0]);
        } else if (it == l->hydro[DMCH_COL].end()) {
            return(l->hydro[DMCC_COL].last());
        } else {
            int lowerIndex = std::distance(l->hydro[DMCH_COL].begin(), it - 1);
            int upperIndex = std::distance(l->hydro[DMCH_COL].begin(), it);

            double dmccl = l->hydro[DMCC_COL][lowerIndex];
            double dmchl = l->hydro[DMCH_COL][lowerIndex];
            double dmccu = l->hydro[DMCC_COL][upperIndex];
            double dmchu = l->hydro[DMCH_COL][upperIndex];

            return (dmccl + (head-dmchl)*(dmccu-dmccl)/(dmchu-dmchl));
        }
    }else {
        // DMCC depends directly on H
        auto it = std::lower_bound(l->hydro[H_COL].begin(), l->hydro[H_COL].end(), head);
        if (it == l->hydro[H_COL].begin()) {
            return(l->hydro[DMCC_COL][0]);
        } else if (it == l->hydro[H_COL].end()) {
            return(l->hydro[DMCC_COL].last());
        } else {
            double lH = *(it - 1);
            double uH = *it;
            double f = (head-lH)/(uH-lH);

            int lowerIndex = std::distance(l->hydro[H_COL].begin(), it - 1); // Index of lower bound
            int upperIndex = std::distance(l->hydro[H_COL].begin(), it);     // Index of upper bound

            double lC = l->hydro[DMCC_COL][lowerIndex];
            double uC = l->hydro[DMCC_COL][upperIndex];

            return (lC+f*(uC-lC));
        }
    }
}
//-----------------------------------------------------------------------------------

/*** swatre 2009
 * double DmcNode(
    double head,
    const  HORIZON *hor)
{
    int i;         // index in LUT where dmch[i] <= head <= dmch[i+1]
    const LUT *l;  // lut of this horizon

     //dit gaat niet goed als profiel van verzadigd naar onversazigd swithched:
    //if (head >= 0) return 0;

    if (head >= -1.0E-2)
            return LUT_LinIntPol(hor->lut, DMCC_COL, head, DMCH_COL);


    l = hor->lut;
    i = LUT_Index_LE(l, head, DMCH_COL);
    i = MIN(LUT_nrRows(l)-2, i);
    i = MAX(i, 0);

    return LUT_ValueAt(l, DMCC_COL, i) +
            (head - LUT_ValueAt(l, DMCH_COL, i)) *
         (LUT_ValueAt(l,DMCC_COL, i+1)-LUT_ValueAt(l,DMCC_COL, i))/
         (LUT_ValueAt(l,DMCH_COL, i+1)-LUT_ValueAt(l,DMCH_COL, i));
}
***/
