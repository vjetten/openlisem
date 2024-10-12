
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
double TWorld::HcoNode(double head,const HORIZON *hor,double calib)
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
/// Differential Moisture Capacity from head
double TWorld::DmcNode(
        double head,           // current head value of this node
        const  HORIZON *hor)   // parameters of horizon this node belongs to
{
    LUT *l = hor->lut;

    if (head >= 0) {
        return l->hydro[DMCC_COL].last();
    }

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
