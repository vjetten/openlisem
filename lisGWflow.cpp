/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/


#include "model.h"


// results in a new Qbin (Qbase in for channel flow)
void TWorld::GroundwaterFlow(void)
{
    cTMap *pore;
    cTMap *ksat;
    cTMap *SoilDepthinit;
    cTMap *SoilDepth;
    if (SwitchTwoLayer) {
        pore = ThetaS2;
        ksat = Ksat2;
        SoilDepthinit = SoilDepth2init;
        SoilDepth = SoilDepth2;
    } else {
        pore = Poreeff;
        ksat = Ksateff;
        SoilDepthinit = SoilDepth1init;
        SoilDepth = SoilDepth1;
    }

    // recharge and deep percolation
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        Perc->Drc = cell_Percolation(r, c, GW_recharge); // in m
        GWrecharge->Drc = Perc->Drc * CellArea->Drc; // m3

        //GWdeep->Drc_ = 0;//GW_deep * _dt/(86400000)*CellArea_; //*qSqrt(GWWH->Drc)
        // percolation from GW to deeper level, to cause decline in dry periods

        GWVol->Drc += GWrecharge->Drc;// - GWdeep->Drc;
        GWWH->Drc = GWVol->Drc/CellArea->Drc/pore->Drc;
        GWout->Drc = 0;
    }}


    if (SwitchSWATGWflow)
        GWFlowLDD();
    // GW contribution to baseflow according to SWAT
    // do this always

    if (SwitchExplicitGWflow)
       GWFlow2D();
    // 2D eplicit flow method based on H+Z differences


    // change the soil depth with GWWH
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        GWWH->Drc = GWVol->Drc/CellArea->Drc/pore->Drc;

        // change soildepth2 with GW changes
        if (SwitchImpermeable && GWWH->Drc > 0) {
            SoilDepth->Drc = SoilDepthinit->Drc - GWWH->Drc;
        }

        GWWHmax->Drc = std::max(GWWHmax->Drc, GWWH->Drc);
    }}


    // do the baseflow
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        // fast component from LDD
        Qbase->Drc = Qbin->Drc * ChannelWidth->Drc/_dx;//m3 added per timestep, for MB

        //add the 2D flow
        if (SwitchExplicitGWflow) {
            double GWWH_ = GWWH->Drc;
            Qbase->Drc += 2 * GW_inflow * ksat->Drc * _dx*GWWH_;// * (1-exp(-4*GWWH_));
            // assumption is a DARCY pressure gradient of hH/dL = 1.0
            // not this: makes it deped on cell size which is not what we want* (GWWH_/(ChannelAdj->Drc/2));
        }
        // sort of Darcy with the pressure term as the GW height over half the distance to the channel
        // do not add pore because Ksat is already the permeability of the matrix, not just the pores in the matrix
        // 2 is form two sides into the channel in the middle

        if (GWVol->Drc*0.9 - Qbase->Drc < 0)
            Qbase->Drc = GWVol->Drc*0.9;

        //double qb = (1-GW_lag)*Qbase->Drc + GW_lag*Qbaseprev->Drc;

        GWVol->Drc -= Qbase->Drc;

        ChannelWaterVol->Drc += Qbase->Drc;

        GWWH->Drc = GWVol->Drc/CellArea->Drc/pore->Drc;

        //Qbaseprev->Drc = Qbase->Drc;
    }}
}

// flow according to SWAT 2009, page 174 manual, eq 2.4.2.8
void TWorld::GWFlowLDD(void)
{
    bool doit = false;

    cTMap *pore;
    cTMap *ksat;
    if (SwitchTwoLayer) {
        pore = ThetaS2;
        ksat = Ksat2;
    } else {
        pore = Poreeff;
        ksat = Ksateff;
    }
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double CellArea_ = CellArea->Drc;
        // between 0 and soildepth - 0.1m
        double maxvol = CellArea_ * (SwitchTwoLayer ? (SoilDepth2->Drc-SoilDepth1->Drc)-0.1 : SoilDepth1->Drc-0.1);
        double GWWH_ = GWWH->Drc;

        double GWout_ = GW_flow * CellArea_ * ksat->Drc * BaseflowL->Drc; // m3 volume out from every cell
 //       double GWout_ = GW_flow * _dx * GWWH->Drc * ksat->Drc * BaseflowL->Drc;
        //m3:  GW_flow* ksat*dt * ((dx/L)^b) *crosssection of flow dh*dx;
        //NOTE cross section changed to cellarea!?

        // DO NOT include pore, ksat is already a flux from a porous soil and includes dt

        GWout_ = GWWH_ > GW_threshold ?  GWout_ * (GWWH_ - GW_threshold) * (1-exp(-GW_threshold*GWWH_)) : 0.0;
        // stop outflow when some minimum GW level, 2.4.2.10 in SWAT
        // apply a smooth threshold with exponential function

        // GWout_ *= (1+Grad->Drc);  // ???? add effect of slope

        if (GWout_ > 0) {
            if (GWVol->Drc - GWout_ < 0)
                GWout_ = GWVol->Drc;

            tmb->Drc = GWout_; // used in accufluwGW
            GWout->Drc = GWout_;

            GWVol->Drc -= GWout_; // subtract from volume
            GWVol->Drc = std::max(GWVol->Drc,0.0);
            GWVol->Drc = std::min(GWVol->Drc, maxvol);

            GWWH->Drc = GWVol->Drc/CellArea_/pore->Drc;

            doit = true;
        }

        Qbin->Drc = 0;
    }}

    if (doit)
        AccufluxGW(crlinkedlddbase_, tmb, Qbin, ChannelWidth);
    // Qbin now has the fast component
}


void TWorld::GWFlow2D(void)
{
    cTMap *pore;
    cTMap *ksat;
    cTMap *SoilDepthinit;
    cTMap *SoilDepth;
    cTMap *z = DEM;
    cTMap *h = GWWH;
    cTMap *vol = GWVol;
    if (SwitchTwoLayer) {
        pore = ThetaS2;
        ksat = Ksat2;
        SoilDepthinit = SoilDepth2init;
        SoilDepth = SoilDepth2;
    } else {
        pore = Poreeff;
        ksat = Ksateff;
        SoilDepthinit = SoilDepth1init;
        SoilDepth = SoilDepth1;
    }

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        z->Drc = DEM->Drc - SoilDepth->Drc + 100;
        tma->Drc = 0;
        tmb->Drc = 0;
    }}

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (GWWH->Drc > GW_threshold) {
            tma->Drc = 1;
            if (c > 0 && !MV(r,c-1)        ) tma->data[r][c-1] = 1;
            if (c < _nrCols-1 && !MV(r,c+1)) tma->data[r][c+1] = 1;
            if (r > 0 && !MV(r-1,c)        ) tma->data[r-1][c] = 1;
            if (r < _nrRows-1 && !MV(r+1,c)) tma->data[r+1][c] = 1;

            if (c > 0 && r > 0 && !MV(r-1,c-1)                ) tma->data[r-1][c-1]=1;
            if (c < _nrCols-1 && r < _nrRows-1 && !MV(r+1,c+1)) tma->data[r+1][c+1]=1;
            if (r > 0 && c < _nrCols-1 && !MV(r-1,c+1)        ) tma->data[r-1][c+1]=1;
            if (c > 0 && r < _nrRows-1 && !MV(r+1,c-1)        ) tma->data[r+1][c-1]=1;
        }
    }}

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(tma->Drc > 0) {
            double H = GWWH->Drc;
            double Z = z->Drc;
            double V = vol->Drc;

            bool bc1 = c > 0 && !MV(r,c-1)        ;
            bool bc2 = c < _nrCols-1 && !MV(r,c+1);
            bool br1 = r > 0 && !MV(r-1,c)        ;
            bool br2 = r < _nrRows-1 && !MV(r+1,c);

            double z_x1 =  bc1 ? z->data[r][c-1] : Z;
            double z_x2 =  bc2 ? z->data[r][c+1] : Z;
            double z_y1 =  br1 ? z->data[r-1][c] : Z;
            double z_y2 =  br2 ? z->data[r+1][c] : Z;

            double h_x1 =  bc1 ? h->data[r][c-1] : H;
            double h_x2 =  bc2 ? h->data[r][c+1] : H;
            double h_y1 =  br1 ? h->data[r-1][c] : H;
            double h_y2 =  br2 ? h->data[r+1][c] : H;

            double v_x1 =  bc1 ?vol->data[r][c-1] : V;
            double v_x2 =  bc2 ?vol->data[r][c+1] : V;
            double v_y1 =  br1 ?vol->data[r-1][c] : V;
            double v_y2 =  br2 ?vol->data[r+1][c] : V;

            double dh_x1 = (h_x1 + z_x1) - (H+Z);
            double dh_x2 = (h_x2 + z_x2) - (H+Z);
            double dh_y1 = (h_y1 + z_y1) - (H+Z);
            double dh_y2 = (h_y2 + z_y2) - (H+Z);

            double f = 1.0/3.0;
            dh_x1 = (1-f)*(z_x1 - Z)+f*dh_x1;
            dh_x2 = (1-f)*(z_x2 - Z)+f*dh_y1;
            dh_y1 = (1-f)*(z_y1 - Z)+f*dh_x2;
            dh_y2 = (1-f)*(z_y2 - Z)+f*dh_y2;

            double df_x1 = GW_flow * ksat->Drc * (h_x1 * _dx) * dh_x1/_dx * pore->Drc;  // ksat has already dt
            double df_y1 = GW_flow * ksat->Drc * (h_y1 * _dx) * dh_y1/_dx * pore->Drc;
            double df_x2 = GW_flow * ksat->Drc * (h_x2 * _dx) * dh_x2/_dx * pore->Drc;
            double df_y2 = GW_flow * ksat->Drc * (h_y2 * _dx) * dh_y2/_dx * pore->Drc;
            // m3 = m/s * s * (m*m) * m/m

            double sign_x1 = df_x1 < 0 ? -1.0 : 1.0;
            double sign_y1 = df_y1 < 0 ? -1.0 : 1.0;
            double sign_x2 = df_x2 < 0 ? -1.0 : 1.0;
            double sign_y2 = df_y2 < 0 ? -1.0 : 1.0;

            f = 1.0;
            df_x1 = std::min(v_x1*f,abs(df_x1)) * sign_x1;
            df_y1 = std::min(v_y1*f,abs(df_y1)) * sign_y1;
            df_x2 = std::min(v_x2*f,abs(df_x2)) * sign_x2;
            df_y2 = std::min(v_y2*f,abs(df_y2)) * sign_y2;

            double dflux = (df_x1 + df_x2 + df_y1 + df_y2);
            double maxvol = CellArea->Drc * (SwitchTwoLayer ? (SoilDepth2->Drc-SoilDepth1->Drc)-0.1 : SoilDepth1->Drc-0.1);

            if (dflux < 0)
                dflux = std::max(-vol->Drc, dflux);
            if (vol->Drc + dflux > maxvol)
                dflux = std::min(dflux, vol->Drc);

            GWout->Drc = dflux;
            vol->Drc += dflux;
            vol->Drc = std::max(0.0, vol->Drc);
            h->Drc = vol->Drc/CellArea->Drc/pore->Drc;

        }}
    }
    //GWout now has the flow but is not used further
}




    /* upstream ovvver baseLDD, fdoes not work very well
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                double CellArea_ = CellArea->Drc;

                //=== GW recharge
                Perc->Drc = cell_Percolation(r, c, GW_recharge); // in m
                GWrecharge->Drc = Perc->Drc * CellArea_; // m3
                // GW recharge same principle as percolation, in m3        //=== lateral GW outflow

                GWdeep->Drc = 0;//GW_deep * _dt/(86400000)*CellArea_; //*qSqrt(GWWH->Drc)
                // percolation from GW to deeper level, to cause decline in dry periods

                double GWVol_ = GWVol->Drc;//outflow m3
                double wh = GWVol_/CellArea_/pore->Drc;
                double GWout_ = GW_flow * _dx * wh * ksat->Drc * (1+Grad->Drc + 0.01); //* pore->Drc
                // m3 volume out from every cell, grad is the sinus is seen as the pressure
                GWout_ = wh > GW_threshold ?  GWout_ * (wh - GW_threshold) * (1-exp(-GW_threshold*wh)) : 0.0;
                // apply a smooth threshold
                GWout_ = std::min(GWout_, GWVol_+ GWrecharge->Drc - GWdeep->Drc);

                GWout->Drc = GWout_;
            }}

            UpstreamGW(crlinkedlddbase_, GWout, Qbin);
            // explicit GW flow from cell to cell

            // ==== update GW level
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                // cannot be more than there is
                GWVol->Drc = std::max(0.0, GWVol->Drc  + GWrecharge->Drc + Qbin->Drc - GWout->Drc - GWdeep->Drc); //m3
                //update GW volume

                GWWH->Drc = GWVol->Drc/CellArea->Drc/pore->Drc;  //for display

                // change soildepth2 with GW changes
                if (GWWH->Drc > 0) {
                    double dh = std::max(0.1,SoilDepthinit->Drc - GWWH->Drc);
                    SoilDepth->Drc = dh;
                    GWWH->Drc = SoilDepthinit->Drc - dh;
                    GWVol->Drc = pore->Drc*GWWH->Drc * CellArea->Drc;
                }
                GWWHmax->Drc = std::max(GWWHmax->Drc, GWWH->Drc);

            }}

    */
