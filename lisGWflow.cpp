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
//    cTMap *SoilDepthinit;
//    cTMap *SoilDepth;
    if (SwitchTwoLayer) {
        pore = ThetaS2;
        ksat = Ksat2;
    //    SoilDepthinit = SoilDepth2init;
     //   SoilDepth = SoilDepth2;
    } else {
        pore = Poreeff;
        ksat = Ksateff;
     //   SoilDepthinit = SoilDepth1init;
     //   SoilDepth = SoilDepth1;
    }

    if (!SwitchExplicitGWflow) {

        // SWAT method, no flow but ditrect contribution to baseflow from perpendicular network

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            double CellArea_ = CellArea->Drc;

            //=== GW recharge
            Perc->Drc = cell_Percolation(r, c, GW_recharge); // in m
            double GWrecharge = Perc->Drc * CellArea_; // m3
            // GW recharge same principle as percolation, in m3        //=== lateral GW outflow
            // ksat is already in m per timestep

            double GWdeep_ = 0;//GW_deep * _dt/(86400000)*CellArea_; //*qSqrt(GWWH->Drc)
            // percolation from GW to deeper level, to cause decline in dry periods

            double GWVol_ = GWVol->Drc;//outflow m3
            double wh = GWVol_/CellArea_/pore->Drc;

            double GWout_ = GW_flow * CellArea_ * ksat->Drc * BaseflowL->Drc * pore->Drc; // m3 volume out from every cell
            //m3:  GW_flow*ksat*dt * ((dx/L)^b) *crosssection of flow dh*dx; //*porosity ??
            GWout_ *= (1+Grad->Drc);  // ???? add effect of slope
            GWout_ = wh > GW_threshold ?  GWout_ * (wh - GW_threshold) * (1-exp(-GW_threshold*wh)) : 0.0;
            // stop outflow when some minimum GW level, 2.4.2.10 in SWAT
            // apply a smooth threshold

            // ==== update GW level

            GWout_ = std::min(GWout_, GWVol_+ GWrecharge - GWdeep_);
            // cannot be more than there is
            GWVol_ = GWVol_ + GWrecharge - GWout_ - GWdeep_; //m3
            GWVol_ = std::max(GWVol_, 0.0);
            double sd = SwitchTwoLayer ? (SoilDepth2->Drc-SoilDepth1->Drc)-0.1 : SoilDepth1->Drc-0.1;
            GWVol_ = std::min(GWVol_, sd*CellArea_);

            GWWH->Drc = GWVol->Drc/CellArea_/pore->Drc;
            //update GW volume
            GWout->Drc = GWout_;
            GWVol->Drc = GWVol_;

            Qbin->Drc = 0;

            // change soildepth2 with GW changes
//            if (GWWH->Drc > 0) {
//                if (GWWH->Drc > 0) {
//                    double dh = std::max(0.1,SoilDepthinit->Drc - GWWH->Drc);
//                    SoilDepth->Drc = dh;
//                    GWWH->Drc = SoilDepthinit->Drc - dh;
//                    GWVol->Drc = pore->Drc*GWWH->Drc * CellArea->Drc;
//                }

//                GWVol->Drc = pore->Drc*GWWH->Drc * CellArea_;
//            }
            GWWHmax->Drc = std::max(GWWHmax->Drc, GWWH->Drc);

        }}

        AccufluxGW(crlinkedlddbase_, GWout, Qbin, ChannelWidth);
        // LDDbase, Qin, Qout, chanwidth used as flag, move the gw flow to the channel,
        // Qbin is inflow to the channel from the surrounding cells in m3 per timestep

        //double factor = exp(-GW_lag);
         #pragma omp parallel for num_threads(userCores)
         FOR_ROW_COL_MV_CHL {
             Qbase->Drc = Qbin->Drc * ChannelWidth->Drc/_dx;//m3 added per timestep, for MB
             // do this or not? for very small channel a lot of water is added but what haoppens to the rest
             ChannelWaterVol->Drc += Qbase->Drc;
             // flow according to SWAT 2009, page 174 manual, eq 2.4.2.8

             GWVol->Drc -= Qbase->Drc;
             GWWH->Drc = GWVol->Drc/CellArea->Drc/ThetaS2->Drc;
         }}

    } else {

        // 2D eplicit flow method based on H+Z differences

        GWFlow2D();

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {

            Qbin->Drc = 2 * ksat->Drc * _dx*GWWH->Drc; // 2 is form two sides into the channel in the middle

            Qbase->Drc = Qbin->Drc;//m3 added per timestep, for MB
            // do this or not? for very small channel a lot of water is added but what haoppens to the rest
            ChannelWaterVol->Drc += Qbase->Drc;

            GWVol->Drc -= Qbin->Drc;
            GWWH->Drc = GWVol->Drc/CellArea->Drc/pore->Drc;
        }}


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
    }
}

void TWorld::GWFlow2D(void)
{
    cTMap *pore;
    cTMap *ksat;
//    cTMap *SoilDepthinit;
    cTMap *SoilDepth;
    cTMap *z = DEM;
    cTMap *h = GWWH;
    cTMap *vol = GWVol;
    if (SwitchTwoLayer) {
        pore = ThetaS2;
        ksat = Ksat2;
//        SoilDepthinit = SoilDepth2init;
        SoilDepth = SoilDepth2;
    } else {
        pore = Poreeff;
        ksat = Ksateff;
//        SoilDepthinit = SoilDepth1init;
        SoilDepth = SoilDepth1;
    }

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        z->Drc = DEM->Drc - SoilDepth->Drc + 100;

        Perc->Drc = cell_Percolation(r, c, GW_recharge); // in m
        GWrecharge->Drc = Perc->Drc * CellArea->Drc; // m3

        GWVol->Drc += GWrecharge->Drc;// - GWdeep->Drc;
        GWWH->Drc = GWVol->Drc/CellArea->Drc/pore->Drc;
        GWout->Drc = 0;
        tma->Drc = 0;
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

            double df_x1 = GW_flow * ksat->Drc * (h_x1 * _dx) * dh_x1/_dx * pore->Drc;  // ksat has already dt
            double df_y1 = GW_flow * ksat->Drc * (h_y1 * _dx) * dh_y1/_dx * pore->Drc;
            double df_x2 = GW_flow * ksat->Drc * (h_x2 * _dx) * dh_x2/_dx * pore->Drc;
            double df_y2 = GW_flow * ksat->Drc * (h_y2 * _dx) * dh_y2/_dx * pore->Drc;

            double sign_x1 = df_x1 < 0 ? -1.0 : 1.0;
            double sign_y1 = df_y1 < 0 ? -1.0 : 1.0;
            double sign_x2 = df_x2 < 0 ? -1.0 : 1.0;
            double sign_y2 = df_y2 < 0 ? -1.0 : 1.0;

            df_x1 = std::min(v_x1,abs(df_x1))*sign_x1;
            df_y1 = std::min(v_y1,abs(df_y1))*sign_y1;
            df_x2 = std::min(v_x2,abs(df_x2))*sign_x2;
            df_y2 = std::min(v_y2,abs(df_y2))*sign_y2;

            double dflux = (df_x1 + df_x2 + df_y1 + df_y2);
            // m3 = m/s * s * (m*m) * m/m
            GWout->Drc = dflux;

        }}
    }

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double CellArea_ = CellArea->Drc;

        //update GW volume
        double sd = SwitchTwoLayer ? (SoilDepth2->Drc-SoilDepth1->Drc)-0.1 : SoilDepth1->Drc-0.1;
        GWVol->Drc = std::max(GWVol->Drc + GWout->Drc,0.0);
        GWVol->Drc = std::min(GWVol->Drc + GWout->Drc, sd*CellArea_);

        GWWH->Drc = GWVol->Drc/CellArea_/pore->Drc;

//        // change soildepth2 with GW changes
//        if (GWWH->Drc > 0) {
//            double dh = std::max(0.1,SoilDepthinit->Drc - GWWH->Drc);
//            SoilDepth->Drc = dh;
//            GWWH->Drc = SoilDepthinit->Drc - dh;
//            GWVol->Drc = pore->Drc*GWWH->Drc * CellArea_;
//        }
        GWWHmax->Drc = std::max(GWWHmax->Drc, GWWH->Drc);
    }}

}
