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
#include "operation.h"

#define MaxGWDepthfrac 0.95

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

    // add recharge and subtract deep percolation
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double dxa = CellArea->Drc;// CHAdjDX->Drc;
        Perc->Drc = cell_Percolation(r, c, GW_recharge); // in m
        GWrecharge->Drc = Perc->Drc * dxa; // m3

        GWdeep->Drc = GW_deep*CellArea->Drc;
        // percolation from GW to deeper level, to cause decline in dry periods

        double maxvol = SoilDepth2init->Drc * CellArea->Drc * pore->Drc;
        if (GWVol->Drc + GWrecharge->Drc - GWdeep->Drc > maxvol)
            GWrecharge->Drc = maxvol - GWVol->Drc + GWdeep->Drc;
        GWVol->Drc += GWrecharge->Drc - GWdeep->Drc;
        GWWH->Drc = GWVol->Drc/(dxa*pore->Drc);

        // GWout->Drc = 0;
        // given a value in GWFlowLDD
    }}

    if (SwitchExplicitGWflow)
        GWFlowLDDKsat();
    // 2D eplicit flow method based on H+Z differences
    // results in Qbin which is the flow in m3

    // change the soil depth with GWWH
    if (SwitchGWChangeSD) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            GWWH->Drc = GWVol->Drc/CellArea->Drc/pore->Drc;

            // change soildepth2 with GW changes
            // when gw flow always impermeable !?
            if (SwitchImpermeable && GWWH->Drc > 0) {
                SoilDepth->Drc = SoilDepthinit->Drc - GWWH->Drc;
            }

            GWWHmax->Drc = std::max(GWWHmax->Drc, GWWH->Drc);
        }}
    }
}

void TWorld::GWFlowLDDKsat(void)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};
    cTMap *pore;
    cTMap *ksat;
    cTMap *SoilDepthinit;
    cTMap *z = DEM;
    if (SwitchTwoLayer) {
        pore = ThetaS2;
        ksat = Ksat2;
        SoilDepthinit = SoilDepth2init;
    } else {
        pore = Poreeff;
        ksat = Ksateff;
        SoilDepthinit = SoilDepth1init;
    }

    // calculate GW flow angle along network
    fill(*tmb, 0.0);
    for(long i_ =  0; i_ < crlinkedlddbase_.size(); i_++)
    {
        int r = crlinkedlddbase_.at(i_).r;
        int c = crlinkedlddbase_.at(i_).c;
        double Zup = 0;

        if (crlinkedlddbase_.at(i_).nr > 0) {
            double cnt = 0;
            for(int j = 0; j < crlinkedlddbase_.at(i_).nr; j++) {
                int rr = crlinkedlddbase_.at(i_).inn[j].r;
                int cr = crlinkedlddbase_.at(i_).inn[j].c;
                Zup += (z->Drcr - SoilDepthinit->Drcr);// + GWWH->Drcr;
                cnt+=1.0;
            }
            Zup /= cnt;
        }
        double Z = z->Drc - SoilDepthinit->Drc;// + GWWH->Drc;

        tmb->Drc = fabs(Zup - Z)/_dx + 0.001;
    }

    fill(*tmc, 0.0);
    int step = 5;
    for (int j = 0; j < step; j++) {
        // loop step times for explicit GW flow
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            tmc->Drc = 1/(double)step * GW_flow * ksat->Drc * GWWH->Drc *_dx * tmb->Drc;//Grad->Drc;
            // flow is ksat over terrain gradient
        }}

        for(long i_ =  0; i_ < crlinkedlddbase_.size(); i_++)
        {
            int r = crlinkedlddbase_.at(i_).r;
            int c = crlinkedlddbase_.at(i_).c;
            double Qin = 0;

            if (crlinkedlddbase_.at(i_).nr > 0) {
                for(int j = 0; j < crlinkedlddbase_.at(i_).nr; j++) {
                    int rr = crlinkedlddbase_.at(i_).inn[j].r;
                    int cr = crlinkedlddbase_.at(i_).inn[j].c;
                    Qin += tmc->Drcr;
                }
            }

            double vol = GWVol->Drc + Qin - tmc->Drc;
            GWVol->Drc = std::max(0.0, vol);
            GWWH->Drc = vol/CellArea->Drc/pore->Drc;
        }
    }

    //GWFlow2D();

    Average3x3(*GWWH, *LDDbaseflow);

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        GWVol->Drc = GWWH->Drc * pore->Drc * CellArea->Drc;
        GWout->Drc = GW_flow * ksat->Drc * GWWH->Drc * _dx * tmb->Drc;
    }}


}


    // DOES not work very well, very unstable
void TWorld::GWFlow2D(void)
{
    cTMap *pore;
    cTMap *ksat;
    cTMap *SoilDepthinit;
    cTMap *z = DEM;
    cTMap *h = GWWH;
    cTMap *vol = GWVol;
    if (SwitchTwoLayer) {
        pore = ThetaS2;
        ksat = Ksat2;
        SoilDepthinit = SoilDepth2init;
    } else {
        pore = Poreeff;
        ksat = Ksateff;
        SoilDepthinit = SoilDepth1init;
    }

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        z->Drc = DEM->Drc - SoilDepthinit->Drc + 10;
        tma->Drc = 0;
    }}

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
       if (GWWH->Drc > HMIN) {
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

//            double f = 0.33; // emphasize terrain slope
//            dh_x1 = (1-f)*(z_x1 - Z)+f*dh_x1;
//            dh_x2 = (1-f)*(z_x2 - Z)+f*dh_y1;
//            dh_y1 = (1-f)*(z_y1 - Z)+f*dh_x2;
//            dh_y2 = (1-f)*(z_y2 - Z)+f*dh_y2;

            // flow = Ksat * cross section * hydraulic gradient
            // ksat has already dt
            double ff = GW_flow;
            double df_x1 = ff* ksat->Drc * (h_x1 * _dx) * dh_x1/_dx;
            double df_y1 = ff* ksat->Drc * (h_y1 * _dx) * dh_y1/_dx;
            double df_x2 = ff* ksat->Drc * (h_x2 * _dx) * dh_x2/_dx;
            double df_y2 = ff* ksat->Drc * (h_y2 * _dx) * dh_y2/_dx;
            // m3 = m/s * s * (m*m) * m/m

            double f = 0.5;//MaxGWDepthfrac;
            df_x1 = std::min(v_x1*f,abs(df_x1)) * df_x1 < 0 ? -1.0 : 1.0;
            df_y1 = std::min(v_y1*f,abs(df_y1)) * df_y1 < 0 ? -1.0 : 1.0;
            df_x2 = std::min(v_x2*f,abs(df_x2)) * df_x2 < 0 ? -1.0 : 1.0;
            df_y2 = std::min(v_y2*f,abs(df_y2)) * df_y2 < 0 ? -1.0 : 1.0;

            double dflux = (df_x1 + df_x2 + df_y1 + df_y2);
            double maxvol = CellArea->Drc * SoilDepthinit->Drc;

            if (vol->Drc + dflux < 0)
                dflux = -vol->Drc * MaxGWDepthfrac;
            if (vol->Drc + dflux > maxvol)
                dflux = (maxvol - vol->Drc) * MaxGWDepthfrac;

            vol->Drc += dflux;
            vol->Drc = std::max(0.0, vol->Drc);
            h->Drc = vol->Drc/CellArea->Drc/pore->Drc;
        }}
    }

}


// flow according to SWAT 2009, page 174 manual, eq 2.4.2.8
//OBSOLETE
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
        double maxvol = CellArea_ * (SwitchTwoLayer ? (SoilDepth2->Drc/*+SoilDepth1->Drc*/) : SoilDepth1->Drc)*MaxGWDepthfrac;
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

