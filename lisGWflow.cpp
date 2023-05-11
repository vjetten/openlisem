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
#define GWpow 1.0

//---------------------------------------------------------------------------
void TWorld::GroundwaterFlow(void)
{
    cTMap *pore;
    cTMap *SoilDepthinit;
    cTMap *SoilDepth;
    if (SwitchTwoLayer) {
        pore = ThetaS2;
        SoilDepthinit = SoilDepth2init;
        SoilDepth = SoilDepth2;
    } else {
        pore = Poreeff;
        SoilDepthinit = SoilDepth1init;
        SoilDepth = SoilDepth1;
    }
//double tot = 0;
//double totr = 0;
//GWdeeptot = 0;

    // add recharge and subtract deep percolation
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        Perc->Drc = cell_Percolation(r, c, GW_recharge); // in m
        GWrecharge->Drc = Perc->Drc * CHAdjDX->Drc; // m3

        GWdeep->Drc = GW_deep * CHAdjDX->Drc;
        // percolation from GW to deeper level, to cause decline in dry periods

        double maxvol = SoilDepthinit->Drc * CellArea->Drc * pore->Drc;

        if (GWVol->Drc + GWrecharge->Drc - GWdeep->Drc > maxvol)
            GWrecharge->Drc = maxvol - GWVol->Drc + GWdeep->Drc;
        if (GWVol->Drc + GWrecharge->Drc - GWdeep->Drc < 0)
            GWdeep->Drc = GWVol->Drc + GWrecharge->Drc;
//        GWdeeptot += GWdeep->Drc;
//        tot += GWVol->Drc;
//        totr += GWrecharge->Drc;
        GWVol->Drc += GWrecharge->Drc - GWdeep->Drc;
        GWWH->Drc = GWVol->Drc/(CellArea->Drc*pore->Drc);
        GWout->Drc = 0;
    }}
//qDebug() << GWdeeptot << totr << tot;

    // results in GWout flux between cells based on pressure differences
    if (SwitchGWflow)
        GWFlow2D();
    // flow with pressure differences
    if (SwitchLDDGWflow)
        GWFlowLDDKsat();
    if (SwitchSWATGWflow)
        GWFlowSWAT();

    // change the soil depth with GWWH
    if (SwitchGWChangeSD) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            GWWH->Drc = GWVol->Drc/CellArea->Drc/pore->Drc;
            // change soildepth2 with GW changes
            if (GWWH->Drc > 0) {
                SoilDepth->Drc = SoilDepthinit->Drc - GWWH->Drc;
            }

            GWWHmax->Drc = std::max(GWWHmax->Drc, GWWH->Drc);
        }}
    }
}
//---------------------------------------------------------------------------
void TWorld::GWFlowLDDKsat(void)
{
    cTMap *pore;
    cTMap *ksat;
    cTMap *SD;
    cTMap *h = GWWH;
    if (SwitchTwoLayer) {
        pore = ThetaS2;
        ksat = Ksat2;
        SD = SoilDepth2init;
    } else {
        pore = Poreeff;
        ksat = Ksateff;
        SD = SoilDepth1init;
    }

    // adjust for threshold
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        h->Drc = std::max(0.0, h->Drc - GW_threshold);
    }}

    // calculate GW flow angle along network
    //Fill(*tma, 0.0);
    Fill(*tmb, 0.0);
    Fill(*tmc, 0.0);
    for(long i_ =  0; i_ < crlinkedlddbase_.size(); i_++)
    {
        int r = crlinkedlddbase_.at(i_).r;
        int c = crlinkedlddbase_.at(i_).c;

        double Hup = 0;
        if (crlinkedlddbase_.at(i_).nr > 0) {
            double cnt = 0;
            for(int j = 0; j < crlinkedlddbase_.at(i_).nr; j++) {
                int rr = crlinkedlddbase_.at(i_).inn[j].r;
                int cr = crlinkedlddbase_.at(i_).inn[j].c;
                Hup += GWz->Drc + h->Drcr;
                cnt+=1.0;
            }
            Hup /= cnt;
        }
        double H = GWz->Drc + h->Drc;

        tmb->Drc = cos(atan(fabs(Hup - H)/_dx)); // hydraulic gradient angle
        //tmb->Drc = fabs(Zup - Z)/_dx;
    }

    int step = 1;
    while (_dt/(_dx*(double)step) > 0.3)
        step++;

    for (int j = 0; j < step; j++) {
        // loop step times for explicit GW flow

        // calculate all fluxes
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            tmc->Drc = 1/(double)step * GW_flow * ksat->Drc * (h->Drc*_dx) * tmb->Drc;
            tmc->Drc = std::min(tmc->Drc, GWVol->Drc*MaxGWDepthfrac);
            // flow is ksat over terrain gradient in m3, cannot be more than volume present
        }}

        //sum all fluxes over network
        for(long i_ =  0; i_ < crlinkedlddbase_.size(); i_++)
        {
            int r = crlinkedlddbase_.at(i_).r;
            int c = crlinkedlddbase_.at(i_).c;

            double Qin = 0;
            // sum fluxes in from incoming branches
            if (crlinkedlddbase_.at(i_).nr > 0) {
                for(int j = 0; j < crlinkedlddbase_.at(i_).nr; j++) {
                    int rr = crlinkedlddbase_.at(i_).inn[j].r;
                    int cr = crlinkedlddbase_.at(i_).inn[j].c;
                    Qin += tmc->Drcr;
                }
            }
            //tma->Drc = Qin;
            double flux = Qin - tmc->Drc;
            double maxvol = CellArea->Drc * SD->Drc * pore->Drc;
            double vol = GWVol->Drc;
            if (vol + flux > maxvol) {
                //flux = maxvol - vol;
                Qin = maxvol - tmc->Drc;
                flux = Qin - tmc->Drc;
            }

            if (vol + flux < 0)
                flux = -vol;
            GWVol->Drc += flux;
            GWWH->Drc = GWVol->Drc/CellArea->Drc/pore->Drc;
            GWout->Drc = flux;// Qin;
        }
    }

    Average3x3(*GWWH, *LDDbaseflow);
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        GWVol->Drc = GWWH->Drc*CellArea->Drc*pore->Drc;
    }}

}
//---------------------------------------------------------------------------
void TWorld::GWFlow2D(void)
{
    cTMap *pore;
    cTMap *ksat;
    cTMap *SD;
    cTMap *z = GWz;
    cTMap *h = GWWH;
    cTMap *vol = GWVol;
    if (SwitchTwoLayer) {
        pore = ThetaS2;
        ksat = Ksat2;
        SD = SoilDepth2init;
    } else {
        pore = Poreeff;
        ksat = Ksateff;
        SD = SoilDepth1init;
    }

    // adjust for threshold
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        h->Drc = std::max(0.0, h->Drc - GW_threshold);
        //tma->Drc = 0;
    }}

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double H = h->Drc;
        double Z = z->Drc;
        double V = vol->Drc;

        bool bc1 = (c > 0 && !MV(r,c-1)        );
        bool bc2 = (c < _nrCols-1 && !MV(r,c+1));
        bool br1 = (r > 0 && !MV(r-1,c)        );
        bool br2 = (r < _nrRows-1 && !MV(r+1,c));

        double z_x1 =  bc1 ? z->data[r][c-1] : Z;
        double z_x2 =  bc2 ? z->data[r][c+1] : Z;
        double z_y1 =  br1 ? z->data[r-1][c] : Z;
        double z_y2 =  br2 ? z->data[r+1][c] : Z;

        double h_x1 =  bc1 ? h->data[r][c-1] : H;
        double h_x2 =  bc2 ? h->data[r][c+1] : H;
        double h_y1 =  br1 ? h->data[r-1][c] : H;
        double h_y2 =  br2 ? h->data[r+1][c] : H;

        double v_x1 =  bc1 ? vol->data[r][c-1] : V;
        double v_x2 =  bc2 ? vol->data[r][c+1] : V;
        double v_y1 =  br1 ? vol->data[r-1][c] : V;
        double v_y2 =  br2 ? vol->data[r+1][c] : V;

        double dh_x1 = (h_x1 + z_x1) - (H+Z);
        double dh_x2 = (h_x2 + z_x2) - (H+Z);
        double dh_y1 = (h_y1 + z_y1) - (H+Z);
        double dh_y2 = (h_y2 + z_y2) - (H+Z);

        double dz_x1 = z_x1 -Z;
        double dz_x2 = z_x2 -Z;
        double dz_y1 = z_y1 -Z;
        double dz_y2 = z_y2 -Z;

        // flow = Ksat * cross section * hydraulic gradient Ks * A * dH/dL
        double ff = GW_flow;

        double df_x1 = ff * ksat->Drc * (h_x1*_dx) * dh_x1/cos(atan(dz_x1/_dx)); // (dh_x1/_dx);
        double df_x2 = ff * ksat->Drc * (h_x2*_dx) * dh_x2/cos(atan(dz_x2/_dx)); // (dh_x2/_dx);
        double df_y1 = ff * ksat->Drc * (h_y1*_dx) * dh_y1/cos(atan(dz_y1/_dx)); // (dh_y1/_dx);
        double df_y2 = ff * ksat->Drc * (h_y2*_dx) * dh_y2/cos(atan(dz_y2/_dx)); // (dh_y2/_dx);
        // m3 = m/s * s * (m*m) * m/m

        // limit flow to a frcation of the volume present
        double f = MaxGWDepthfrac;
        df_x1 = std::min(v_x1*f, fabs(df_x1)) * (df_x1 < 0 ? -1.0 : 1.0);
        df_x2 = std::min(v_x2*f, fabs(df_x2)) * (df_x2 < 0 ? -1.0 : 1.0);
        df_y1 = std::min(v_y1*f, fabs(df_y1)) * (df_y1 < 0 ? -1.0 : 1.0);
        df_y2 = std::min(v_y2*f, fabs(df_y2)) * (df_y2 < 0 ? -1.0 : 1.0);

        //avoid single pixels with MV on 3 sides that fill up
        if( df_x1 < 0 && bc1) df_x1 = 0.0;
        if( df_x2 > 0 && bc2) df_x2 = 0.0;
        if( df_y1 < 0 && br1) df_y1 = 0.0;
        if( df_y2 > 0 && br2) df_y2 = 0.0;

        // sum and correct all fluxes
        double dflux = (df_x1 + df_x2 + df_y1 + df_y2);
        double maxvol = CellArea->Drc * SD->Drc * pore->Drc;
        if (V + dflux > maxvol)
            dflux =  maxvol - V;
        if (V + dflux < 0)
            dflux = -V;
        //fill tma with the resulting flux of a cell
        GWout->Drc = dflux;
    }}

    // adjust the vol
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        GWVol->Drc += GWout->Drc;
        // update gwvol with flux
        GWWH->Drc = GWVol->Drc/CellArea->Drc/pore->Drc;
        // recalc gwwh

       // GWout->Drc = tma->Drc;
        // flux is gwout for baseflow
    }}

}

//---------------------------------------------------------------------------

// flow according to SWAT 2009, page 174 manual, eq 2.4.2.8
void TWorld::GWFlowSWAT(void)
{
    bool doit = false;

    cTMap *ksat = Ksat2;
    cTMap *pore = ThetaS2;
    if (!SwitchTwoLayer) {
        ksat = Ksateff;
        pore = Thetaeff;
    }

    // calculated lateral flow
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        //        double GWout_ = GW_flow *  CellArea->Drc * ksat->Drc * BaseflowL->Drc; // m3 volume out from every cell
        double GWout_ = GW_flow * ksat->Drc * _dx * std::max(0.0, GWWH->Drc-GW_threshold) * BaseflowL->Drc;
        //m3:  ksat*dt  * dh*dx * ((dx/L)^b);  ksat * cross section * distance factor
        // stop outflow when some minimum GW level, 2.4.2.10 in SWAT

        //GWout_ = GWout_ * std::max(0.0, GWWH_ - GW_threshold) * (1-exp(-GW_threshold*GWWH_));
        // apply a smooth threshold with exponential function
        GWout_ = std::min(GWVol->Drc*MaxGWDepthfrac, GWout_);
        tmb->Drc = GWout_;

        // adjust volume with outflow
        if (GWout_ > 0) {
           GWVol->Drc -= GWout_; // subtract from volume
           GWWH->Drc = GWVol->Drc/ CellArea->Drc/pore->Drc;
           doit = true; // is someqwhere GWout > 0 do accuflux
        }
    }}

    if (doit)
        AccufluxGW(crlinkedlddbase_, tmb, GWout, ChannelWidth);
    // GWout has not the accumulated flow pattern, do NOT use this anymore for the mass balance
}

