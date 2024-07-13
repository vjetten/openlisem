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
#define GWS 3.0

//---------------------------------------------------------------------------
void TWorld::GroundwaterFlow(void)
{
    cTMap *pore;
    cTMap *SoilDepthinit;
    cTMap *SoilDepth;
    if (SwitchTwoLayer) {
        pore = ThetaS2;
        SoilDepthinit = SoilDepth2init;
//        FOR_ROW_COL_MV_L {
//            SoilDepthinit->Drc = SoilDepthinit->Drc - SoilDepth1init->Drc;
//        }}
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

        double maxvol = SoilDepthinit->Drc * CHAdjDX->Drc * pore->Drc;

        if (GWVol->Drc + GWrecharge->Drc - GWdeep->Drc > maxvol)
            GWrecharge->Drc = maxvol - GWVol->Drc + GWdeep->Drc;
        if (GWVol->Drc + GWrecharge->Drc - GWdeep->Drc < 0)
            GWdeep->Drc = GWVol->Drc + GWrecharge->Drc;

        GWVol->Drc += GWrecharge->Drc - GWdeep->Drc;
        GWVol->Drc = std::min(maxvol, GWVol->Drc);
        GWWH->Drc = GWVol->Drc/(CHAdjDX->Drc*pore->Drc);
        GWout->Drc = 0;

    }}

    //flow with pressure differences
    if (SwitchGW2Dflow) {
        Fill(*tma, 0);
        for (int j_ = 0; j_ < 5; j_++) {
            GWFlow2D(0.2);
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                tma->Drc += GWout->Drc;
            }}
        }
        Copy(*tma,*GWout);

         //GWFlow2D(1.0);
    }

    if (SwitchGWSWOFflow) {
        double er = fullSWOF2GW(GWWH, GWU, GWV, GWz);
        //Fill(*tma,0);
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            GWout->Drc = sqrt(GWU->Drc*GWU->Drc + GWV->Drc*GWV->Drc) * _dx * GWWH->Drc * _dt;
            GWVol->Drc = GWWH->Drc*pore->Drc*CHAdjDX->Drc;
            //tma->Drc = sqrt(GWU->Drc*GWU->Drc + GWV->Drc*GWV->Drc)*1000*3600/_dt;
        }}
        //report(*tma,"gwv");
    }

    if (SwitchLDDGWflow)
        GWFlowLDDKsat(); // ldd with ksat based flow

    if (SwitchSWATGWflow)
        GWFlowSWAT();    // swat based flow using ldd and accuflux

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double maxvol = SoilDepthinit->Drc * CHAdjDX->Drc * pore->Drc;
        GWVol->Drc = std::min(maxvol, GWVol->Drc);
        GWWH->Drc = GWVol->Drc/CHAdjDX->Drc/pore->Drc;
        // change soildepth2 with GW changes
        if (GWWH->Drc > 0) {
            SoilDepth->Drc = SoilDepthinit->Drc - GWWH->Drc;
        }

        GWWHmax->Drc = std::max(GWWHmax->Drc, GWWH->Drc);
    }}
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
        tma->Drc = 0;
        tmb->Drc = 0;
        tmc->Drc = 0;
    }}

    // calculate GW flow angle along network
    for(long i_ =  0; i_ < crlinkedlddbase_.size(); i_++)
    {
        int r = crlinkedlddbase_.at(i_).r;
        int c = crlinkedlddbase_.at(i_).c;

        double Hup = 0;
        //double Zup = 0;
        if (crlinkedlddbase_.at(i_).nr > 0) {
            double cnt = 0;
            for(int j = 0; j < crlinkedlddbase_.at(i_).nr; j++) {
                int rr = crlinkedlddbase_.at(i_).inn[j].r;
                int cr = crlinkedlddbase_.at(i_).inn[j].c;
                //Hup += (GWS*GWz->Drcr + h->Drcr)/(GWS+1.0);  // GWS is a weight to emphasize the Z!
                Hup += (GWz->Drcr + h->Drcr);
                //Zup += GWz->Drcr;
                cnt+=1.0;
            }
            Hup /= cnt;  // average gradient angle for all upstream ceels
            //Zup /= cnt;
        }
        //double H = (GWS*GWz->Drc + h->Drc)/(GWS+1.0);
        double H = (GWz->Drc + h->Drc);

        tmb->Drc = cos(atan(fabs(Hup - H)/_dx)); // hydraulic gradient angle, always positive in the dircetion of the LDD
       // tmb->Drc = cos(atan(fabs(Zup - GWz->Drc)/_dx));
    }

    // calculate all incoming fluxes
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        tmc->Drc = GW_flow * ksat->Drc * (h->Drc*_dx) * tmb->Drc;
        tmc->Drc = std::min(tmc->Drc, GWVol->Drc*MaxGWDepthfrac);
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

        double Q1 = GW_flow * ksat->Drc * (h->Drc*_dx) * tmb->Drc; // before Qin
        double Q2 = GW_flow * ksat->Drc * ((h->Drc+Qin/CHAdjDX->Drc)*_dx) * tmb->Drc; // with Qin
        double Qn = 0.5*(Q1+Q2); // average flow out

        double dflux = Qin - Qn;
        double maxvol = CHAdjDX->Drc * SD->Drc * pore->Drc;
        double vol = GWVol->Drc;
        if (vol + Qin - Qn > maxvol) {
            Qin = maxvol + Qn;
        }
        if (vol + Qin - Qn < 0) {
            Qn = vol + Qin;
        }
        dflux = Qin - Qn;

        GWVol->Drc += dflux;
        GWWH->Drc = GWVol->Drc/CHAdjDX->Drc/pore->Drc;
        GWout->Drc = Qn;
   }

//    Average3x3(*GWWH, *LDDbaseflow, true);
//    #pragma omp parallel for num_threads(userCores)
//    FOR_ROW_COL_MV_L {
//        GWVol->Drc = GWWH->Drc*CHAdjDX->Drc*pore->Drc;
//    }}

}
//---------------------------------------------------------------------------
void TWorld::GWFlow2D(double factor)
{
    cTMap *pore;
    cTMap *ksat;
    cTMap *SD;
    cTMap *z = GWz;
    cTMap *h = GWWH;
    cTMap *vol = GWVol;
    if (SwitchTwoLayer) {
        pore = ThetaS2;
        ksat = Ksat2;   // has the unit m, *_dt time is included
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

        double dh_x1 = ((h_x1 + GWS*z_x1) - (H+GWS*Z))/(GWS+1.0);
        double dh_x2 = ((h_x2 + GWS*z_x2) - (H+GWS*Z))/(GWS+1.0);
        double dh_y1 = ((h_y1 + GWS*z_y1) - (H+GWS*Z))/(GWS+1.0);
        double dh_y2 = ((h_y2 + GWS*z_y2) - (H+GWS*Z))/(GWS+1.0);

        double dz_x1 = z_x1 -Z;
        double dz_x2 = z_x2 -Z;
        double dz_y1 = z_y1 -Z;
        double dz_y2 = z_y2 -Z;

        // flow = Ksat * cross section * hydraulic gradient Ks * A * dH/dL
        double ff = GW_flow * factor;

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
        double maxvol = CHAdjDX->Drc * SD->Drc * pore->Drc;
        if (V + dflux > maxvol)
            dflux =  maxvol - V;
        if (V + dflux < 0)
            dflux = -V;
        //fill with the resulting flux of a cell
        GWout->Drc = dflux;//std::max(0.0,dflux);
    }}

    // adjust the vol
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        GWVol->Drc += GWout->Drc;
        // update gwvol with flux
        GWWH->Drc = GWVol->Drc/CHAdjDX->Drc/pore->Drc;
        // recalc gwwh
        //GWout->Drc = fabs(GWout->Drc);
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
        tmb->Drc = 0;
        tmc->Drc = 0;
        double GWout_ = GW_flow *  CHAdjDX->Drc * std::max(0.0, GWWH->Drc-GW_threshold) * ksat->Drc * BaseflowL->Drc; // m3 volume out from every cell
        GWout_ *= (1+Grad->Drc);

        //  GWout_ *= (1-exp(-GW_threshold*GWWH->Drc));
        //m3:  ksat*dt  * dh*dx * ((dx/L)^b);  ksat * cross section * distance factor
        // stop outflow when some minimum GW level, 2.4.2.10 in SWAT
        // apply a smooth threshold with exponential function
        GWout_ = std::min(GWVol->Drc*MaxGWDepthfrac, GWout_);
        GWout_ = ChannelWidth->Drc > 0 ? 0.0 : GWout_; // set GWout in channel cell to zero else accumulation to the outlet
        tmb->Drc = GWout_;

        // adjust volume with outflow
        if (GWout_ > 0) {
           GWVol->Drc -= GWout_; // subtract from volume
           GWWH->Drc = GWVol->Drc/ CHAdjDX->Drc/pore->Drc;
           doit = true; // is someqwhere GWout > 0 do accuflux
        }
    }}

    if (doit)
        AccufluxGW(crlinkedlddbase_, tmb, tmc, ChannelWidth);
    // tmc has now the accumulated flow pattern, do NOT use this anymore for the mass balance
    // channelwidth is flag: stop accumulating when you reach channel

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        GWout->Drc = tmc->Drc;
    }}

}

double TWorld::fullSWOF2GW(cTMap *h, cTMap *u, cTMap *v, cTMap *z)
{
    double timesum = 0;
    double dt_max = std::min(_dt, _dx/2);
    int count = 0;
    double sumh = 0;
    bool stop;
    double dt_req_min = dt_max;
    int step = 0;

    cTMap *pore;
    cTMap *ksat;
    cTMap *SD;
    if (SwitchTwoLayer) {
        pore = ThetaS2;
        ksat = Ksat2;   // has the unit m, *_dt time is included
        SD = SoilDepth2init;
    } else {
        pore = Poreeff;
        ksat = Ksateff;
        SD = SoilDepth1init;
    }

    sumh = getMass(h, 0);

    if (sumh == 0)
        return 0;

    do {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            hs->Drc = h->Drc;
            tmc->Drc = dt_max;
            tmd->Drc = 0;

            if (hs->Drc > 0.001) {
                tmd->Drc = 1;
                if (c > 0 && !MV(r,c-1)        ) tmd->data[r][c-1] = 1;
                if (c < _nrCols-1 && !MV(r,c+1)) tmd->data[r][c+1] = 1;
                if (r > 0 && !MV(r-1,c)        ) tmd->data[r-1][c] = 1;
                if (r < _nrRows-1 && !MV(r+1,c)) tmd->data[r+1][c] = 1;

                if (c > 0 && r > 0 && !MV(r-1,c-1)                ) tmd->data[r-1][c-1]=1;
                if (c < _nrCols-1 && r < _nrRows-1 && !MV(r+1,c+1)) tmd->data[r+1][c+1]=1;
                if (r > 0 && c < _nrCols-1 && !MV(r-1,c+1)        ) tmd->data[r-1][c+1]=1;
                if (c > 0 && r < _nrRows-1 && !MV(r+1,c-1)        ) tmd->data[r+1][c-1]=1;
            }
        }}

        //do all flow and state calculations
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {

            if (tmd->Drc > 0) {
                    //double dt = FloodDT->Drc; //dt_req_min;
                double dt = dt_req_min;
                double vxn, vyn;
                double Ks = GW_flow*ksat->Drc/_dt;

                //FloodT->Drc += FloodDT->Drc;

                vec4 hll_x1;
                vec4 hll_x2;
                vec4 hll_y1;
                vec4 hll_y2;

                double dx = _dx;
                double dy = _dx;

                double H = hs->Drc;
                double n = GWN->Drc;
                double Z = z->Drc;
                double Vx = GW_flow*u->Drc;
                double Vy = GW_flow*v->Drc;

                bool bc1 = c > 0 && !MV(r,c-1)        ;
                bool bc2 = c < _nrCols-1 && !MV(r,c+1);
                bool br1 = r > 0 && !MV(r-1,c)        ;
                bool br2 = r < _nrRows-1 && !MV(r+1,c);

                double z_x1 =  bc1 ? z->data[r][c-1] : Z;
                double z_x2 =  bc2 ? z->data[r][c+1] : Z;
                double z_y1 =  br1 ? z->data[r-1][c] : Z;
                double z_y2 =  br2 ? z->data[r+1][c] : Z;

                double h_x1 =  bc1 ? hs->data[r][c-1] : H;
                double h_x2 =  bc2 ? hs->data[r][c+1] : H;
                double h_y1 =  br1 ? hs->data[r-1][c] : H;
                double h_y2 =  br2 ? hs->data[r+1][c] : H;

                double vx_x1 = bc1 ? GW_flow*u->data[r][c-1] : Vx;
                double vx_x2 = bc2 ? GW_flow*u->data[r][c+1] : Vx;
                double vx_y1 = br1 ? GW_flow*u->data[r-1][c] : Vx;
                double vx_y2 = br2 ? GW_flow*u->data[r+1][c] : Vx;

                double vy_x1 = bc1 ? GW_flow*v->data[r][c-1] : Vy;
                double vy_x2 = bc2 ? GW_flow*v->data[r][c+1] : Vy;
                double vy_y1 = br1 ? GW_flow*v->data[r-1][c] : Vy;
                double vy_y2 = br2 ? GW_flow*v->data[r+1][c] : Vy;

                double dz_x1 = (Z - z_x1);
                double dz_x2 = (z_x2 - Z);
                double dz_y1 = (Z - z_y1);
                double dz_y2 = (z_y2 - Z);

                 vx_x1 = std::min(Ks,fabs(vx_x1))*(vx_x1 < 0 ? -1.0 : 1.0);
                 vx_x2 = std::min(Ks,fabs(vx_x2))*(vx_x2 < 0 ? -1.0 : 1.0);
                 vx_y1 = std::min(Ks,fabs(vx_y1))*(vx_y1 < 0 ? -1.0 : 1.0);
                 vx_y2 = std::min(Ks,fabs(vx_y2))*(vx_y2 < 0 ? -1.0 : 1.0);
                 vy_x1 = std::min(Ks,fabs(vy_x1))*(vy_x1 < 0 ? -1.0 : 1.0);
                 vy_x2 = std::min(Ks,fabs(vy_x2))*(vy_x2 < 0 ? -1.0 : 1.0);
                 vy_y1 = std::min(Ks,fabs(vy_y1))*(vy_y1 < 0 ? -1.0 : 1.0);
                 vy_y2 = std::min(Ks,fabs(vy_y2))*(vy_y2 < 0 ? -1.0 : 1.0);
                 Vx = std::min(Ks,fabs(Vx))*(Vx < 0 ? -1.0 : 1.0);
                 Vy = std::min(Ks,fabs(Vy))*(Vy < 0 ? -1.0 : 1.0);


                // z is blocking to prevent flow when water is flat and Z is not flat, described in article SWOF
                double h_x1r = std::max(0.0, h_x1 - std::max(0.0,  dz_x1));
                double H_l   = std::max(0.0, H    - std::max(0.0, -dz_x1));
                if(bc1)  // if inside
                    hll_x1 = F_Rusanov(h_x1r,vx_x1,vy_x1, H_l,Vx,Vy); // c-1 and c  //
                else
                    hll_x1 = F_Rusanov(0,0,0, H_l,Vx,Vy);

                double H_r   = std::max(0.0, H    - std::max(0.0,  dz_x2));
                double h_x2l = std::max(0.0, h_x2 - std::max(0.0, -dz_x2));
                if(bc2)
                    hll_x2 = F_Rusanov(H_r,Vx,Vy, h_x2l,vx_x2,vy_x2); // c and c+1
                else
                    hll_x2 = F_Rusanov(H_r,Vx,Vy, 0,0,0);

                double h_y1d = std::max(0.0, h_y1 - std::max(0.0,  dz_y1));
                double H_u   = std::max(0.0, H    - std::max(0.0, -dz_y1));
                if (br1)
                    hll_y1 = F_Rusanov(h_y1d,vy_y1,vx_y1, H_u,Vy,Vx); // r-1 and r
                else
                    hll_y1 = F_Rusanov(0,0,0, H_u,Vy,Vx);

                double H_d   = std::max(0.0, H    - std::max(0.0,  dz_y2));
                double h_y2u = std::max(0.0, h_y2 - std::max(0.0, -dz_y2));
                if(br2)
                    hll_y2 = F_Rusanov(H_d,Vy,Vx, h_y2u,vy_y2,vx_y2); // r and r+1
                else
                    hll_y2 = F_Rusanov(H_d,Vy,Vx, 0,0,0);

                // determine smallest dt in x and y for each cell
                double dtx = dx/std::max(hll_x1.v[3],hll_x2.v[3]);
                double dty = dy/std::max(hll_y1.v[3],hll_y2.v[3]); // v[3] is max U and V in x and y

                double dt_req = std::max(TimestepfloodMin, std::min(dt_max, courant_factor*std::min(dtx, dty)));
                tmc->Drc = dt_req; // dt does not need to be a map, left over from earlier code
                // if step = 0 do not calculate new fluxes and states yet because the first dt is always dt_max
                // find a smallest dt of the flow domain first

                //########### after finding the smallest dt, do the st venant eq)
                if (step > 0) {

                    double tx = dt/dx;
                    double ty = dt/dy;

                    double hn = std::max(0.0, H + dt/_dx*(hll_x1.v[0]-hll_x2.v[0] + hll_y1.v[0]-hll_y2.v[0]));
                    // mass balance, hll_....v[0] is the height

                    // momentum balance for cells with water
                    if(hn > 1e-6) {
                        double gflow_x = GRAV*0.5*( (H_l-H)*(H_l+H)+(H-H_r)*(H+H_r));
                        double gflow_y = GRAV*0.5*( (H_u-H)*(H_u+H)+(H-H_d)*(H+H_d));
                        // graviy term: gh

                        //double qxn = H * Vx - tx*(gflow_x);
                        //double qyn = H * Vy - ty*(gflow_y);
                        double qxn = H * Vx - tx*(hll_x2.v[1] - hll_x1.v[1] + gflow_x) - ty*(hll_y2.v[2] - hll_y1.v[2]);
                        double qyn = H * Vy - tx*(hll_x2.v[2] - hll_x1.v[2]) - ty*(hll_y2.v[1] - hll_y1.v[1] + gflow_y);

                        double vsq = sqrt(Vx * Vx + Vy * Vy);
                        double nsq1 = n*n*GRAV/pow(hn,4.0/3.0);//
                        double nsq = nsq1*vsq*dt;

                        vxn = (qxn/(1.0+nsq))/hn;
                        vyn = (qyn/(1.0+nsq))/hn;
                    } else { // hn < ha
                        hn = H; // if no fluxes then also no change in h
                        vxn = 0;
                        vyn = 0;
                    }

                    // dan maar even met geweld!
                    if (std::isnan(vxn) || std::isnan(vyn)  )
                    {
                        vxn = 0;
                        vyn = 0;
                    }
                    if (FlowBoundaryType == 0 || (FlowBoundaryType == 2 && FlowBoundary->Drc == 0)) {

                        if (DomainEdge->Drc == 4 && vxn < 0) {
                            vxn = 0;
                        }
                        if (DomainEdge->Drc == 6 && vxn > 0) {
                            vxn = 0;
                        }
                        if (DomainEdge->Drc == 2 && vyn > 0) {
                            vyn = 0;
                        }
                        if (DomainEdge->Drc == 8 && vyn < 0) {
                            vyn = 0;
                        }

                    }
                    if (vyn == 0 && vxn == 0)
                        hn = H;

                    vxn = std::min(Ks,fabs(vxn))*(vxn < 0 ? -1.0 : 1.0);
                    vyn = std::min(Ks,fabs(vyn))*(vyn < 0 ? -1.0 : 1.0);

                    h->Drc = hn;
                    u->Drc = vxn;
                    v->Drc = vyn;

                } // step > 0
            } // tmd > 0, active cells + 1
        }}

        // find smallest domain dt
        #pragma omp parallel for reduction(min:dt_req_min) num_threads(userCores)
        FOR_ROW_COL_MV_L {
            dt_req_min = std::min(dt_req_min, tmc->Drc);
        }}

        dt_req_min = std::min(dt_req_min, _dt-timesum);

        if (step > 0) {
            timesum += dt_req_min;
            count++; // nr loops
        }

        step += 1; // now we have a good dt min, do the real calculations

        stop = timesum > _dt-0.001;
        if(count > F_MaxIter)
        stop = true;

    } while (!stop);

    correctMassBalance(sumh, h, 0);

    //qDebug() << _dt/count << count << dt_req_min;
    iter_n = std::max(1,count);
    return(count > 0 ? _dt/count : _dt);
}
