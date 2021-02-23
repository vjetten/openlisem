#include "lisemqt.h"
#include "model.h"



void TWorld::do_SplashDetachment(int r, int c, double WH)
{
    DETSplash->Drc = 0;
    if(WH > 0.0001)// || hmx->Drc > 0.0001)
    {
        double DetDT1 = 0, DetDT2 = 0, DetLD1, DetLD2;
        double g_to_kg = 0.001;
        double Lc = Litter->Drc;
        double Cv = Cover->Drc;
        double strength = AggrStab->Drc;
        double b = 0;//splashb->Drc;
        double Int = Rain->Drc * 3600/_dt * 1000; // intensity in mm/h, Rain is in m
        double KE_DT = 0.0;
        double DETSplash_;

        switch (KEequationType)
        {
        case KE_EXPFUNCTION: KE_DT = KEParamater_a1*(1-(KEParamater_b1*exp(-KEParamater_c1*Int))); break;
        case KE_LOGFUNCTION: KE_DT = (Int > 1 ? KEParamater_a2 + KEParamater_b2*log10(Int) : 0); break;
        case KE_POWERFUNCTION: KE_DT = KEParamater_a3*pow(Int, KEParamater_b3); break;
            // kin energy in J/m2/mm
        }
        //VJ 110706  KE equations

        double directrain = (1-Lc) * (1-Cv)*Rainc->Drc * 1000;
        // Added litter also to directrain, assme it covers the entire cell, not only under the plant
        // rainfall between plants in mm

        double KE_LD = std::max(15.3*sqrt(PlantHeight->Drc)-5.87, 0.0);
        // kin energy in J/m2/mm
        double throughfall = (1-Lc) * Cv * LeafDrain->Drc * 1000;
        // leaf drip in mm, is calculated as plant leaf drip in interception function so mult cover
        // VJ 110206 stemflow is also accounted for

        double WH0 = exp(-1.48*hmxWH->Drc*1000);
        // water buffer effect on surface, WH in mm in this empirical equation from Torri ?

        if(SwitchUseMaterialDepth)
        {
            double depdepth = std::max((StorageDep->Drc / BulkDens)/(_dx * DX->Drc),0.0);
            double fac1 = std::max(0.0,1.0 - depdepth/(SedimentMixingDepth->Drc+0.01));
            double fac2 = 1.0 - fac1;

            strength = strength * fac2 + (0.1033/DepositedCohesion) * fac1;
            b = b * fac2 + 3.58 * fac1;
        }

        double FPA = 1.0;
        if (RR->Drc > 0.1)
            FPA =  1-exp(-1.875*(WH/(0.01*RR->Drc)));

        // Between plants, directrain is already with 1-cover
        DetDT1 = g_to_kg * FPA*(strength*KE_DT+b)*WH0 * directrain;
        //ponded areas, kg/m2/mm * mm = kg/m2
        DetDT2 = hmxWH->Drc > 0 ? g_to_kg * (1-FPA)*(strength*KE_DT+b) * directrain * SplashDelivery : 0.0;
        //dry areas, kg/m2/mm * mm = kg/m2

        if (SwitchKETimebased)
        {
            if (directrain > 0)
            {
                DetDT1 = g_to_kg * FPA*(strength*KE_DT+b)*WH0 * _dt/3600;
                //ponded areas, kg/m2/sec * sec = kg/m2
                DetDT2 = g_to_kg * (1-FPA)*(strength*KE_DT+b) * _dt/3600 * SplashDelivery;
                //dry areas, kg/m2/sec * sec = kg/m2
            }
        }
        //based on work by Juan Sanchez

        // Under plants, throughfall is already with cover
        DetLD1 = g_to_kg * FPA*(strength*KE_LD+b)*WH0 * throughfall;
        //ponded areas, kg/m2/mm * mm = kg/m2
        DetLD2 = g_to_kg * (1-FPA)*(strength*KE_LD+b) * throughfall * SplashDelivery;
        //dry areas, kg/m2/mm * mm = kg/m2

        DETSplash_ = DetLD1 + DetLD2 + DetDT1 + DetDT2;
        // Total splash kg/m2

        // Deal with all exceptions:

        DETSplash_ *= (SoilWidthDX->Drc*DX->Drc);
        // kg/cell, only splash over soilwidth, not roads and channels
        // FROM KG/M2 TO KG/CELL

        DETSplash_ = (1-StoneFraction->Drc) * DETSplash_;
        // no splash on stone surfaces

        if (SwitchGrassStrip)
            DETSplash_ = (1-GrassFraction->Drc) * DETSplash_;

        //      if(SwitchSedtrap)
        //          DETSplash->Drc = (1-SedimentFilter->Drc) * DETSplash->Drc;

        if (SwitchHardsurface)
            DETSplash_ = (1-HardSurface->Drc)*DETSplash_;
        // no splash on hard surfaces

        if (SwitchHouses)
            DETSplash_ = (1-HouseCover->Drc)*DETSplash_;
        //is already contained in soilwidth
        // no splash from house roofs

        if (SwitchSnowmelt)
            DETSplash_ = (1-Snowcover->Drc)*DETSplash_;
        // no splash on snow deck

        if(SwitchUseMaterialDepth)
        {
            //check wat we can detach from the top and bottom layer of present material
            double dleft = DETSplash_;
            double deptake = 0;
            double mattake = 0;
            double detachment = 0;

            deptake = std::min(dleft,StorageDep->Drc);
            StorageDep->Drc -= deptake;
            // det not more than storage
            // decrease store depth

            detachment += deptake;
            // detachment is now taken material


            if(!(Storage->Drc < 0))
            {
                mattake = std::min(dleft,Storage->Drc);
                Storage->Drc -= mattake;

                detachment += mattake;
            }else
            {
                detachment += dleft;
            }
            DETSplash_ = detachment;
        }


        if (hmx->Drc > 0) {
            SSFlood->Drc += DETSplash_;
            SSCFlood->Drc = MaxConcentration(CHAdjDX->Drc * hmx->Drc, &SSFlood->Drc, &DepFlood->Drc);

        } else {
            Sed->Drc += DETSplash_;
            Conc->Drc = MaxConcentration(WaterVolall->Drc, &Sed->Drc, &DEP->Drc);
        }

        DETSplash->Drc = DETSplash_;
        // IN KG/CELL
    }
}

//---------------------------------------------------------------------------
void TWorld::do_Interception(int r, int c)
{
    // all variables are in m
    double Cv = Cover->Drc;
    double Rainc_ = Rainc->Drc;
    double RainCum_ = RainCum->Drc;
    double AreaSoil = SoilWidthDX->Drc * DX->Drc;
    double RainNet_ = Rainc_;

    if (Cv > 0)
    {
        double CS = CStor->Drc;
        //actual canopy storage in m
        double Smax = CanopyStorage->Drc;
        //max canopy storage in m

        if (Smax > 0) {
            CS = Smax*(1-exp(-kLAI->Drc*RainCum_/Smax));
        }
        //else
        //  CS = 0;

        LeafDrain->Drc = std::max(0.0, Cv*(Rainc_ - (CS - CStor->Drc)));
        // diff between new and old strage is subtracted from rainfall
        // rest reaches the soil surface. ASSUMPTION: with the same intensity as the rainfall!
        // note: cover already implicit in LAI and Smax, part falling on LAI is cover*rainfall

        CStor->Drc = CS;
        // put new storage back in map
        Interc->Drc =  Cv * CS * AreaSoil;
        // Interc->Drc =  Cv * CS * _dx * DX->Drc;
        // WHY: cvover already takes care of this, trees can be above a road or channel

        RainNet_ = LeafDrain->Drc + (1-Cv)*Rainc_;
        // net rainfall is direct rainfall + drainage
        // rainfall that falls on the soil, used in infiltration
    }

    if (SwitchLitter) {
        double CvL = Litter->Drc;
        if (hmx->Drc == 0 && WH->Drc == 0 && CvL > 0 && RainNet_ > 0)
        {

            double Smax = LitterSmax/1000.0;
            // assume simply that the cover linearly scales between 0 and LtterSmax of storage

            double LCS = LCStor->Drc;
            //actual canopy storage in m

            LCS = std::min(LCS + RainNet_, Smax);
            // add water to the storage, not more than max

            double drain = std::max(0.0, CvL*(RainNet_ - (LCS - LCStor->Drc)));
            // diff between new and old storage is subtracted from leafdrip

            LCStor->Drc = LCS;
            // put new storage back in map for next dt

            LInterc->Drc =  CvL * LCS * AreaSoil;
            // only on soil surface, not channels or roads, in m3

            RainNet_ = drain + (1-CvL)*RainNet_;
            //recalc
        }
    }

    // all variables are in m
    if (SwitchHouses)
    {
        double CvH = HouseCover->Drc;
        if (CvH > 0 &&  RainNet_ > 0)
        {
            //house on cell in m2
            double HS = HStor->Drc;
            //actual roof storage in m

            double Hmax = RoofStore->Drc;
            //max roof storage in m

            HS = std::min(HS + RainNet_, Hmax);

            double housedrain = std::max(0.0, CvH * (RainNet_ - (HS - HStor->Drc)));
            // overflow in m3/m2 of house

            HStor->Drc = HS;
            // put new storage back in maps in m

            double roofsurface = (AreaSoil * CvH); // m2
            // double roofsurface = (_dx * DX->Drc * CvH); // m2
            // user should assure housecover is correct with respect to channel and roads
            IntercHouse->Drc =  roofsurface * HS;
            // total interception in m3,exclude roads, channels

            RainNet_ = housedrain + (1-CvH)*RainNet_;
            // net rainfall is direct rainfall + drainage
            // rainfall that falls on the soil, used in infiltration

            // filling raindrums with surplus drainage from roofs
            // drum is recalculated to m based on roof surface
            double DS = 0;
            if (SwitchRaindrum && DrumStore->Drc > 0)
            {
                double Dmax = DrumStore->Drc/roofsurface;
                //max drum storage in m as if roof storage is more

                DS = DStor->Drc;
                //actual drum storage in m

                DS = std::min(DS + RainNet_, Dmax);
                // fill tank to max
                double drumdrain = std::max(0.0, HouseCover->Drc * (RainNet_ - (DS - DStor->Drc)));

                DStor->Drc = DS;
                // put new drum storage back in maps in m3

                IntercHouse->Drc += roofsurface * DS;
                // total interception in m3,exclude roads, channels

                RainNet_ = drumdrain + (1-CvH)*RainNet_;
            }
        }
    }
    RainNet->Drc = RainNet_;
}


void TWorld::do_InfiltrationGA(int r, int c, double fwh, double SW, double flooddomain)
{
    // default vars are first layer vars
    double fact_ = fact->Drc;
    double Ks = Ksateff->Drc*_dt/3600000.0;  //in m
    double Psi = Psi1->Drc/100; // in m
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;
    double SoilDep2;
    double Ks2;
    double fpot = 0;

    //calculate potential insiltration rate fpot
    if (SwitchTwoLayer) {
        SoilDep2 = SoilDepth2->Drc;
        Ks2 = Ksat2->Drc*_dt/3600000.0;
        // if wetting front in second layer set those vars
        if (Lw_ > SoilDep1)
        {
            Ks = std::min(Ksateff->Drc, Ks2); // !!! was wrong because _dt/3600000.0 was for both
            // if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2
            Psi = Psi2->Drc/100;
        }
    }

    if (InfilMethod == INFIL_GREENAMPT || InfilMethod == INFIL_GREENAMPT2)
        fpot = Ks*(1.0+(Psi+fwh)/std::max(1e-4, Lw_));
    else {
        double space = SwitchTwoLayer ? std::max(ThetaS2->Drc-ThetaI2->Drc, 0.0) : std::max(Poreeff->Drc-Thetaeff->Drc, 0.0);
        double B = (fwh + Psi)*space;
        if (B > 0.01) {
            fpot = Ks*exp(Fcum->Drc/B)/(exp(Fcum->Drc/B)-1);
        } else
            fpot = Ks;
    }

    fact_ = std::min(fpot, fwh);
    // actual infil in m, cannot have more infil than water on the surface

    fact_ = IncreaseInfiltrationDepthNew(fact_, Lw_, r, c);
    // adjust fact and increase Lw, for twolayer, impermeable etc
    Lw_ = Lw->Drc;

    if (fwh < fact_)
    {
        fact_ = fwh;
        fwh = 0;
    }
    else
        fwh -= fact_;

    // adjust the WH in the correct domain with new fact
    if(flooddomain == 0)
        WH->Drc = fwh;
    else
        hmx->Drc = fwh;

    Fcum->Drc += fact_;  // increase cumulative infil in m for Smith and Parlange

    InfilVol->Drc = fact_*SW*DX->Drc;
    // calc infiltrated volume for mass balance

    // calc surplus infiltration (negative in m) for kin wave
    if(SwitchKinematic2D != K2D_METHOD_DYN) {
        double space = 0;
        if (Lw_ < SoilDep1)
            space = (SoilDep1 - Lw->Drc)*(Poreeff->Drc-Thetaeff->Drc);
        if (SwitchTwoLayer) {
            if (Lw_ > SoilDep1 && Lw_ < SoilDep2)
                space = (SoilDep2 - Lw_)*(ThetaS2->Drc-ThetaI2->Drc);
        }

        FSurplus->Drc = -1.0 * std::min(space, fact_);//std::max(0.0, fpot_-fact_));
        // negative and smallest of space or fpot-fact ???
    }
}

void TWorld::do_Percolation(int r, int c)
{
    double Percolation, dL, pore, theta, thetar, theta_E;
    Percolation = 0;
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;

    if(SwitchTwoLayer) {
        pore = ThetaS2->Drc;
        thetar = 0.025 * pore;
        theta = ThetaI2->Drc;
        double SoilDep2 = SoilDepth2->Drc;

        if(theta > thetar) {
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = Ksat2->Drc * pow(theta_E, bca->Drc);
            // percolation in m

            if (Lw->Drc > SoilDepth1->Drc)
                dL = SoilDep2 - Lw_;
            else
                dL = SoilDep2 - SoilDep1;
            // if Wet Fr still in first layer percolation only make 2nd drier

            double moisture = dL*(theta-thetar);

            if (moisture > Percolation) {
                // decrease thetaeff because of percolation
                moisture -= Percolation;
                theta = moisture/dL+thetar;
            } else {
                // wetting front = soildepth1, dL = 0, moisture = 0
                // assume theta goes back to 0.7 pore and decrease the wetting fornt
                theta = 0.7*(pore - thetar);
                Lw_ -= std::max(0.0, Percolation/(pore - theta));
            }
            ThetaI2->Drc = theta;
        }
    } else {
        // one layer
        double pore = Poreeff->Drc;
        thetar = 0.025 * pore;
        double theta = Thetaeff->Drc;

        if(theta > thetar) {
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = Ksateff->Drc*_dt/3600000.0 * pow(theta_E, bca->Drc);
        }

        if (Percolation > 0) {
            dL = std::max(0.0, SoilDep1 - Lw_);
            double moisture = dL*(theta-thetar);
            if (moisture > Percolation) {
                // wetting front has not reached bottom, make soil drier
                // decrease thetaeff because of percolation
                moisture -= Percolation;
                theta = moisture/dL+thetar;
            } else {
                // wetting front = soildepth1, dL = 0, moisture = 0
                // assume tehta goes back to half pore and decrease the wetting fornt
                theta = 0.7*(pore - thetar);
                Lw_ -= std::max(0.0, Percolation/(pore - theta));
            }
            Thetaeff->Drc = theta;
        }
    }

    if (Percolation > 0) {
        double moisture = dL*(theta-thetar);
        if (moisture > Percolation) {
            // wetting front has not reached bottom, make soil drier
            // decrease thetaeff because of percolation
            moisture -= Percolation;
            theta = moisture/dL+thetar;
        } else {
            // wetting front = soildepth1, dL = 0, moisture = 0
            // assume tehta goes back to half pore and decrease the wetting fornt
            theta = 0.7*(pore - thetar);
            Lw_ -= std::max(0.0, Percolation/(pore - theta));
        }
    }

    Lw->Drc = Lw_;
    Perc->Drc = Percolation;
}


void TWorld::do_CellProcesses()
{
      RainfallMap();         // get rainfall from table or mpas


#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {

        if (SwitchRainfall) {
//            double timeminprev = (time-_dt) / 60; //prev time in minutes
//            int  rainplace;
//            double tt = 3600000.0;

//            for (rainplace = 0; rainplace < nrRainfallseries; rainplace++)
//                if (timeminprev < RainfallSeriesM[rainplace].time)
//                    break;

//            if (RainfallSeriesM[rainplace].isMap)
//            {
//                auto _M = std::unique_ptr<cTMap>(new cTMap(readRaster(RainfallSeriesM[rainplace].name)));
//                if (pcr::isMV(_M->Drc))
//                {
//                    QString sr, sc;
//                    sr.setNum(r); sc.setNum(c);
//                    ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+RainfallSeriesM[rainplace].name;
//                    throw 1;
//                }
//                else
//                    Rain->Drc = _M->Drc *_dt/tt;
//            }
//            else
//            {
//                Rain->Drc = RainfallSeriesM[rainplace].intensity[(int) RainZone->Drc-1]*_dt/tt;
//            }

//            if (!rainStarted) {
//                if(Rain->Drc > 0)
//                    rainStarted = true;
//            }
//            if (rainStarted && RainstartTime == -1)
//                RainstartTime = time;

//            Rainc->Drc = Rain->Drc * _dx/DX->Drc;
//            // correction for slope dx/DX, water spreads out over larger area
//            RainCumFlat->Drc += Rain->Drc;
//            // cumulative rainfall
//            RainCum->Drc += Rainc->Drc;
//            // cumulative rainfall corrected for slope, used in interception
//            RainNet->Drc = Rainc->Drc;
//            // net rainfall in case of interception

            if (Rainc->Drc > 0)
                do_Interception(r,c);
        }

        //========================

        double FloodDomain_ = FloodDomain->Drc;
        if (FloodDomain_ == 0) {
            WH->Drc += RainNet->Drc + Snowmeltc->Drc;
        } else {
            hmx->Drc += RainNet->Drc + Snowmeltc->Drc;
        }
        double RW = RoadWidthHSDX->Drc;
        if (SwitchRoadsystem && RW > 0) {
                WHroad->Drc += Rainc->Drc + Snowmeltc->Drc;
        }

        //========================

        if (InfilMethod != INFIL_NONE && InfilMethod != INFIL_SWATRE) {

          //  do_InfiltrationGA(r,c, fwh, SW, FloodDomain_);
            InfilMethodsNew(r, c);
            if (!SwitchImpermeable)
                do_Percolation(r, c);
        }

        //===== SURFACE STORAGE =====

        double SW = SoilWidthDX->Drc;
        double WHr = WHroad->Drc;
        double WHs; //WHstore
        double WH_ = WH->Drc;

        //### surface storage on rough surfaces
        WHs = std::min(WH_, MDS->Drc*(1-exp(-1.875*(WH_/std::max(0.01,0.01*RR->Drc)))));
        // non-linear release fo water from depression storage
        // resembles curves from GIS surface tests, unpublished

        double FW = std::min(ChannelAdj->Drc, SW + RW);
        // calculate flowwidth by fpa*surface + road, excludes channel already

        WHrunoff->Drc = ((WH_ - WHs)*SW + WHr*RW)/FW;
        FlowWidth->Drc = FW;

        WaterVolall->Drc = DX->Drc*(WH_*SW + WHr*RW);
        WHstore->Drc = WHs;

        //========================

        if (SwitchErosion) {
            double wh = FloodDomain_ == 0 ? WH->Drc : hmx->Drc;
            do_SplashDetachment(r,c,wh);
        }




    }}
}

