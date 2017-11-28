
/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/


/*!
 \file lisUnifiedFlow.cpp
 \brief

functions: \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"

//connection
void TWorld::UF2D1D_Connection(int thread,cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                               cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                               cTMap * _fu1D,cTMap * _s1D,
                               cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                               cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                               cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                               cTMap * _su2D,cTMap * _sv2D)
{
    FOR_ROW_COL_UF2DMT
    {
        if(!UF_OUTORMV(_ldd,r,c))
        {
            double h = (_f1D->Drc + _s1D->Drc)/(_dx*_lddw->Drc);

            double fV = sqrt(_fu2D->Drc * _fu2D->Drc + _fv2D->Drc * _fv2D->Drc);
            double sV = sqrt(_su2D->Drc * _su2D->Drc + _sv2D->Drc * _sv2D->Drc);

            if(UF_CHANNELFLOOD)
            {
                //inflow
                if(h < _lddh->Drc)
                {

                    double maxvol= _lddw->Drc * _dx * (_lddh->Drc - h);


                    double fractionf = std::max(0.0,std::min(0.95, dt->Drc*fV/std::max(0.01*_dx,0.5*(_dx - _lddw->Drc))));

                    double fractions = std::max(0.0,std::min(0.95, dt->Drc*sV/std::max(0.01*_dx,0.5*(_dx - _lddw->Drc))));

                    double volf = _f2D->Drc * (fractionf);
                    double vols = _s2D->Drc * (fractions);

                    if(volf+vols > maxvol)
                    {
                        double totalvol = (volf+vols);
                        volf = volf * maxvol/totalvol;
                        vols = vols * maxvol/totalvol;
                    }

                    double vf = _f1D->Drc + volf;
                    double vs = _s1D->Drc + vols;

                    //_fu1D->Drc = vf > UF_VERY_SMALL? (_fu1D->Drc *_f1D->Drc + fV * _f2D->Drc * (fractionf))/vf: 0.0;
                    //_su1D->Drc = vs > UF_VERY_SMALL? (_su1D->Drc *_s1D->Drc + sV * _s2D->Drc * (fractions))/vs: 0.0;
                    //_visc1D->Drc = vf > UF_VERY_SMALL? (_visc1D->Drc *_f1D->Drc + _visc2D->Drc * _f2D->Drc * (fractionf))/vf: 0.0;
                    //_d1D->Drc = vs > UF_VERY_SMALL? (_d1D->Drc *_s1D->Drc + _d2D->Drc * vols)/vs: 0.0;
                    _ifa1D->Drc = vs > UF_VERY_SMALL? (_ifa1D->Drc *_s1D->Drc + _ifa2D->Drc * vols)/vs: 0.0;
                    _rocksize1D->Drc = vs > UF_VERY_SMALL? (_rocksize1D->Drc *_s1D->Drc + _rocksize2D->Drc * vols)/vs: 0.0;

                    _f1D->Drc = vf;
                    _s1D->Drc = vs;

                    _f2D->Drc = _f2D->Drc - volf;
                    _s2D->Drc = _s2D->Drc - vols;

                    if(SwitchErosion)
                    {
                        double blconc = _f2D > 0? UF2D_blm->Drc/ _f2D->Drc : 0.0;
                        double ssconc = _f2D > 0? UF2D_ssm->Drc/ _f2D->Drc : 0.0;
                        double qbls = std::min(UF2D_blm->Drc,blconc * volf);
                        double qsss = std::min(UF2D_ssm->Drc,ssconc * volf);
                        UF1D_blm->Drc += qbls;
                        UF2D_blm->Drc -= qbls;
                        UF1D_ssm->Drc += qsss;
                        UF2D_ssm->Drc -= qsss;

                        if(SwitchUseGrainSizeDistribution)
                        {
                            FOR_GRAIN_CLASSES
                            {
                                double blconc2 = _f2D > 0? UF2D_blm_D.Drcd / _f2D->Drc :0.0;
                                double ssconc2 = _f2D > 0? UF2D_ssm_D.Drcd / _f2D->Drc : 0.0;
                                double qbls2 = std::min(UF2D_blm_D.Drcd,blconc2 * volf);
                                double qsss2 = std::min(UF2D_ssm_D.Drcd,ssconc2 * volf);
                                UF1D_blm_D.Drcd += qbls2;
                                UF2D_blm_D.Drcd -= qbls2;
                                UF1D_ssm_D.Drcd += qsss2;
                                UF2D_ssm_D.Drcd -= qsss2;
                            }
                        }
                    }

                //outflow
                }else
                {

                    double vol = h * _lddw->Drc * _dx;
                    double extravol = _lddw->Drc * _dx * (h-_lddh->Drc);
                    double fractionf = (extravol/vol) ;
                    double fractions = (extravol/vol) ;

                    double volf = fractionf * _f1D->Drc;
                    double vols = fractions * _s1D->Drc;

                    double vfn = _f2D->Drc + volf;
                    double vsn = _s2D->Drc + vols;

                    _d2D->Drc = vsn > UF_VERY_SMALL? (_d2D->Drc *_s2D->Drc + _d1D->Drc * vols)/vsn: 0.0;
                    _ifa2D->Drc = vsn > UF_VERY_SMALL? (_ifa2D->Drc *_s2D->Drc + _ifa1D->Drc * vols)/vsn: 0.0;
                    _rocksize2D->Drc = vsn > UF_VERY_SMALL? (_rocksize2D->Drc *_s2D->Drc + _rocksize1D->Drc * vols)/vsn: 0.0;

                    _f1D->Drc = _f1D->Drc - volf;
                    _s1D->Drc = _s1D->Drc - vols;

                    _f2D->Drc = vfn;
                    _s2D->Drc = vsn;

                    if(SwitchErosion)
                    {
                        double blconc = _f1D > 0? UF1D_blm->Drc/ _f1D->Drc : 0.0;
                        double ssconc = _f1D > 0? UF1D_ssm->Drc/ _f1D->Drc : 0.0;
                        double qbls = std::min(UF1D_blm->Drc,blconc * volf);
                        double qsss = std::min(UF1D_ssm->Drc,ssconc * volf);
                        UF1D_blm->Drc -= qbls;
                        UF2D_blm->Drc += qbls;
                        UF1D_ssm->Drc -= qsss;
                        UF2D_ssm->Drc += qsss;

                        if(SwitchUseGrainSizeDistribution)
                        {
                            FOR_GRAIN_CLASSES
                            {
                                double blconc2 = _f1D > 0? UF1D_blm_D.Drcd / _f1D->Drc :0.0;
                                double ssconc2 = _f1D > 0? UF1D_ssm_D.Drcd / _f1D->Drc : 0.0;
                                double qbls2 = std::min(UF1D_blm_D.Drcd,blconc2 * volf);
                                double qsss2 = std::min(UF1D_ssm_D.Drcd,ssconc2 * volf);
                                UF1D_blm_D.Drcd -= qbls2;
                                UF2D_blm_D.Drcd += qbls2;
                                UF1D_ssm_D.Drcd -= qsss2;
                                UF2D_ssm_D.Drcd += qsss2;
                            }
                        }
                    }
                }
            }else
            {
                double fractionf = std::max(0.0,std::min(0.95, dt->Drc*fV/std::max(0.01*_dx,0.5*(_dx - _lddw->Drc))));
                double fractions = std::max(0.0,std::min(0.95, dt->Drc*sV/std::max(0.01*_dx,0.5*(_dx - _lddw->Drc))));

                double volf = _f2D->Drc * (fractionf);
                double vols = _s2D->Drc * (fractions);

                double vf = _f1D->Drc + volf;
                double vs = _s1D->Drc + vols;

                //_fu1D->Drc = vf > UF_VERY_SMALL? (_fu1D->Drc *_f1D->Drc + fV * _f2D->Drc * (fractionf))/vf: 0.0;
                //_su1D->Drc = vs > UF_VERY_SMALL? (_su1D->Drc *_s1D->Drc + sV * _s2D->Drc * (fractions))/vs: 0.0;
                //_visc1D->Drc = vf > UF_VERY_SMALL? (_visc1D->Drc *_f1D->Drc + _visc2D->Drc * _f2D->Drc * (fractionf))/vf: 0.0;
                //_d1D->Drc = vs > UF_VERY_SMALL? (_d1D->Drc *_s1D->Drc + _d2D->Drc * vols)/vs: 0.0;
                _ifa1D->Drc = vs > UF_VERY_SMALL? (_ifa1D->Drc *_s1D->Drc + _ifa2D->Drc * vols)/vs: 0.0;
                _rocksize1D->Drc = vs > UF_VERY_SMALL? (_rocksize1D->Drc *_s1D->Drc + _rocksize2D->Drc * vols)/vs: 0.0;

                _f1D->Drc = vf;
                _s1D->Drc = vs;

                _f2D->Drc = _f2D->Drc - volf;
                _s2D->Drc = _s2D->Drc - vols;

                if(SwitchErosion)
                {
                    if(UF1D_blm->Drc < 0 || UF1D_ssm->Drc < 0)
                    {
                        qDebug() << "neg 1" << r << c << UF1D_blm->Drc << UF1D_ssm->Drc << UF2D_blm->Drc << UF2D_ssm->Drc;
                    }
                    double blconc = _f2D > 0? UF2D_blm->Drc/ _f2D->Drc : 0.0;
                    double ssconc = _f2D > 0? UF2D_ssm->Drc/ _f2D->Drc : 0.0;
                    double qbls = std::min(UF2D_blm->Drc,blconc * volf);
                    double qsss = std::min(UF2D_ssm->Drc,ssconc * volf);
                    UF1D_blm->Drc += qbls;
                    UF2D_blm->Drc -= qbls;
                    UF1D_ssm->Drc += qsss;
                    UF2D_ssm->Drc -= qsss;

                    if(UF1D_blm->Drc < 0 || UF1D_ssm->Drc < 0)
                    {
                        qDebug() << "neg 2" << r << c << UF1D_blm->Drc << UF1D_ssm->Drc << UF2D_blm->Drc << UF2D_ssm->Drc;
                    }

                    if(SwitchUseGrainSizeDistribution)
                    {
                        FOR_GRAIN_CLASSES
                        {
                            double blconc2 = _f2D > 0? UF2D_blm_D.Drcd / _f2D->Drc :0.0;
                            double ssconc2 = _f2D > 0? UF2D_ssm_D.Drcd / _f2D->Drc : 0.0;
                            double qbls2 = std::min(UF2D_blm_D.Drcd,blconc2 * volf);
                            double qsss2 = std::min(UF2D_ssm_D.Drcd,ssconc2 * volf);
                            UF1D_blm_D.Drcd += qbls2;
                            UF2D_blm_D.Drcd -= qbls2;
                            UF1D_ssm_D.Drcd += qsss2;
                            UF2D_ssm_D.Drcd -= qsss2;
                        }
                    }
                }


            }

        }
    }}}

    //could be done with manning type velocity based on channel/overland flow height difference.
    //This simulates sub-gridcell flow velocity into the channel
    /*double NN = N->Drc;
    double hf = (_f2D->Drc)/(_dx*_lddw->Drc);
    double hs = (_s2D->Drc)/(_dx*_lddw->Drc);

    double Perim = 2*hf+_dx;
    double R = 0;
    if (Perim > 0)
        R = hf*_dx/Perim;
    else
        R = 0;

    double fVkin = pow(R,3.0/2.0) * sqrt(hf/_dx)/NN;

    Perim = 2*hs+_dx;
    R = 0;
    if (Perim > 0)
        R = hs*_dx/Perim;
    else
        R = 0;
    double sVkin = pow(R,3.0/2.0) * sqrt(hs/_dx)/NN;*/

}

void TWorld::UF2D1D_Infiltration(int thread,cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                                 cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                                 cTMap * _fu1D,cTMap * _s1D,
                                 cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                                 cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                                 cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                                 cTMap * _su2D,cTMap * _sv2D)
{

    FOR_ROW_COL_UF2DMT_DT
    {
        double cdx = DX->Drc;
        double cdy = ChannelAdj->Drc;

        //calculate infiltartion in time step
        double infil = std::min((FSurplus->Drc *SoilWidthDX->Drc*cdx*dt->Drc/_dt + UF2D_Infiltration->Drc),0.0);
        if(_f2D->Drc < fabs(infil))
        {
            infil = -_f2D->Drc;

        }
        //keep track of infiltration
        UF2D_Infiltration->Drc -= (infil);
        _f2D->Drc = std::max(_f2D->Drc + infil,0.0);

    }}}
    FOR_ROW_COL_UF1DMT_DT
    {
        //channel infiltration
        if(SwitchChannelInfil)
        {
            double infilvol = std::min(((ChannelKsat->Drc*_lddw->Drc + 2.0 * _f1D->Drc/(_dx * _lddw->Drc))* _dx * dt->Drc/ _dt + UF1D_Infiltration->Drc), _f1D->Drc);
            UF1D_Infiltration->Drc += infilvol;
            _f1D->Drc -= infilvol;
        }
    }}}


}

void TWorld::UF2D1D_ChannelWater(int thread,cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                                 cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                                 cTMap * _fu1D,cTMap * _s1D,
                                 cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                                 cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                                 cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                                 cTMap * _su2D,cTMap * _sv2D)
{

    FOR_ROW_COL_UF1DMT_DT
    {
        // add rainfall in m3, no interception
        _f1D->Drc += (dt->Drc/_dt) * Rainc->Drc*ChannelWidthUpDX->Drc *ChannelDX->Drc;


        if(SwitchInflow)
        {


            if(InflowID->Drc > 0)
            {

                int col = InflowID->Drc;
                int i = 0;
                bool end = false;
                while(!end)
                {
                    if(i <= IFTime.length())
                    {
                        if(IFTime.at(i) < time/60)
                        {
                            i++;
                        }else
                        {

                            end = true;
                        }
                    }else
                    {
                        end = true;
                    }
                }
                i = i - 1;

                double q = 0;

                if(i < 0)
                {
                    q = 0;
                }else
                {
                    q = IFQ.at(col-1)->at(std::min(IFQ.at(col-1)->length(),i));
                }


                _f1D->Drc += dt->Drc * q;

            }
        }

        //add baseflow
        if(SwitchChannelBaseflow)
        {
            if(!addedbaseflow)
            {
                addedbaseflow = true;
                _f1D->Drc += BaseFlowInitialVolume->Drc;
                BaseFlow += BaseFlowInitialVolume->Drc;

            }
            _f1D->Drc += BaseFlowInflow->Drc * _dt;
            BaseFlow += BaseFlowInflow->Drc * _dt;

        }
    }}}


}


void TWorld::UFDEMLDD_Connection(int thread)
{
    /*cTMap * _dem = UF2D_DEM;
    FOR_ROW_COL_UF2DMT
    {
        UF2D_DEM->Drc += DEMChange->Drc;
        DEMChange->Drc = 0;
    }}}*/

}
