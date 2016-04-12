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

#include "model.h"

void TWorld::DrainSewer(double dt)
{
    if(!SwitchSewer)
    {
        return;
    }

    FOR_ROW_COL_MV
    {
        if(SewerID->Drc > 0)
        {
            double surfrac = SewerArea->Drc * _dx * _dx / DX->Drc * _dx;
            double frac = std::min(1.0,pow(surfrac,dt/SewerDrainageTime->Drc));

            SewerQ->Drc = std::min((ChannelAdj->Drc * DX->Drc * WHrunoff->Drc), ChannelAdj->Drc * DX->Drc * WHrunoff->Drc * frac);
            SewerQt->Drc += SewerQ->Drc;

            if(ChannelAdj->Drc > 0)
            {
                WHrunoff->Drc = WHrunoff->Drc - SewerQ->Drc/(ChannelAdj->Drc * DX->Drc);
            }

            if(SwitchErosion)
            {

                SewerQs->Drc = std::min(Sed->Drc,SewerQ->Drc* Conc->Drc * frac);
                SewerQts->Drc += SewerQs->Drc;
                Sed->Drc -=SewerQs->Drc;

                if(SwitchUseGrainSizeDistribution)
                {
                    Sed->Drc = 0;
                    SewerQs->Drc = 0;

                    FOR_GRAIN_CLASSES
                    {
                        SewerQs_D.Drcd = std::min(Sed_D.Drcd,SewerQ->Drc* Conc_D.Drcd * frac);
                        SewerQts_D.Drcd += SewerQs_D.Drcd;
                        Sed_D.Drcd -= SewerQs_D.Drcd;
                        Sed->Drc += Sed_D.Drcd;
                        SewerQs->Drc += SewerQs_D.Drcd;
                    }

                    SewerQts->Drc += SewerQs->Drc;
                }
            }

            SewerQTotal += SewerQ->Drc;
            SewerQsTotal += SewerQs->Drc;
        }
    }

}
