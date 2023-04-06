/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
 \file version.h
 \brief version information, add your name here if you contribute to the development

  Auhtors maintaining code: \n
 - VJ = Victor Jetten, Dept of Earth Systems Analysis, ITC, Twente University, v.g.jetten@utwente.nl\n
 - BB = Bastian van den Bout, Dept of Earth Systems Analysis, ITC, Twente University, b.vandenbout@utwente.nl\n
 - XX = <name, affiliation, email>\n
 */

#ifndef VERSION_H_
#define VERSION_H_

#define VERSIONNR "6.885"
#define DATE "2023/04/05"
#define VERSIONNR_pest "0.62"


#define VERSION QString("openLISEM version %1 - %2 pesticides %3").arg(VERSIONNR).arg(DATE).arg(VERSIONNR_pest)


#endif /* VERSION_H_ */
