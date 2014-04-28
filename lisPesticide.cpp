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

// check if cell From flows to To
#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )



/*!
  \file lisPesticide.cpp
  \brief Transport and partitionning of pesticides

functions: \n
- void TWorld::TransPesticide() \n
- double TWorld:: cmx=(double t,double k1,double k2,double A, double B,double C, double D,double Cm0,double Cs0)
'diff(Cm(t),t,1)+k1*Cm(t)-A*C-B*Cs(t)=0;
Solution from Maxima
- double TWorld:: csx=(double t,double k1,double k2,double A,double B,double C,double D,double Cm0,double Cs0)
'diff(Cs(t),t,1)+k2*Cs(t)-D*Cm(t)=0; Solution from Maxima
- double TWorld:: factorize(solMax,A,info)

*/


//---------------------------------------------------------------------------
//Mixing layer and infiltration
void TWorld :: Pestmobilisation(void)
{

    if (!SwitchPesticide)
        return;

    FOR_ROW_COL_MV
    {
        pestiinf->Drc = (InfilVol->Drc)/(_dt*_dx*_dx); //DX->Drc);
        //if (WaterVolall->Drc < dx*DX->Drc*1e-6)
        if (Q->Drc < 1e-6)//< 1e-6) //reflexion a avoir sur la def de ruissellement
        {
            C->Drc=0.0;
            C_N->Drc=0.0;
            Pest->Drc = 0.0;
            // gaat niet goed, pest moet toch ergens blijven, kan niet zomaar verdwijnen?
            PCinfilt->Drc = 0.0;
            PCfilmexit->Drc = 0.0;
            Pdetach->Drc = 0.0;

            if (pestiinf->Drc <= 1e-6)
            {
      //          qDebug() << "no inf";
                // mise à l'équilibre: pas de ruissellement ni d'infiltration
                CM->Drc=(poro->Drc*CM_N->Drc+rhob->Drc*CS_N->Drc)/(rhob->Drc*KD->Drc+poro->Drc);
                CS->Drc=KD->Drc*CM_N->Drc;
            }
            else
            {
        //        qDebug() << "no runoff but inf";
                // calcul de CM et CS avec Kfilm=0
                CM->Drc=cmx_analytique(_dt,0.0,pestiinf->Drc,epsil->Drc,rhob->Drc,kr->Drc,KD->Drc,poro->Drc,
                                       CM_N->Drc,CS_N->Drc,C_N->Drc);
                CS->Drc=csx_analytique(_dt,0.0,pestiinf->Drc,epsil->Drc,rhob->Drc,kr->Drc,KD->Drc,poro->Drc,
                                       CM_N->Drc,CS_N->Drc,C_N->Drc);
            }
        }
        else // ruissellement
        {
            //qDebug() << "yes runoff";

            CM->Drc=cmx_analytique(_dt,Kfilm->Drc,pestiinf->Drc,epsil->Drc,rhob->Drc,kr->Drc,KD->Drc,poro->Drc,
                                   CM_N->Drc,CS_N->Drc,C_N->Drc);
            CS->Drc=csx_analytique(_dt,Kfilm->Drc,pestiinf->Drc,epsil->Drc,rhob->Drc,kr->Drc,KD->Drc,poro->Drc,
                                   CM_N->Drc,CS_N->Drc,C_N->Drc);

            PCinfilt->Drc = pestiinf->Drc*C_N->Drc*_dx*_dx*_dt*1000*1000*1000;
            PCfilmexit->Drc = Kfilm->Drc*C_N->Drc*_dx*_dx*_dt*1000*1000*1000;
            Pdetach->Drc = Kfilm->Drc*CM->Drc*_dx*_dx*_dt*1000*1000*1000;

        }

        Pest->Drc = Pest->Drc+ (Kfilm->Drc*(CM->Drc-C_N->Drc)-pestiinf->Drc*C_N->Drc)*_dx*_dx*_dt;
        Pest->Drc = max(0, Pest->Drc);

        C->Drc =  ConcentrationP(WaterVolall->Drc, Pest->Drc);

        //conditions initiales du pas de temps suivant
        CM_N->Drc=CM->Drc;
        CS_N->Drc=CS->Drc;
        C_N->Drc=C->Drc;

    }
    CM_N->report("CM_N");
    CS_N->report("CS_N");
    C_N->report("C_N");
    Pest->report("pest");
//    WaterVolall->report("wva");
   // double hoi1 = Pest->mapTotal();
    //qDebug() << hoi << hoi1;
    //      qDebug()<< "Q<:C"<< C->DrcOutlet <<"cm"<< CM->DrcOutlet << "CS"<< CS->DrcOutlet;
}

//---------------------------------------------------------------------------
double TWorld::Implicitscheme(double Qj1i1, double Qj1i, double Qji1,double Pj1i, double Pji1, double alpha, double dt,double dx, double Kfilm, double CMi1j1)
{
    double Pj1i1, Cavg, Qavg, aQb, abQb_1, A, B, C;
    const double beta = 0.6;
    double s=Kfilm*CMi1j1;
    if (Qj1i1 < MIN_FLUX)
        return (0);

    Qavg = 0.5*(Qji1+Qj1i);
    if (Qavg <= MIN_FLUX)
        return (0);

    Cavg = (Pj1i+Pji1)/(Qj1i+Qji1);
    aQb = alpha*pow(Qavg,beta);
    abQb_1 = alpha*beta*pow(Qavg,beta-1);

    A = dt*Pj1i;
    B = -dx*Cavg*abQb_1*(Qj1i1-Qji1);
    C = (Qji1 <= MIN_FLUX ? 0 : dx*aQb*Pji1/Qji1);
    if (Qj1i1 > MIN_FLUX)
        Pj1i1 = (dx*dt*s+A+C+B)/(dt+dx*aQb/Qj1i1);
    else
        Pj1i1 = 0;

    return max(0,Pj1i1);
}
//---------------------------------------------------------------------------
double TWorld::ConcentrationP(double watvol, double pest)
{
    double conc = Maxsolubility;

    if (watvol > 0)
    {
        conc = pest/watvol;
        // 1e-6 is 1 ml/m2
        //    qDebug()<<"concentration ok ";
    }
    else
    {
        conc = 0;
        //qDebug()<< "forcer concentration nulle" << pest;
        //     qDebug()<< "watvol reajust oui";
    }

    //   if (conc > Maxsolubility)
    //      {
    //       conc = Maxsolubility;
    //       qDebug()<< "max xolubility";
    //      }

    //   pest = conc * watvol;

    return conc;
}
//---------------------------------------------------------------------------
//n*eps*diff(CM(t),t)-K*(C-CM(t))-f*(C-CM(t))+eps*rho*kr*(kd*CM(t)-CS(t));
//Solution analytique from Maxima : Philippe Ackerer

double TWorld::cmx_analytique(double t,double dKfi,double dpestiinf,double depsil,double drhob,double dkr,double dKD,double dn, double CM0,double CS0, double Cr)

{
    double a = -(dKfi+dpestiinf+depsil*drhob*dkr*dKD)/(dn*depsil);
    double b = drhob*dkr/dn;
    double c = dkr*dKD;
    double d = -dkr;
    double e = (dKfi+dpestiinf)*Cr/(dn*depsil);

    double CmX=0.0;

    //    qDebug()<< "time"<< time/60;
    //    qDebug()<< "a" << a;
    //    qDebug()<< "b" << b;
    //    qDebug()<< "c" << c;
    //    qDebug()<< "d" << d;
    //    qDebug()<< "e" << e;
    //    qDebug()<< "Cr" << Cr;
    //    qDebug()<< "dkfi+dpestiinf" << dKfi+dpestiinf;
    //    qDebug()<< "dn" << dn;
    //    qDebug()<< "depsil" << depsil;


    double c1 = (-2*a*b*d*CS0+2*pow(b,2)*c*CS0+(d*(a*sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))-b*c-pow(a,2))
                                                +b*c*(a-sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2)))+a*pow(d,2))*CM0+(d*(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))-a)+pow(d,2)+2*b*c)*e)
            /(2*a*d*sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))-2*b*c*sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2)));

    //  qDebug()<< "C1" << c1;

    double c2 = -(-2*a*b*d*CS0+2*pow(b,2)*c*CS0+(d*(-a*sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))-b*c-pow(a,2))
                                                 +b*c*(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))+a)+a*pow(d,2))*CM0+(d*(-sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))-a)+pow(d,2)+2*b*c)*e)
            /(2*a*d*sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))-2*b*c*sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2)));

    //qDebug()<< "C2" << c2;

    // solutions particulières constante
    CmX = exp(-(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t/2))*((a*c2*d-b*c*c2)*exp(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t+d*t/2+a*t/2)+(a*c1*d-b*c*c1)*exp(d*t/2+a*t/2)
                                                           -d*e*exp(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t/2))/(a*d-b*c);

    // qDebug()<< "ad-bc" << a*d-b*c;

    // solutions particulieres phi.int(phi-1.b)
    //  CmX = -exp(-(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t/2))*(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*((pow(d,2)-a*d+2*b*c)*e*exp(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))
    //      *t+d*t/2+a*t/2)+(-pow(d,2)+a*d-2*b*c)*e*exp(d*t/2+a*t/2))+((-pow(d,3)+2*a*pow(d,2)+(-4*b*c-pow(a,2))*d)*e-2*a*c2*pow(d,3)+(2*b*c+4*pow(a,2))*c2*pow(d,2)
    //      +(-12*a*b*c-2*pow(a,3))*c2*d+(8*pow(b,2)*pow(c,2)+2*pow(a,2)*b*c)*c2)*exp(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t+d*t/2+a*t/2)
    //     +((-pow(d,3)+2*a*pow(d,2)+(-4*b*c-pow(a,2))*d)*e-2*a*c1*pow(d,3)+(2*b*c+4*pow(a,2))*c1*pow(d,2)+(-12*a*b*c-2*pow(a,3))*c1*d
    //     +(8*pow(b,2)*pow(c,2)+2*pow(a,2)*b*c)*c1)*exp(d*t/2+a*t/2)+(2*pow(d,3)-4*a*pow(d,2)+(8*b*c+2*pow(a,2))*d)*e*exp(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t/2))
    //    /(2*a*pow(d,3)+(-2*b*c-4*pow(a,2))*pow(d,2)+(12*a*b*c+2*pow(a,3))*d-8*pow(b,2)*pow(c,2)-2*pow(a,2)*b*c);

    return CmX;

}

//--------------------------------------------------------------------------------------------------------------------------------------
//diff(CS(t),t)-kr*(kd*CM(t)-CS(t));
//Solution analytique from Maxima : Philippe Ackerer

double TWorld::csx_analytique(double t, double dKfi,double dpestiinf,double depsil,double drhob,double dkr,double dKD,double dn, double CM0,double CS0,double Cr)
{
    double a = -(dKfi+dpestiinf+depsil*drhob*dkr*dKD)/(dn*depsil);
    double b = drhob*dkr/dn;
    double c = dkr*dKD;
    double d = -dkr;
    double e =(dKfi+dpestiinf)*Cr/(dn*depsil);

    double CsX=0.0;

    double c1 = (-2*a*b*d*CS0+2*pow(b,2)*c*CS0+(d*(a*sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))-b*c-pow(a,2))
                                                +b*c*(a-sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2)))+a*pow(d,2))*CM0+(d*(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))-a)+pow(d,2)+2*b*c)*e)
            /(2*a*d*sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))-2*b*c*sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2)));

    double c2 = -(-2*a*b*d*CS0+2*pow(b,2)*c*CS0+(d*(-a*sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))-b*c-pow(a,2))
                                                 +b*c*(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))+a)+a*pow(d,2))*CM0+(d*(-sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))-a)+pow(d,2)+2*b*c)*e)
            /(2*a*d*sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))-2*b*c*sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2)));


    // solutions particulières constante
    CsX =exp(-(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t/2))*(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*((a*c2*d-b*c*c2)*exp(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t+d*t/2+a*t/2)+(b*c*c1-a*c1*d)
                                                                                               *exp(d*t/2+a*t/2))+(a*c2*pow(d,2)+(-b*c-pow(a,2))*c2*d+a*b*c*c2)*exp(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t+d*t/2+a*t/2)+(a*c1*pow(d,2)+(-b*c-pow(a,2))*c1*d+a*b*c*c1)
                                                          *exp(d*t/2+a*t/2)+2*b*c*e*exp(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t/2))/(2*a*b*d-2*pow(b,2)*c);

    // solutions particulieres phi.int(phi-1.b)
    //     CsX =exp(-(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t/2))*(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*(((b*c*d+a*b*c)*e+a*c2*pow(d,3)+(-b*c-2*pow(a,2))*c2*pow(d,2)
    //           +(6*a*b*c+pow(a,3))*c2*d+(-4*pow(b,2)*pow(c,2)-pow(a,2)*b*c)*c2)*exp(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t+d*t/2+a*t/2)+((-b*c*d-a*b*c)*e
    //           -a*c1*pow(d,3)+(b*c+2*pow(a,2))*c1*pow(d,2)+(-6*a*b*c-pow(a,3))*c1*d+(4*pow(b,2)*pow(c,2)+pow(a,2)*b*c)*c1)*exp(d*t/2+a*t/2))+((-b*c*pow(d,2)+2*a*b*c*d-4*pow(b,2)*pow(c,2)
    //           -pow(a,2)*b*c)*e+a*c2*pow(d,4)+(-b*c-3*pow(a,2))*c2*pow(d,3)+(7*a*b*c+3*pow(a,3))*c2*pow(d,2)+(-4*pow(b,2)*pow(c,2)-7*pow(a,2)*b*c-pow(a,4))*c2*d+(4*a*pow(b,2)*pow(c,2)+pow(a,3)*b*c)*c2)
    //           *exp(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t+d*t/2+a*t/2)+((-b*c*pow(d,2)+2*a*b*c*d-4*pow(b,2)*pow(c,2)-pow(a,2)*b*c)*e+a*c1*pow(d,4)+(-b*c-3*pow(a,2))*c1*pow(d,3)
    //           +(7*a*b*c+3*pow(a,3))*c1*pow(d,2)+(-4*pow(b,2)*pow(c,2)-7*pow(a,2)*b*c-pow(a,4))*c1*d+(4*a*pow(b,2)*pow(c,2)+pow(a,3)*b*c)*c1)*exp(d*t/2+a*t/2)+(2*b*c*pow(d,2)-4*a*b*c*d+8*pow(b,2)*pow(c,2)
    //           +2*pow(a,2)*b*c)*e*exp(sqrt(pow(d,2)-2*a*d+4*b*c+pow(a,2))*t/2))/(2*a*b*pow(d,3)+(-2*pow(b,2)*c-4*pow(a,2)*b)*pow(d,2)+(12*a*pow(b,2)*c+2*pow(a,3)*b)*d-8*pow(b,3)*pow(c,2)
    //           -2*pow(a,2)*pow(b,2)*c);

    return CsX;
}


