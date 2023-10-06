
(************************************)
{** SOAP version 2.1               **)
(**                                **)
(** UNIT SOAPSoil                  **)
(**                                **)
(** Author: V.G. Jetten            **)
(**    Dept. of Physical Geography **)
(**    Utrecht University          **)
(** Last revision: 24/8/94         **)
(**                                **)
(************************************)


{$F+}
UNIT SOAPSoil;
{----------------------------------------------------------------}
{ UNIT SoapSoil:                                                 }
{ Calculates one dimensional soil water balance based on         }
{ Richards equation, using a fully implicit, mass lumped,        }
{ Galerkin scheme Finite Element Method.                         }
{ Based on FORTRAN program WORM by Van Genuchten (1987),         }
{ Research report no.121 of USDA Salinity Laboratory and         }
{ modified by J.B.Kool and G.H. de Rooij (1989) LUW.             }
{                                                                }
{ Adapted for rainfall, evapotranspiration and groundwater level }
{ change. Conductivity and retention curves according to         }
{ Van Genuchten functions (1980).                                }
{ Important variables:                                           }
{ NN : number of nodes N                                         }
{ H, Hnew : pressure head begin and end of timestep dt           }
{ W, Wnew : water content begin and end of timestep dt           }
{ K1, K2 : hydraulic conductivity                                }
{ C1, C1 & C2 : differential water capacity, 2 dif. methods  }
{ F, A, D : Galerkin scheme tri-diagonal matrix coefficients     }
{ S : sink term                                                  }
{                                                                }
{ NOTE: all fluxes are in CM/MIN and heads in CM !!!             }
{                                                                }
{----------------------------------------------------------------}

INTERFACE

USES Dos, crt, graph, SOAPVar;

PROCEDURE SoilMoisture;


IMPLEMENTATION

PROCEDURE CalcSinkTerm;
VAR i : integer;

   PROCEDURE Transpiration;
   {---------------------------------------------------------------}
   { calculate sink term S with transpiration reduction term       }
   { and canopy wetness fraction and time of day, inpu Evapday []  }
   {---------------------------------------------------------------}
   VAR redterm : real;
       i : integer;
   BEGIN
     with profile^ do
     for i := 1 to NN do
     Begin
      { calculate reduction term to multiply Tp with }
      { PET reduction term used by Koorevaar }

    {  redterm := 1/(1+power((-200+Hx^[i])/-1200,1.6));  }

      { original van Genuchten, set 50%transp on pF 3.5 }
      with veg^ do
      redterm :=  1/(1+power(H^[i]/H50,Hbeta));
      { linear reduction when very wet }
      if H^[i] > -10 then redterm := H^[i]/-10;
      { no uptake when node is saturated }
      if H^[i] >= 0 then redterm := 0;
      { no uptake when dryer then wilting point }
      if H^[i] < Hwilting then redterm := 0;
      { calc sink term in cm/min with relative root distr in layer }
      S^[i] := redterm * VEG^.RD^[i]/roottot * ATM^.Trans
     End;
   END;


   PROCEDURE Evaporation;
   {---------------------------------------------------------------}
   { note that flux is evaluated with average node 1and 2 and 0.5  }
   { of thickness of first layer dz, potential of air above soil   }
   { with ln(RH)*RT/Mg                                             }
   {---------------------------------------------------------------}
   CONST R = 8.3143; { universal Gas constant, J/mol/kg }
         M = 0.018;  { H2O kg/m3 }
         g = 9.87;   { gravity }
   VAR HairAbove, Hsoil : real;
   BEGIN
    with profile^, meteoact[VEG^.nol+1], VEG^ do
    if not ponded then
    Begin
     { minimum H allowed evaluated for 0.5 first layer }
     Hsoil := (H^[0]+H^[1])/2;
     HairAbove := ln(RH/100) * (R*(Temp+273))/(M*g);
     Sevap := -sqrt(fKWH(mat[0], 2, HairAbove, true) *
                    fKWH(mat[0], 2, Hsoil, true))
               * ((HairAbove-Hsoil)/(0.5*dz^[0])+1);

     { soil evaopration must be smaller then
       demand atmosphere at lowest canopy level }
     If Sevap > ATM^.SoilEvapMax then Sevap := ATM^.SoilEvapMax;
     if Sevap < 0 then Sevap := 0;
    End
    else
     Sevap := ATM^.SoilEvapMax;
   END;


 { calc sinkterm }
BEGIN
 with Profile^ do
 Begin
  { transpiration }
  fillchar(S^,sizeof(S^),0);
  if useuptake and (ATM^.trans > 0) then Transpiration;

  { soil evaporation }
  Sevap := 0;
  If usesoilevap and not ponded
   and (ATM^.SoilEvapMax > 0) then Evaporation;

  { calc total uptake with trapezium rule }
  if useuptake or usesoilevap then
  Begin
   { let op niet dz als dx maar 1 omdat S een flux is dus S/dz }
   Stot := S^[0]+S^[NN];
   for i := 1 to NN-1 do Stot := Stot + S^[i]*2;
   Stot := Stot/2;

   { include soil evap in top nodes }
   S^[0] := S^[0] + Sevap;
   for i := 0 to NN do S^[i] := S^[i]/dz^[i];
  End;

 End;
END;


PROCEDURE CalcNewNodalValues;
{---------------------------------------------------------------}
{ calculate nodal values of conductivity and diff. moist. cap.  }
{ nodes are evaluated at 0.5 new and old H to ensure stability  }
{ NOTE: dz/dt included here in capacity calculations            }
{---------------------------------------------------------------}
VAR i : integer;
    dw, dw1, dh, dh1 : double;
    top : boolean;
BEGIN
 with profile^ do
 Begin
  for i := 0 to NN do
  with soil^[mat[i]] do
  Begin
   top := i = 0;
   { calc average H }
   Hx^[i] := 0.5 * (Hnew^[i] + H^[i]);
   { calc average K }
   K^[i] := fKWH(mat[i], 2, Hx^[i], top);
   { calc wnew for hnew }
   Wnew^[i] := fKWH(mat[i], 1, Hnew^[i], top);

   { calc capacity C for Hnew }
   { in original WORM always analytical:}
   { if ILCP = 0 then C1^[i] := fKWH(mat[i], 3, Hx^[i], top) * dz/dt; }
   { Fastest way of calculating capacity, Karvonen :}

   if abs(Hnew^[i]-H^[i]) <= 3*tol2 then
   Begin
    { if dh is too small then dé/dh impossible and analytical solution }
    C1^[i] := fKWH(mat[i], 3, Hx^[i], top) * dz^[i]/dt
   End
   else
    { else C literally from dé/dh }
    C1^[i] := (Wnew^[i]-W^[i]) / (Hnew^[i]-H^[i]) * dz^[i]/dt;

  End;
  { mass conservative C solution }
  if ILCP = 2 then
  Begin
   for i := 0 to NN-1 do
   Begin
    dh := Hnew^[i]-H^[i];
    dh1 := Hnew^[i+1]-H^[i+1];
    dw := Wnew^[i]-W^[i];
    dw1 := Wnew^[i+1]-W^[i+1];
    if abs(2*dh+dh1) <= 3*tol2 then
     C1^[i] := fKWH(mat[i], 3, Hx^[i], top) * dz^[i]/dt
    else
     C1^[i] := 3*dw/(2*dh+dh1) * dz^[i]/dt;

    if abs(dh+2*dh1) <= 3*tol2 then
     C2^[i] := fKWH(mat[i], 3, Hx^[i+1], top) * dz^[i+1]/dt
    else
     C2^[i] := 3*dw1/(dh+2*dh1) * dz^[i+1]/dt;
   End;
  End;

 End;
END;


PROCEDURE SolveFiniteElement;
{---------------------------------------------------------------------}
{ Galerkin scheme finite element solution, matrix dioganals F,A and D }
{ Steps: 1-conductivity calc; 2 set up matrix; 3 include K & sink term}
{ 4-include boundary conditions in matrix ; 5-Gaussian elimination &  }
{ and back substitution; 6-calc boundary fluxes                       }
{---------------------------------------------------------------------}
VAR i, U1, U2 : integer;
    ANE, FNN1, FNN2, A1, A2 : double;
BEGIN
 with profile^ do
 Begin

  { CONDUCTIVITIES }

  for i := 0 to NN do K^[i] := K^[i]/dz^[i];
  Case UPW of
  0 : Begin
       { arithmetic weighing }

       for i := 0 to NN-1 do K2^[i] := 0.5*(K^[i]+K^[i+1]);
       for i := 1 to NN   do K1^[i] := 0.5*(K^[i-1]+K^[i]);
       K1^[0] := K2^[0];
       K2^[NN] := K1^[NN];
      End;
  1 : Begin
       { no weighing: K = K of wettest node if dH <= dz or
                      K = K of dryest node if dH > dz or   }
       for i := 1 to NN do
        if Hx^[i] > Hx^[i-1]+dz^[i] then K1^[i] := K^[i]
                                else K1^[i] := K^[i-1];
       for i := 0 to NN-1 do
        if Hx^[i+1]-dz^[i] > Hx^[i] then K2^[i] := K^[i+1]
                                else K2^[i] := K^[i];
       K1^[0] := K2^[0];
       K2^[NN] := K1^[NN];
      End;
  2 : Begin
       { weighting with K = geometric mean }
       for i := 0 to NN-1 do K2^[i] := sqrt(K^[i]*K^[i+1]);
       for i := 1 to NN   do K1^[i] := sqrt(K^[i-1]*K^[i]);
       K1^[0] := K2^[0];
       K2^[NN] := K1^[NN];
      End;
  End;
  for i := 0 to NN do K^[i] := K^[i]*dz^[i];


  { TRIDIAGONAL MATRIX for Galerkin scheme }

  { top node }
  case ILCP of
   0 : Begin { original setup, looks like Simpson num. integral }
        F^[0] := (2*C1^[0] + C1^[1])/6;
        A^[0] := -K1^[0];
        D^[0] := -A^[0] + F^[0];
        for i := 1 to NN-1 do
        Begin
         F^[i] := (C1^[i-1] + 4*C1^[i] + C1^[i+1])/6;
         A^[i] := -K2^[i];
         D^[i] := K1^[i] + K2^[i] + F^[i];
        End;
        ANE := A^[NN-1];
        F^[NN] := (C1^[NN-1] + 2*C1^[NN])/6;
        D^[NN] := -A^[NN-1] + F^[NN];
       End;

   1 : Begin { alternative, simple trapezium solution }
        F^[0] := 0.5*C1^[0];
        A^[0] := -K1^[0];
        D^[0] := -A^[0] + F^[0];
        for i := 1 to NN-1 do
        Begin
         F^[i] := C1^[i];
         A^[i] := -K2^[i];
         D^[i] := K1^[i] + K2^[i] + F^[i];
        End;
        ANE := A^[NN-1];
        F^[NN] := 0.5*C1^[NN];
        D^[NN] := -A^[NN-1] + F^[NN];
       End;

   2 : Begin { alternative, mass conservative solution }
        F^[0] := (2*C1^[0] + C2^[0])/6;
        A^[0] := -K1^[0];
        D^[0] := -A^[0] + F^[0];
        for i := 1 to NN-1 do
        Begin
         F^[i] := (C1^[i-1] + 2*C2^[i-1] + 2*C1^[i] + C2^[i])/6;
         A^[i] := -K2^[i];
         D^[i] := K1^[i] + K2^[i] + F^[i];
        End;
        ANE := A^[NN-1];
        F^[NN] := (C1^[NN-1] + 2*C2^[NN-1])/6;
        D^[NN] := -A^[NN-1] + F^[NN];
       End;

  End;
  { save F[nn] for lower boundary }
  FNN1 := F^[NN];

  { include gravity term (i.e. k only) and sink term }
  { NOTE in WORM for first term originally -K1 + S1 but in paper
    by Van Genuchten -K1 -S1 !!! }

  F^[0] := F^[0]*H^[0] - dz^[0] * (K1^[0] - (2*S^[0]+S^[1])/6);
  for i := 1 to NN-1 do
    F^[i] := F^[i]*H^[i] + 0.5*(dz^[i]+dz^[i-1]) * K1^[i]
                         - 0.5*(dz^[i]+dz^[i+1]) * K2^[i]
                         - dz^[i] * (S^[i-1]+4*S^[i]+S^[i+1])/6;
  F^[NN] := F^[NN]*H^[NN] + dz^[NN] * (K2^[NN] - (S^[NN-1]+2*S^[NN])/6);


  { BOUNDARY CONDITIONS }

  { Upper boundary condition }
  if ponded then
  Begin
   Hnew^[0] := H^[0] + dt * (throughfall-influx);

   influx := D^[0]*Hnew^[0] - F^[0];

   F^[0] := Hnew^[0];
   A1 := A^[0];
   A^[0] := 0;
   D^[0] := 1;
   F^[1] := F^[1] - A1 * Hnew^[0];
  End
  else
  Begin
   { no ponding so influx equals flux from above }
   influx := throughfall;
   F^[0] := F^[0] + influx;
  End;

  { lower boundary condition }
  Case lowerBC of
   0 : Begin { free drainage }
        D^[NN] := -A^[NN-1];
        FNN2 := F^[NN];
        F^[NN] := 0;
       End;
   1 : Begin { fixed potential in lower node }
        Drain := F^[NN] - D^[NN]*Hnew^[NN];
        F^[NN] := Hnew^[NN];
        D^[NN] := 1.0;
        F^[NN-1] := F^[NN-1] - ANE*Hnew^[NN];
        A^[NN-1] := 0;
       End;
   2 : Begin { fixed flux: drain spec by user, e.g. when groundwater }
        Drain := GroundwaterFlux/1440;
        F^[NN] := F^[NN] - Drain;
       End;
  End;

  { calc Hnew with Gaussian elimination and back substitution }
  for i := 1 to NN do
  Begin
   A2 := A^[i-1]/D^[i-1];
   D^[i] := D^[i] - A2 * A^[i-1];
   F^[i] := F^[i] - A2 * F^[i-1];
  End;

  Hnew^[NN] := F^[NN]/D^[NN];
  for i := NN-1 downto 0 do
   Hnew^[i] := (F^[i] - A^[i]*Hnew^[i+1])/D^[i];

  for i := 1 to NN do
  Begin
   if (lowerBC = 2) and (z^[i] >= GWlevel) then break;
   if (Hnew^[i] >= Hair) then Hnew^[i] := Hair;
  End;

  { calc boundary fluxes }

  if ponded and (Hnew^[1] < 0) then influx := influx + A1 * Hnew^[1];
    { de influx is nu gelijk aan:
      K1*((hnew1-hnew2)/dz+1) + 0.5*dé/dh * (hnew1-h1) }

  Case LowerBC of
   { free drainage }
   0 : Begin
        Drain := 0.5*(fKWH(mat[NN], 2, H^[NN], false) +
                      fKWH(mat[NN], 2, Hnew^[NN], false));
        if Hnew^[NN] < 0 then
         Drain := 0.5*(drain + FNN2 - Hnew^[NN]*FNN1)
        else
         Drain := 0.5*(drain + FNN2);
       End;
   { fixed potential }
   1 : Drain := Drain - ANE*Hnew^[NN-1];
  End;
 End;
END;



PROCEDURE Corrections;
{-------------------------------------------------------------------}
{ Corrections needed if situation changes from ponded to not ponded }
{ or vice versa within iteration: surface storage, first node and   }
{ influx have to be adjusted                                        }
{ Calculation of Wnew and K for output screen                       }
{-------------------------------------------------------------------}
VAR i, j : integer;
    influxmax, balerr, oldsurfstor : double;
BEGIN
 with profile^ do
 Begin

  if lowerBC = 2 then
  Begin
   for i := NN downto 1 do
    if (Hnew^[i] < Hair) or (Hnew^[i] < 0) then break;
   GWlevel := z^[i+1] - Hnew^[i+1]/(Hnew^[i+1]-Hnew^[i])*dz^[i+1];
  End;


  oldsurfstor := Surfstor;

  Surfstor := 0;

  if ponded then
  Begin
   { correct influx if change from non-ponded to ponded at
     last iteration step }
   if (influx = throughfall) then
   Begin
    influx := (soil^[mat[0]].Wsat - Wnew^[0])*dz^[0]/dt;
    Hnew^[0] := dt*(throughfall-influx);
   End;

   { calc overland flow }
   OLflow := (Hnew^[0] - SstorageMax);
   if OLflow < 0 then OLflow := 0;
   if Hnew^[0] > Sstoragemax then Hnew^[0] := Sstoragemax;

   { set storage to H first node }
   Surfstor := Hnew^[0];
   if Surfstor < 0 then Surfstor := 0;

   balerr := (throughfall-influx)*dt - (Surfstor-oldsurfstor+OLflow);
   if abs(balerr) > 1e-6 then
   Begin
    if surfstor < Sstoragemax then
    Begin
     Surfstor :=  Surfstor + balerr;
     Hnew^[0] := Surfstor;
    End
     else OLflow := OLflow + balerr;
   End;
  End;

 End;
END;


PROCEDURE SoilMoisture;
{---------------------------------------------------------------------}
{ iteration loop until dt too small or error in H less than tolerance }
{---------------------------------------------------------------------}
VAR tol : real;
    stopit : boolean;
    i, j : integer;
BEGIN
 with profile^ do
 Begin
  { Calc sink term from transpiration and evaporation }
  CalcSinkTerm;

  NIT := 0;
  Repeat
   { check ponded }
   ponded := Hnew^[0] > 0;
   { increase iteration step counter }
   NIT := NIT + 1;
   { copy previous iter step potentials }
   move(Hnew^, Hold^, (NN+1) * sizeof(double));
   { Calc new K, C, Wnew from Hnew }
   CalcNewNodalValues;
   { solve matrix with FEM }
   SolveFiniteElement;
   { check ponded }
   ponded := Hnew^[0] > 0;
   { check for convergence by comparing new with old potentials }
   stopit := true;
   for i := 0 to NN do
   Begin
    { tol = 0.2 + 1% of potential }
    tol := tol1*abs(Hold^[i]) + tol2;
    if ponded then tol := tol/2;
    if abs(Hnew^[i] - Hold^[i]) > tol then stopit := false;
    if not stopit then break;
   End;

   { if too many iter then 0.5 dt and redo }
   If not stopit and (NIT > nitmax) then
   Begin
    dt := dt/2;
    time := time - dt;
    for i := 0 to NN do Hnew^[i] := (H^[i]+Hnew^[i])/2;
   End;

  Until stopit or (dt < dtmin);

  if dt < dtmin then dt := dtorig;
  Corrections;

 End;
END;


END.

