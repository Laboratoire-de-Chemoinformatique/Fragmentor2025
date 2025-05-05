{ Fragmentor of the ISIDA Project

  Copyright (C) 2022 Laboratoire de Chemoinformatique, UMR 7140 CNRS (http://complex-matter.unistra.fr/equipes-de-recherche/laboratoire-de-chemoinformatique/home/)
  contact: Gilles Marcou g.marcou@unistra.fr

  This library is free software; you can redistribute it and/or modify it
  under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 3 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
  for more details.

  You should have received a copy of the GNU Library General Public License
  along with this library; if not, write to the Free Software Foundation,
  Inc., 51 Franklin Street - Fifth Floor, Boston, MA 02110-1335, USA.
}
unit UnitShortestPath;
//Unit to fragment based on shortest path sequences
{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,math,unitfragmentbase,unitsequences, UnitAtmPrpWeight;

type

    TShrtPthException = class(Exception)
    end;


    { TShortestPath }

    TShortestPath = class(TSequences)
    private
           PrpWA: TAtmPropWeight;

    protected
             procedure WayStep(waS, waF, waP: Integer; UseMarkAtom: boolean); override;
             function DoAtomS(UseMarkAtom: boolean): integer;
    public
          constructor Create;
          destructor Destroy; override;
          procedure SDFToFrgLst(MolSDF: TStringList); override;
          procedure NBSCR; virtual; abstract; // Sequences
    end;

implementation

{ TShortestPath }

constructor TShortestPath.Create;
begin
     inherited Create;
     PrpWA:=TAtmPropWeight.Create;
     PrpWA.init_atmPrTable;
end;

destructor TShortestPath.Destroy;
begin
     inherited Destroy;
     PrpWA.Free;
end;

procedure TShortestPath.SDFToFrgLst(MolSDF: TStringList);
begin
     TrailSDFInput(MolSDF);
     NBSCR;
end;

procedure TShortestPath.WayStep(waS, waF, waP: Integer; UseMarkAtom: boolean);
var
   iwm, iwc: integer;
   I,J,Idx, delta, c1, c2: integer;
   TEMPd: double;
   debug: boolean;
begin
     if (waS = waF) then // points s and f are equal!
     begin
          //debug
          debug:=false;
          writeln;
          if (waP-1 = 4) then
          begin
               debug:=true;
               write('('+IntToStr(Road[1])+','+IntToStr(Road[waP-1])+')');
               for iwm := 1 to waP - 1 do
                   write(IntToStr(G[Road[iwm],Road[iwm]])+' ');//+':'+FloatToStr(PrpWA.GetAtmPrp(G[Road[iwm],Road[iwm]]))+' ');
               writeln;
          end;
          //Estimate the weight of the current Road.
          TEMPd:=0;
          c1:=0;
          c2:=waP;
          for iwm := 1 to waP - 1 do
          begin
               inc(c1);
               dec(c2);
               delta:=min(c1,c2);
               TEMPd:=TEMPd+delta*PrpWA.GetAtmPrp(G[Road[iwm],Road[iwm]]);
               write(IntToStr(delta)+':');
          end;
          writeln;
          //If it is a better Road than an existing one
          if ((waP - 2) < W[Start,Finish]) or
             (((waP - 2) = W[Start,Finish]) and (TEMPd > Wd[Start,Finish])) then
          begin
               //Calculate the index of the Road or identify it as a new Road
               if (W[Start,Finish]=NOTCONNECTED) then
               begin
                    inc(NAS);
                    Z[Start,Finish]:=NAS;
                    Z[Finish,Start]:=NAS;
               end;
               Idx:=Z[Start,Finish];
               //Update length and weight of the Road
               Wd[Start,Finish]:=TEMPd;
               Wd[Finish,Start]:=TEMPd;
               W[Start,Finish]:=waP-2;
               W[Finish,Start]:=waP-2;
               //Store the Road
               //writeln(IntToStr(Idx)+'**'+IntToStr(waP-1)+'**'+IntToStr(Start)+'**'+IntToStr(Finish));
               for iwm := 1 to waP - 1 do
               //begin
                   AtomS[Idx,iwm]:=Road[iwm];
               //    write(IntToStr(Road[iwm])+'//');
               //end;
               //writeln;
               //debug
               if (waP-1 = 4) then
               begin
                    write('*('+IntToStr(Road[1])+','+IntToStr(Road[waP-1])+')');
                    for iwm := 1 to waP - 1 do
                        write(IntToStr(G[Road[iwm],Road[iwm]])+' ');
                    writeln;
               end;
          end;
          //debug
          writeln('W:');
          for I:=1 to N do
          begin
              for J:=1 to N do
                  write(' '+IntToStr(W[I,J]));
              writeln;
          end;
          writeln;
          writeln('Wd:');
          for I:=1 to N do
          begin
              for J:=1 to N do
                  write(' '+FloatToStr(Wd[I,J]));
              writeln;
          end;
          writeln;
          writeln('Z:');
          for I:=1 to N do
          begin
              for J:=1 to N do
                  write(' '+IntToStr(Z[I,J]));
              writeln;
          end;
          writeln;
          writeln('Paths:');
          for I:=1 to N do
          begin
              for J:=I+1 to N do
                  if (Z[I,J] <> 0) then
                  begin
                       Idx:=Z[I,J];
                       for iwm:=1 to W[I,J]+1 do
                           write(' '+IntToStr(AtomS[Idx,iwm]));
                       writeln;
                  end;
          end;
          writeln;
     end else                      // selecting next poit for step
     begin
          for iwc := 1 to N do // try all points (vertexes)
          begin
               // this point is connected with current point, but it is not included in the road (way)
               if (G[waS, iwc] <> 0) and (not Incl[iwc]) then
               begin
                    Road[waP] := iwc;          // add this point into way
                    Incl[iwc] := TRUE;         // mark this point as included in the way
                    WayStep(iwc, waF, waP + 1,UseMarkAtom);
                    Incl[iwc] := FALSE;
                    Road[waP] := 0;
               end;
          end;
     end;
end;

function TShortestPath.DoAtomS(UseMarkAtom: boolean): integer;
var
   //NAS: integer;
   I,J,K,L,M,II,MM: integer;
   TEMP,IGARE: integer;
   T,IW,JW:T1;
   TEMPd: double;
   cpt: integer;
   TIGARE: array of integer;
begin
     { Define what to do }
     NAS:=0;
     for I:=0 to N do
         for J:=0 to N do
             Z[I,J]:=0;
     AllWays(UseMarkAtom);
     //    calculation of W and Z matrix with using Floid algorithm
     {for K := 1 to N do
         for I := 1 to N do
             for J := 1 to N do
                 if (I<>J) and (K<>I) and (K<>J) and (W[I,K]<>NOTCONNECTED) and (W[K,J]<>NOTCONNECTED) then
                 begin
                      TEMPd := Wd[I, K] + Wd[K, J];
                      if TEMPd < Wd[I, J] then
                      begin
                           Z[I, J] := Z[I, K];
                           Wd[I, J] := TEMPd;
                           W[I, J] := W[I, K] + W[K, J];
                      end;
                 end;}
     //    creation from matrix W vector T,
     //    keeping indexes (IW, JW) of matrix
     //for I := 1 to N do
     {for I := 1 to N-1 do
     begin //      for orgraph (NOR=1) J will be from 1 to N
           //for J := 1 to I do
           for J := I+1 to N do
           begin
                if (I<>J) then
                begin
                     //          conditions for sequence length of atom-bond-atom...
                     //          so trails shoter, than LenMin, and longer, than LenMax
                     //          are not take into account
                     if ((W[I, J] >= LenMin - 1) and (W[I, J] <= LenMax - 1)) then
                     //if ((W[I, J] >= LenMin) and (W[I, J] <= LenMax)) then
                     begin
                          if UseMarkAtom then
                          begin
                               if (Stereo[I] = 7) or (Stereo[J] = 7) then
                                  Inc(NAS);
                          end else
                              Inc(NAS);
                          if NAS > MaxSMF then
                             TShrtPthException.Create('ERROR ShrtPthAtm: MaxSMF overflow');
                          T[NAS] := W[I, J];
                          IW[NAS] := I;
                          JW[NAS] := J;
                     end;                          // LenMinMax
                end;                            // I<>J
           end;                              // J=1..I
     end;                                // I=1..N}
     for I:=1 to N-1 do
         for J:=I+1 to N do
         begin
              II:=Z[I,J];
              T[II]:=W[I,J];
              IW[II]:=I;
              JW[II]:=J;
         end;
     //    sort array T and simultanously IW and JW arrays
     {for II := 1 to NAS do
     begin
          I := IW[II];
          J := JW[II];
          write('Sequence '+IntToStr(II)+': '+IntToStr(I)+' ');
          MM := 1;
          repeat
                MM := MM + 1;
                I := Z[I, J];
                write(IntToStr(I)+' ');
          until I = J;
          writeln;
     end;}
     //
     {SetLength(TIGARE,N+1);
     for I := 1 to NAS - 1 do
         for J := I+1 to NAS do
             if (T[J] < T[I]) then
              begin
                   for II:=1 to T[I]+1 do
                       TIGARE[II]:=AtomS[I,II];
                   for II:=1 to T[J]+1 do
                       AtomS[I,II]:=AtomS[J,II];
                   for II:=1 to T[I]+1 do
                       AtomS[J,II]:=TIGARE[II];
                   IGARE := T[I];
                   T[I] := T[M];
                   T[M] := IGARE;
                   IGARE := IW[I];
                   IW[I] := IW[M];
                   IW[M] := IGARE;
                   IGARE := JW[I];
                   JW[I] := JW[M];
                   JW[M] := IGARE;
              end;}
     //    calculation of atom sequences along trail from atom I
     //    to amtom J for each pair I and J.
     //    atom sequences - matrix AS, NAS - number atom sequences,
     //    where eatch element AS - number of atoms in matrix G.
     for II := 1 to NAS do
     {begin
          I := IW[II];
          J := JW[II];
          AtomS[II, 1] := I;
          MM := 1;
          repeat
                MM := MM + 1;
                I := Z[I, J];
                AtomS[II, MM] := I;
          until I = J;}
          LEN[II] := T[II]+1;
     //end;
     //
     Result:=NAS;
end;

end.

