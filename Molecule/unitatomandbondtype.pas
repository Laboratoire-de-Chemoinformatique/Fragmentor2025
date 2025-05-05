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
unit unitAtomAndBondType;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, typinfo;

const
     ZMin=1;
     ZMax=105;
     BMin=1;
     BMax=99;//BMax=98;
     SMin=0;
     SMax=51;
     dZmin=0;
     dZmax=15;
     AtomSymbol: array [ZMin..ZMax] of string[2]=('H','He','Li','Be','B','C','N','O',
     'F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
     'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y',
     'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
     'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',
     'Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
     'At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es',
     'Fm','Md','No','Lr','*','XX');//'CD','CT','CB','CA','CO','CN','NI','NA','*','XX');
     DynAtomSymbol: array [dZmin..dZmax] of string [3] = ('','c+1','c-1','c+2','c-2','r+1','r-1','r+2','r-2','i+1','i-1','i+2','i-2','i+3','i-3','dXX');
                                                         //0, 1,    2,    3,    4,    5,    6,    7,    8,    9,    10,   11,   12,   13,   14,   15
     StereoSymbol: array [SMin..SMax] of string[3] = ('','n','u','m','c','e','z','s','a','n>u','n>c','n>e','n>z','n>s','n>a','u>n','u>c','u>e','u>z','u>s','u>a','c>n','c>u','c>e','c>z','c>s','c>a','e>n','e>u','e>c','e>z','e>s','e>a','z>n','z>u','z>c','z>e','z>s','z>a','s>n','s>u','s>c','s>e','s>z','s>a','a>n','a>u','a>c','a>e','a>z','a>s','ss');
                                                     //0, 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10  , 11  , 12  , 13  , 14  , 15  , 16  , 17  , 18  , 19  , 20  , 21  , 22  , 23  , 24  , 25  , 26  , 27  , 28  , 29  , 30  , 31  , 32  , 33  , 34  , 35  , 36  , 37  , 38  , 39  , 40  , 41  , 42  , 43  , 44  , 45  , 46  , 47  , 48  , 49  , 50  , 51  , 52
     AtomWeight: array [ZMin..ZMax] of double=(1.00794,4.00260,6.941,9.01218,10.81,
     12.011,14.0067,15.9994,18.998403,20.179,22.98977,24.30,26.98154,28.0855,
     30.97376,32.06,35.453,39.948,39.0983,40.08,44.9559,47.88,50.9415,51.996,
     54.9380,55.847,58.9332,55.69,63.546,65.38,69.72,72.59,74.9216,78.96,79.904,
     83.80,85.4678,87.62,88.9069,91.22,92.9064,95.94,98.00,101.07,102.9055,
     106.42,104.8682,112.41,114.82,118.69,121.75,127.60,126.9045,131.29,
     132.9054,137.33,138.9055,140.12,140.9077,144.24,145.00,150.36,151.96,
     157.25,158.9254,162.50,164.9304,167.26,168.9342,173.04,174.967,178.49,
     180.9479,183.85,186.207,190.2,192.22,195.08,196.9665,200.59,204.383,207.2,
     208.9804,209.00,210.00,222.00,223.00,226.0254,227.0278,232.0381,231.0359,
     238.0289,237.0482,244.00,243.00,247.00,247.00,251.00,252.00,257.00,258.00,
     259.00,260.00,-1.00,-10.00);
     {BondSymbol: array [BMin..BMax] of string[2]=('-','=','+','*','.',':','"',
     '~','&','1','2','3','4','5','6','7','8','9','D','T','YY');}
     {BondSymbol: array [BMin..BMax] of string[2]=('-','=','+','*','5','6','7','?','_','.',':','#','~','81','82','83','84','18','28','38','48','12','13','14','21','23','24','31','32','34','41','42','43','YY');}
                                                 //1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10, 11, 12, 13, 14 , 15 , 16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 , 32 , 33 , 34
     BondSymbol: array [BMin..BMax] of string[3]=('-','=','+','*','5','6','7','?','_','.',':','#','h','c','n','u','m','c','e','z','s','a','0>1','0>2','0>3','0>4','0>c','0>h','1>0','2>0','3>0','4>0','c>0','h>0','1>2','1>3','1>4','1>c','1>h','2>1','2>3','2>4','2>c','3>1','3>2','3>4','3>c','4>1','4>2','4>3','c>1','c>2','c>3','c>4','h>1','n>u','n>c','n>e','n>z','n>s','n>a','u>n','u>c','u>e','u>z','u>s','u>a','c>n','c>u','c>e','c>z','c>s','c>a','e>n','e>u','e>c','e>z','e>s','e>a','z>n','z>u','z>c','z>e','z>s','z>a','s>n','s>u','s>c','s>e','s>z','s>a','a>n','a>u','a>c','a>e','a>z','a>s', '0', 'YY');
                                                 //1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23  , 24  , 25  , 26  , 27  , 28  , 29  , 30  , 31  , 32  , 33  , 34  , 35  , 36  , 37  , 38  , 39  , 40  , 41  , 42  , 43  , 44  , 45  , 46  , 47  , 48  , 49  , 50  , 51  , 52  , 53  , 54  , 55  , 56  , 57  , 58  , 59  , 60  , 61  , 62  , 63  , 64  , 65  , 66  , 67  , 68  , 69  , 70  , 71  , 72  , 73  , 74  , 75  , 76  , 77  , 78  , 79  , 80  , 81  , 82  , 83  , 84  , 85  , 86  , 87  , 88  , 89  , 90  , 91  , 92  , 93  , 94  , 95  , 96  , 97  , 98  , 99


type
    APrp=array of double;
    AByt=array of byte;
    TZ=ZMin..ZMax;
    TdZ=dZmin..dZmax;
    TB=BMin..BMax;
    TSt=SMin..SMax;
    TS=string[6];

    EDynWrd=(extrabond,dynbond,dynbondstereo,bondstereo,atomstereo,atomcharge,atomradical,atomisotope,dyncharge,dynradical,dynisotope,dynatom,UNK);

    PRAtom=^RAtom;
    RAtom = record
            Z: TZ;  // atom symbol to int
            S: TS; // atom symbol
            W: double;
            P: APrp;
            I: AByt;
    end;

    PRBond=^RBond;
    RBond = record
            B: TB;//type
            S: TS;//character
            h,t: integer;//head, tail
            P: APrp;//double properties - your responsible to keep track of size and index
            I: AByt;//byte properties
    end;

function StringToEDynWrd(s: string): EDynWrd;
function EDynWrdToString(e: EDynWrd): string;
function AtomSymbolToInt(str: string): TZ;
function IntToAtomSymbol(Z: TZ): string;
function AtomSymbolToWeight(str: string): double;
function IntToWeight(Z: TZ): double;
function AtomStereoToInt(str: string): TSt;
function IntToAtomStereo(St: TSt): string;
function BondSymbolToInt(str: string): TB;
function IntToBondSymbol(B: TB): string;
function IntToTB(col3, col7: integer): TB;
function IntToS(col3, col7: integer): TS;
procedure TBtoInt(Bnd:TB; var col3, col7:integer);
function DynAtmAnnotation(wrd1,wrd2: string): string;
function IntToDynAtomSymbol(dZ:TdZ): string;
//function IntToDynAtomSymbol(dZ:byte): string;
function DynAtomSymbolToInt(str: string): TdZ;

implementation

function StringToEDynWrd(s: string): EDynWrd;
begin
  try
    Result:=EDynWrd(GetEnumValue(TypeInfo(EDynWrd),s));
  except
    Result:=UNK;
    //Exception.Create('WARNING: unable to parse keyword '+s+'.');
  end;
  {if (s='extrabond') then Result:=extrabond
  else if (s='dynbond') then Result:=dynbond
  else if (s='dynbondstereo') then Result:=dynbondstereo
  else if (s='bondstereo') then Result:=bondstereo
  else if (s='atomstereo') then Result:=atomstereo
  else if (s='dyncharge') then Result:=dyncharge
  else if (s='dynradical') then Result:=dynradical
  else if (s='dynisotope') then Result:=dynisotope
  else if (s='dynatom') then Result:=dynatom
  else Result:=UNK;}
  //if Result=UNK then writeln('WARNING: unable to parse keyword '+s+'.');
end;

function EDynWrdToString(e: EDynWrd): string;
begin
  try
    Result:=GetEnumName(TypeInfo(EDynWrd),Ord(e));
  except
    Result:='UNK';
    Exception.Create('WARNING: an unknown keyword was found.');
  end;
  {case e of
    extrabond : Result:='extrabond';
    dynbond   : Result:='dynbond';
    dynbondstereo : Result:='dynbondstereo';
    bondstereo : Result:='bondstereo';
    atomstereo : Result:='atomstereo';
    dyncharge : Result:='dyncharge';
    dynradical: Result:='dynradical';
    dynisotope: Result:='dynisotope';
    dynatom: Result:='dynatom';
  else
    Result:='UNK'
  end;
  if e=UNK then writeln('WARNING: an unknown keyword was found.');}
end;

function AtomSymbolToInt(str: string): TZ;
begin
  Result := Low(AtomSymbol);
  While ((Result<High(AtomSymbol)) and (AtomSymbol[Result]<>str)) do
        inc(Result);
     //Result := Low(AtomSymbol){%H-}-1;
     //repeat
     //      inc(Result);
     //until ((Result>=High(AtomSymbol)) or (AtomSymbol[Result]=str));
     if (AtomSymbol[Result]<>str) then Result:=ZMax;
     if (Result=ZMax) then writeln('WARNING: unable to parse atom '+str+'.');
end;

function IntToAtomSymbol(Z: TZ): string;
begin
     Result:=AtomSymbol[Z];
     if Z=ZMax then writeln('WARNING: an unknown atom was found.');
end;

function AtomSymbolToWeight(str: string): double;
begin
     Result:=AtomWeight[AtomSymbolToInt(str)];
end;

function IntToWeight(Z: TZ): double;
begin
     Result:=AtomWeight[Z];
end;

function AtomStereoToInt(str: string): TSt;
begin
  Result := Low(StereoSymbol);
  While ((Result<High(StereoSymbol)) and (StereoSymbol[Result]<>str)) do
        inc(Result);
  {Result := Low(StereoSymbol){%H-}-1;
  repeat
        inc(Result);
  until ((Result>=High(StereoSymbol)) or (StereoSymbol[Result]=str));}
  if (StereoSymbol[Result]<>str) then Result:=SMax;
  if Result=SMax then writeln('WARNING: unable to parse stereo key word '+str+'.');
end;

function IntToAtomStereo(St: TSt): string;
begin
  Result:=StereoSymbol[St];
  if St=SMax then writeln('WARNING: an unknown stereo keyword was found.');
end;

function BondSymbolToInt(str: string): TB;
begin
  Result := Low(BondSymbol);
  while ((Result<High(BondSymbol)) and (BondSymbol[Result]<>str)) do
        inc(Result);
     {Result := Low(BondSymbol){%H-}-1;
     repeat
           inc(Result);
     until ((Result>=High(BondSymbol)) or (BondSymbol[Result]=str));}
     if (BondSymbol[Result]<>str) then Result:=BMax;
     if (Result=BMax) then writeln('WARNING: impossible to parse bond '+str+'.');
end;

function IntToBondSymbol(B: TB): string;
begin
  if (B=BMax) then writeln('WARNING: an unknown bond was found');
  Result:=BondSymbol[B];
end;

function IntToTB(col3, col7: integer): TB;
begin
  if (col7 = 0) then
  begin
    Result := col3; //for case 1 to 9
    case col3 of
      50: Result := 10;  //bond .
      60: Result := 11;  //bond :
      70: Result := 12;  //bond #
      80: Result := 13;  //bond ~
    end;
  end
  else if (col3 = 1) then
  begin
    case col7 of
      8: Result := 14;   //bond 81
      -1: Result := 18;  //bond 18
      4: Result := 25;   //bond 21
      12: Result := 28;  //bond 31
      1: Result := 31;   //bond 41
    end;
  end
  else if (col3 = 2) then
  begin
    case col7 of
      4: Result := 15;   //bond 82
      -1: Result := 19;  //bond 28
      8: Result := 22;   //bond 12
      12: Result := 29;  //bond 32
      1: Result := 32;   //bond 42
    end;
  end
  else if (col3 = 3) then
  begin
    case col7 of
      12: Result := 16;  //bond 83
      -1: Result := 20;  //bond 38
      8: Result := 23;   //bond 13
      4: Result := 26;   //bond 23
      1: Result := 33;   //bond 43
    end;
  end
  else if (col3 = 4) then
  begin
    case col7 of
      1: Result := 17;   //bond 84
      -1: Result := 21;  //bond 48
      8: Result := 24;   //bond 14
      4: Result := 27;   //bond 24
      12: Result := 30;  //bond 34
    end;
  end
  else
  begin
    Result := BMax;
  end;
  writeln('WARNING: problem parsing column '+IntToStr(col3)+' and '+IntToStr(col7)+'.');
end;

procedure TBtoInt(Bnd:TB; var col3,col7:integer);
begin
  if(Bnd=10) then
  begin
    col3:=50;
    col7:=0;
  end
  else if(Bnd=11) then
  begin
    col3:=60;
    col7:=0;
  end
  else if(Bnd=12) then
  begin
    col3:=70;
    col7:=0;
  end
  else if(Bnd=13) then
  begin
    col3:=80;
    col7:=0;
  end
  else if(Bnd=14) then
  begin
    col3:=1;
    col7:=8;
  end
  else if(Bnd=15) then
  begin
    col3:=2;
    col7:=4;
  end
  else if(Bnd=16) then
  begin
    col3:=3;
    col7:=12;
  end
  else if(Bnd=17) then
  begin
    col3:=4;
    col7:=1;
  end
  else if(Bnd=18) then
  begin
    col3:=1;
    col7:=-1;
  end
  else if(Bnd=19) then
  begin
    col3:=2;
    col7:=-1;
  end
  else if(Bnd=20) then
  begin
    col3:=3;
    col7:=-1;
  end
  else if(Bnd=21) then
  begin
    col3:=4;
    col7:=-1;
  end
  else if(Bnd=22) then
  begin
    col3:=2;
    col7:=8;
  end
  else if(Bnd=23) then
  begin
    col3:=3;
    col7:=8;
  end
  else if(Bnd=24) then
  begin
    col3:=4;
    col7:=8;
  end
  else if(Bnd=25) then
  begin
    col3:=1;
    col7:=4;
  end
  else if(Bnd=26) then
  begin
    col3:=3;
    col7:=4;
  end
  else if(Bnd=27) then
  begin
    col3:=4;
    col7:=4;
  end
  else if(Bnd=28) then
  begin
    col3:=1;
    col7:=12;
  end
  else if(Bnd=29) then
  begin
    col3:=2;
    col7:=12;
  end
  else if(Bnd=30) then
  begin
    col3:=4;
    col7:=12;
  end
  else if(Bnd=31) then
  begin
    col3:=1;
    col7:=1;
  end
  else if(Bnd=32) then
  begin
    col3:=2;
    col7:=1;
  end
  else if(Bnd=33) then
  begin
    col3:=2;
    col7:=1;
  end;
  if Bnd=BMax then writeln('WARNING: an unknown bound was found.');
end;

function DynAtmAnnotation(wrd1, wrd2: string): string;
begin
  Result:='';
  if (wrd1='dyncharge') and (LeftStr(wrd2,1)='c') then Result:=wrd2;
  if (wrd1='dynradical') and (LeftStr(wrd2,1)='r') then Result:=wrd2;
  if (wrd1='dynisotope') and (LeftStr(wrd2,1)='i') then Result:=wrd2;
end;

function IntToDynAtomSymbol(dZ: TdZ): string;
begin
  Result:=DynAtomSymbol[dz];
end;

//function IntToDynAtomSymbol(dZ: byte): string;
//begin
//  Result:=DynAtomSymbol[dz];
//end;

function DynAtomSymbolToInt(str: string): TdZ;
var
  i: TdZ;
  bFound: Boolean;
begin
  i:=0;
  bFound:=False;
  while (i<dZmax) and (bFound=False) do
  begin
    if DynAtomSymbol[i]=str then bFound:=True
    else i:=i+1;
  end;
  if (bFound=True) then
    Result:=i
  else
    Result:=dZmax;
end;

function IntToS(col3, col7: integer): TS;
begin
  Result := BondSymbol[IntToTB(col3, col7)];
end;

end.

