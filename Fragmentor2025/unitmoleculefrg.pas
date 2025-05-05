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
unit UnitMoleculeFrg;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, contnrs, U_TYPE_GRAPHES, UnitMoleculeBase,
  unitAtomAndBondType, UnitAtomBase;

type

  recordmol = record
    AdynA, AdynB, Adynr, Adync, Adyni, AChg, AIso, ARad: TList;
  end;
  { TMoleculeFrg }

  TMoleculeFrg = class(TMoleculeBase)
  private
    fOffSetMarkAt: byte;
    fOffSetMarkBd: byte;
    fOffSetStereoAt: byte;
    fOffSetStereoBd: byte;
    fOffSetTopolBd: byte;
    fOffFormalCharge: byte;
    fOffRadicalAt: byte;
    fOffIsotopeAt: byte;
    fOffDynAtm: byte;
    fOffDynCharge: byte;
    fOffDynRadical: byte;
    fOffDynIsotope: byte;
    fOffCircuitBd: byte;
    fOffCircuitAt: byte;
    fOffAtomMap: byte;
    fAtISze: byte;
    fBdISze: byte;

  public
    constructor Create;
    destructor Destroy; override;
    function IsMarkedAt(Id: AtomID): boolean;
    function IsMarkedAt(A: PRAtom): boolean;
    function IsChargedAt(Id: AtomID): boolean;
    function IsChargedAt(A: PRAtom): boolean;
    function GetFormalCharge(Id: AtomID): byte;
    function GetFormalCharge(A: PRAtom): byte;
    function GetRadicalAt(Id: AtomID): byte;
    function GetRadicalAt(A: PRAtom): byte;
    function GetIsotopeAt(Id: AtomID): byte;
    function GetIsotopeAt(A: PRAtom): byte;
    function GetAtomMap(Id: AtomID): byte;    // Needed for the output
    function GetAtomMap(A: PRAtom): byte;
    function GetX(Id: AtomID): double;
    function GetX(A: PRAtom): double;
    function GetY(Id: AtomID): double;
    function GetY(A: PRAtom): double;
    function GetZ(Id: AtomID): double;
    function GetZ(A: PRAtom): double;
    function IsMarkedBd(Id: BondID): boolean;
    function IsMarkedBd(B: PRBond): boolean;
    function AtStereoParity(Id: AtomId): byte;
    function AtStereoParity(A: PRAtom): byte;
    function AtInCircuit(A: PRAtom): boolean;
    function BdInCricuit(B: PRBond): boolean;

    function AtDynAtm(Id: AtomId): byte;
    function AtDynAtm(A: PRAtom): byte;
    function IsDynAtom(Id: AtomId): boolean;
    function IsDynAtom(A: PRAtom): boolean;

    function AtDynCharge(Id: AtomId): byte;
    function AtDynCharge(A: PRAtom): byte;
    function AtDynRadical(Id: AtomId): byte;
    function AtDynRadical(A: PRAtom): byte;
    function AtDynIsotope(Id: AtomId): byte;
    function AtDynIsotope(A: PRAtom): byte;

    function AtDynAtmS(Id: AtomId): string;
    function AtDynAtmS(A: PRAtom): string;

    function AtDynChargeS(Id: AtomId): string;
    function AtDynChargeS(A: PRAtom): string;
    function AtDynRadicalS(Id: AtomId): string;
    function AtDynRadicalS(A: PRAtom): string;
    function AtDynIsotopeS(Id: AtomId): string;
    function AtDynIsotopeS(A: PRAtom): string;
    function BdStereo(Id: BondId): byte;
    function BdStereo(B: PRBond): byte;
    function BdTopol(Id: BondId): byte;
    function BdTopol(B: PRBond): byte;
    procedure LoadSDF(sdfstr: TStringList);// override;
    //procedure LoadSDF(sdfstr: TStringList; bMarkCircuit: Boolean=False);
    function {%H-}LoadSDFTStringList(sdfstr: TStringList; prpstr: string;
      msstr: TStringList): TFPStringHashTable;
    procedure LoadSDFAGroup(sdfstr:TStringList);
    function LoadSDFField(sdfstr: TStringList; prpstr: string;
      msstr: TStringList): TFPStringHashTable;
    function LoadSDFField(sdfstr: TStringList; prpstr: TStringList;
      msstr: TStringList): TFPStringHashTable;
    procedure LoadSDFField(sdfstr: TStringList; prpstr: string;
      msstr: TStringList; strhash: TFPStringHashTable);
    procedure LoadSDFField(sdfstr: TStringList; prpstr: TStringList;
      msstr: TStringList; strhash: TFPStringHashTable);
    procedure LoadSDFFieldBdColor(sdfstr: TStringList; prpstr: string);
    procedure AnnotateCycleAtm(bReset: boolean);
    procedure AnnotateCycleBnd(bReset: boolean);
    procedure SaveSDF(sdfout: string);
    function ReturnSdfTstringList(): TStringList;
    function InitRecord(): recordmol;
    procedure disposeRecord(recordm: recordmol);
  end;

implementation

{ TMoleculeFrg }

constructor TMoleculeFrg.Create;
begin
  //Atoms
  fOffSetMarkAt := 0;
  fOffSetStereoAt := 1;
  fOffFormalCharge := 2;
  fOffRadicalAt:=3;
  fOffIsotopeAt:=4;
  fOffDynAtm:=5;//3;
  fOffDynCharge:=6;
  fOffDynRadical:=7;
  fOffDynIsotope:=8;
  fOffCircuitAt:=9;
  fOffAtomMap:=10;
  //Bonds
  fOffSetMarkBd := 0;
  fOffSetStereoBd := 1;
  fOffSetTopolBd := 2;
  fOffCircuitBd:=3;
  fAtISze := 11;
  fBdISze := 4;

  inherited Create;
end;

destructor TMoleculeFrg.Destroy;
begin
  inherited Destroy;
end;

function TMoleculeFrg.IsMarkedAt(Id: AtomID): boolean;
begin
  Result := IsMarkedAt(AtmSet[Id]);
end;

function TMoleculeFrg.IsMarkedAt(A: PRAtom): boolean;
begin
  Result := (A^.I[fOffSetMarkAt] = 1);
end;

function TMoleculeFrg.IsChargedAt(Id: AtomID): boolean;
begin
  Result := IsChargedAt(AtmSet[Id]);
end;

function TMoleculeFrg.IsChargedAt(A: PRAtom): boolean;
begin
  if (A^.I[fOffFormalCharge] <> 0) then
    Result := True
  else
    Result := False;
end;

function TMoleculeFrg.GetFormalCharge(Id: AtomID): byte;
begin
  Result := GetFormalCharge(AtmSet[Id]);
end;

function TMoleculeFrg.GetFormalCharge(A: PRAtom): byte;
begin
  Result := A^.I[fOffFormalCharge];
end;

function TMoleculeFrg.GetRadicalAt(Id: AtomID): byte;
begin
  Result:=GetRadicalAt(AtmSet[Id]);
end;

function TMoleculeFrg.GetRadicalAt(A: PRAtom): byte;
begin
  Result:=A^.I[fOffRadicalAt];
end;

function TMoleculeFrg.GetIsotopeAt(Id: AtomID): byte;
begin
  Result:=GetIsotopeAt(AtmSet[Id]);
end;

function TMoleculeFrg.GetIsotopeAt(A: PRAtom): byte;
begin
  Result:=A^.I[fOffIsotopeAt];
end;

function TMoleculeFrg.GetAtomMap(Id: AtomID): byte;                                                        /// atom mapping
begin
  Result := GetAtomMap(AtmSet[Id]);
end;

function TMoleculeFrg.GetAtomMap(A: PRAtom): byte;
begin
  Result := A^.I[fOffAtomMap];
end;

function TMoleculeFrg.GetX(Id: AtomID): double;
begin
  Result := GetX(AtmSet[Id]);
end;

function TMoleculeFrg.GetX(A: PRAtom): double;
begin
  Result := A^.P[0];
end;

function TMoleculeFrg.GetY(Id: AtomID): double;
begin
  Result := GetY(AtmSet[Id]);
end;

function TMoleculeFrg.GetY(A: PRAtom): double;
begin
  Result := A^.P[1];
end;

function TMoleculeFrg.GetZ(Id: AtomID): double;
begin
  Result := GetZ(AtmSet[Id]);
end;

function TMoleculeFrg.GetZ(A: PRAtom): double;
begin
  Result := A^.P[2];
end;

function TMoleculeFrg.IsMarkedBd(Id: BondID): boolean;
begin
  Result := IsMarkedBd(BndSet[Id]);
end;

function TMoleculeFrg.IsMarkedBd(B: PRBond): boolean;
begin
  Result := (B^.I[fOffSetMarkBd] = 1);
end;

function TMoleculeFrg.AtStereoParity(Id: AtomId): byte;
begin
  Result := AtStereoParity(AtmSet[Id]);
end;

function TMoleculeFrg.AtStereoParity(A: PRAtom): byte;
begin
  Result := A^.I[fOffSetStereoAt];
end;

function TMoleculeFrg.AtInCircuit(A: PRAtom): boolean;
begin
  Result:=(A^.I[fOffCircuitAt]=1);
end;

function TMoleculeFrg.BdInCricuit(B: PRBond): boolean;
begin
  Result:=(B^.I[fOffCircuitBd]=1);
end;

function TMoleculeFrg.AtDynAtm(Id: AtomId): byte;
begin
  Result := AtDynAtm(AtmSet[Id]);
end;

function TMoleculeFrg.AtDynAtm(A: PRAtom): byte;
begin
  Result := A^.I[fOffDynAtm];
end;

function TMoleculeFrg.IsDynAtom(Id: AtomId): boolean;
begin
  Result:=IsDynAtom(AtmSet[Id]);
end;

function TMoleculeFrg.IsDynAtom(A: PRAtom): boolean;
begin
  If A^.I[fOffDynAtm]=0 then Result:=False
  else Result:=True;
end;

function TMoleculeFrg.AtDynCharge(Id: AtomId): byte;
begin
  Result := AtDynCharge(AtmSet[Id]);
end;

function TMoleculeFrg.AtDynCharge(A: PRAtom): byte;
begin
  Result := A^.I[fOffDynCharge];
end;

function TMoleculeFrg.AtDynRadical(Id: AtomId): byte;
begin
  Result := AtDynRadical(AtmSet[Id]);
end;

function TMoleculeFrg.AtDynRadical(A: PRAtom): byte;
begin
  Result := A^.I[fOffDynRadical];
end;

function TMoleculeFrg.AtDynIsotope(Id: AtomId): byte;
begin
  Result := AtDynIsotope(AtmSet[Id]);
end;

function TMoleculeFrg.AtDynIsotope(A: PRAtom): byte;
begin
  Result := A^.I[fOffDynIsotope];                           //LP
end;

function TMoleculeFrg.AtDynAtmS(Id: AtomId): string;
begin
  Result := AtDynAtmS(AtmSet[Id]);
end;

function TMoleculeFrg.AtDynAtmS(A: PRAtom): string;
var
  i: byte;
begin
  i := A^.I[fOffDynAtm];
  Result:=IntToDynAtomSymbol(i);
end;

function TMoleculeFrg.AtDynChargeS(Id: AtomId): string;
begin
  Result := AtDynChargeS(AtmSet[Id]);
end;

function TMoleculeFrg.AtDynChargeS(A: PRAtom): string;
var
  i: byte;
begin
  i := A^.I[fOffDynCharge];
  if (i=0) then Result:=''
  else if (i mod 2 <> 0) then Result:='c+'+IntToStr((i div 2)+1) //odd are positive
  else Result:='c'+IntToStr(-(i div 2));                         //even are negative
end;

function TMoleculeFrg.AtDynRadicalS(Id: AtomId): string;
begin
  Result := AtDynRadicalS(AtmSet[Id]);
end;

function TMoleculeFrg.AtDynRadicalS(A: PRAtom): string;
var
  i: byte;
begin
  i := A^.I[fOffDynRadical];
  if (i=0) then Result:=''
  else if (i mod 2 <> 0) then Result:='r+'+IntToStr((i div 2)+1) //odd are positive
  else Result:='r'+IntToStr(-(i div 2));                         //even are negative
end;

function TMoleculeFrg.AtDynIsotopeS(Id: AtomId): string;
begin
  Result := AtDynIsotopeS(AtmSet[Id]);
end;

function TMoleculeFrg.AtDynIsotopeS(A: PRAtom): string;
var
  i: integer;
begin
  i := A^.I[fOffDynIsotope];
  //if i>0 then writeln('BUG!'); //GM
  if (i=0) then Result:=''
  else if (i mod 2 <> 0) then Result:='i+'+IntToStr((i div 2)+1) //odd are positive
  else Result:='i'+IntToStr(-(i div 2));                         //even are negative
end;

function TMoleculeFrg.BdStereo(Id: BondId): byte;
begin
  Result := BdStereo(BndSet[Id]);
end;

function TMoleculeFrg.BdStereo(B: PRBond): byte;
begin
  Result := B^.I[fOffSetStereoBd];
end;

function TMoleculeFrg.BdTopol(Id: BondId): byte;
begin
  Result := BdTopol(BndSet[Id]);
end;

function TMoleculeFrg.BdTopol(B: PRBond): byte;
begin
  Result := B^.I[fOffSetTopolBd];
end;

{procedure TMoleculeFrg.LoadSDF(sdfstr: TStringList);
const
  maxbyte=255;
var
  i, j:   integer;
  LineNo: integer;
  PAt:    PRAtom;
  PBo:    PRBond;
  atmp:   CostMatrix;
  wtmp:   array of PRBond;
  s, t:   Node;
  M:      ArcNum;
  alist:  TList;
  pitem:  PRSAL;
  stmp :string;
  b : TB;
  e: EDynWrd;
  ival: integer;
  iival: byte;
  ErrCode: Word;
begin
  //User defined
  APrpSze := 0;
  ABytSze := fAtISze;
  BPrpSze := 0;
  BBytSze := fBdISze;
  //Read each line and create the base molecule
  LineNo  := 0;
  //---Line 1 : Molecule name---
  MolName := sdfstr[LineNo];
  //---Line 2 of MOLfile---
  //  User's first and last initials (UserI), program name (ProgN), date/time:
  //  M/D/Y,H:m (DatTim), dimensional codes (DimCo), scaling factors (ScaI, ScaR),
  //  energy (Ene) if modeling program input, internal registry number (RegNo)
  //  if input through MDL format:
  Inc(LineNo);
  //---Line 3 of MOLfile---
  // A line for comments. If no comment is entered, a blank line must be present.
  Inc(LineNo);
  //---Line 4 of MOLfile : The Counts Line---
  Inc(LineNo);
  p_NX := 0;
  p_NY := 0; //Molecular graphs are not bipartite a priori
  p_M  := 0;
  p_NX := int_readpos(sdfstr[LineNo], 1, 3);     // number of atoms
  p_M  := int_readpos(sdfstr[LineNo], 4, 3);      // number of bonds
  //Position 31 : Number of lines of additional properties
  //AddPrp:=int_readpos(sdfstr[LineNo],31,3);
  //---Lines 5-... of MOLfile : The Atom Block---
  AtmSet[0] := nil;
  for i := 1 to p_NX do
  begin
    Inc(LineNo);
    new(PAt);
    SetLength(PAt^.P, APrpSze);
    SetLength(PAt^.I, ABytSze);
    for j := Low(PAt^.P) to High(PAt^.P) do
      PAt^.P[j] := 0;
    for j := Low(PAt^.I) to High(PAt^.I) do
      PAt^.I[j] := 0;
    PAt^.S      := Trim(Copy(sdfstr[LineNo], 32, 3));        // atom symbol
    PAt^.Z      := AtomSymbolToInt(PAt^.S);
    //Only these infos are stored. Add whichever you would like to add.
    PAt^.I[fOffFormalCharge] := int_readpos(sdfstr[LineNo], 37, 3);
    // formal charge -> PAt^.I[0]
    PAt^.I[fOffSetMarkAt] := int_readpos(sdfstr[LineNo], 55, 3);
    //Unused field in V2000 used to mark atoms
    PAt^.I[fOffSetStereoAt] := int_readpos(sdfstr[LineNo], 40, 3); //Atom Stereo Parity

    AtmSet[i] := PAt;
  end;  // 1..p_NX
  //---Lines of MOLfile : The Bond Block---
  //Express the molecular graph as a simple matrix graph
  for s := 0 to MaxAtom + 1 do
    for t := 0 to MaxAtom + 1 do
      atmp[s, t] := 0;
  SetLength(wtmp, p_M + 1);
  BndSet[0] := nil;
  for i := 1 to p_M do
  begin
    Inc(LineNo);
    new(PBo);
    SetLength(PBo^.P, BPrpSze);
    SetLength(PBo^.I, BBytSze);
    for j := Low(PBo^.P) to High(PBo^.P) do
      PBo^.P[j] := 0;
    for j := Low(PBo^.I) to High(PBo^.I) do
      PBo^.I[j] := 0;
    PBo^.t      := int_readpos(sdfstr[LineNo], 1, 3);//Tail of bond
    PBo^.h      := int_readpos(sdfstr[LineNo], 4, 3);//Head of bond
    PBo^.B      := IntToTB(int_readpos(sdfstr[LineNo], 7, 3),int_readpos(sdfstr[LineNo], 19, 3));//Type of bond
    PBo^.S      := BondSymbol[PBo^.B];
    //Only these infos are stored. Add whichever you would like to add.
    PBo^.I[fOffSetMarkBd] := int_readpos(sdfstr[LineNo], 13, 3);
    //Unused field in V2000 used to mark bonds
    PBo^.I[fOffSetStereoBd] := int_readpos(sdfstr[LineNo], 10, 3); // Bond Stereo
    PBo^.I[fOffSetTopolBd] := int_readpos(sdfstr[LineNo], 16, 3); //Topology
    atmp[PBo^.t, PBo^.h] := i;
    atmp[PBo^.h, PBo^.t] := i;
    wtmp[i]     := PBo;
  end;
  //Store the graph in a packed format
  M := 0;
  for s := 1 to p_NX do
  begin
    p_HEAD[s] := M + 1;
    for t := 1 to p_NX do
      if (atmp[s, t] <> 0) then
      begin
        Inc(M);
        p_SUCC[M] := t;
        BndSet[M] := wtmp[atmp[s, t]];
      end;
  end;
  p_HEAD[p_NX + p_NY + 1] := M + 1;
  p_M := M;
  //Annotate bonds and atoms using groups for dynamical bonds
  alist:=LoadSDFSGroup(sdfstr);
  for i:=0 to alist.Count-1 do
  begin
    e:=UNK;
    pitem:=alist[i];
    if (pitem<>nil) then
    begin
      if (pitem^.asize=2) then
      begin
        s:=pitem^.alist[1];
        t:=pitem^.alist[2];
        b:=BondSymbolToInt(pitem^.aword); //if the input file is wrong a default symbol will be used
        PBo:=FindBond(s,t);
        if PBo<>nil then //if nil this is an attempt to break the connectivity table and it must be ignored
        begin
          PBo^.B:=b;
          PBo^.S:=IntToBondSymbol(b);
        end;
      end else if (pitem^.asize=1) then
      begin
        s:=pitem^.alist[1];
        PAt:=AtmSet[s];
        if (PAt<>nil) then
          if (upcase(pitem^.aword[1])<>pitem^.aword[1]) then //test that the dynamic atom coding is lower case else unkown coding
          begin
            //writeln(PAt^.S);
            e:=pitem^.etype;
            if e=atomstereo then
            begin
              PAt^.I[fOffSetStereoAt]:=AtomStereoToInt(pitem^.aword);
            end else if e=dynatom then
            begin
              PAt^.I[fOffDynAtm]:=DynAtomSymbolToInt(pitem^.aword);
              PAt^.S:=PAt^.S+IntToDynAtomSymbol(PAt^.I[fOffDynAtm]);
              //if PAt^.I[fOffDynAtm]=dZmax then WriteLn('WARNING: Unknown dynatom description');
            end else
            begin
              stmp:=Copy(pitem^.aword,2,4);
              Val(stmp,ival,ErrCode);
              if (ival>0) then iival:=2*(ival-1)+1 //odd for positive value
              else if (ival<0) then iival:=-2*ival //even for negative
              else iival:=0;
              if ErrCode<>0 then
              begin
                //Writeln('ERROR: reading sdf entry ' + MolName);
                Exception.Create('ERROR: reading sdf entry ' + MolName);
                halt(1);
              end;
              case e of
                dyncharge: PAt^.I[fOffDynCharge]:=iival;
                dynradical: PAt^.I[fOffDynRadical]:=iival;
                dynisotope: PAt^.I[fOffDynIsotope]:=iival;
                //dynatom: PAt^.S:=PAt^.S+pitem^.aword;
              end;
            end;
          end else
          begin //Creating abberations in the output fragments if the CGR syntax is wrong
            //WriteLn('WARNING: Illegal atom id into a dynatom description');
            e:=pitem^.etype;
            case e of
              atomstereo: PAt^.I[fOffSetStereoAt]:=maxbyte;
              dyncharge: PAt^.I[fOffDynCharge]:=maxbyte;
              dynradical: PAt^.I[fOffDynRadical]:=maxbyte;
              dynisotope: PAt^.I[fOffDynIsotope]:=maxbyte;
              dynatom: PAt^.S:=AtomSymbol[ZMax];
            end;
          end;
      end;
    end;
  end;
  //
  for i:=0 to alist.Count-1 do
  begin
    pitem:=alist[i];
    if pitem<>nil then
    begin
      Dispose(pitem);
      pitem:=nil;
    end;
  end;
  alist.Clear;
  FreeAndNil(alist);
end;}

procedure TMoleculeFrg.LoadSDF(sdfstr: TStringList);//; bMarkCircuit: Boolean=False);
const
  maxbyte=255;
var
  i, j:   integer;
  LineNo: integer;
  PAt:    PRAtom;
  PBo:    PRBond;
  atmp:   CostMatrix;
  wtmp:   array of PRBond;
  s, t:   Node;
  M:      ArcNum;
  alist:  TList;
  pitem:  PRSAL;
  stmp, strtest :string;
  b : TB;
  e: EDynWrd;
  ival: integer;
  iival: byte;
  ErrCode: Word;
  //
begin
  //User defined
  APrpSze := 3;    // Coordinates
  ABytSze := fAtISze;
  BPrpSze := 0;
  BBytSze := fBdISze;
  //Read each line and create the base molecule
  LineNo  := 0;
  //---Line 1 : Molecule name---
  MolName := sdfstr[LineNo];
  //---Line 2 of MOLfile---
  //  User's first and last initials (UserI), program name (ProgN), date/time:
  //  M/D/Y,H:m (DatTim), dimensional codes (DimCo), scaling factors (ScaI, ScaR),
  //  energy (Ene) if modeling program input, internal registry number (RegNo)
  //  if input through MDL format:
  Inc(LineNo);
  //---Line 3 of MOLfile---
  // A line for comments. If no comment is entered, a blank line must be present.
  Inc(LineNo);
  //---Line 4 of MOLfile : The Counts Line---
  Inc(LineNo);
  p_NX := 0;
  p_NY := 0; //Molecular graphs are not bipartite a priori
  p_M  := 0;
  p_NX := int_readpos(sdfstr[LineNo], 1, 3);     // number of atoms
  p_M  := int_readpos(sdfstr[LineNo], 4, 3);      // number of bonds
  //Position 31 : Number of lines of additional properties
  //AddPrp:=int_readpos(sdfstr[LineNo],31,3);
  //---Lines 5-... of MOLfile : The Atom Block---
  AtmSet[0] := nil;
  for i := 1 to p_NX do
  begin
    Inc(LineNo);
    new(PAt);
    SetLength(PAt^.P, APrpSze);
    SetLength(PAt^.I, ABytSze);
    for j := Low(PAt^.P) to High(PAt^.P) do
      PAt^.P[j] := 0;
    for j := Low(PAt^.I) to High(PAt^.I) do
      PAt^.I[j] := 0;

    // Coordinates
    j := 0;
    begin
      PAt^.P[j] := dbl_readpos(sdfstr[LineNo], 1, 10); //X
      Inc(j);
    end;
    begin
      PAt^.P[j] := dbl_readpos(sdfstr[LineNo], 11, 10);//Y
      Inc(j);
    end;
    begin
      PAt^.P[j] := dbl_readpos(sdfstr[LineNo], 21, 10);//Z
      Inc(j);
    end;

    //PAt^.P[X] := dbl_readpos(sdfstr[LineNo], 1, 10); //X
    PAt^.S      := Trim(Copy(sdfstr[LineNo], 32, 3));        // atom symbol
    PAt^.Z      := AtomSymbolToInt(PAt^.S);
    //Only these infos are stored. Add whichever you would like to add.
    PAt^.I[fOffFormalCharge] := int_readpos(sdfstr[LineNo], 37, 3);
    // formal charge -> PAt^.I[0]
    PAt^.I[fOffSetMarkAt] := int_readpos(sdfstr[LineNo], 55, 3);
    //Unused field in V2000 used to mark atoms
    PAt^.I[fOffSetStereoAt] := int_readpos(sdfstr[LineNo], 40, 3); //Atom Stereo Parity
    PAt^.I[fOffAtomMap] := int_readpos(sdfstr[LineNo], 62, 3); //Atom-atom mapping number
    strtest:=Copy(sdfstr[LineNo], 62,3);
    AtmSet[i] := PAt;
  end;  // 1..p_NX
  //---Lines of MOLfile : The Bond Block---
  //Express the molecular graph as a simple matrix graph
  for s := 0 to MaxAtom + 1 do
    for t := 0 to MaxAtom + 1 do
      atmp[s, t] := 0;
  SetLength(wtmp, p_M + 1);
  BndSet[0] := nil;
  for i := 1 to p_M do
  begin
    Inc(LineNo);
    new(PBo);
    SetLength(PBo^.P, BPrpSze);
    SetLength(PBo^.I, BBytSze);
    for j := Low(PBo^.P) to High(PBo^.P) do
      PBo^.P[j] := 0;
    for j := Low(PBo^.I) to High(PBo^.I) do
      PBo^.I[j] := 0;

    PBo^.t      := int_readpos(sdfstr[LineNo], 1, 3);//Tail of bond
    PBo^.h      := int_readpos(sdfstr[LineNo], 4, 3);//Head of bond
    PBo^.B      := IntToTB(int_readpos(sdfstr[LineNo], 7, 3),int_readpos(sdfstr[LineNo], 19, 3));//Type of bond
    PBo^.S      := BondSymbol[PBo^.B];
    //Only these infos are stored. Add whichever you would like to add.
    PBo^.I[fOffSetMarkBd] := int_readpos(sdfstr[LineNo], 13, 3);
    //Unused field in V2000 used to mark bonds
    PBo^.I[fOffSetStereoBd] := int_readpos(sdfstr[LineNo], 10, 3); // Bond Stereo
    PBo^.I[fOffSetTopolBd] := int_readpos(sdfstr[LineNo], 16, 3); //Topology
    atmp[PBo^.t, PBo^.h] := i;
    atmp[PBo^.h, PBo^.t] := i;
    wtmp[i]     := PBo;
  end;
  //Store the graph in a packed format
  M := 0;
  for s := 1 to p_NX do
  begin
    p_HEAD[s] := M + 1;
    for t := 1 to p_NX do
      if (atmp[s, t] <> 0) then
      begin
        Inc(M);
        p_SUCC[M] := t;
        BndSet[M] := wtmp[atmp[s, t]];
      end;
  end;
  p_HEAD[p_NX + p_NY + 1] := M + 1;
  p_M := M;
  //Annotate bonds and atoms using groups for dynamical bonds
  LoadSDFAGroup(sdfstr);
  alist:=LoadSDFSGroup(sdfstr);
  for i:=0 to alist.Count-1 do
  begin
    e:=UNK;
    pitem:=alist[i];
    if (pitem<>nil) then
    begin
      if (pitem^.asize=2) then
      begin
        s:=pitem^.alist[1];
        t:=pitem^.alist[2];
        b:=BondSymbolToInt(pitem^.aword); //if the input file is wrong a default symbol will be used
        PBo:=FindBond(s,t);
        if PBo<>nil then //if nil this is an attempt to break the connectivity table and it must be ignored
        begin
          PBo^.B:=b;
          PBo^.S:=IntToBondSymbol(b);
        end;
      end else if (pitem^.asize=1) then
      begin
        s:=pitem^.alist[1];
        PAt:=AtmSet[s];
        if (PAt<>nil) then
          //ADD HERE CHG RAD AND ISO
          if (upcase(pitem^.aword[1])<>pitem^.aword[1]) then //test that the dynamic atom coding is lower case else unkown coding
          begin
            //writeln(PAt^.S);
            e:=pitem^.etype;
            if e=atomstereo then
            begin
              PAt^.I[fOffSetStereoAt]:=AtomStereoToInt(pitem^.aword);
            end else if e=dynatom then
            begin
              PAt^.I[fOffDynAtm]:=DynAtomSymbolToInt(pitem^.aword);
              PAt^.S:=PAt^.S+IntToDynAtomSymbol(PAt^.I[fOffDynAtm]);
              //if PAt^.I[fOffDynAtm]=dZmax then WriteLn('WARNING: Unknown dynatom description');
            end else
            begin
              stmp:=Copy(pitem^.aword,2,4);
              Val(stmp,ival,ErrCode);
              if (ival>0) then iival:=2*(ival-1)+1 //odd for positive value
              else if (ival<0) then iival:=-2*ival //even for negative
              else iival:=0;
              if ErrCode<>0 then
              begin
                //Writeln('ERROR: reading sdf entry ' + MolName);
                Exception.Create('ERROR: reading sdf entry ' + MolName);
                halt(1);
              end;
              case e of
                dyncharge: PAt^.I[fOffDynCharge]:=iival;
                dynradical: PAt^.I[fOffDynRadical]:=iival;
                dynisotope: PAt^.I[fOffDynIsotope]:=iival;
                //dynatom: PAt^.S:=PAt^.S+pitem^.aword;
              end;
            end;
          end else
          begin //Creating abberations in the output fragments if the CGR syntax is wrong
            //WriteLn('WARNING: Illegal atom id into a dynatom description');
            e:=pitem^.etype;
            case e of
              atomstereo: PAt^.I[fOffSetStereoAt]:=maxbyte;
              dyncharge: PAt^.I[fOffDynCharge]:=maxbyte;
              dynradical: PAt^.I[fOffDynRadical]:=maxbyte;
              dynisotope: PAt^.I[fOffDynIsotope]:=maxbyte;
              dynatom: PAt^.S:=AtomSymbol[ZMax];
            end;
          end;
      end;
    end;
  end;
  //
  for i:=0 to alist.Count-1 do
  begin
    pitem:=alist[i];
    if pitem<>nil then
    begin
      Dispose(pitem);
      pitem:=nil;
    end;
  end;
  alist.Clear;
  FreeAndNil(alist);
end;

function TMoleculeFrg.LoadSDFTStringList(sdfstr: TStringList;
  prpstr: string; msstr: TStringList): TFPStringHashTable;
var
  strhash: TFPStringHashTable;
begin
  LoadSDF(sdfstr);
  strhash := TFPStringHashTable.Create;
  LoadSDFField(sdfstr, prpstr, msstr, strhash);
  Result := strhash;
end;

procedure TMoleculeFrg.LoadSDFAGroup(sdfstr: TStringList);
var
  i: integer;
  lne: string;
  function NextS(offset: integer): integer;
  var
    i1: integer;
    bStop: boolean;
  begin
    i1:=offset-1;
    bStop:=False;
    repeat
      Inc(i1);
      if (Pos('M  CHG',sdfstr[i1])>0) then bStop:=True;
      if (Pos('M  RAD',sdfstr[i1])>0) then bStop:=True;
      if (Pos('M  ISO',sdfstr[i1])>0) then bStop:=True;
      if (i1>=sdfstr.Count) then bStop:=True;
    until (bStop) or (i1>=(sdfstr.Count-1)) or (Pos('M  END',sdfstr[i1])>0);
    if (bStop) then NextS:=i1 else NextS:=-1;
  end;
  procedure ParseGeneric(lne: string; etype:EDynWrd);
  var
    i1: integer;
    loffset,i1max: integer;
    atId: integer;
    aInt: integer;
    aword: string;
    pitem: PRSAL;
    pAtm: PRAtom;

  begin
    i1max:=int_readpos(lne,7,3);
    loffset:=10;
    for i1:=1 to i1max do
    begin
      try
        atId:=int_readpos(lne,loffset+1,3);
        aInt:=int_readpos(lne,loffset+5,3);// aByte:=Byte(int_readpos(lne,loffset+4,3)) => SDF with another offset?
        pAtm:=AtmSet[atId];
        case etype of
          atomcharge: pAtm^.I[fOffFormalCharge]:=UnConvertFormalCharge(aInt);
          atomisotope: pAtm^.I[fOffIsotopeAt]:=Byte(aInt);
          atomradical: pAtm^.I[fOffRadicalAt]:=Byte(aInt);
        end;
        //
        loffset:=loffset+8;
      except
        raise Exception.Create('ERROR LoadSDFAGroup/ParseGeneric: improperly formated line in the SDF: '+lne);
        halt;
      end;
    end;
  end;
begin
  i:=0;
  i:=NextS(0);//index i points to the first line of interest
  if (i>0) and (i<sdfstr.Count) then begin
    repeat
      lne:=sdfstr[i];
      if (Pos('M  CHG',lne)>0) then
        ParseGeneric(sdfstr[i],atomcharge)
      else if (Pos('M  RAD',lne)>0) then
        ParseGeneric(sdfstr[i],atomradical)
      else if (Pos('M  ISO',lne)>0) then
        ParseGeneric(sdfstr[i],atomisotope);
      i:=NextS(i+1);
    until (i<0) or (i>=sdfstr.Count);
  end;
end;

function TMoleculeFrg.LoadSDFField(sdfstr: TStringList; prpstr: string;
  msstr: TStringList): TFPStringHashTable;
var
  strhash: TFPStringHashTable;
begin
  strhash := TFPStringHashTable.Create;
  LoadSDFField(sdfstr, prpstr, msstr, strhash);
  Result := strhash;
end;

function TMoleculeFrg.LoadSDFField(sdfstr: TStringList; prpstr: TStringList;
  msstr: TStringList): TFPStringHashTable;
var
  strhash: TFPStringHashTable;

begin
  strhash := TFPStringHashTable.Create;
  LoadSDFField(sdfstr, prpstr, msstr, strhash);
  Result := strhash;
end;

procedure TMoleculeFrg.LoadSDFField(sdfstr: TStringList; prpstr: TStringList;
  msstr: TStringList; strhash: TFPStringHashTable);
var
  i, indexof: integer;
  tmp:      string;
  ms_nb:    integer;//counts the number of microspecies and gives each of them a new key
  {debugstr: string;}
begin
  strhash.Clear;
  msstr.Clear;
  //For each line of the dictionary, retrieve the prop in the tstringlist of the pdf
  for i := 0 to prpstr.Count - 1 do
  begin //Test search
    indexof := sdfstr.IndexOf('>  <' + prpstr[i] + '>');
    tmp     := '';
    ms_nb   := 0;
    if (prpstr[i] = 'Default') then
      msstr.Add(prpstr[i])
    else
    begin
      if (indexof > 0) then
      begin
        while (Copy(sdfstr[indexof + 1], 0, 4) <> '>  <') and
          (Copy(sdfstr[indexof + 1], 0, 4) <> '$$$$') do
          // Loop over all the lines in the property field
        begin
          {debugstr := Trim(sdfstr[indexof + 1]);}
          if (Trim(sdfstr[indexof + 1]) <> '') then
          begin
            ms_nb := ms_nb + 1;
            tmp   := prpstr[i] + '_' + IntToStr(ms_nb);
            msstr.Add(tmp);
            strhash.Add(tmp, Trim(sdfstr[indexof + 1]));
          end;
          Inc(indexof);
        end;
      end;
    end;
    if (msstr.Count = 0) then
    begin
      Writeln('ERROR - atom colour indicated could not be found as an SDF field');
      Exception.Create('ERROR - atom colour indicated could not be found as an SDF field');
      Halt;
    end;
  end;
end;

procedure TMoleculeFrg.LoadSDFFieldBdColor(sdfstr: TStringList; prpstr:string);
var
  nAt, nBd: integer;
  i, BdId,PrpLne,LineNo: integer;
  t, h: integer;
  atmp: array of array of integer;
  SLtmp,SLTuple: TStringList;
  stmp: string;
  PB: PRBond;
begin
  if (prpstr<>'') and (prpstr<>'Default') then
  begin
    SLtmp:=TStringList.Create;
    SLTuple:=TStringList.Create;
    SLTuple.Delimiter:=':';
    nAt := int_readpos(sdfstr[3], 1, 3);     // number of atoms
    nBd := int_readpos(sdfstr[3], 4, 3);      // number of bonds
    PrpLne:=sdfstr.IndexOf('>  <'+prpstr+'>');
    If PrpLne>=0 then
    begin
      stmp:=sdfstr[PrpLne+1];
      SLtmp.DelimitedText:=stmp;
      for i:=1 to SLtmp.Count-1 do //ignore 1st column
      begin
        SLTuple.DelimitedText:=SLtmp[i];
        BdId:=StrToInt(SLTuple[0]);
        LineNo:=3+nAt+BdId; //locate and interpret the bond
        t := int_readpos(sdfstr[LineNo], 1, 3);//Tail of bond
        h := int_readpos(sdfstr[LineNo], 4, 3);//Head of bond
        PB:=FindBond(t,h);
        PB^.S:=SLTuple[1]; //substitute the bond symbol
      end;
    end else
    begin
      Writeln('ERROR - bond colour indicated could not be found as an SDF field');
      Exception.Create('ERROR - bond colour indicated could not be found as an SDF field');
      Halt;
    end;
    FreeAndNil(SLTuple);
    FreeAndNil(SLtmp);
  end;
end;

procedure TMoleculeFrg.AnnotateCycleAtm(bReset: boolean);
var
  CrcN: TNodeInfo;
  SzeN: Node;
  i: integer;
  PAt: PRAtom;
begin
  //Detect cycles and annotate atoms
  GetCircuitsNodes(CrcN,SzeN);//outputs an array of atom IDs
  for i:=1 to SzeN do
  begin
    PAt:=AtmSet[CrcN[i]];
    PAt^.I[fOffCircuitAt]:=1;
    if bReset then //If reset, restore original atom annotation
      PAt^.S:=IntToAtomSymbol(PAt^.Z)
    else //else annotate the atom as inside a cycle
      PAt^.S:=PAt^.S+'.';
  end;
end;

procedure TMoleculeFrg.AnnotateCycleBnd(bReset: boolean);
var
  CrcE: TArcBool;
  i: integer;
  PBo: PRBond;
begin
  //Detect cycles and annotate bonds
  GetCircuitsEdges(CrcE);//outputs a boolean array, one cell per arc.
  for i:=1 to p_M do
  begin
    PBo:=BndSet[i];
    if(CrcE[i]) then
    begin
      //Beware that bonds are bidrectional and for this reason are processed twice
      if PBo^.I[fOffCircuitBd]=0 then//Forbid to reprocess the same bond
      begin
        PBo^.I[fOffCircuitBd]:=1;
        if bReset then //If reset, restore original bond annotation
          PBo^.S:=IntToBondSymbol(PBo^.B)
        else //else annotate the bond as inside a cycle
          PBo^.S:=PBo^.S+'.'; //Beware that PBo^.S is large enough to accomodate extra character
      end;
      //Shall we use old notation from Vitaly or add a "." to the bond symbol
      {if (PBo^.B=1) then
      begin
        PBo^.B:=10;
        PBo^.S:=IntToBondSymbol(PBo^.B);
      end else if (PBo^.B=2) then
      begin
        PBo^.B:=11;
        PBo^.S:=IntToBondSymbol(PBo^.B);
      end else if (PBo^.B=3) then
      begin
        PBo^.B:=12;
        PBo^.S:=IntToBondSymbol(PBo^.B);
      end else
      begin
        if (PBo^.B<10) or (PBo^.B>12) then //Beware that bonds are bidrectional and for this reason are processed two times
          PBo^.S:=PBo^.S+IntToBondSymbol(13);
      end;}
    end else
      PBo^.I[fOffCircuitBd]:=0;
  end;
end;

procedure TMoleculeFrg.LoadSDFField(sdfstr: TStringList; prpstr: string;
  msstr: TStringList; strhash: TFPStringHashTable);
var
  strlist: TStringList;
begin
  strhash.Clear;
  strlist := convert_dict_list(prpstr);
  LoadSDFField(sdfstr, strlist, msstr, strhash);
end;

procedure TMoleculeFrg.SaveSDF(sdfout: string);
var
  SLtmp, SLtmp_2: TStringList;
  stmp: string;
  PAt: PRAtom;
  PBd: PRBond;
  btmp: array of array of boolean;
  i,j,LineNo, ibondord: integer;
  recordM: recordmol;
  Xbnd,Ybnd: double;

  function arrayProp(Lst: TStringList; typeprop: String; Aprop: TList): TStringList;
  var
    SList: string;
    i, j: integer;
    prop: byte;
  begin
    j := 0;
    if (Aprop.Count > 8) then      // Can't have info on more than 8 atoms in a line
    begin
      SList:= typeprop + Format('%0:3D', [8]);
    end else SList := typeprop + Format('%0:3D', [Aprop.Count]);

    for i:=0 to (Aprop.Count-1) do
    begin
      if ((i mod 7 = 0) and (i<>0)) then
      begin
        Lst.Add(SList);
        if (Aprop.Count - j  > 8) then
        begin
          SList := typeprop + Format('%0:3D', [8]);
        end else SList := typeprop + Format('%0:3D', [Aprop.Count-j]);
      end;
      case typeprop of
        'M  ISO': prop := GetIsotopeAt(Aprop.Items[i]);
        'M  CHG': prop := GetFormalCharge(Aprop.Items[i]);
        'M  RAD': prop := GetRadicalAt(Aprop.Items[i]);
      end;
      SList := SList + Format('%0:4D%1:4D', [FindAtom(Aprop.Items[i]), prop]);
      Inc(j);
    end;
    Lst.Add(SList);
    Result:= Lst;
  end;
begin
  SLtmp := TStringList.Create;
  //------Line 1 of MOLfile------------------------------------------------------
  LineNo := 0;
  SLtmp.Add(MolName);
  //------Line 2 of MOLfile------------------------------------------------------
  //  User's first and last initials (UserI), program name (ProgN), date/time:
  //  M/D/Y,H:m (DatTim), dimensional codes (DimCo), scaling factors (ScaI, ScaR),
  //  energy (Ene) if modeling program input, internal registry number (RegNo)
  //  if input through MDL format:
  Inc(LineNo);
  SLtmp.Add('LPISIDA'+DateToStr(Date)+'111NA'+MolName);
  //------Line 3 of MOLfile------------------------------------------------------
  // A line for comments. If no comment is entered, a blank line must be present.
  Inc(LineNo);
  SLtmp.Add('');
  //------Line 4 of MOLfile : The Counts Line------------------------------------
  Inc(LineNo);
  stmp:=Format('%0:3d%1:3d',[nAtom,round(nBonds/2)])+'  0  0  0  0             999 V2000';
  SLtmp.Add(stmp);
  j := 0;
  //Position 7 : Number of atom lists
  //begin atmL:=int_readpos(sdfstr[LineNo],7,3); inc(j); end;   // number of atom lists
  //position 10 is obsolete
  //Position 13 : Chiral flag
  //begin cf:=int_readpos(sdfstr[LineNo],13,3); inc(j); end;  //chiral flag

  //Position 16 : Number of stext entries (ISIS)
  //begin stxt:=int_readpos(sdfstr[LineNo],16,3); inc(j); end;  //number of stext entries (ISIS)
  //position 19 is obsolete
  //position 22 is obsolete
  //position 25 is obsolete
  //position 28 is obsolete

  //Position 31 : Number of lines of additional properties
  //begin ap:=int_readpos(sdfstr[LineNo],31,3); inc(j); end;  // number of lines of additional properties

  //------Lines 5-... of MOLfile : The Atom Block--------------------------------

  // Initiate the arrays
  recordM := InitRecord();
  for i := 1 to p_NX do
  begin
    PAt:=AtmSet[i];
    Inc(LineNo);
    stmp:=Format('%0:10.4F%1:10.4F%2:10.4F',[GetX(i),GetY(i),GetZ(i)]);
    stmp:=stmp+' '+Format('%-3S',[PAt^.S]);        // atom symbol
    stmp:=stmp+'  0  0  0  0  0  0  0  0  0';
    stmp:=stmp+Format('%3D',[PAt^.I[fOffAtomMap]]);//Atom mapping
    stmp:=stmp+'  0  0';
    SLtmp.Add(stmp);
  end;  // 1..p_NX
  //------Lines of MOLfile : The Bond Block--------------------------------------
  //ii:=natom;
  SetLength(btmp,nAtom+1,nAtom+1);
  for i:=1 to nAtom do
    for j:=1 to nAtom do
      btmp[i,j]:=False;
  for i:=1 to nBonds do
  begin
    PBd:=BndSet[i];
    if PBd<>nil then
      if btmp[PBd^.h,PBd^.t]=False then
      begin
        if (PBd^.B > 22) and (PBd^.B < 98) then
        begin
          ibondord:=8;         //DynBnd
        end else ibondord:=PBd^.B;
        stmp:=Format('%0:3D%1:3D%2:3D%3:3D',[PBd^.t,PBd^.h,ibondord, BdStereo(PBd)])+ '  0'+
        Format('%0:3D', [BdTopol(PBd)])+'  0';
        SLtmp.Add(stmp);
        btmp[PBd^.h,PBd^.t]:=True;
        btmp[PBd^.t,PBd^.h]:=True;
      end;
  end;

  // M CHG, M ISO, M RAD
  if (recordM.AChg.Count > 0) then
  begin
    SLtmp_2 := TStringList.Create;
    SLtmp_2 := arrayProp(SLtmp_2, 'M  CHG', recordM.AChg);
    for i:=0 to SLtmp_2.Count -1 do
      begin
        SLtmp.Add(SLtmp_2[i]);
      end;
    FreeAndNil(SLtmp_2);
  end;

  if (recordM.AIso.Count >0) then
  begin
    SLtmp_2 := TStringList.Create;
    SLtmp_2 := arrayProp(SLtmp_2, 'M  ISO', recordM.AIso);
    for i:=0 to (SLtmp_2.Count -1) do
      begin
        SLtmp.Add(SLtmp_2[i]);
      end;
    FreeAndNil(SLtmp_2);
  end;

  if (recordM.ARad.Count > 0) then
  begin
    SLtmp_2 := TStringList.Create;
    SLtmp_2 := arrayProp(SLtmp_2, 'M  RAD', recordM.ARad);
    for i:=0 to SLtmp_2.Count -1 do
      begin
        SLtmp.Add(SLtmp_2[i]);
      end;
    FreeAndNil(SLtmp_2);
  end;

  // Dyn atom and bonds
  if ((recordM.AdynA.Count>0) or (recordM.AdynB.Count>0)) then
  begin
    // First line
    stmp:='M  STY' + Format('%0:3D', [recordM.AdynA.Count+recordM.AdynB.Count]);
    for j:=1 to (recordM.AdynA.Count+recordM.AdynB.Count) do
    begin
      stmp:= stmp + Format('%0:4D', [j]) + ' DAT';
    end;
    SLtmp.Add(stmp);

  //DynAtom
    j:=1;
    //DynIso
    for i:=0 to (recordM.Adyni.Count-1) do
    begin
      stmp := 'M  SAL' + Format('%0:4D', [j]) + '  1' + Format('%0:4D', [FindAtom(recordM.Adyni[i])]);
      SLtmp.Add(stmp);
      stmp := 'M  SDT' + Format('%0:4D', [j]) + ' dynatom';
      SLtmp.Add(stmp);
      stmp := 'M  SDD' + Format('%0:4D%1:1S%2:10.4F%3:10.4F', [j, ' ',GetX(FindAtom(recordM.Adyni[i]))+0.5, GetY(FindAtom(recordM.Adyni[i]))+0.5])+ '    DAU   ALL  0       0';
      SLtmp.Add(stmp);
      stmp := 'M  SED' + Format('%0:4D%1:4S', [j, AtDynIsotopeS(recordM.Adyni.Items[i])]);
      SLtmp.Add(stmp);
      Inc(j);
    end;
    //DynChg
    for i:=0 to (recordM.Adync.Count-1) do
    begin
      stmp := 'M  SAL' + Format('%0:4D', [j]) + '  1' + Format('%0:4D', [FindAtom(recordM.Adync[i])]);
      SLtmp.Add(stmp);
      stmp := 'M  SDT' + Format('%0:4D', [j]) + ' dynatom';
      SLtmp.Add(stmp);
      stmp := 'M  SDD' + Format('%0:4D%1:1S%2:10.4F%3:10.4F', [j, ' ',GetX(FindAtom(recordM.Adync[i]))+0.5, GetY(FindAtom(recordM.Adync[i]))+0.5])+ '    DAU   ALL  0       0';
      SLtmp.Add(stmp);
      stmp := 'M  SED' + Format('%0:4D%1:4S', [j, AtDynChargeS(recordM.Adync.Items[i])]);
      SLtmp.Add(stmp);
      Inc(j);
    end;
    //DynRad
    for i:=0 to (recordM.Adynr.Count-1) do
    begin
      stmp := 'M  SAL' + Format('%0:4D', [j]) + '  1' + Format('%0:4D', [FindAtom(recordM.Adynr[i])]);
      SLtmp.Add(stmp);
      stmp := 'M  SDT' + Format('%0:4D', [j]) + ' dynatom';
      SLtmp.Add(stmp);
      stmp := 'M  SDD' + Format('%0:4D%1:1S%2:10.4F%3:10.4F', [j, ' ',GetX(FindAtom(recordM.Adynr[i]))+0.5, GetY(FindAtom(recordM.Adynr[i]))+0.5])+ '    DAU   ALL  0       0';
      SLtmp.Add(stmp);
      stmp := 'M  SED' + Format('%0:4D%1:4S', [j, AtDynRadicalS(recordM.Adynr.Items[i])]);
      SLtmp.Add(stmp);
      Inc(j);
    end;

  //DynBond
    for i:=0 to (recordM.AdynB.Count-1) do
    begin
      stmp := 'M  SAL' + Format('%0:4D', [j]) + '  2' + Format('%0:4D%1:4D', [BndSet[FindBond(recordM.AdynB.Items[i])]^.t, BndSet[FindBond(recordM.AdynB.Items[i])]^.h]);
      SLtmp.Add(stmp);
      stmp := 'M  SDT' + Format('%0:4D', [j]) + ' dynbond';
      SLtmp.Add(stmp);
      Xbnd := (GetX(BndSet[FindBond(recordM.AdynB.Items[i])]^.t) + GetX(BndSet[FindBond(recordM.AdynB.Items[i])]^.h))/2;
      Ybnd := (GetY(BndSet[FindBond(recordM.AdynB.Items[i])]^.t) + GetY(BndSet[FindBond(recordM.AdynB.Items[i])]^.h))/2;
      stmp := 'M  SDD' + Format('%0:4D%1:1S%2:10.4F%3:10.4F', [j, ' ',Xbnd, Ybnd]) + '    DAU   ALL  0       0';
      SLtmp.Add(stmp);
      stmp := 'M  SED' + Format('%0:4D%1:4S', [j, BndSet[FindBond(recordM.AdynB.Items[i])]^.S]);
      SLtmp.Add(stmp);
      Inc(j);
    end;
  end;

  //End molecule
  SLtmp.Add('M  END');
  SLtmp.SaveToFile(sdfout);
  FreeAndNil(SLtmp);
  //disposeRecord(recordM);
  FreeAndNil(recordm.AChg);
  FreeAndNil(recordm.AdynA);
  FreeAndNil(recordm.AdynB);
  FreeAndNil(recordm.Adync);
  FreeAndNil(recordm.Adyni);
  FreeAndNil(recordm.Adynr);
  FreeAndNil(recordm.AIso);
  FreeAndNil(recordm.ARad);
end;

{Output a TStringList, you will need to clean it.}
function TMoleculeFrg.ReturnSdfTstringList(): TStringList;
var
  SLtmp_2: TStringList;
  stmp: string;
  PAt: PRAtom;
  PBd: PRBond;
  btmp: array of array of boolean;
  i,j,LineNo, ibondord: integer;
  recordM: recordmol;
  Xbnd,Ybnd: double;

  function arrayProp(Lst: TStringList; typeprop: String; Aprop: TList): TStringList;
  var
    SList: string;
    i, j: integer;
    prop: byte;
  begin
    j := 0;
    if (Aprop.Count > 8) then      // Can't have info on more than 8 atoms in a line
    begin
      SList:= typeprop + Format('%0:3D', [8]);
    end else SList := typeprop + Format('%0:3D', [Aprop.Count]);

    for i:=0 to (Aprop.Count-1) do
    begin
      if ((i mod 7 = 0) and (i<>0)) then
      begin
        Lst.Add(SList);
        if (Aprop.Count - j  > 8) then
        begin
          SList := typeprop + Format('%0:3D', [8]);
        end else SList := typeprop + Format('%0:3D', [Aprop.Count-j]);
      end;
      case typeprop of
        'M  ISO': prop := GetIsotopeAt(Aprop.Items[i]);
        'M  CHG': prop := GetFormalCharge(Aprop.Items[i]);
        'M  RAD': prop := GetRadicalAt(Aprop.Items[i]);
      end;
      SList := SList + Format('%0:4D%1:4D', [FindAtom(Aprop.Items[i]), prop]);
      Inc(j);
    end;
    Lst.Add(SList);
    Result:= Lst;
  end;
begin
  Result := TStringList.Create;
  //------Line 1 of MOLfile------------------------------------------------------
  LineNo := 0;
  Result.Add(MolName);
  //------Line 2 of MOLfile------------------------------------------------------
  //  User's first and last initials (UserI), program name (ProgN), date/time:
  //  M/D/Y,H:m (DatTim), dimensional codes (DimCo), scaling factors (ScaI, ScaR),
  //  energy (Ene) if modeling program input, internal registry number (RegNo)
  //  if input through MDL format:
  Inc(LineNo);
  Result.Add('LPISIDA'+DateToStr(Date)+'111NA'+MolName);
  //------Line 3 of MOLfile------------------------------------------------------
  // A line for comments. If no comment is entered, a blank line must be present.
  Inc(LineNo);
  Result.Add('');
  //------Line 4 of MOLfile : The Counts Line------------------------------------
  Inc(LineNo);
  stmp:=Format('%0:3d%1:3d',[nAtom,round(nBonds/2)])+'  0  0  0  0             999 V2000';
  Result.Add(stmp);
  j := 0;
  //Position 7 : Number of atom lists
  //begin atmL:=int_readpos(sdfstr[LineNo],7,3); inc(j); end;   // number of atom lists
  //position 10 is obsolete
  //Position 13 : Chiral flag
  //begin cf:=int_readpos(sdfstr[LineNo],13,3); inc(j); end;  //chiral flag

  //Position 16 : Number of stext entries (ISIS)
  //begin stxt:=int_readpos(sdfstr[LineNo],16,3); inc(j); end;  //number of stext entries (ISIS)
  //position 19 is obsolete
  //position 22 is obsolete
  //position 25 is obsolete
  //position 28 is obsolete

  //Position 31 : Number of lines of additional properties
  //begin ap:=int_readpos(sdfstr[LineNo],31,3); inc(j); end;  // number of lines of additional properties

  //------Lines 5-... of MOLfile : The Atom Block--------------------------------

  // Initiate the arrays
  recordM := InitRecord();
  for i := 1 to p_NX do
  begin
    PAt:=AtmSet[i];
    Inc(LineNo);
    stmp:=Format('%0:10.4F%1:10.4F%2:10.4F',[GetX(i),GetY(i),GetZ(i)]);
    stmp:=stmp+' '+Format('%-3S',[PAt^.S]);        // atom symbol
    stmp:=stmp+'  0  0  0  0  0  0  0  0  0';
    stmp:=stmp+Format('%3D',[PAt^.I[fOffAtomMap]]);//Atom mapping
    stmp:=stmp+'  0  0';
    Result.Add(stmp);
  end;  // 1..p_NX
  //------Lines of MOLfile : The Bond Block--------------------------------------
  //ii:=natom;
  SetLength(btmp,nAtom+1,nAtom+1);
  for i:=1 to nAtom do
    for j:=1 to nAtom do
      btmp[i,j]:=False;
  for i:=1 to nBonds do
  begin
    PBd:=BndSet[i];
    if PBd<>nil then
      if btmp[PBd^.h,PBd^.t]=False then
      begin
        if (PBd^.B > 22) and (PBd^.B < 98) then
        begin
          ibondord:=8;         //DynBnd
        end else ibondord:=PBd^.B;
        stmp:=Format('%0:3D%1:3D%2:3D%3:3D',[PBd^.t,PBd^.h,ibondord, BdStereo(PBd)])+ '  0'+
        Format('%0:3D', [BdTopol(PBd)])+'  0';
        Result.Add(stmp);
        btmp[PBd^.h,PBd^.t]:=True;
        btmp[PBd^.t,PBd^.h]:=True;
      end;
  end;

  // M CHG, M ISO, M RAD
  if (recordM.AChg.Count > 0) then
  begin
    SLtmp_2 := TStringList.Create;
    SLtmp_2 := arrayProp(SLtmp_2, 'M  CHG', recordM.AChg);
    for i:=0 to SLtmp_2.Count -1 do
      begin
        Result.Add(SLtmp_2[i]);
      end;
    FreeAndNil(SLtmp_2);
  end;

  if (recordM.AIso.Count >0) then
  begin
    SLtmp_2 := TStringList.Create;
    SLtmp_2 := arrayProp(SLtmp_2, 'M  ISO', recordM.AIso);
    for i:=0 to (SLtmp_2.Count -1) do
      begin
        Result.Add(SLtmp_2[i]);
      end;
    FreeAndNil(SLtmp_2);
  end;

  if (recordM.ARad.Count > 0) then
  begin
    SLtmp_2 := TStringList.Create;
    SLtmp_2 := arrayProp(SLtmp_2, 'M  RAD', recordM.ARad);
    for i:=0 to SLtmp_2.Count -1 do
      begin
        Result.Add(SLtmp_2[i]);
      end;
    FreeAndNil(SLtmp_2);
  end;

  // Dyn atom and bonds
  if ((recordM.AdynA.Count>0) or (recordM.AdynB.Count>0)) then
  begin
    // First line
    stmp:='M  STY' + Format('%0:3D', [recordM.AdynA.Count+recordM.AdynB.Count]);
    for j:=1 to (recordM.AdynA.Count+recordM.AdynB.Count) do
    begin
      stmp:= stmp + Format('%0:4D', [j]) + ' DAT';
    end;
    Result.Add(stmp);

  //DynAtom
    j:=1;
    //DynIso
    for i:=0 to (recordM.Adyni.Count-1) do
    begin
      stmp := 'M  SAL' + Format('%0:4D', [j]) + '  1' + Format('%0:4D', [FindAtom(recordM.Adyni[i])]);
      Result.Add(stmp);
      stmp := 'M  SDT' + Format('%0:4D', [j]) + ' dynatom';
      Result.Add(stmp);
      stmp := 'M  SDD' + Format('%0:4D%1:1S%2:10.4F%3:10.4F', [j, ' ',GetX(FindAtom(recordM.Adyni[i]))+0.5, GetY(FindAtom(recordM.Adyni[i]))+0.5])+ '    DAU   ALL  0       0';
      Result.Add(stmp);
      stmp := 'M  SED' + Format('%0:4D%1:4S', [j, AtDynIsotopeS(recordM.Adyni.Items[i])]);
      Result.Add(stmp);
      Inc(j);
    end;
    //DynChg
    for i:=0 to (recordM.Adync.Count-1) do
    begin
      stmp := 'M  SAL' + Format('%0:4D', [j]) + '  1' + Format('%0:4D', [FindAtom(recordM.Adync[i])]);
      Result.Add(stmp);
      stmp := 'M  SDT' + Format('%0:4D', [j]) + ' dynatom';
      Result.Add(stmp);
      stmp := 'M  SDD' + Format('%0:4D%1:1S%2:10.4F%3:10.4F', [j, ' ',GetX(FindAtom(recordM.Adync[i]))+0.5, GetY(FindAtom(recordM.Adync[i]))+0.5])+ '    DAU   ALL  0       0';
      Result.Add(stmp);
      stmp := 'M  SED' + Format('%0:4D%1:4S', [j, AtDynChargeS(recordM.Adync.Items[i])]);
      Result.Add(stmp);
      Inc(j);
    end;
    //DynRad
    for i:=0 to (recordM.Adynr.Count-1) do
    begin
      stmp := 'M  SAL' + Format('%0:4D', [j]) + '  1' + Format('%0:4D', [FindAtom(recordM.Adynr[i])]);
      Result.Add(stmp);
      stmp := 'M  SDT' + Format('%0:4D', [j]) + ' dynatom';
      Result.Add(stmp);
      stmp := 'M  SDD' + Format('%0:4D%1:1S%2:10.4F%3:10.4F', [j, ' ',GetX(FindAtom(recordM.Adynr[i]))+0.5, GetY(FindAtom(recordM.Adynr[i]))+0.5])+ '    DAU   ALL  0       0';
      Result.Add(stmp);
      stmp := 'M  SED' + Format('%0:4D%1:4S', [j, AtDynRadicalS(recordM.Adynr.Items[i])]);
      Result.Add(stmp);
      Inc(j);
    end;

  //DynBond
    for i:=0 to (recordM.AdynB.Count-1) do
    begin
      stmp := 'M  SAL' + Format('%0:4D', [j]) + '  2' + Format('%0:4D%1:4D', [BndSet[FindBond(recordM.AdynB.Items[i])]^.t, BndSet[FindBond(recordM.AdynB.Items[i])]^.h]);
      Result.Add(stmp);
      stmp := 'M  SDT' + Format('%0:4D', [j]) + ' dynbond';
      Result.Add(stmp);
      Xbnd := (GetX(BndSet[FindBond(recordM.AdynB.Items[i])]^.t) + GetX(BndSet[FindBond(recordM.AdynB.Items[i])]^.h))/2;
      Ybnd := (GetY(BndSet[FindBond(recordM.AdynB.Items[i])]^.t) + GetY(BndSet[FindBond(recordM.AdynB.Items[i])]^.h))/2;
      stmp := 'M  SDD' + Format('%0:4D%1:1S%2:10.4F%3:10.4F', [j, ' ',Xbnd, Ybnd]) + '    DAU   ALL  0       0';
      Result.Add(stmp);
      stmp := 'M  SED' + Format('%0:4D%1:4S', [j, BndSet[FindBond(recordM.AdynB.Items[i])]^.S]);
      Result.Add(stmp);
      Inc(j);
    end;
  end;

  //End molecule
  Result.Add('M  END');
  //Result:=SLtmp;

  //disposeRecord(recordM);
  FreeAndNil(recordm.AChg);
  FreeAndNil(recordm.AdynA);
  FreeAndNil(recordm.AdynB);
  FreeAndNil(recordm.Adync);
  FreeAndNil(recordm.Adyni);
  FreeAndNil(recordm.Adynr);
  FreeAndNil(recordm.AIso);
  FreeAndNil(recordm.ARad);
end;


function TMoleculeFrg.InitRecord(): recordmol;
// Gives back the id of the atoms witch contains dynAtoms, dynBond, Charges, Isotopoes and Radicals.
var
  recordM: recordmol;
  i, j: integer;
  PAt: PRAtom;
  PBd: PRBond;
  btmp: array of array of boolean;

begin
  recordM.AdynA  := TList.Create;
  recordM.AdynB  := TList.Create;
  recordM.Adynr  := TList.Create;
  recordM.Adyni  := TList.Create;
  recordM.Adync  := TList.Create;
  recordM.ARad   := TList.Create;
  recordM.AIso   := TList.Create;
  recordM.AChg   := TList.Create;

  // Atom record
  for i := 1 to p_NX do
    begin
      PAt:=AtmSet[i];
      if (IsDynAtom(PAt))         then recordM.AdynA.Add(PAt);
      if (AtDynCharge(Pat) <> 0)  then recordM.Adync.Add(PAt);
      if (AtDynIsotope(Pat) <> 0) then recordM.Adyni.Add(PAt);
      if (AtDynRadical(Pat) <> 0) then recordM.Adynr.Add(PAt);
      if (IsChargedAt(PAt))       then recordM.AChg.Add(PAt);
      if (GetIsotopeAt(PAt) <> 0) then recordM.AIso.Add(PAt);
      if (GetRadicalAt(PAt) <> 0) then recordM.ARad.Add(PAt);
    end;

  // Bond record
  SetLength(btmp,nAtom+1,nAtom+1);
  for i:=1 to nAtom do
    for j:=1 to nAtom do
      btmp[i,j]:=False;

  for i:=1 to nBonds do
  begin
    PBd:=BndSet[i];
    if PBd<>nil then
      if btmp[PBd^.h,PBd^.t]=False then
      begin
        if (PBd^.B > 22) and (PBd^.B < 98) then
        begin
          btmp[PBd^.h,PBd^.t]:=True;
          btmp[PBd^.t,PBd^.h]:=True;
          recordM.AdynB.Add(PBd);
        end;
      end;
  end;
  Result:=recordM;
end;

procedure TMoleculeFrg.disposeRecord(recordm: recordmol);
var
  i:integer;
begin
  FreeAndNil(recordm.AChg);
  FreeAndNil(recordm.AdynA);
  FreeAndNil(recordm.AdynB);
  FreeAndNil(recordm.Adync);
  FreeAndNil(recordm.Adyni);
  FreeAndNil(recordm.Adynr);
  FreeAndNil(recordm.AIso);
  FreeAndNil(recordm.ARad);
end;

end.
















