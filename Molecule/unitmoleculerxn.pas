{ Fragmentor of the ISIDA Project

  Copyright (C) 2023 Laboratoire de Chemoinformatique, UMR 7140 CNRS (http://complex-matter.unistra.fr/equipes-de-recherche/laboratoire-de-chemoinformatique/home/)
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
unit UnitMoleculeRXN;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, contnrs, U_TYPE_GRAPHES, UnitMoleculeBase,
  unitAtomAndBondType, UnitAtomBase;

type

  { TMoleculeRXN }

  TMoleculeRXN = class(TMoleculeBase)
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
    function GetAtomMap(Id: AtomID): byte;
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
    procedure CreateCGR(const reacts, prods: TFPObjectList);//TObjectList);
    procedure CreateUncGraph(const mols: TFPObjectList);//TObjectList);
    procedure ModifyPAtCGR(PAtr, PAtp: PRAtom);
    procedure ModifyPBoCGR(PBor, PBop: PRBond);
    procedure InitiatePBo(PBo: PRBond; at, ah: Node);
    procedure CreateCGRtest(RxnFile: TStringList);
  end;

implementation

{ TMoleculeRXN }

constructor TMoleculeRXN.Create;
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

destructor TMoleculeRXN.Destroy;
begin
  inherited Destroy;
end;

function TMoleculeRXN.IsMarkedAt(Id: AtomID): boolean;
begin
  Result := IsMarkedAt(AtmSet[Id]);
end;

function TMoleculeRXN.IsMarkedAt(A: PRAtom): boolean;
begin
  Result := (A^.I[fOffSetMarkAt] = 1);
end;

function TMoleculeRXN.IsChargedAt(Id: AtomID): boolean;
begin
  Result := IsChargedAt(AtmSet[Id]);
end;

function TMoleculeRXN.IsChargedAt(A: PRAtom): boolean;
begin
  if (A^.I[fOffFormalCharge] <> 0) then
    Result := True
  else
    Result := False;
end;

function TMoleculeRXN.GetFormalCharge(Id: AtomID): byte;
begin
  Result := GetFormalCharge(AtmSet[Id]);
end;

function TMoleculeRXN.GetFormalCharge(A: PRAtom): byte;
begin
  Result := A^.I[fOffFormalCharge];
end;

function TMoleculeRXN.GetRadicalAt(Id: AtomID): byte;
begin
  Result:=GetRadicalAt(AtmSet[Id]);
end;

function TMoleculeRXN.GetRadicalAt(A: PRAtom): byte;
begin
  Result:=A^.I[fOffRadicalAt];
end;

function TMoleculeRXN.GetIsotopeAt(Id: AtomID): byte;
begin
  Result:=GetIsotopeAt(AtmSet[Id]);
end;

function TMoleculeRXN.GetIsotopeAt(A: PRAtom): byte;
begin
  Result:=A^.I[fOffIsotopeAt];
end;

function TMoleculeRXN.GetAtomMap(Id: AtomID): byte;                                                        /// atom mapping
begin
  Result := GetAtomMap(AtmSet[Id]);
end;

function TMoleculeRXN.GetAtomMap(A: PRAtom): byte;
begin
  Result := A^.I[fOffAtomMap];
end;

function TMoleculeRXN.GetX(Id: AtomID): double;
begin
  Result := GetX(AtmSet[Id]);
end;

function TMoleculeRXN.GetX(A: PRAtom): double;
begin
  Result := A^.P[0];
end;

function TMoleculeRXN.GetY(Id: AtomID): double;
begin
  Result := GetY(AtmSet[Id]);
end;

function TMoleculeRXN.GetY(A: PRAtom): double;
begin
  Result := A^.P[1];
end;

function TMoleculeRXN.GetZ(Id: AtomID): double;
begin
  Result := GetZ(AtmSet[Id]);
end;

function TMoleculeRXN.GetZ(A: PRAtom): double;
begin
  Result := A^.P[2];
end;

function TMoleculeRXN.IsMarkedBd(Id: BondID): boolean;
begin
  Result := IsMarkedBd(BndSet[Id]);
end;

function TMoleculeRXN.IsMarkedBd(B: PRBond): boolean;
begin
  Result := (B^.I[fOffSetMarkBd] = 1);
end;

function TMoleculeRXN.AtStereoParity(Id: AtomId): byte;
begin
  Result := AtStereoParity(AtmSet[Id]);
end;

function TMoleculeRXN.AtStereoParity(A: PRAtom): byte;
begin
  Result := A^.I[fOffSetStereoAt];
end;

function TMoleculeRXN.AtInCircuit(A: PRAtom): boolean;
begin
  Result:=(A^.I[fOffCircuitAt]=1);
end;

function TMoleculeRXN.BdInCricuit(B: PRBond): boolean;
begin
  Result:=(B^.I[fOffCircuitBd]=1);
end;

function TMoleculeRXN.AtDynAtm(Id: AtomId): byte;
begin
  Result := AtDynAtm(AtmSet[Id]);
end;

function TMoleculeRXN.AtDynAtm(A: PRAtom): byte;
begin
  Result := A^.I[fOffDynAtm];
end;

function TMoleculeRXN.IsDynAtom(Id: AtomId): boolean;
begin
  Result:=IsDynAtom(AtmSet[Id]);
end;

function TMoleculeRXN.IsDynAtom(A: PRAtom): boolean;
begin
  If A^.I[fOffDynAtm]=0 then Result:=False
  else Result:=True;
end;

function TMoleculeRXN.AtDynCharge(Id: AtomId): byte;
begin
  Result := AtDynCharge(AtmSet[Id]);
end;

function TMoleculeRXN.AtDynCharge(A: PRAtom): byte;
begin
  Result := A^.I[fOffDynCharge];
end;

function TMoleculeRXN.AtDynRadical(Id: AtomId): byte;
begin
  Result := AtDynRadical(AtmSet[Id]);
end;

function TMoleculeRXN.AtDynRadical(A: PRAtom): byte;
begin
  Result := A^.I[fOffDynRadical];
end;

function TMoleculeRXN.AtDynIsotope(Id: AtomId): byte;
begin
  Result := AtDynIsotope(AtmSet[Id]);
end;

function TMoleculeRXN.AtDynIsotope(A: PRAtom): byte;
begin
  Result := A^.I[fOffDynIsotope];         //LP
end;

function TMoleculeRXN.AtDynAtmS(Id: AtomId): string;
begin
  Result := AtDynAtmS(AtmSet[Id]);
end;

function TMoleculeRXN.AtDynAtmS(A: PRAtom): string;
var
  i: byte;
begin
  i := A^.I[fOffDynAtm];
  Result:=IntToDynAtomSymbol(i);
end;

function TMoleculeRXN.AtDynChargeS(Id: AtomId): string;
begin
  Result := AtDynChargeS(AtmSet[Id]);
end;

function TMoleculeRXN.AtDynChargeS(A: PRAtom): string;
var
  i: byte;
begin
  i := A^.I[fOffDynCharge];
  if (i=0) then Result:=''
  else if (i mod 2 <> 0) then Result:='c+'+IntToStr((i div 2)+1) //odd are positive
  else Result:='c'+IntToStr(-(i div 2));                         //even are negative
end;

function TMoleculeRXN.AtDynRadicalS(Id: AtomId): string;
begin
  Result := AtDynRadicalS(AtmSet[Id]);
end;

function TMoleculeRXN.AtDynRadicalS(A: PRAtom): string;
var
  i: byte;
begin
  i := A^.I[fOffDynRadical];
  if (i=0) then Result:=''
  else if (i mod 2 <> 0) then Result:='r+'+IntToStr((i div 2)+1) //odd are positive
  else Result:='r'+IntToStr(-(i div 2));                         //even are negative
end;

function TMoleculeRXN.AtDynIsotopeS(Id: AtomId): string;
begin
  Result := AtDynIsotopeS(AtmSet[Id]);
end;

function TMoleculeRXN.AtDynIsotopeS(A: PRAtom): string;
var
  i: integer;
begin
  i := A^.I[fOffDynIsotope];
  //if i>0 then writeln('BUG!'); //GM
  if (i=0) then Result:=''
  else if (i mod 2 <> 0) then Result:='i+'+IntToStr((i div 2)+1) //odd are positive
  else Result:='i'+IntToStr(-(i div 2));                         //even are negative
end;

function TMoleculeRXN.BdStereo(Id: BondId): byte;
begin
  Result := BdStereo(BndSet[Id]);
end;

function TMoleculeRXN.BdStereo(B: PRBond): byte;
begin
  Result := B^.I[fOffSetStereoBd];
end;

function TMoleculeRXN.BdTopol(Id: BondId): byte;
begin
  Result := BdTopol(BndSet[Id]);
end;

function TMoleculeRXN.BdTopol(B: PRBond): byte;
begin
  Result := B^.I[fOffSetTopolBd];
end;

procedure TMoleculeRXN.LoadSDF(sdfstr: TStringList);//; bMarkCircuit: Boolean=False);
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
  //
begin
  //User defined
  APrpSze := 3;       //Coordinates
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

    PAt^.S      := Trim(Copy(sdfstr[LineNo], 32, 3));        // atom symbol
    PAt^.Z      := AtomSymbolToInt(PAt^.S);
    //Only these infos are stored. Add whichever you would like to add.
    PAt^.I[fOffFormalCharge] := int_readpos(sdfstr[LineNo], 37, 3);
    // formal charge -> PAt^.I[0]
    PAt^.I[fOffSetMarkAt] := int_readpos(sdfstr[LineNo], 55, 3);
    //Unused field in V2000 used to mark atoms
    PAt^.I[fOffSetStereoAt] := int_readpos(sdfstr[LineNo], 40, 3); //Atom Stereo Parity
    PAt^.I[fOffAtomMap] := int_readpos(sdfstr[LineNo], 62, 3); //Atom-atom mapping number

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

function TMoleculeRXN.LoadSDFTStringList(sdfstr: TStringList;
  prpstr: string; msstr: TStringList): TFPStringHashTable;
var
  strhash: TFPStringHashTable;
begin
  LoadSDF(sdfstr);
  strhash := TFPStringHashTable.Create;
  LoadSDFField(sdfstr, prpstr, msstr, strhash);
  Result := strhash;
end;

procedure TMoleculeRXN.LoadSDFAGroup(sdfstr: TStringList);
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
    aInt: Integer;
    aword: string;
    pitem: PRSAL;
    pAtm: PRAtom;

  begin
    i1max:=int_readpos(lne,7,3);
    loffset:=10;
    for i1:=1 to i1max do
    begin
      atId:=int_readpos(lne,loffset+1,3);
      aInt:=int_readpos(lne,loffset+5,3);
      pAtm:=AtmSet[atId];
      case etype of
        atomcharge: pAtm^.I[fOffFormalCharge]:= UnConvertFormalCharge(aInt);
        atomisotope: pAtm^.I[fOffIsotopeAt]:= Byte(aInt);
        atomradical: pAtm^.I[fOffRadicalAt]:= Byte(aInt);
      end;
      //
      loffset:=loffset+8;
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

function TMoleculeRXN.LoadSDFField(sdfstr: TStringList; prpstr: string;
  msstr: TStringList): TFPStringHashTable;
var
  strhash: TFPStringHashTable;
begin
  strhash := TFPStringHashTable.Create;
  LoadSDFField(sdfstr, prpstr, msstr, strhash);
  Result := strhash;
end;

function TMoleculeRXN.LoadSDFField(sdfstr: TStringList; prpstr: TStringList;
  msstr: TStringList): TFPStringHashTable;
var
  strhash: TFPStringHashTable;

begin
  strhash := TFPStringHashTable.Create;
  LoadSDFField(sdfstr, prpstr, msstr, strhash);
  Result := strhash;
end;

procedure TMoleculeRXN.LoadSDFField(sdfstr: TStringList; prpstr: TStringList;
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

procedure TMoleculeRXN.LoadSDFFieldBdColor(sdfstr: TStringList; prpstr:string);
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

procedure TMoleculeRXN.AnnotateCycleAtm(bReset: boolean);
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

procedure TMoleculeRXN.AnnotateCycleBnd(bReset: boolean);
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

procedure TMoleculeRXN.LoadSDFField(sdfstr: TStringList; prpstr: string;
  msstr: TStringList; strhash: TFPStringHashTable);
var
  strlist: TStringList;
begin
  strhash.Clear;
  strlist := convert_dict_list(prpstr);
  LoadSDFField(sdfstr, strlist, msstr, strhash);
end;

procedure TMoleculeRXN.CreateCGR(const reacts, prods: TFPObjectList);//TObjectList);
var
  reactUncGraph, prodUncGraph: TMoleculeRXN;
  i, nbBonds: Integer;
  PAtr, PAt: PRAtom;
  PBor, PBo: PRBond;
  AOrdAtmReact, AOrdAtmProd: array of integer;
  s, t, ah, at: Node;
  atmp: CostMatrix;
  wtmp: array of PRBond;
  M: ArcNum;
begin
  p_NX := 0;
  p_M := 0;
  // Create the unconnected graph for both products and reactifs.
  //WriteLn(IntToStr(GetTickCount64) + ' list of react created');
  reactUncGraph := TMoleculeRXN.Create;
  reactUncGraph.CreateUncGraph(reacts);
  //WriteLn(IntToStr(GetTickCount64) + ' unc graph 1 created');
  prodUncGraph  := TMoleculeRXN.Create;
  prodUncGraph.CreateUncGraph(prods);
  //WriteLn(IntToStr(GetTickCount64) + ' unc graph 2 created');
  // Create two arrays containing indexes of atoms, ordered by atom mapping number.
  SetLength(AOrdAtmReact, reactUncGraph.nAtom +1);
  SetLength(AOrdAtmProd, prodUncGraph.nAtom +1);

  for i := 1 to reactUncGraph.nAtom do
  begin
    AOrdAtmReact[GetAtomMap(reactUncGraph.AtmSet[i])] := i; // Giving the atom mapping gives the index of the atom
  end;

  for i := 1 to prodUncGraph.nAtom do
  begin
    AOrdAtmProd[GetAtomMap(prodUncGraph.AtmSet[i])] := i ;
  end;

  // Loop threw the AtmSet of prodUncGraph and create the CGR AtmSet using the OrdAtm sets

  for i := 1 to prodUncGraph.nAtom do
  begin
    new(PAt);
    PAt^ := prodUncGraph.AtmSet[i]^;
    PAtr := reactUncGraph.AtmSet[AOrdAtmReact[GetAtomMap(PAt)]]; // Atom with the same atom mapping number in the reactUncGraph.Atmset
    ModifyPAtCGR(Patr, PAt);
    AtmSet[i] := PAt;
    p_NX := p_NX + 1;
  end;

  // Loop threw the BndSet of prodUncGraph and create the CGR BndSet using the OrdAtm sets

  i := 0;
  for s := 0 to MaxAtom + 1 do
    for t := 0 to MaxAtom + 1 do
      atmp[s, t] := 0;
  // We don't know the number of bonds, overshoot ?
  SetLength(wtmp, prodUncGraph.nBonds*2);

  nbBonds:=0;
  BndSet[0] := nil;

  for i := 1 to prodUncGraph.nBonds do      // We check every bond from the product, and compare with the bonds from the reactif
  begin
    if (atmp[prodUncGraph.BndSet[i]^.t, prodUncGraph.BndSet[i]^.h]=0) then
    begin
      new(PBo);
      PBo^ := prodUncGraph.BndSet[i]^;
      ah := AOrdAtmReact[GetAtomMap(prodUncGraph.AtmSet[PBo^.h])];
      at := AOrdAtmReact[GetAtomMap(prodUncGraph.AtmSet[PBo^.t])];
      PBor := reactUncGraph.FindBond(at, ah);        // If none, PBor := nil : created bond
      if (PBor=nil) then                            // Then we initiate it
      begin
        new(PBor);
        InitiatePBo(PBor, at, ah);
        ModifyPBoCGR(PBor, PBo);
        dispose(PBor);
        PBor := nil;
      end else
      begin
        ModifyPBoCGR(PBor, PBo);
      end;
      atmp[PBo^.t, PBo^.h] := i;
      atmp[PBo^.h, PBo^.t] := i;
      wtmp[i]              := PBo;
      nbBonds              := nbBonds + 1;
    end;
  end;

  // Check for bonds in react but not in product, use nbBonds instead of i in: atmp[PBo^.t, PBo^.h] := i;

  for i := 1 to reactUncGraph.nBonds do
  begin
    ah := AOrdAtmProd[GetAtomMap(reactUncGraph.AtmSet[reactUncGraph.BndSet[i]^.h])];
    at := AOrdAtmProd[GetAtomMap(reactUncGraph.AtmSet[reactUncGraph.BndSet[i]^.t])];

    // if the bond is not present in atmp, it's new, then we add it to the CGR
    if (atmp[at, ah] = 0) then
    begin
      // As the bond is not present in atmp, it's not present in prodUncGraph, then we initiate it.
      new(PBo);
      InitiatePBo(PBo, at, ah);
      PBor := reactUncGraph.BndSet[i];
      ModifyPBoCGR(PBor, PBo);
      atmp[PBo^.t, PBo^.h] := i + nbBonds;
      atmp[PBo^.h, PBo^.t] := i + nbBonds;
      wtmp[i + nbBonds]    := PBo;
    end;
  end;

  // Store the graph in a packed format
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

  FreeAndNil(reactUncGraph);
  FreeAndNil(prodUncGraph);
end;

procedure TMoleculeRXN.CreateUncGraph(Const mols: TFPObjectList);
// Create an unconnected graph from a tobjectlist Tobjectlist of molecules
var
  i, j, k, nbAtom, nbBonds: Integer;
  mol: TMoleculeRXN;
  PBo: PRBond;
  PAt: PRAtom;
  atmp: CostMatrix;
  s, t, ah, at: Node;
  wtmp: array of PRBond;
  M: ArcNum;
begin
  i := 0;
  nbAtom:= 0;
  nbBonds:=0;
  M:=0;
  p_NX := 0;
  p_M := 0;

  // Add the atoms to the UncGraph
  for k:=0 to mols.Count-1 do
  begin
    mol := mols.Items[k] as TMoleculeRXN;
    if (i = 0) then
    begin
      AtmSet[0] := nil;
    end;
    Inc(i);

    for j := 1 to mol.nAtom do
      begin
        new(PAt);
        PAt^ := mol.AtmSet[j]^;
        AtmSet[(j + nbAtom)] := PAt;
        p_NX := p_NX + 1;
      end;
    nbAtom := nAtom;
    p_M := p_M + mol.p_M;
  end;

  // Add the bonds to the UncGraph
  i := 0;
  nbAtom:=0;
  for s := 0 to p_NX + 1 do
    for t := 0 to p_NX + 1 do
      atmp[s, t] := 0;
  SetLength(wtmp, p_M + 1);

  for k:=0 to mols.Count-1 do
  begin
    mol := mols.Items[k] as TMoleculeRXN;
    if (i = 0) then
    begin
      BndSet[0] := nil;
    end;
    Inc(i);

    for j := 1 to mol.nBonds do
    begin
      ah := mol.BndSet[j]^.h + nbAtom;
      at := mol.BndSet[j]^.t + nbAtom;
      if (atmp[at, ah] = 0) then   // Bonds are represented in both ways, the second time we write on top of it
      begin
        new(PBo);
        PBo^ := mol.BndSet[j]^;
        PBo^.h := (mol.BndSet[j]^.h + nbAtom); // ah ?
        PBo^.t := (mol.BndSet[j]^.t + nbAtom); // at ?
        atmp[PBo^.t, PBo^.h] := j + nbBonds;
        atmp[PBo^.h, PBo^.t] := j + nbBonds;
        wtmp[j + nbBonds] := PBo;
      end;
    end;
    nbAtom  := nbAtom + mol.nAtom;
    nbBonds := nbBonds + mol.p_M;
  end;

  // Store the graph in a packed format
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
end;

procedure TMoleculeRXN.ModifyPAtCGR(PAtr, PAtp: PRAtom);
var
  ichr, ichp, irar, irap, iisr, iisp: Integer;

  function chargeUncoding(bcha: byte): Integer;
  begin
    case bcha of
      0:       Result := 0;
      1:       Result := 3;
      2:       Result := 2;
      3:       Result := 1;
      5:       Result := -1;
      6:       Result := -2;
      7:       Result := -3;
    end;
  end;

  function compareproperty(iprpr, iprpp: Integer): byte;
  var
  iprpcgr: Integer;

  begin
    if (iprpr <> iprpp) then
      begin
        iprpcgr := iprpp - iprpr;

        // To get the byte value
        if (iprpcgr>0) then Result := 2*(iprpcgr-1)+1 //odd for positive value
        else if (iprpcgr<0) then Result := -2*iprpcgr; //even for negative
      end;
  end;

begin
  // Check charge I[fOffDynCharge] | I[fOffDynAtm] ?
  if (IsChargedAt(PAtr)) then ichr := chargeUncoding(GetFormalCharge(PAtr)) else ichr := 0;
  if (IsChargedAt(PAtp)) then ichp := chargeUncoding(GetFormalCharge(PAtp)) else ichp := 0;

  if (ichr <> ichp) then
    begin
      Patp^.I[fOffDynCharge] := compareproperty(ichr, ichp);
    end;

  // Check radical I[fOffDynRadical]

  if (GetRadicalAt(PAtr) <> 0) then irar := trunc(GetRadicalAt(PAtr).ToDouble) else irar := 0;
  if (GetRadicalAt(PAtp) <> 0) then irap := trunc(GetRadicalAt(PAtp).ToDouble) else irap := 0;

  if (irar <> irap) then
    begin
      Patp^.I[fOffDynRadical] := compareproperty(irar, irap);
    end;

  // Check isotope I[fOffDynIsotope]

  if (GetIsotopeAt(PAtr) <> 0) then iisr := trunc(GetIsotopeAt(PAtr).ToDouble) else iisr := trunc(AtomSymbolToWeight(Patr^.S));
  if (GetIsotopeAt(PAtp) <> 0) then iisp := trunc(GetIsotopeAt(PAtp).ToDouble) else iisp := trunc(AtomSymbolToWeight(Patp^.S));

  if (iisr <> iisp) then
    begin
      Patp^.I[fOffDynIsotope] := compareproperty(iisr, iisp);
    end;

  // Check si un des trois champs est utilisés, si c'est le cas [fOffDynAtm] := 1

  if (AtDynCharge(PAtp) <> 0) or (AtDynIsotope(PAtp) <> 0) or (AtDynRadical(PAtp) <> 0) then Patp^.I[fOffDynAtm] := 1;
end;

procedure TMoleculeRXN.ModifyPBoCGR(PBor, PBop: PRBond);

  // If PBor or PBop are = nil then we initiate them
  function modifybondorder(TBr, TBp: TB): TB;
  begin
    // e>Z is 65 and 77 : 65 = u>z ?

    if (TBr=TBp) then
      begin
        Result := TBp;
      end
    else if (TBr = 98) then
      begin
        case TBp of
          1: Result := 23;  //bond 0>1
          2: Result := 24;  //bond 0>2
          3: Result := 25;  //bond 0>3
          4: Result := 26;  //bond 0>4
          14: Result := 27;  //bond 0>h
          13: Result := 28;  //bond 0>c
        end;
      end

    else if (TBp = 98) then
      begin
        case TBr of
          1: Result := 29;  //bond 1>0
          2: Result := 30;  //bond 2>0
          3: Result := 31;  //bond 3>0
          4: Result := 32;  //bond 4>0
          14: Result := 33;  //bond c>0
          13: Result := 34;  //bond h>0
        end;
      end

    else if (TBr = 1) then
      begin
        case TBp of
          2: Result := 35;  //bond 1>2
          3: Result := 36;  //bond 1>3
          4: Result := 37;  //bond 1>4
         14: Result := 38;  //bond 1>c
         13: Result := 39;  //bond 1>h
        end;
      end
    else if (TBr = 2) then
      begin
        case TBp of
          1: Result := 40;  //bond 2>1
          3: Result := 41;  //bond 2>3
          4: Result := 42;  //bond 2>4
         14: Result := 43;  //bond 2>c
        end;
      end
    else if (TBr = 3) then
      begin
        case TBp of
          1: Result := 44;  //bond 3>1
          2: Result := 45;  //bond 3>2
          4: Result := 46;  //bond 3>4
         14: Result := 47;  //bond 3>c
        end;
      end
    else if (TBr = 4) then
      begin
        case TBp of
          1: Result := 48;  //bond 4>1
          2: Result := 49;  //bond 4>2
          3: Result := 50;  //bond 4>3
        end;
      end
    else if (TBr = 18) then
      begin
        case TBp of
          1: Result := 51;  //bond c>1
          2: Result := 52;  //bond c>2
          3: Result := 53;  //bond c>3
          4: Result := 54;  //bond c>4
        end;
      end
    else if (TBr = 13) then
      begin
        case TBp of
          1: Result := 55;  //bond h>1
        end;
      end
    else if (TBr = 15) then
      begin
        case TBp of
         16: Result := 56;  //bond n>u
         14: Result := 57;  //bond n>c
         19: Result := 58;  //bond n>e
         20: Result := 59;  //bond n>z
         21: Result := 60;  //bond n>s
         22: Result := 61;  //bond n>a
        end;
      end
    else if (TBr = 16) then
      begin
        case TBp of
         15: Result := 62;  //bond u>n
         14: Result := 63;  //bond u>c
         19: Result := 64;  //bond u>e
         20: Result := 65;  //bond u>z
         21: Result := 66;  //bond u>s
         22: Result := 67;  //bond u>a
        end;
      end
    else if (TBr = 14) then
      begin
        case TBp of
         15: Result := 68;  //bond c>n
         16: Result := 69;  //bond c>u
         19: Result := 70;  //bond c>e
         20: Result := 71;  //bond c>z
         21: Result := 72;  //bond c>s
         22: Result := 73;  //bond c>a
        end;
      end
    else if (TBr = 19) then
      begin
        case TBp of
         15: Result := 74;  //bond e>n
         16: Result := 75;  //bond e>u
         18: Result := 76;  //bond e>c
         20: Result := 77;  //bond e>z
         21: Result := 78;  //bond e>s
         22: Result := 79;  //bond e>a
        end;
      end
    else if (TBr = 20) then
      begin
        case TBp of
         15: Result := 80;  //bond z>n
         16: Result := 81;  //bond z>u
         18: Result := 82;  //bond z>c
         19: Result := 83;  //bond z>e
         21: Result := 84;  //bond z>s
         22: Result := 85;  //bond z>a
        end;
      end
    else if (TBr = 21) then
      begin
        case TBp of
         15: Result := 86;  //bond s>n
         16: Result := 87;  //bond s>u
         18: Result := 88;  //bond s>c
         19: Result := 89;  //bond s>e
         20: Result := 90;  //bond s>z
         22: Result := 91;  //bond s>a
        end;
      end
    else if (TBr = 21) then
      begin
        case TBp of
         15: Result := 92;  //bond a>n
         16: Result := 93;  //bond a>u
         18: Result := 94;  //bond a>c
         19: Result := 95;  //bond a>e
         20: Result := 96;  //bond a>z
         21: Result := 97;  //bond a>s
        end;
      end
    else
    begin
      Result := BMax;
    end;
  end;

begin
  {if (PBor^.B = 98) then
    begin
      PBop^.B := modifybondorder(98, PBop^.B);
      PBop^.S      := BondSymbol[PBop^.B];
    end
  else if (PBop^.B = 98) then
    begin
      PBop^.B := modifybondorder(PBor^.B, 98);
      PBop^.S      := BondSymbol[PBop^.B];
    end
  else
    begin
      PBop^.B := modifybondorder(PBor^.B, PBop^.B);
      PBop^.S      := BondSymbol[PBop^.B];
    end;}
  // Modify PBo^.S, PBo^.B
  PBop^.B   := modifybondorder(PBor^.B, PBop^.B);
  PBop^.S   := BondSymbol[PBop^.B];
  // PBo^.S Modified


  // OffSetMarkBd

  // OffSetStereoBd

  // OffSetTopolBd


end;

procedure TMoleculeRXN.InitiatePBo(PBo: PRBond; at, ah: Node);
var
j: Integer;
begin
  BPrpSze := 0;
  BBytSze := fBdISze;
  SetLength(PBo^.P, BPrpSze);
  SetLength(PBo^.I, BBytSze);
  for j := Low(PBo^.P) to High(PBo^.P) do
    PBo^.P[j] := 0;
  for j := Low(PBo^.I) to High(PBo^.I) do
    PBo^.I[j] := 0;
  PBo^.t      := at;//Tail of bond
  PBo^.h      := ah;//Head of bond
  PBo^.S := BondSymbol[98];
  PBo^.B := BondSymbolToInt(PBo^.S);
  PBo^.I[fOffSetMarkBd]   := 0;
  PBo^.I[fOffSetStereoBd] := 0;
  PBo^.I[fOffSetTopolBd]  := 0;
end;


procedure TMoleculeRXN.CreateCGRtest(RxnFile: TStringList);
var
  reactUncGraph, prodUncGraph: TMoleculeRXN;
  LineNo, nbProd, nbReact, i, nbBonds, testdeb, testdeb2: Integer;
  reacts, prods: TFPObjectList;
  mol: TStringList;
  count_line : array of string;
  tmol: TMoleculeRXN;
  PAtr, PAt: PRAtom;
  PBor, PBo: PRBond;
  AOrdAtmReact, AOrdAtmProd: array of integer;
  s, t, ah, at: Node;
  atmp: CostMatrix;
  wtmp: array of PRBond;
  M: ArcNum;
begin
  //---Line 1 : $RXN in the first position on this line identifies the file as a reaction file---
  LineNo  := 0;
  //---Line 2 : Reaction name---
  Inc(LineNo);
  //---Line 3 of RxnFile---
  //  User's first and last initials (UserI), program name (ProgN), date/time:
  //  M/D/Y,H:m (DatTim), dimensional codes (DimCo), scaling factors (ScaI, ScaR),
  //  energy (Ene) if modeling program input, internal registry number (RegNo)
  //  if input through MDL format:
  Inc(LineNo);
  //---Line 4 of RxnFile---
  // A line for comments. If no comment is entered, a blank line must be present.
  Inc(LineNo);
  //---Line 5 of RxnFile : A line identifying the number of reactants and products, in that order.---
  //  We get the number of React and Prod
  Inc(LineNo);
  count_line := RxnFile[LineNo].Split([' '], TStringSplitOptions.ExcludeEmpty);
  nbReact := StrToInt(count_line[0]);
  nbProd := StrToInt(count_line[1]);

  reacts := TFPObjectList.Create(True);//TObjectList.Create(True);
  prods := TFPObjectList.Create(True);//TObjectList.Create(True);

  //WriteLn(IntToStr(GetTickCount64) + ' start to create reactobjlst');

 // Reacts
  for i := 1 to nbReact do
  begin
    // $MOL
    Inc(LineNo);
    mol := TStringList.Create;
    // Molfile
    while (RxnFile[LineNo] <> 'M  END') do
    begin
      Inc(LineNo);
      mol.Add(RxnFile[LineNo]);
    end;
    // Create TMol
    tmol :=  TMoleculeRXN.Create;
    tmol.LoadSDF(mol);
    reacts.Add(tmol);
    FreeAndNil(mol);
  end;
  //WriteLn(IntToStr(GetTickCount64) + ' start to create prodobjlst');
  // Products
  for i := 1 to nbProd do
  begin
    // $MOL
    Inc(LineNo);
    mol := TStringList.Create;
    // Molfile.
    while (RxnFile[LineNo] <> 'M  END') do
    begin
      Inc(LineNo);
      mol.Add(RxnFile[LineNo]);
    end;
    // Create TMol
    tmol :=  TMoleculeRXN.Create;
    tmol.LoadSDF(mol);
    prods.Add(tmol);
    FreeAndNil(mol);
  end;


  p_NX := 0;
  p_M := 0;

  // Create the unconnected graph for both products and reactifs.
  //WriteLn(IntToStr(GetTickCount64) + ' list of react created');
  reactUncGraph := TMoleculeRXN.Create;
  reactUncGraph.CreateUncGraph(reacts);
  //WriteLn(IntToStr(GetTickCount64) + ' unc graph 1 created');
  prodUncGraph  := TMoleculeRXN.Create;
  prodUncGraph.CreateUncGraph(prods);
  //WriteLn(IntToStr(GetTickCount64) + ' unc graph 2 created');
  // Create two arrays containing indexes of atoms, ordered by atom mapping number.
  SetLength(AOrdAtmReact, reactUncGraph.nAtom +1);
  SetLength(AOrdAtmProd, prodUncGraph.nAtom +1);

  for i := 1 to reactUncGraph.nAtom do
  begin
    AOrdAtmReact[GetAtomMap(reactUncGraph.AtmSet[i])] := i; // Giving the atom mapping gives the index of the atom
  end;

  for i := 1 to prodUncGraph.nAtom do
  begin
    AOrdAtmProd[GetAtomMap(prodUncGraph.AtmSet[i])] := i ;
  end;

  // Loop threw the AtmSet of prodUncGraph and create the CGR AtmSet using the OrdAtm sets

  for i := 1 to prodUncGraph.nAtom do
  begin
    new(PAt);
    PAt^ := prodUncGraph.AtmSet[i]^;
    PAtr := reactUncGraph.AtmSet[AOrdAtmReact[GetAtomMap(PAt)]]; // Atom with the same atom mapping number in the reactUncGraph.Atmset
    ModifyPAtCGR(Patr, PAt);   // We use information about both reactif and product state of the atom to create the "dynamic" atom
    AtmSet[i] := PAt;
    p_NX := p_NX + 1;
  end;

  // Loop threw the BndSet of prodUncGraph and create the CGR BndSet using the OrdAtm sets

  i := 0;
  for s := 0 to MaxAtom + 1 do
    for t := 0 to MaxAtom + 1 do
      atmp[s, t] := 0;
  // We don't know the number of bonds, overshoot ?
  SetLength(wtmp, prodUncGraph.nBonds*2);

  nbBonds:=0;
  BndSet[0] := nil;

  for i := 1 to prodUncGraph.nBonds do      // We check every bond from the product, and compare with the bonds from the reactif
  begin
    testdeb := prodUncGraph.BndSet[i]^.t;
    testdeb2 := prodUncGraph.BndSet[i]^.h;
    if (atmp[prodUncGraph.BndSet[i]^.t, prodUncGraph.BndSet[i]^.h]=0) then
    begin                                                                     //TODO Prb here ça saute la liaison entre 4 5
      new(PBo);
      PBo^ := prodUncGraph.BndSet[i]^;
      ah := AOrdAtmReact[GetAtomMap(prodUncGraph.AtmSet[PBo^.h])];
      at := AOrdAtmReact[GetAtomMap(prodUncGraph.AtmSet[PBo^.t])];
      PBor := reactUncGraph.FindBond(at, ah);        // If none, PBor := nil
      if (PBor=nil) then                            // Then we initiate it
      begin
        new(PBor);
        InitiatePBo(PBor, at, ah);
        ModifyPBoCGR(PBor, PBo);
        dispose(PBor);
        PBor := nil;
      end else
      begin
        ModifyPBoCGR(PBor, PBo);
      end;
      atmp[PBo^.t, PBo^.h] := i;
      atmp[PBo^.h, PBo^.t] := i;
      wtmp[i]              := PBo;
      nbBonds              := nbBonds + 1;
    end;
  end;

  // Check for bonds in react but not in product, use nbBonds instead of i in: atmp[PBo^.t, PBo^.h] := i;

  for i := 1 to reactUncGraph.nBonds do
  begin
    ah := AOrdAtmProd[GetAtomMap(reactUncGraph.AtmSet[reactUncGraph.BndSet[i]^.h])];
    at := AOrdAtmProd[GetAtomMap(reactUncGraph.AtmSet[reactUncGraph.BndSet[i]^.t])];

    // if the bond is not present in atmp, it's new, then we add it to the CGR
    if (atmp[at, ah] = 0) then
    begin
      // As the bond is not present in atmp, it's not present in prodUncGraph, then we initiate it.
      new(PBo);
      InitiatePBo(PBo, at, ah);
      PBor := reactUncGraph.BndSet[i];
      ModifyPBoCGR(PBor, PBo);
      atmp[PBo^.t, PBo^.h] := i + nbBonds;
      atmp[PBo^.h, PBo^.t] := i + nbBonds;
      wtmp[i + nbBonds]    := PBo;
    end;
  end;

  // Store the graph in a packed format
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

  FreeAndNil(reactUncGraph);
  FreeAndNil(prodUncGraph);

  for i:=reacts.Count-1 downto 0 do
  begin
    tmol := reacts.Extract(reacts[i]) as TMoleculeRXN;
    FreeAndNil(tmol);
  end;
  for i:=prods.Count-1 downto 0 do
  begin
    tmol := prods.Extract(prods[i]) as TMoleculeRXN;
    FreeAndNil(tmol);
  end;

  reacts.Destroy;
  prods.Destroy;
end;
end.


