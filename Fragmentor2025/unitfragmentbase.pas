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
unit UnitFragmentBase;
//Unit for base definition of fragmentation
{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, UnitAtmPrpWeight, U_GRAPHES, U_TYPE_GRAPHES,
  UnitMoleculeFrg, unitAtomAndBondType, UnitAtomBase, UnitBondBase, contnrs;

//Keep for historical reasons. Old nomencalture of Vitaly
     {AtomSymbol: array [1..111] of string[2]=('H','He','Li','Be','B','C','N','O',
     'F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
     'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y',
     'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
     'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',
     'Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
     'At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es',
     'Fm','Md','No','Lr','CD','CT','CB','CA','CO','CN','NI','NA');}
type

  TAoList = array of TList;
  PTAoList = ^TAoList;

  RAtmFrg = record
    idx, len: integer;
  end;
  PRAtmFrg = ^RAtmFrg;

  { TFrgBase }

  TFrgBase = class(TObject)
  private
    fN: integer; // Number of atoms
    fColorTerms: TStringList;
    fColorKeys: TStringList;
    fColorHash: TFPStringHashTable;
    fColorBondSDField: string; //Store the keyword of the SDF Field containing bond colors. Should contain only the field value?
    fStrictFrg: boolean;
    fNewFrg: boolean;
    fDynBnd: byte;
    fUseFormalCharge: boolean;
    fUseRadical: boolean;
    fUseIsotope: boolean;
    fIsPair: boolean;
    fbCycle: boolean;
    fMarkAtom: byte;
    function AtomSymbolToInt(str: string): integer;
  protected
    fFrgLst: TStringList;
    fColorPerAtom: TObjectList; //For each atom a list of possible colors
    fFrgPerAtom: TObjectList;
    fGetFrgPerAtom: boolean;
    Stereo: array of integer;           // Stereo indicator array
  public
    property N: integer read fN;
    property FrgLst: TStringList read fFrgLst write fFrgLst;
    property StrictFrg: boolean read fStrictFrg write fStrictFrg;
    property NewFrg: boolean read fNewFrg write fNewFrg;
    property IsPair: boolean read fIsPair write fIsPair;
    property DynBnd: byte read fDynBnd write fDynBnd;
    property MarkAtom: byte read fMarkAtom write fMarkAtom;
    property UseFormalCharge: boolean read fUseFormalCharge write fUseFormalCharge;
    property UseRadical: boolean read fUseRadical write fUseRadical;
    property UseIsotope: boolean read fUseIsotope write fUseIsotope;
    property ColorTerms: TStringList read fColorTerms write fColorTerms;
    property ColorKeys: TStringList read fColorKeys write fColorKeys;
    property ColorHash: TFPStringHashTable read fColorHash write fColorHash;
    property ColorBondSDField: string read fColorBondSDField write fColorBondSDField;
    property FrgPerAtom: TObjectList read fFrgPerAtom write fFrgPerAtom;
    property GetFrgPerAtom: boolean read fGetFrgPerAtom write fGetFrgPerAtom;
    property bCycle: boolean read fbCycle write fbCycle;
    constructor Create;
    constructor Create(s: TStringList);
    destructor Destroy; override;
    procedure Clear;
    function AddFragment(frg: string): integer;
    function AddFragment(frg: string; iw: integer): integer;
    procedure SDFToFrgLst(MolSDF: TStringList); virtual;
    procedure MolToFrgLst(Mol: TMoleculeFrg); virtual; abstract;
    procedure SetFrgLst(s: TStringList);
    function GetFrgLst: TStringList;
    procedure InitAtomString(sdline: string); virtual; abstract;
    procedure InitMF(mfstr: string; mf: TStringList); virtual;
    procedure InitMF(ColorKey: string; out iw: integer; var RepA: TAtomBase;
      Mol: TMoleculeFrg); virtual;
    procedure InitColorPerAtom(Mol: TMoleculeFrg);
    function DynBondInFrg(Pth: NodeMatrix; LneIndx, LneLen: Node;
      Mol: TMoleculeFrg): boolean;
    function AllDynBondInFrg(Pth: NodeMatrix; LneIndx, LneLen: Node;
      Mol: TMoleculeFrg): boolean;
    function MrkAtmInFrg(Pth: NodeMatrix; LneIndx, LneLen: Node;
      Mol: TMoleculeFrg): boolean;
    //Pth: 2D array of sequences; LneIndx: line of interest in Pth;
    //LneLen: sequence length;Mol: the molecule.
  end;

  function TSLCompareStr( List: TStringList; Index1: Integer; Index2: Integer): Integer;

implementation

function TSLCompareStr(List: TStringList; Index1: Integer; Index2: Integer
  ): Integer;
//Use this to avoid unstable results of AnsiCompareStr using into StringList.Sort
begin
  Result:=CompareStr(List.Strings[Index1],List.Strings[Index2]);
end;



{ FrgBase }

function TFrgBase.AtomSymbolToInt(str: string): integer;
  //Définition nécessaire pour dissocier les fragments de l'interprétation des SDF
begin
  Result := Low(AtomSymbol) - 1;
  repeat
    Inc(Result);
  until ((Result >= High(AtomSymbol)) or (AtomSymbol[Result] = str));
  if (AtomSymbol[Result] <> str) then
    Result := 0;
end;

constructor TFrgBase.Create;
begin
  inherited Create;
  fFrgLst := nil;
  fFrgPerAtom := nil;
  fColorPerAtom := TObjectList.Create;
  fColorPerAtom.OwnsObjects := True;
  fColorTerms := TStringList.Create;
  fColorTerms.Delimiter := ';';
  fColorKeys := TStringList.Create;
  fColorHash := TFPStringHashTable.Create;
  fColorBondSDField:='Default';
  fNewFrg := False;
  fIsPair := False;
  fDynBnd := 0;
  fUseFormalCharge := False;
  fUseRadical:=False;
  fUseIsotope:=False;
  fMarkAtom := 0;
  fbCycle:=False;
end;

constructor TFrgBase.Create(s: TStringList);
begin
  Create;
  fFrgLst := s;
  fFrgLst.CaseSensitive:=True;
end;

destructor TFrgBase.Destroy;
begin
  Clear;
  FreeAndNil(fColorTerms);
  FreeAndNil(fColorKeys);
  FreeAndNil(fColorHash);
  FreeAndNil(fColorPerAtom);
  inherited Destroy;
end;

procedure TFrgBase.Clear;
var
  i: integer;
begin
  fNewFrg := False;
  fFrgLst := nil;
  fFrgPerAtom := nil;
  fColorPerAtom.Clear;
  fColorTerms.Clear;
  fColorKeys.Clear;
  fColorHash.Clear;
  fColorBondSDField:='Default';
end;

function TFrgBase.AddFragment(frg: string): integer;
var
  idx: integer;
  cnt: Pinteger;
begin
  Result := -1;
  idx := fFrgLst.IndexOf(frg);
  if ((idx < 0) and not (fStrictFrg)) then
  begin
    new(cnt);
    cnt^ := 1;
    Result := fFrgLst.AddObject(frg, TObject(cnt));
  end
  else if (idx >= 0) then
  begin
    Inc(PInteger(fFrgLst.Objects[idx])^);
    Result := idx;
  end;
  if (idx < 0) then
    fNewFrg := True;
end;

function TFrgBase.AddFragment(frg: string; iw: integer): integer;
var
  idx: integer;
  cnt: Pinteger;
begin
  Result := -1;
  idx := fFrgLst.IndexOf(frg);
  //writeln('---> '+frg+' :: '+IntToStr(idx)+' :: '+IntToStr(iw));
  if ((idx < 0) and not (fStrictFrg)) then
  begin
    new(cnt);
    cnt^ := iw;
    Result := fFrgLst.AddObject(frg, TObject(cnt));
  end
  else if (idx >= 0) then
  begin
    PInteger(fFrgLst.Objects[idx])^ := PInteger(fFrgLst.Objects[idx])^ + iw;
    Result := idx;
  end;
  if (idx < 0) then
    fNewFrg := True;
end;

procedure TFrgBase.SDFToFrgLst(MolSDF: TStringList);
var
  Mol: TMoleculeFrg;
begin
  Mol := TMoleculeFrg.Create;
  Mol.LoadSDF(MolSDF);
  MolToFrgLst(Mol);
  FreeAndNil(Mol);
end;

procedure TFrgBase.SetFrgLst(s: TStringList);
var
  i: integer;
begin
  fFrgLst := s;
  fFrgLst.CaseSensitive:=True;
end;

function TFrgBase.GetFrgLst: TStringList;
begin
  Result := fFrgLst;
end;

procedure TFrgBase.InitMF(mfstr: string; mf: TStringList);
begin
  mf.Clear;
  mf.StrictDelimiter := True;
  mf.Delimiter := '/';
  mf.DelimitedText := mfstr;
end;

procedure TFrgBase.InitMF(ColorKey: string; out iw: integer;
  var RepA: TAtomBase; Mol: TMoleculeFrg);
var
  ii, u: integer;
  mf: TStringList;
begin
  for u := 1 to Mol.nAtom do
  begin
    mf := fColorPerAtom[u] as TStringList;
    mf.Clear;
  end;
  iw := 1;
  if (ColorKey = 'Default') then
    for u := 1 to Mol.nAtom do
      RepA.AtomString[u] := Mol.S_[u]
  else
    iw := RepA.InitAtomStringWeight(ColorHash.Items[ColorKey]);
  for u := 1 to Mol.nAtom do
  begin
    mf := fColorPerAtom[u] as TStringList;
    mf.Clear;
    mf.StrictDelimiter := True;
    mf.Delimiter := '/';
    mf.DelimitedText := RepA.AtomString[u];
    //writeln(); writeln('mf='+RepA.AtomString[u]);
    //for ii:=0 to mf.Count-1 do write(mf[ii]+' ');
    //writeln;
    for ii := 0 to mf.Count - 1 do//Add dyn atom annotations
      mf[ii]:=mf[ii]+Mol.AtDynChargeS(u)+Mol.AtDynRadicalS(u)+Mol.AtDynIsotopeS(u);
    if (UseFormalCharge and Mol.IsChargedAt(u)) then
      for ii := 0 to mf.Count - 1 do
        mf[ii] := mf[ii] + '&FC' + RepA.ConvertFormalCharge(
          Mol.GetFormalCharge(u)) + '&';
    if UseRadical then begin
      for ii := 0 to mf.Count - 1 do
        if Mol.GetRadicalAt(u)<>0 then mf[ii]:=mf[ii]+'&RA'+IntToStr(Mol.GetRadicalAt(u))+'&';
    end;
    if UseIsotope then begin
      for ii := 0 to mf.Count - 1 do
        if Mol.GetIsotopeAt(u)<>0 then mf[ii]:=mf[ii]+'&IS'+IntToStr(Mol.GetIsotopeAt(u))+'&';
    end;
    if ((MarkAtom = 3) and Mol.IsMarkedAt(u)) then
      for ii := 0 to mf.Count - 1 do
        mf[ii] := mf[ii] + '&MA&';
  end;
end;

procedure TFrgBase.InitColorPerAtom(Mol: TMoleculeFrg);
var
  s: integer;
begin
  fColorPerAtom.Add(TStringList.Create);
  for s := 1 to Mol.nAtom do
    fColorPerAtom.Add(TStringList.Create);
end;

function TFrgBase.DynBondInFrg(Pth: NodeMatrix; LneIndx, LneLen: Node;
  Mol: TMoleculeFrg): boolean;
var
  i, j: integer;
  at1, at2: Node;
  prbd: PRBond;
begin
  i := 1;
  Result := False;
  while ((i <= LneLen - 1) and (Result = False)) do
  begin
    at1 := Pth[LneIndx, i];
    at2 := Pth[LneIndx, i + 1];
    prbd := Mol.FindBond(at1, at2);
    //if ((prbd^.B>9) and (prbd^.B<21)) then Result:=True;
    //if ((prbd^.B > 13) and (prbd^.B < 34)) then Result := True;
    if ((prbd^.B > 22) and (prbd^.B < 98)) then Result := True;
    if Mol.IsDynAtom(at1) or Mol.IsDynAtom(at2) then Result:=True;
    Inc(i);
  end;
end;

function TFrgBase.AllDynBondInFrg(Pth: NodeMatrix; LneIndx, LneLen: Node;
  Mol: TMoleculeFrg): boolean;
var
  i: integer;
  at1, at2: Node;
  prbd: PRBond;
begin
  i := 1;
  Result := True;
  while ((i <= LneLen - 1) and (Result = True)) do
  begin
    at1 := Pth[LneIndx, i];
    at2 := Pth[LneIndx, i + 1];
    prbd := Mol.FindBond(at1, at2);
    //if ((prbd^.B<10) or (prbd^.B>20)) then Result:=False;
    //if ((prbd^.B < 14) or (prbd^.B > 33)) then Result := False;
    if ((prbd^.B < 23) or (prbd^.B > 97)) then Result := False;
    Inc(i);
  end;
end;

function TFrgBase.MrkAtmInFrg(Pth: NodeMatrix; LneIndx, LneLen: Node;
  Mol: TMoleculeFrg): boolean;
var
  i: integer;
  at: Node;
begin
  i := 1;
  Result := False;
  while ((i <= LneLen) and (Result = False)) do
  begin
    at := Pth[LneIndx, i];
    if (Mol.IsMarkedAt(at)) then
      Result := True;

    Inc(i);
  end;
end;

end.

