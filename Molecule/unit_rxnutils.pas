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
unit Unit_RXNUtils;

{$mode objfpc}{$H+}

interface

uses
  Classes,
  SysUtils,
  contnrs,
  UnitMoleculeFrg,
  UnitMoleculeRXN,
  UnitMoleculeBase,
  unitAtomAndBondType,
  U_TYPE_GRAPHES;

function RXN_To_FRG(RxnFile: TStringList): TMoleculeFrg;
//function RXN_To_CGR(RxnFile: TStringList): TMoleculeRXN;
function CGR_To_Frg(CGR: TMoleculeRXN): TMoleculeFrg;
procedure RXN_To_Uncgraphs(RxnFile: TStringList; reactUnc, prodUnc: TStringList);

implementation

function RXN_To_FRG(RxnFile: TStringList): TMoleculeFrg;
var
  CGR: TMoleculeRXN;
  num: Integer;
begin
  CGR := TMoleculeRXN.Create;
  num := 0;
  //test := RxnFile[num];
  if (RxnFile[0] <> '$RXN') then
    begin
      CGR.LoadSDF(RxnFile);
    end
  else
    begin
      CGR.CreateCGRtest(RxnFile);
    end;

  //WriteLn(IntToStr(GetTickCount64) + ' CGR creation end of procedure');
  Result := CGR_To_Frg(CGR);
  FreeAndNil(CGR);
end;

function CGR_To_Frg(CGR: TMoleculeRXN): TMoleculeFrg;
var
  i:      Integer;
  PAt:    PRAtom;
  PBo:    PRBond;
  M:      ArcNum;
  s, t:   Node;
  atmp:   CostMatrix;
  wtmp:   array of PRBond;
begin
  Result := TMoleculeFrg.Create;

  M:=0;
  Result.p_NX := 0;
  Result.p_M := 0;

  // Copy of the atom set from the CGR to the Result
  Result.AtmSet[0] := nil;

  for i := 1 to CGR.nAtom do
      begin
        new(PAt);
        PAt^ := CGR.AtmSet[i]^;
        Result.AtmSet[i] := PAt;
        Result.p_NX := Result.p_NX + 1;
      end;
  Result.p_M := CGR.p_M;

 // Copy of the bond set from the CGR to the Result
  i := 0;
    for s := 0 to MaxAtom + 1 do
      for t := 0 to MaxAtom + 1 do
        atmp[s, t] := 0;
    SetLength(wtmp, Result.p_M + 1);

  Result.BndSet[0] := nil;

  for i := 1 to CGR.nBonds do
      begin
        if (atmp[CGR.BndSet[i]^.t, CGR.BndSet[i]^.h] = 0) then
        begin
          new(PBo);
          PBo^ := CGR.BndSet[i]^;
          atmp[PBo^.t, PBo^.h] := i;
          atmp[PBo^.h, PBo^.t] := i;
          wtmp[i] := PBo;
        end;
      end;
  M := 0;
  for s := 1 to Result.p_NX do
  begin
    Result.p_HEAD[s] := M + 1;
    for t := 1 to Result.p_NX do
      if (atmp[s, t] <> 0) then
      begin
        Inc(M);
        Result.p_SUCC[M] := t;
        Result.BndSet[M] := wtmp[atmp[s, t]];
      end;
  end;
  Result.p_HEAD[Result.p_NX + Result.p_NY + 1] := M + 1;
  Result.p_M := M;
end;

{function RXN_To_CGR(RxnFile: TStringList): TMoleculeRXN;
var
  LineNo, nbProd, nbReact, i: integer;
  ReacName: String;
  count_line : array of string;
  reacts, prods: TFPObjectList;//TObjectList; // of TMoleculeRXN;
  mol: TStringList;
  tmol: TMoleculeRXN;
  CGR: TMoleculeRXN;

begin
//Read each line and create the base molecules
  LineNo  := 0;
  //---Line 1 : $RXN in the first position on this line identifies the file as a reaction file---
  Inc(LineNo);
  //---Line 2 : Reaction name---
  ReacName := RxnFile[LineNo];
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
  //WriteLn(IntToStr(GetTickCount64) + ' start to create CGR');

  // Create cgr
  CGR := TMoleculeRXN.Create;
  CGR.CreateCGR(reacts, prods);
  //WriteLn(IntToStr(GetTickCount64) + ' CGR Created');

  // Cleaning

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

  Result := CGR;
end;
}

procedure RXN_To_Uncgraphs(RxnFile: TStringList; reactUnc, prodUnc: TStringList);
var
  LineNo, nbProd, nbReact, i, j: integer;
  ReacName: String;
  count_line : array of string;
  mol, reactuncTemp, produncTemp: TStringList;
  reacts, prods: TFPObjectList;
  tmol, reactUncGraph, prodUncGraph: TMoleculeRXN;
  reactUncGraphFrg, prodUncGraphFrg: TMoleculeFrg;

begin
//Read each line and create the base molecules
  LineNo  := 0;
  //---Line 1 : $RXN in the first position on this line identifies the file as a reaction file---
  Inc(LineNo);
  //---Line 2 : Reaction name---
  ReacName := RxnFile[LineNo];
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

  reactUncGraph := TMoleculeRXN.Create;
  reactUncGraph.CreateUncGraph(reacts);

  reactUncGraphFrg := TMoleculeFrg.Create;
  reactUncGraphFrg := CGR_To_Frg(reactUncGraph);
  reactuncTemp := reactUncGraphFrg.ReturnSdfTstringList();

  for j:=0 to reactuncTemp.Count-1 do
    begin
      reactUnc.Add(reactuncTemp.Strings[j]);
    end;
  FreeAndNil(reactuncTemp);

  prodUncGraph := TMoleculeRXN.Create;
  prodUncGraph.CreateUncGraph(prods);

  prodUncGraphFrg := TMoleculeFrg.Create;
  prodUncGraphFrg := CGR_To_Frg(prodUncGraph);
  produncTemp := prodUncGraphFrg.ReturnSdfTstringList();

  for j:=0 to produncTemp.Count-1 do
    begin
      prodUnc.Add(produncTemp.Strings[j]);
    end;
  FreeAndNil(produncTemp);

  // Cleaning
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
  FreeAndNil(reactUncGraph);
  FreeAndNil(reactUncGraphFrg);
  FreeAndNil(prodUncGraph);
  FreeAndNil(prodUncGraphFrg);
end;
end.


