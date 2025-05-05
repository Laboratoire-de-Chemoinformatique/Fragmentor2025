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
unit UnitACFXPathBnd;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, contnrs, UnitFragmentBase, UnitACFXBase, UnitMoleculeFrg, UnitMoleculeBase,
  unitAtomAndBondType, U_TYPE_GRAPHES, UnitAtomBase; //UnitBondBase,

type

  { TACFXPathBnd }

  TACFXPathBnd = class(TACFXBase)
  private
    //RepB: TBondBase;
  protected
    procedure PathToString(Mol: TMoleculeFrg; TOL: TObjectList); override;
    procedure PathToString(Mol: TMoleculeFrg; LOPath: TIACHolder;
      TSL: TStringList; ACID: integer); override;
  public
    constructor Create;
    constructor Create(s: TStringList);
    destructor Destroy; override;
  end;

implementation

{ TACFXPathBnd }

procedure TACFXPathBnd.PathToString(Mol: TMoleculeFrg; TOL: TObjectList);
var
  u, v, vmax: Node;
  head, tail: AtomID;
  PB:    PRBond;
  FORWD: string;
begin
  for u := 1 to fNP do
  begin
    FORWD := '';
    vmax  := fL[u];
    for v := 1 to fL[u] - 1 do
    begin
      tail  := fPath[u, v];
      head  := fPath[u, v + 1];
      PB    := Mol.FindBond(tail, head);
      //RepB.StereoBond := Mol.BdStereo(PB);
      FORWD := FORWD + PB^.S;
    end;
    (TOL.Items[fV[fPath[u, vmax]] - LenMin + 1] as TStringList).Add(
      '(' + FORWD + ')');
  end;
end;

procedure TACFXPathBnd.PathToString(Mol: TMoleculeFrg; LOPath: TIACHolder;
  TSL: TStringList; ACID: integer);
var
   s, tail, head: Node;
   i, k, ii, kk, kkmax: integer;
   col: integer;
   iw, u: integer;
   piw, pu: PInteger;
   FORWD,sFrg: string;
   TmpSL: TStringList;
   LOPS: TObjectList;
   mf: TStringList;
   aIPthHold: TIPathHolder;
   aACHold: TSACHolder;
   PB: PRBond;
begin
  //LOPath contains all path starting at a given atom
  //LOPS contains all alternatives of given fragment
  TmpSL:=TStringList.Create;
  TmpSL.Duplicates:=dupAccept;
  TmpSL.Delimiter:=',';
  LOPS:=TObjectList.create;
  LOPS.OwnsObjects:=True;
  for col := 0 to ColorKeys.Count - 1 do
  begin
    InitMF(ColorKeys[col],iw,TAtomBase(RepAC),Mol);
    mf:=fColorPerAtom[ACID] as TStringList;
    LOPS.Add(TSACHolder.Create);
    aACHold:=LOPS.Last as TSACHolder;
    for ii:=0 to LOPS.Count-1 do
    begin
      aACHold:=LOPS[ii] as TSACHolder;
      aACHold.AC:=mf[0];
    end;
    for k:=1 to mf.Count-1 do
    begin
      //Add a copy of the current state of the fragment
      kkmax:=LOPS.Count-1;
      for kk:=0 to kkmax do
      begin
        LOPS.Add(TSACHolder.Create);
        aACHold:=LOPS.Last as TSACHolder;
        aACHold.AC:=mf[k];
      end;
    end;
    //Parse the LOPS. Each object is a list containing the representation of a
    //path relevent for the fragment. Each object must be translated to a string,
    //and combined canonicaly for a string representation of the fragment.
    TmpSL.Clear;
    for ii:=0 to LOPS.Count-1 do
    begin
      aACHold:=LOPS[ii] as TSACHolder;
      for i:=0 to LOPath.LIPathHolder.Count-1 do
      begin
        aIPthHold:=LOPath.LIPathHolder[i] as TIPathHolder;
        FORWD:='';
        for s:=1 to High(aIPthHold.IPath)-1 do
        begin
          tail:=aIPthHold.IPath[s];
          head:=aIPthHold.IPath[s+1];
          PB:=Mol.FindBond(tail,head);
          FORWD:=FORWD+PB^.S;
        end;
        FORWD:='('+FORWD+')';
        TmpSL.Add(FORWD);
      end;
      TmpSL.CustomSort(@TSLCompareStr);//TmpSL.Sort;
      if TmpSL.Count>0 then
      begin
        new(piw); piw^:=iw;
        sFrg:=TmpSL.DelimitedText+',x'+aACHold.AC;
        TSL.AddObject(sFrg, TObject(piw));
      end;
    end;
    LOPS.Clear;
  end;
  FreeAndNil(LOPS);
  FreeAndNil(TmpSL);
end;

constructor TACFXPathBnd.Create;
begin
  inherited Create;
  //RepB := TBondBase.Create;
end;

constructor TACFXPathBnd.Create(s: TStringList);
begin
  inherited Create;
  //RepB    := TBondBase.Create;
  fFrgLst := s;
end;

destructor TACFXPathBnd.Destroy;
begin
  //FreeAndNil(RepB);
  inherited Destroy;
end;


end.

