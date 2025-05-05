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
unit UnitACPathBnd;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, strutils, UnitShortestPath, UnitFragmentBase, UnitACBase, UnitMoleculeFrg,
  UnitMoleculeBase, unitAtomAndBondType, UnitAtomBase, U_TYPE_GRAPHES, contnrs; //UnitBondBase,

type

  { TACPathBnd }

  TACPathBnd = class(TACBase)
  private
    //fRepB: TBondBase;
  protected
    procedure PathToString(Mol: TMoleculeFrg; TSL: TStringList); override;
    procedure PathToString(Mol: TMoleculeFrg; LOPath: TIACHolder;
      TSL: TStringList; ACID: integer); override;
  public
    constructor Create;
    constructor Create(s: TStringList);
    destructor Destroy; override;
  end;

implementation

{ TACPathBnd }

constructor TACPathBnd.Create;
begin
  inherited Create;
  //fRepB := TBondBase.Create;
end;

constructor TACPathBnd.Create(s: TStringList);
begin
  inherited Create;
  //fRepB   := TBondBase.Create;
  fFrgLst := s;
end;

destructor TACPathBnd.Destroy;
begin
  //FreeAndNil(fRepB);
  inherited Destroy;
end;

procedure TACPathBnd.PathToString(Mol: TMoleculeFrg; TSL: TStringList);
var
  u, v:  Node;
  head, tail: AtomID;
  PB:    PRBond;
  FORWD: string;
begin
  for u := 1 to fNP do
  begin
    if (fL[u] >= LenMin) then
    begin

      FORWD := '';
      for v := 1 to fL[u] - 1 do
      begin
        tail  := fPath[u, v];
        head  := fPath[u, v + 1];
        PB    := Mol.FindBond(tail, head);
        //fRepB.StereoBond := Mol.BdStereo(PB);
        FORWD := FORWD + PB^.S;
      end;

      TSL.Add('(' + FORWD + ')');
    end;
  end;
end;

procedure TACPathBnd.PathToString(Mol: TMoleculeFrg; LOPath: TIACHolder;
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
    InitMF(ColorKeys[col],iw,TAtomBase(fRepAC),Mol);
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
      new(piw); piw^:=iw;
      TmpSL.CustomSort(@TSLCompareStr);//TmpSL.Sort;
      if TmpSL.Count>0 then
      begin
        sFrg:=TmpSL.DelimitedText+',x'+aACHold.AC;
        TSL.AddObject(sFrg, TObject(piw));
      end;
    end;
    LOPS.Clear;
  end;
  FreeAndNil(LOPS);
  FreeAndNil(TmpSL);
end;

end.

