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
unit UnitAllPathBond;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, strutils, U_TYPE_GRAPHES, UnitMoleculeFrg, UnitMoleculeBase,
  unitAtomAndBondType, unitfragmentbase, unitsequences,
  UnitAllPath;        //UnitBondBase,

type

  { TAllPathBond }

  TAllPathBond = class(TAllPath)
  private
    //Rep: TBondBase;
  protected
    function PathToString(Mol: TMoleculeFrg): TStringList; override;
    procedure PathToString(Mol: TMoleculeFrg; TLSFrg: TStringList); override;
    procedure PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList; PathLst: TList); override;
  public
    constructor Create;
    constructor Create(s: TStringList);
    destructor Destroy; override;
  end;

implementation

{ TAllPathBond }
procedure TAllPathBond.PathToString(Mol: TMoleculeFrg; TLSFrg: TStringList);
var
  u, v:    Node;
  head, tail: AtomID;
  PB:      PRBond;
  FORWD, BCKWD: string;
  piw:     PInteger;
  piwused: boolean;
begin
  TLSFrg.Clear;
  new(piw);
  piw^    := 1;
  piwused := False;
  for u := 1 to fNPath do
  begin
    if fPath[u,0]=1 then
    {if (fLen[u] >= LenMin) then
    begin
      if (((DynBnd = 1) and (DynBondInFrg(fPath, u, fLen[u], Mol))) or
        ((DynBnd = 2) and (AllDynBondInFrg(fPath, u, fLen[u], Mol))) or
        (DynBnd = 0)) then}
      begin
        piwused := True;
        FORWD   := '';
        BCKWD   := '';
        for v := 1 to fLen[u] - 1 do
        begin
          tail  := fPath[u, v];
          head  := fPath[u, v + 1];
          PB    := Mol.FindBond(tail, head);
          //Rep.StereoBond := Mol.BdStereo(PB);
          FORWD := FORWD + PB^.S;
        end;

        BCKWD := ReverseString(FORWD);
        //only usable because a bond is represented by one char!
        if FORWD <= BCKWD then
        begin
          ;
          TLSFrg.AddObject(FORWD, TObject(piw));
        end
        else
        begin
          TLSFrg.AddObject(BCKWD, TObject(piw));
        end;
      end;
    //end;
  end;
  if (piwused = False) then
    dispose(piw);
end;

procedure TAllPathBond.PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList;
  PathLst: TList);
var
  piw: PInteger;
  j, k, iw, ii, jj, iimax, jjmax: integer;
  i, u, v, rv, tail, head: Node;
  FORWD, BCKWD: string;
  THERE, BACK, Code: string;
  PB: PRBond;
  pu: PInteger;
begin
  for j := 0 to ColorKeys.Count - 1 do
  begin
    //initialization of the representation and extraction of weight associated to micro-species
    for u := 1 to fNPath do
      if fPath[u,0]=1 then
      begin
        ffw.Clear;
        fbw.Clear;
        FORWD := '';
        BCKWD := '';
        ffw.Add(FORWD);
        fbw.Add(BCKWD);
        for v := 1 to fLen[u]-1 do
        begin
          tail := fPath[u, v];
          head := fPath[u, v + 1];
          PB := Mol.FindBond(tail, head);
          iimax := ffw.Count - 1;
          for ii := 0 to iimax do
            ffw.Add(ffw[ii] +PB^.S);
          for ii := iimax downto 0 do
            ffw.Delete(ii);
          //Generation of backward sequence - needed to have canonical sequences
          tail := fPath[u, fLen[u] - v + 1];
          head := fPath[u, fLen[u] - v];
          PB := Mol.FindBond(tail, head);
          iimax := fbw.Count - 1;
          for ii := 0 to iimax do
            fbw.Add(fbw[ii] + PB^.S);
          for ii := iimax downto 0 do
            fbw.Delete(ii);
        end;
        //Count fragment
        for ii := 0 to ffw.Count - 1 do
        begin
          new(piw); piw^:=1;//iw;
          if (ffw[ii] <= fbw[ii]) then
            TSLFrg.AddObject(ffw[ii], TObject(piw))
          else
            TSLFrg.AddObject(fbw[ii], TObject(piw));
          new(pu); pu^:=u;
          PathLst.Add(TObject(pu));
        end;
      end;
  end;
end;

function TAllPathBond.PathToString(Mol: TMoleculeFrg): TStringList;
var
  PathLst: TList;
  i: integer;
  pu: PInteger;
begin
  PathLst:=TList.Create;
  Result := TStringList.Create;
  PathToString(Mol, Result,PathLst);
  for i:=0 to PathLst.Count-1 do
  begin
    pu:=PInteger(PathLst[i]);
    if pu<>nil then
       dispose(pu);
    pu:=nil;
  end;
  FreeAndNil(PathLst);
end;

constructor TAllPathBond.Create;
begin
  inherited Create;
  //Rep := TBondBase.Create;
end;

constructor TAllPathBond.Create(s: TStringList);
begin
  inherited Create;
  //Rep     := TBondBase.Create;
  fFrgLst := s;
end;

destructor TAllPathBond.Destroy;
begin
  //FreeAndNil(Rep);
  inherited Destroy;
end;

end.

