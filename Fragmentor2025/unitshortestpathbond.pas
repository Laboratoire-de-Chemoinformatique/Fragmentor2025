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
unit UnitShortestPathBond;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, strutils, unitfragmentbase, unitsequences,
  UnitShortestPath, UnitMoleculeFrg, UnitMoleculeBase, unitAtomAndBondType,
  U_TYPE_GRAPHES; //UnitBondBase

type

  { TShrtPthBnd }

  TShrtPthBnd = class(TShortestPath)
  private
    //Rep: TBondBase;
  protected
    function PathToString(Mol: TMoleculeFrg): TStringList; override;
    procedure PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList); override;
    procedure PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList; PathLst: TList); override;
  public
    constructor Create;
    constructor Create(s: TStringList);
    destructor Destroy; override;
  end;

implementation

{ TShrtPthBnd }
procedure TShrtPthBnd.PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList);
var
  u, v:    Node;
  head, tail: AtomID;
  PB:      PRBond;
  FORWD, BCKWD: string;
  THERE, BACK, Code: string;
  piw:     PInteger;
  piwused: Boolean;
begin
  TSLFrg.Clear;
  new(piw);
  piw^    := 1;
  piwused := False;
  for u := 1 to fNP do
  begin
    if (fL[u] >= LenMin) then
    begin
      if (((DynBnd = 1) and (DynBondInFrg(fPath, u, fL[u], Mol))) or
        ((DynBnd = 2) and (AllDynBondInFrg(fPath, u, fL[u], Mol))) or (DynBnd = 0)) then
      begin
        piwused := True;
        BCKWD   := '';
        FORWD   := '';
        THERE   := '';
        BACK    := '';
        for v := 1 to fL[u] - 1 do
        begin
          tail  := fPath[u, v];
          head  := fPath[u, v + 1];
          PB    := Mol.FindBond(tail, head);
          Str(PB^.B,Code); THERE:=THERE+Code;
          //Rep.StereoBond := Mol.BdStereo(PB);
          //FORWD := FORWD + Rep.GetBondSymbol(PB^.B);
          FORWD := FORWD + PB^.S;
          tail:=fPath[u,fL[u]-v+1];
          head:=fPath[u,fL[u]-v];
          PB:=Mol.FindBond(tail,head);
          Str(PB^.B,Code); BACK:=BACK+Code;
          BCKWD := BCKWD + PB^.S;
          //Rep.StereoBond:=Mol.BdStereo(PB);
          //BCKWD:=BCKWD+Rep.GetBondSymbol(PB^.B);
        end;
        if THERE <= BACK then
          TSLFrg.AddObject(FORWD, TObject(piw))
        else
          TSLFrg.AddObject(BCKWD, TObject(piw));
      end;
    end;
  end;
  if (piwused = False) then
    dispose(piw);
end;

procedure TShrtPthBnd.PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList;
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
    for u := 1 to fNP do
      if fPath[u,0]=1 then
      begin
        ffw.Clear;
        fbw.Clear;
        FORWD := '';
        BCKWD := '';
        ffw.Add(FORWD);
        fbw.Add(BCKWD);
        for v := 1 to fL[u]-1 do
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
          tail := fPath[u, fL[u] - v + 1];
          head := fPath[u, fL[u] - v];
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

function TShrtPthBnd.PathToString(Mol: TMoleculeFrg): TStringList;
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

constructor TShrtPthBnd.Create;
begin
  inherited Create;
  //Rep := TBondBase.Create;
end;

constructor TShrtPthBnd.Create(s: TStringList);
begin
  inherited Create;
  //Rep     := TBondBase.Create;
  fFrgLst := s;
end;

destructor TShrtPthBnd.Destroy;
begin
  //FreeAndNil(Rep);
  inherited Destroy;
end;

end.

