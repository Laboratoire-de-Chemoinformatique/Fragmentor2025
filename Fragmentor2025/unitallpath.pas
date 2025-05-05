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
unit UnitAllPath;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, U_TYPE_GRAPHES,
  UnitMoleculeFrg, UnitMolecule, UnitSequences, UnitFragmentBase;

type

  TAllPathException = class(Exception)
  end;

  { TAllPath }

  TAllPath = class(TSequences)
  private
  protected
    fNPath: Node;       //Number of discovered paths
    fLen:   TNodeInfo;  //Length of each path in atom number
    fPath:  NodeMatrix; //Array of each path
    ffw,fbw,fmf: TStringList;
    function PathToString(Mol: TMoleculeFrg): TStringList; virtual; abstract;
    procedure PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList); virtual; abstract;
    procedure PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList; PathLst: TList); virtual; abstract;
  public
    constructor Create;
    destructor Destroy; override;
    procedure MolToFrgLst(Mol: TMoleculeFrg); override;
    procedure PathFilter(Mol: TMoleculeFrg); override;
  end;

implementation

{ TAllPath }

constructor TAllPath.Create;
begin
  inherited Create;
  NAS := 0;
  ffw:=TStringList.Create;
  fbw:=TStringList.Create;
  fmf:=TStringList.Create;
end;

destructor TAllPath.Destroy;
begin
  inherited Destroy;
  FreeAndNil(ffw);
  FreeAndNil(fbw);
  FreeAndNil(fmf);
end;

procedure TAllPath.MolToFrgLst(Mol: TMoleculeFrg);
var
  j:   AtomID;
  s, t, u, v: Node;
  TSLFrg: TStringList;
  PathList: TList;
  ii, iw, idx: integer;
  ptmpi:  Pinteger;
  pidx: PRAtmFrg;

begin
  //Initialization
  TSLFrg := TStringList.Create;
  PathList:=TList.Create;
  fNPath := 0;
  for s := 1 to Mol.nAtom do
  begin
    fLen[s] := 0;
    for j := 1 to Mol.nAtom do
      fPath[s, j] := 0;
  end;
  InitColorPerAtom(Mol);
  //Process each atom pair
  for s := 1 to Mol.nAtom do
  begin
    for t := s + 1 to Mol.nAtom do
    begin
      Mol.DFSAllPath(s, t, fNPath, fLen, fPath, LenMax);
      PathFilter(Mol);
      if fNPath > 0 then
      begin
        PathToString(Mol, TSLFrg, PathList);
        for ii := 0 to TSLFrg.Count - 1 do
        begin
          iw := PInteger(TSLFrg.Objects[ii])^;
          idx:=AddFragment(TSLFrg[ii], iw);
          if fGetFrgPerAtom then
          begin
            u:=PInteger(PathList[ii])^;
            for v:=1 to fLen[u] do
            begin
              new(pidx); pidx^.idx:=idx; pidx^.len:=fLen[u];
              (fFrgPerAtom[fPath[u,v]] as TList).Add(TObject(pidx));//Attention à la destruction
            end;
          end;
        end;
        for ii:=0 to PathList.Count-1 do
        begin
          ptmpi:=PInteger(PathList[ii]);
          dispose(ptmpi);
        end;
        PathList.Clear;
        for ii := 0 to TSLFrg.Count - 1 do
        begin
          ptmpi := Pinteger(TSLFrg.Objects[ii]);
          dispose(ptmpi);
        end;
        TSLFrg.Clear;
      end;
      for u := 1 to fNPath do
      begin
        for v := 1 to fLen[u] do
          fPath[u, v] := 0;
        fLen[u] := 0;
      end;
      fNPath := 0;
    end;
  end;
  fColorPerAtom.Clear;
  FreeAndNil(PathList);
  FreeAndNil(TSLFrg);
end;

procedure TAllPath.PathFilter(Mol: TMoleculeFrg);
var
  u: Node;
  s,t,l :Node;
begin
  //Filtering of fragments.
  for u := 1 to fNPath do
  begin
    fPath[u,0]:=0;
    l:=fLen[u];
    s:=fPath[u,1];//first atom of the fragment
    t:=fPath[u,l];//last atom of the fragment
    if (fLen[u] >= LenMin) then //minimal length
      if ((MarkAtom = 0) or (MarkAtom = 3) or ((MarkAtom = 2) and //MA0 or MA3
        (MrkAtmInFrg(fPath, u, l, Mol))) or //MA2 and mark atom in fragment
        ((MarkAtom = 1) and (Mol.IsMarkedAt(s) or Mol.IsMarkedAt(t))) or //MA1 and mark atom start or end fragment
        ((MarkAtom = 4) and Mol.IsMarkedAt(s) and Mol.IsMarkedAt(t))) then //MA4 and mark atoms start and end fragment
        if (((DynBnd = 1) and (DynBondInFrg(fPath, u, l, Mol))) or //fragment contains a dynamic bond
          ((DynBnd = 2) and (AllDynBondInFrg(fPath, u, l, Mol))) or //fragment contains only dynamic bonds
          (DynBnd = 0)) then   //dynamic bond is irrelavant
          fPath[u,0]:=1;
  end;
end;

end.

