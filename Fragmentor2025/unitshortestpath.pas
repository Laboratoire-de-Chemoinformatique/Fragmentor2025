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
unit UnitShortestPath;
//Unit to fragment based on shortest path sequences
{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Math, unitfragmentbase, UnitMoleculeFrg, UnitMoleculeBase,
  unitAtomAndBondType, U_TYPE_GRAPHES, unitsequences;

type

  TShrtPthException = class(Exception)
  end;

  { TShortestPath }

  TShortestPath = class(TSequences)
  private
  protected
    fNP: Node; //Number of discovered paths, Number of t terminated path
    fL: TNodeInfo;
    //Cost, Length of each path, Lenght of t terminated path
    fV: TNodeCost;
    fPath: NodeMatrix; //Array of each path
    ffw,fbw,fmf: TStringList;
    function PathToString(Mol: TMoleculeFrg): TStringList; virtual; abstract;
    procedure PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList); virtual; abstract;
    procedure PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList; PathLst: TList); virtual; abstract;
  public
    constructor Create;
    constructor Create(s: TStringList);
    destructor Destroy; override;
    procedure MolToFrgLst(Mol: TMoleculeFrg); override;
    procedure PathFilter(Mol:TMoleculeFrg); override;
  end;

implementation

{ TShortestPath }

constructor TShortestPath.Create;
var
  i,j: Node;
begin
  ffw:=TStringList.Create;
  fbw:=TStringList.Create;
  fmf:=TStringList.Create;
  for i:=Low(fPath) to High(fPath) do
    for j:=Low(fPath[i]) to High(fPath[i]) do
      fPath[i,j]:=0;
  for i:=Low(fL) to High(fL) do
    fL[i]:=0;
  for i:=Low(fV) to High(fV) do
    fV[i]:=0;
  fNP:=0;
  inherited Create;
end;

constructor TShortestPath.Create(s: TStringList);
var
  i,j: Node;
begin
  ffw:=TStringList.Create;
  fbw:=TStringList.Create;
  fmf:=TStringList.Create;
  for i:=Low(fPath) to High(fPath) do
    for j:=Low(fPath[i]) to High(fPath[i]) do
      fPath[i,j]:=0;
  for i:=Low(fL) to High(fL) do
    fL[i]:=0;
  for i:=Low(fV) to High(fV) do
    fV[i]:=0;
  fNP:=0;
  inherited Create(s)
end;

destructor TShortestPath.Destroy;
var
  i,j: Node;
begin
  FreeAndNil(ffw);
  FreeAndNil(fbw);
  FreeAndNil(fmf);
  for i:=Low(fPath) to High(fPath) do
    for j:=Low(fPath[i]) to High(fPath[i]) do
      fPath[i,j]:=0;
  for i:=Low(fL) to High(fL) do
    fL[i]:=0;
  for i:=Low(fV) to High(fV) do
    fV[i]:=0;
  fNP:=0;
  inherited Destroy;
end;

procedure TShortestPath.MolToFrgLst(Mol: TMoleculeFrg);
const
  WDflt=1;
var
  e: BondID;
  W: TArcCost;
  i, j: AtomID;
  s, t, u, v: Node;
  P: TNodeInfo;
  TSLFrg: TStringList;
  PathList: TList;
  ii, iw: integer;
  ptmpi:  PInteger;
  idx: integer;
  pidx: PRAtmFrg;
  tmpsl: TStringList;
  PM: PathMatrix;
begin
  //writeln('Adress ShortestPath='+IntToStr(Byte(Pointer(fFrgPerAtom))));
  TSLFrg := TStringList.Create;
  PathList:=TList.Create;
  //Initialization
  fColorPerAtom.Add(TStringList.Create);
  for s := 1 to Mol.nAtom do
  begin
    fL[s] := 0;
    for t := 1 to Mol.nAtom do
      fPath[s, t] := 0;
    fColorPerAtom.Add(TStringList.Create);
  end;
  for e:=0 to Mol.p_M do W[e]:=WDflt; //Default weights are expected on bonds
  //Process each atom pair
  for s := 1 to Mol.nAtom do
  begin
    for t := s + 1 to Mol.nAtom do
    begin
      //writeln('ICI '+IntToStr(s)+' '+IntToStr(t)+' '+IntToStr(maxint)+' '+IntToStr(maxsmallint));
      //Mol.Yen(W,s,t,fPath,fL,fV,4,LenMax);
      Mol.YenEq(W,s,t,fPath,fL,fV,fNP,LenMax);
      //writeln('ICI1');
      if (fV[1] + 1 >= LenMin) and (fV[1] + 1 <= LenMax) then //Check the cost of the optimal path
      begin
        PathFilter(Mol);//Filter useless path
        {for i:=1 to fNP do
        begin
          write('s='+IntToStr(s)+' t='+IntToStr(t)+' : V='+IntToStr(fV[i])+' : Keep='+IntToStr(fPath[i,0])+' : '+IntToStr(fPath[i,1]));
          for u:=2 to fL[i] do
            write('-'+IntToStr(fPath[i,u]));
          writeln;
        end;}
        PathToString(Mol, TSLFrg, PathList); //check for fNPath and fLen that are now useless
        for ii := 0 to TSLFrg.Count - 1 do
        begin
          iw := PInteger(TSLFrg.Objects[ii])^;
          idx:=AddFragment(TSLFrg[ii], iw); //Add fragments
          //writeln(IntToStr(iw)+' '+IntToStr(idx)+' '+TSLFrg[ii]);
          if fGetFrgPerAtom then  //Manage fragment info per atom if needed
          begin
            u:=PInteger(PathList[ii])^;
            if not (IsPair) then
              for v:=1 to fL[u] do
              begin
                new(pidx); pidx^.idx:=idx; pidx^.len:=fL[u];
                (fFrgPerAtom[fPath[u,v]] as TList).Add(TObject(pidx));//Attention à la destruction
              end
            else
            begin
              new(pidx); pidx^.idx:=idx; pidx^.len:=2;
              (fFrgPerAtom[fPath[u,1]] as TList).Add(TObject(pidx));//Attention à la destruction
              new(pidx); pidx^.idx:=idx; pidx^.len:=2;
              (fFrgPerAtom[fPath[u,fL[u]]] as TList).Add(TObject(pidx));//Attention à la destruction
            end;
          end;
        end;
        if PathList.Count>0 then //Release memory
          for ii:=0 to PathList.Count-1 do
          begin
            ptmpi:=PInteger(PathList[ii]);
            dispose(ptmpi);
          end;
        PathList.Clear;
        if TSLFrg.Count>0 then
          for ii := 0 to TSLFrg.Count - 1 do
          begin
            ptmpi := Pinteger(TSLFrg.Objects[ii]);
            dispose(ptmpi);
          end;
        TSLFrg.Clear;
      end;
    end;
    for u := 1 to fNP do //Reinitialize
    begin
      for v := 0 to fL[u] do
        fPath[u, v] := 0;
      fL[u] := 0;
      fV[u]:=0;
    end;
    fNP := 0;
  end;
  fColorPerAtom.Clear; //Finishing
  FreeAndNil(PathList);
  FreeAndNil(TSLFrg);
end;

procedure TShortestPath.PathFilter(Mol: TMoleculeFrg);
var
  u: Node;
  s,t,l :Node;
begin
  //Filtering of fragments.
  for u := 1 to fNP do
  begin
    fPath[u,0]:=0;
    l:=fL[u];
    s:=fPath[u,1];//first atom of the fragment
    t:=fPath[u,l];//last atom of the fragment
    //writeln('U'+IntToStr(s)+' '+IntToStr(t));
    if (fV[u]+1 >= LenMin) then //minimal length
      if ((MarkAtom = 0) or (MarkAtom = 3) or ((MarkAtom = 2) and //MA0 or MA3
        (MrkAtmInFrg(fPath, u, l, Mol))) or //MA2 and mark atom in fragment
        ((MarkAtom = 1) and (Mol.IsMarkedAt(s) or Mol.IsMarkedAt(t))) or //MA1 and mark atom start or end fragment
        ((MarkAtom = 4) and Mol.IsMarkedAt(s) and Mol.IsMarkedAt(t))) then //MA4 and mark atoms start and end fragment
        if (((DynBnd = 1) and (DynBondInFrg(fPath, u, l, Mol))) or //fragment contains a dynamic bond
          ((DynBnd = 2) and (AllDynBondInFrg(fPath, u, l, Mol))) or //fragment contains only dynamic bonds
          (DynBnd = 0)) then   //dynamic bond is irrelavant
          begin fPath[u,0]:=1; end; //writeln('F'+IntToStr(s)+' '+IntToStr(t)); end;
  end;
end;

end.
