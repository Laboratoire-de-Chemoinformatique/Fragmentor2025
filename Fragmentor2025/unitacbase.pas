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
unit UnitACBase;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, UnitSequences, UnitFragmentBase, UnitAtomSymbol,
  UnitMoleculeFrg, U_TYPE_GRAPHES, contnrs;

type

  PPathPos=^PathPos;
  PathPos=record
    Path, Pos: integer;
  end;

  { TIACHolder }

  TIACHolder = class(TObject)
  private
    fAC: Node;
    fLIPathHolder: TObjectList;
  public
    property AC: Node read fAC write fAC;
    property LIPathHolder: TObjectList read fLIPathHolder write fLIPathHolder;
    procedure Clear;
    constructor Create;
    destructor Destroy; override;
  end;

  { TSACHolder }

  TSACHolder = class(TObject)
  private
    fAC: string;
    fLSPathHolder: TObjectList;
  public
    property AC: string read fAC write fAC;
    property LSPathHolder: TObjectList read fLSPathHolder write fLSPathHolder;
    procedure Clear;
    constructor Create;
    destructor Destroy; override;
  end;

  { TIPathHolder }

  TIPathHolder=Class(TObject)
    ICst: integer;
    IPath: array of Node;// index: position in path
  end;

  { TSPathHolder }

  TSPathHolder=Class(TObject)
    SPath: array of string;// index: position in path
  end;


  { TACBase }

  TACBase = class(TSequences)
  private
  protected
    fRepAC: TAtomSymbol;
    fNP: Node;
    //Number of discovered paths, Number of t terminated path
    fL: TNodeInfo;
    //Cost, Length of each path, Lenght of t terminated path
    fV: TNodeCost;
    fPath: NodeMatrix; //Array of each path
    fatinfrg: TList;//atom id in list if in fragment
    procedure PathToString(Mol: TMoleculeFrg; TSL: TStringList);
      virtual; abstract;
    procedure PathToString(Mol: TMoleculeFrg; LOPath: TIACHolder;
      TSL: TStringList; ACID: integer); virtual; abstract;
  public
    constructor Create;
    destructor Destroy; override;
    procedure MolToFrgLst(Mol: TMoleculeFrg); override;
    procedure PathFilter(Mol:TMoleculeFrg) override;
  end;

implementation

{ TSACHolder }

procedure TSACHolder.Clear;
begin
  fAC:='§';
  fLSPathHolder.Clear;
end;

constructor TSACHolder.Create;
begin
  fAC:='§';
  fLSPathHolder:=TObjectList.create;
  fLSPathHolder.OwnsObjects:=True;
end;

destructor TSACHolder.Destroy;
begin
  FreeAndNil(fLSPathHolder);
  inherited Destroy;
end;

{ TIACHolder }

procedure TIACHolder.Clear;
begin
  fAC:=0;
  fLIPathHolder.Clear;
end;

constructor TIACHolder.Create;
begin
  fAC:=0;
  fLIPathHolder:=TObjectList.create;
  fLIPathHolder.OwnsObjects:=True;
end;

destructor TIACHolder.Destroy;
begin
  FreeAndNil(fLIPathHolder);
  inherited Destroy;
end;

{ TACBase }

constructor TACBase.Create;
begin
  inherited Create;
  fRepAC := TAtomSymbol.Create;
  fatinfrg:=TList.Create;
end;

destructor TACBase.Destroy;
var
  pat: PInteger;
  i: integer;
begin
  FreeAndNil(fRepAC);
  for i:=0 to fatinfrg.Count-1 do
  begin
    pat:=PInteger(fatinfrg[i]);
    if (pat<>nil) then
    begin
     dispose(pat);
     pat:=nil;
    end;
  end;
  FreeAndNil(fatinfrg);
  inherited Destroy;
end;

procedure TACBase.MolToFrgLst(Mol: TMoleculeFrg);
var
  P: TNodeInfo;
  ACFrgLst: TStringList;
  aACHlder: TIACHolder;
  aPthHlder: TIPathHolder;
  iw, idx, u, v: integer;
  piw, pat, pu: PInteger;
  pidx: PRAtmFrg;
  ii: integer;
  s,t : Node;
  DoAdd, DoAddF: Boolean;
  W: TArcCost;
  e: ArcNum;
begin
  ACFrgLst:=TStringList.Create;
  InitColorPerAtom(Mol);
  aACHlder:=TIACHolder.Create;
  for e:=0 to Mol.p_M do W[e]:=1;
  for s := 1 to Mol.nAtom do
  begin
    aACHlder.AC:=s;
    DoAddF:=False;
    for t := 1 to Mol.nAtom do
    begin
      Mol.YenEq(W,s,t,fPath,fL,fV,fNP,LenMax);//fPath is indexed at 1
      if (fV[1] + 1 >= LenMin) and (fV[1] + 1 <= LenMax) then //Check the length of the shortest of k-path
      begin
        PathFilter(Mol);
          for u:=1 to fNP do
          begin
              //store reference to path u starting at s
              aACHlder.LIPathHolder.Add(TIPathHolder.Create);
              aPthHlder:=aACHlder.LIPathHolder.Last as TIPathHolder;
              aPthHlder.ICst:=fV[u];
              SetLength(aPthHlder.IPath,fL[u]+1);
              aPthHlder.IPath[0]:=fPath[u,0];//Index 0 is a keep/not keep indicator variable
              for v := 1 to fL[u] do
              begin
                aPthHlder.IPath[v]:=fPath[u,v];
                ii:=0;
                DoAdd:=True;
                if IsPair and ((v<>1) and (v<>fL[u])) then
                  DoAdd:=False;
                while((ii<fatinfrg.Count) and (DoAdd)) do
                begin
                  if PInteger(fatinfrg[ii])^=fPath[u,v] then
                     DoAdd:=False;
                  Inc(ii);
                end;
                if DoAdd then
                begin
                  new(pat); pat^:=fPath[u,v];
                  fatinfrg.Add(pat);
                end;
              end;
          end;
      end;
    end;
    //Scan all stored path. If the fragment contains no valid path, it is skipped
    DoAddF:=False;
    for ii:=0 to aACHlder.LIPathHolder.Count-1 do
    begin
      aPthHlder:=aACHlder.LIPathHolder[ii] as TIPathHolder;
      if (aPthHlder.IPath[0]=1) then DoAddF:=True;
    end;
    //Here all paths in fragment centered at s are recorded in aACHlder
    if DoAddF then
    begin
      PathToString(Mol,aACHlder,ACFrgLst,s);
      for ii := 0 to ACFrgLst.Count - 1 do
      begin
        piw:= PInteger(ACFrgLst.Objects[ii]);
        iw := piw^;
        idx:= AddFragment(ACFrgLst[ii], iw);
        if fGetFrgPerAtom then
        begin
          for v:=0 to fatinfrg.Count-1 do
          begin
            pat:=PInteger(fatinfrg[v]);
            new(pidx); pidx^.idx:=idx; pidx^.len:=fatinfrg.Count;
            (fFrgPerAtom[pat^] as TList).Add(TObject(pidx));//Attention à la destruction
          end;
        end;
      end;
    end;
    aACHlder.Clear;
    for ii:=0 to ACFrgLst.Count-1 do
    begin
      piw:= PInteger(ACFrgLst.Objects[ii]);
      if piw<>nil then dispose(piw);
    end;
    for ii:=0 to fatinfrg.Count-1 do
    begin
      pat:= PInteger(fatinfrg[ii]);
      if pat<>nil then dispose(pat);
    end;
    ACFrgLst.Clear;
    fatinfrg.Clear;
  end;
  FreeAndNil(ACFrgLst);
  FreeAndNil(aACHlder);
end;

procedure TACBase.PathFilter(Mol: TMoleculeFrg);
var
  u: Node;
  s,t,l :Node;
  uu: Node;
begin
  //Filtering of fragments.
  //Return True if at least one fragment is not filtered out
  for u := 1 to fNP do
  begin
    fPath[u,0]:=0;
    l:=fL[u];
    s:=fPath[u,1];//first atom of the fragment
    t:=fPath[u,l];//last atom of the fragment
//    for uu:=1 to l do
//      write(IntToStr(fPath[u,uu])+' ');
    if (fL[u] >= LenMin) then //minimal length
      if ((MarkAtom = 0) or (MarkAtom = 3) or ((MarkAtom = 2) and //MA0 or MA3
        (MrkAtmInFrg(fPath, u, l, Mol))) or //MA2 and mark atom in fragment
        ((MarkAtom = 1) and (Mol.IsMarkedAt(s) or Mol.IsMarkedAt(t))) or //MA1 and mark atom start or end fragment
        ((MarkAtom = 4) and Mol.IsMarkedAt(s) and Mol.IsMarkedAt(t))) then //MA4 and mark atoms start and end fragment
        if (((DynBnd = 1) and (DynBondInFrg(fPath, u, l, Mol))) or //fragment contains a dynamic bond
          ((DynBnd = 2) and (AllDynBondInFrg(fPath, u, l, Mol))) or //fragment contains only dynamic bonds
          (DynBnd = 0)) then   //dynamic bond is irrelavant
          begin
            fPath[u,0]:=1;
          end;
//    writeln(' ->'+IntToStr(fPath[u,0]));
  end;
//  writeln('+++++++++');
end;

end.
