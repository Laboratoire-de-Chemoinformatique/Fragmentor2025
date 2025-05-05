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
unit UnitACFXBase;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, contnrs, UnitSequences, UnitFragmentBase, UnitAtomSymbol,
  UnitMoleculeFrg, U_TYPE_GRAPHES, UnitMoleculeBase;

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
    fatinfrg: TList;
  public
    property AC: Node read fAC write fAC;
    property LIPathHolder: TObjectList read fLIPathHolder write fLIPathHolder;
    property atinfrg: TList read fatinfrg write fatinfrg;
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

  { TACFXBase }

  TACFXBase = class(TSequences)
  private
  protected
    RepAC: TAtomSymbol;
    fNP: Node;
    //Number of discovered paths, Number of t terminated path
    fL: TNodeInfo;
    //Cost, Length of each path, Lenght of t terminated path
    fV: TNodeCost;
    fPath: NodeMatrix; //Array of each path
    procedure PathToString(Mol: TMoleculeFrg; TOL: TObjectList);
      virtual; abstract;
    procedure PathToString(Mol: TMoleculeFrg; LOPath: TIACHolder;
      TSL: TStringList; ACID: integer); virtual; abstract;
  public
    constructor Create;
    destructor Destroy; override;
    procedure MolToFrgLst(Mol: TMoleculeFrg); override;
    procedure PathFilter(Mol: TMoleculeFrg);
    //procedure InitMF(mfstr: string; mf: TStringList); override;
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
var
  i: integer;
  pidx: PInteger;
begin
  fAC:=0;
  fLIPathHolder.Clear;
  for i:=0 to fatinfrg.Count-1 do
  begin
    pidx:=PInteger(fatinfrg[i]);
    if pidx<>nil then
    begin
      dispose(pidx);
      pidx:=nil;
    end;
  end;
  fatinfrg.Clear;
end;

constructor TIACHolder.Create;
begin
  fAC:=0;
  fLIPathHolder:=TObjectList.create;
  fLIPathHolder.OwnsObjects:=True;
  fatinfrg:=TList.Create;
end;

destructor TIACHolder.Destroy;
begin
  FreeAndNil(fLIPathHolder);
  FreeAndNil(fatinfrg);
  inherited Destroy;
end;

{ TACFXBase }

procedure TACFXBase.PathFilter(Mol: TMoleculeFrg);
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
    if (fL[u] >= LenMin) then //minimal length
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

constructor TACFXBase.Create;
var
  u,v: Node;
begin
  inherited Create;
  RepAC := TAtomSymbol.Create;
  for u:=0 to MaxNode do
    for v:=0 to MaxNode do
      fPath[u,v]:=0;
end;

destructor TACFXBase.Destroy;
begin
  FreeAndNil(RepAC);
  inherited Destroy;
end;

procedure TACFXBase.MolToFrgLst(Mol: TMoleculeFrg);
var
  P: TNodeInfo;
  ACFrgLst: TStringList;
  LOPI: TObjectList;
  aACHlder: TIACHolder;
  aPthHlder: TIPathHolder;
  iw, idx, u, v: integer;
  piw, pat, pu: PInteger;
  pidx: PRAtmFrg;
  ii: integer;
  s,t : Node;
  DoAdd, DoAddF: Boolean;
  lf: integer;
  W: TArcCost;
  e: BondID;
begin
  ACFrgLst:=TStringList.Create;
  InitColorPerAtom(Mol);
  LOPI:=TObjectList.create;
  LOPI.OwnsObjects:=True;
  for lf:=LenMin to LenMax do
    LOPI.Add(TIACHolder.Create);
  for e:=0 to Mol.nBonds do W[e]:=1;
  for s := 1 to Mol.nAtom do
  begin
    for lf:=0 to LOPI.Count-1 do
    begin
      aACHlder:=LOPI[lf] as TIACHolder;
      aACHlder.AC:=s;
    end;
    for t := 1 to Mol.nAtom do
    begin
      Mol.YenEq(W,s,t,fPath,fL,fV,fNP,LenMax);//fPath is indexed at 1
      if (fV[1] + 1 >= LenMin) and (fV[1] + 1 <= LenMax) then //Check for smallest size of the k-path
      begin
        PathFilter(Mol);
        lf:=fV[1]+1-LenMin;//find index for the fragment holder containg fragments of a given size
        aACHlder:=LOPI[lf] as TIACHolder;
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
              //write(IntToStr(fPath[u,v])+'-');//GM
              DoAdd:=True;
              if IsPair and ((v<>1) and (v<>fL[u]) )then
                DoAdd:=False;
              while((ii<aACHlder.fatinfrg.Count) and (DoAdd)) do
              begin
                if PInteger(aACHlder.fatinfrg[ii])^=fPath[u,v] then
                   DoAdd:=False;
                Inc(ii);
              end;
              if DoAdd then
              begin
                new(pat); pat^:=fPath[u,v];
                aACHlder.fatinfrg.Add(pat);
              end;
            end;
            //writeln;//GM
        end;
      end;
    end;
    //Here all paths in fragment centered at s are recorded in aACHlder
    for lf:=0 to LOPI.Count-1 do
    begin
      aACHlder:=LOPI[lf] as TIACHolder;
      //Scan all stored path. If the fragment contains no valid path, it is skipped
      DoAddF:=False;
      for ii:=0 to aACHlder.LIPathHolder.Count-1 do
      begin
        aPthHlder:=aACHlder.LIPathHolder[ii] as TIPathHolder;
        if (aPthHlder.IPath[0]=1) then DoAddF:=True;
      end;
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
            //writeln;//GM
            for v:=0 to aACHlder.fatinfrg.Count-1 do
            begin
              pat:=PInteger(aACHlder.fatinfrg[v]);
              //write(IntToStr(pat^)+',');//GM
              new(pidx); pidx^.idx:=idx; pidx^.len:=aACHlder.fatinfrg.Count;
              (fFrgPerAtom[pat^] as TList).Add(TObject(pidx));//Attention à la destruction
            end;
          end;
        end;
      end;
      aACHlder.Clear;
      if ACFrgLst.Count>0 then
      begin
        for ii:=0 to ACFrgLst.Count-1 do
        begin
          piw:= PInteger(ACFrgLst.Objects[ii]);
          if piw<>nil then dispose(piw);
        end;
      end;
      ACFrgLst.Clear;
    end;
    for lf:=0 to LOPI.Count-1 do
    begin
      aACHlder:=LOPI[lf] as TIACHolder;
      aACHlder.Clear;
    end;
  end;
  FreeAndNil(ACFrgLst);
  FreeAndNil(LOPI);
end;

end.
