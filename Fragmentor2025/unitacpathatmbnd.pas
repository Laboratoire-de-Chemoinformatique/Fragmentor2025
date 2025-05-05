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
unit UnitACPathAtmBnd;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, strutils, UnitFragmentBase, UnitACBase, UnitAtomSymbol, UnitMoleculeBase,
  UnitMoleculeFrg, unitAtomAndBondType, U_TYPE_GRAPHES, UnitAtomBase, contnrs;// cwstring; //UnitBondBase,
//Remark cwstring might be needed for fragment sorting into a centered fragment
//containing non alpha-numeric characters (=,#,etc). If needed, requires
//Lazarus1.2 or later.

type

  { TACPathBnd }

  { TACPathAtmBnd }

  TACPathAtmBnd = class(TACBase)
  private
    //RepB: TBondBase;
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

{ TACPathAtmBnd }

procedure TACPathAtmBnd.PathToString(Mol: TMoleculeFrg; TSL: TStringList);
var
  u, v:  Node;
  head, tail: AtomID;
  PA:    PRAtom;
  PB:    PRBond;
  FORWD: string;
begin
  for u := 1 to fNP do
  begin
    if (fL[u] >= LenMin) then
    begin
      FORWD := '';
      for v := 1 to fL[u]-1 do
      begin
        tail  := fPath[u, v];
        head  := fPath[u, v + 1];
        PB    := Mol.FindBond(tail, head);
        PA    := Mol.AtmSet[tail];
        fRepAC.StereoParity := Mol.AtStereoParity(PA);
        //RepB.StereoBond := Mol.BdStereo(PB);
        FORWD := FORWD + fRepAC.AtomString[tail] + PB^.S;
      end;
      PA    := Mol.AtmSet[fPath[u, fL[u]]];
      fRepAC.StereoParity := Mol.AtStereoParity(PA);
      FORWD := FORWD + fRepAC.AtomString[fPath[u, fL[u]]];
      TSL.Add('(' + FORWD + ')');
    end;
  end;
end;

procedure TACPathAtmBnd.PathToString(Mol: TMoleculeFrg; LOPath: TIACHolder;
  TSL: TStringList; ACID: integer);
var
   s, tail, head: Node;
   i, j, k, ii, jj, kk, kkmax: integer;
   col, sze: integer;
   iat, at: integer;
   iw, u: integer;
   pat, piw, pu: PInteger;
   FORWD,sFrg: string;
   TmpSL: TStringList;
   LOPS: TObjectList;
   mf: TStringList;
   aACHold, aACHold0: TSACHolder;
   aIPthHold: TIPathHolder;
   aSPthHold, aSPthHold0: TSPathHolder;
   LPathPos: TList;
   PPthPs: PPathPos;
   PB: PRBond;
begin
  //LOPath contains all path starting at a given atom
  //LOPS contains all alternatives of given fragment
  TmpSL:=TStringList.Create;
  TmpSL.Duplicates:=dupAccept;
  TmpSL.Delimiter:=',';
  LOPS:=TObjectList.create;
  LOPS.OwnsObjects:=True;
  LPathPos:=TList.Create;
  for col := 0 to ColorKeys.Count - 1 do
  begin
    InitMF(ColorKeys[col],iw,TAtomBase(fRepAC),Mol);
    LOPS.Add(TSACHolder.Create);
    aACHold:=LOPS.Last as TSACHolder;
    //For each path reserve space to store atomic representation
    for ii:=0 to LOPath.LIPathHolder.Count-1 do
    begin
      aIPthHold:=LOPath.LIPathHolder[ii] as TIPathHolder;
      sze:=Length(aIPthHold.IPath);
      aACHold.LSPathHolder.Add(TSPathHolder.Create);
      aSPthHold:=aACHold.LSPathHolder.Last as TSPathHolder;
      SetLength(aSPthHold.SPath,sze);
      for jj:=0 to High(aSPthHold.SPath) do
        aSPthHold.SPath[jj]:='§';
    end;
    for iat:=0 to fatinfrg.Count-1 do
    begin
      pat:=PInteger(fatinfrg[iat]);
      at:=pat^;
      mf:=fColorPerAtom[at] as TStringList;
      if at=ACID then
      begin
        //Manage central atom
        for ii:=0 to LOPS.Count-1 do
        begin
          aACHold:=LOPS[ii] as TSACHolder;
          aACHold.AC:=mf[0];
        end;
      end;
      for j:=0 to LOPath.LIPathHolder.Count-1 do
      begin
        aIPthHold:=LOPath.LIPathHolder[j] as TIPathHolder;
        for s:=1 to High(aIPthHold.IPath) do
        begin
          if (aIPthHold.IPath[s]=at) then
          begin
            //Beware to add alternative colorations to all existing element of LOPS
            for ii:=0 to LOPS.Count-1 do
            begin
              aACHold:=LOPS[ii] as TSACHolder;
              aSPthHold:=aACHold.LSPathHolder[j] as TSPathHolder;
              aSPthHold.SPath[s]:=mf[0];
            end;
            new(PPthPs);
            PPthPs^.Path:=j;
            PPthPs^.Pos:=s;
            LPathPos.Add(PPthPs);
          end;
        end;
      end;
      for k:=1 to mf.Count-1 do
      begin
        //Add a copy of the current state of the fragment
        kkmax:=LOPS.Count-1;
        for kk:=0 to kkmax do
        begin
          LOPS.Add(TSACHolder.Create);
          aACHold:=LOPS.Last as TSACHolder;
          aACHold0:=LOPS[0] as TSACHolder;
          for ii:=0 to aACHold0.LSPathHolder.Count-1 do
          begin
            aSPthHold0:=aACHold0.LSPathHolder[ii] as TSPathHolder;
            sze:=Length(aSPthHold0.SPath);
            aACHold.LSPathHolder.Add(TSPathHolder.Create);
            aSPthHold:=aACHold.LSPathHolder.Last as TSPathHolder;
            SetLength(aSPthHold.SPath,sze);
            for jj:=0 to High(aSPthHold.SPath) do
              aSPthHold.SPath[jj]:=aSPthHold0.SPath[jj];
          end;
          //update representation of current atom
          if at=ACID then aACHold.AC:=mf[k]
          else aACHold.AC:=aACHold0.AC;
          for j:=0 to LPathPos.Count-1 do
          begin
            PPthPs:=PPathPos(LPathPos[j]);
            aSPthHold:=aACHold.LSPathHolder[PPthPs^.Path] as TSPathHolder;
            aSPthHold.SPath[PPthPs^.Pos]:=mf[k];
          end;
        end;
      end;
      for j:=0 to LPathPos.Count-1 do
      begin
        PPthPs:=PPathPos(LPathPos[j]);
        dispose(PPthPs);
      end;
      LPathPos.Clear;
    end;
    //Parse the LOPS. Each object is a list containing the representation of a
    //path relevent for the fragment. Each object must be translated to a string,
    //and combined canonicaly for a string representation of the fragment.
    TmpSL.Clear;
    for ii:=0 to LOPS.Count-1 do
    begin
      aACHold:=LOPS[ii] as TSACHolder;
      for i:=0 to aACHold.LSPathHolder.Count-1 do
      begin
        aIPthHold:=LOPath.LIPathHolder[i] as TIPathHolder;
        aSPthHold:=aACHold.LSPathHolder[i] as TSPathHolder;
        FORWD:='';
        for s:=1 to High(aSPthHold.SPath)-1 do
        begin
          tail:=aIPthHold.IPath[s];
          head:=aIPthHold.IPath[s+1];
          PB:=Mol.FindBond(tail,head);
          FORWD:=FORWD+aSPthHold.SPath[s]+PB^.S;
        end;
        FORWD:='('+FORWD+aSPthHold.SPath[High(aSPthHold.SPath)]+')';
        TmpSL.Add(FORWD);
      end;
      new(piw); piw^:=iw;
      //writeln('***');
      //writeln(TmpSL.DelimitedText);
      //writeln(TmpSL[0]+' '+TmpSL[3]+' '+IntToStr(WideCompareStr(TmpSL[0],TmpSL[3])));
      TmpSL.CustomSort(@TSLCompareStr);//TmpSL.Sort;
      if TmpSL.Count>0 then
      begin
        sFrg:=TmpSL.DelimitedText+',x'+aACHold.AC;
        //writeln('Frg='+sFrg);
        TSL.AddObject(sFrg, TObject(piw));
      end;
    end;
    LOPS.Clear;
  end;
  FreeAndNil(LPathPos);
  FreeAndNil(LOPS);
  FreeAndNil(TmpSL);
end;

constructor TACPathAtmBnd.Create;
begin
  inherited Create;
end;

constructor TACPathAtmBnd.Create(s: TStringList);
begin
  inherited Create;
  fFrgLst := s;
end;

destructor TACPathAtmBnd.Destroy;
begin
  inherited Destroy;
end;

end.

