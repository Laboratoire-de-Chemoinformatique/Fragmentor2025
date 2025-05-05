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
unit UnitACFXPathAtmBnd;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, contnrs, UnitFragmentBase, UnitACFXBase, UnitMoleculeFrg, UnitMoleculeBase,
  unitAtomAndBondType, U_TYPE_GRAPHES, UnitAtomBase; //UnitBondBase,

type

  { TACFXPathAtmBnd }

  TACFXPathAtmBnd = class(TACFXBase)
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

{ TACPathAtmBnd }

procedure TACFXPathAtmBnd.PathToString(Mol: TMoleculeFrg; TOL: TObjectList);
var
  u, v, vmax: Node;
  head, tail: AtomID;
  PA:    PRAtom;
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
      PA    := Mol.AtmSet[tail];
      RepAC.StereoParity := Mol.AtStereoParity(PA);
      //RepB.StereoBond := Mol.BdStereo(PB);
      FORWD := FORWD + RepAC.AtomString[tail] + PB^.S;
    end;
    PA    := Mol.AtmSet[fPath[u, fL[u]]];
    RepAC.StereoParity := Mol.AtStereoParity(PA);
    FORWD := FORWD + RepAC.AtomString[fPath[u, fL[u]]];
    (TOL.Items[fV[fPath[u, vmax]] - LenMin + 1] as TStringList).Add(
      '(' + FORWD + ')');
  end;
end;

procedure TACFXPathAtmBnd.PathToString(Mol: TMoleculeFrg; LOPath: TIACHolder;
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
  TmpSL.CaseSensitive:=True;
  LOPS:=TObjectList.create;
  LOPS.OwnsObjects:=True;
  LPathPos:=TList.Create;
  for col := 0 to ColorKeys.Count - 1 do
  begin
    InitMF(ColorKeys[col],iw,TAtomBase(RepAC),Mol);
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
    for iat:=0 to LOPath.atinfrg.Count-1 do
    begin
      pat:=PInteger(LOPath.atinfrg[iat]);
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
      {if TmpSL.Count>1 then
      begin
        writeln('. '+TmpSL[0]+' '+TmpSL[1]+' AnsiCompareStr='+IntToStr(AnsiCompareStr(TmpSL[0],TmpSL[1]))+' AnsiCompareText='+IntToStr(AnsiCompareText(TmpSL[0],TmpSL[1])));
        writeln(', '+TmpSL[0]+' '+TmpSL[1]+' CompareStr='+IntToStr(CompareStr(TmpSL[0],TmpSL[1]))+' CompareText='+IntToStr(CompareText(TmpSL[0],TmpSL[1])));
      end;}
      TmpSL.CustomSort(@TSLCompareStr);//TmpSL.Sort;
      {if TmpSL.Count>1 then
      begin
        writeln('.. '+TmpSL[0]+' '+TmpSL[1]+' AnsiCompareStr='+IntToStr(AnsiCompareStr(TmpSL[0],TmpSL[1]))+' AnsiCompareText='+IntToStr(AnsiCompareText(TmpSL[0],TmpSL[1])));
        writeln(',, '+TmpSL[0]+' '+TmpSL[1]+' CompareStr='+IntToStr(CompareStr(TmpSL[0],TmpSL[1]))+' CompareText='+IntToStr(CompareText(TmpSL[0],TmpSL[1])));
      end;}
      if TmpSL.Count>0 then
      begin
        new(piw); piw^:=iw;
        sFrg:=TmpSL.DelimitedText+',x'+aACHold.AC;
        //writeln('ACFXPath '+sFrg);
        TSL.AddObject(sFrg, TObject(piw));
      end;
    end;
    LOPS.Clear;
  end;
  FreeAndNil(LPathPos);
  FreeAndNil(LOPS);
  FreeAndNil(TmpSL);
end;

constructor TACFXPathAtmBnd.Create;
begin
  inherited Create;
  //RepB := TBondBase.Create;
end;

constructor TACFXPathAtmBnd.Create(s: TStringList);
begin
  inherited Create;
  //RepB    := TBondBase.Create;
  fFrgLst := s;
end;

destructor TACFXPathAtmBnd.Destroy;
begin
  //FreeAndNil(RepB);
  inherited Destroy;
end;

end.
