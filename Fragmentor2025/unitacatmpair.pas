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
unit UnitACAtmPair;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, strutils, UnitFragmentBase, UnitACBase, UnitSequences, U_TYPE_GRAPHES,
  UnitMoleculeFrg, UnitAtomSymbol, UnitMoleculeBase, unitAtomAndBondType,
  UnitAtomBase, contnrs;

type
{ TACAtmPair }
  TACAtmPair = class(TACBase)
  private
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

{ TACAtmPair }

procedure TACAtmPair.PathToString(Mol: TMoleculeFrg; TSL: TStringList);
var
  u, v:  Node;
  shead, stail, thead, ttail: AtomID;
  PA:    PRAtom;
  FORWD: string;
  LocST: TStringList;
begin
  LocST := TStringList.Create;
  LocST.Sorted := True; //Ignore identical fragments between the same atom-bond pairs
  LocST.Duplicates := dupIgnore;
  for u := 1 to fNP do
  begin

    stail := fPath[u, 1];
    ttail := fPath[u, fL[u]];
    PA    := Mol.AtmSet[stail];
    fRepAC.StereoParity := Mol.AtStereoParity(PA);
    FORWD := fRepAC.AtomString[stail] + IntToStr(fV[u]);
    PA    := Mol.AtmSet[ttail];
    fRepAC.StereoParity := Mol.AtStereoParity(PA);;
    FORWD := FORWD + fRepAC.AtomString[ttail];

    LocST.Add('(' + FORWD + ')');
  end;
  TSL.AddStrings(LocST);
  FreeAndNil(LocST);
end;

procedure TACAtmPair.PathToString(Mol: TMoleculeFrg; LOPath: TIACHolder;
  TSL: TStringList; ACID: integer);
var
   s: Node;
   first, last: Node;
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
        first:=1;
        last:=High(aIPthHold.IPath);
        if (aIPthHold.IPath[first]=at) then
        begin
          //Beware to add alternative colorations to all existing element of LOPS
          for ii:=0 to LOPS.Count-1 do
          begin
            aACHold:=LOPS[ii] as TSACHolder;
            aSPthHold:=aACHold.LSPathHolder[j] as TSPathHolder;
            aSPthHold.SPath[first]:=mf[0];
          end;
          new(PPthPs);
          PPthPs^.Path:=j;
          PPthPs^.Pos:=first;
          LPathPos.Add(PPthPs);
        end;
        if (aIPthHold.IPath[last]=at) then
        begin
          //Beware to add alternative colorations to all existing element of LOPS
          for ii:=0 to LOPS.Count-1 do
          begin
            aACHold:=LOPS[ii] as TSACHolder;
            aSPthHold:=aACHold.LSPathHolder[j] as TSPathHolder;
            aSPthHold.SPath[last]:=mf[0];
          end;
          new(PPthPs);
          PPthPs^.Path:=j;
          PPthPs^.Pos:=last;
          LPathPos.Add(PPthPs);
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
          else aACHold0.AC:=aACHold0.AC;
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
        aSPthHold:=aACHold.LSPathHolder[i] as TSPathHolder;
        aIPthHold:=LOPath.LIPathHolder[i] as TIPathHolder;
        first:=1;
        last:=High(aSPthHold.SPath);
        FORWD:=aSPthHold.SPath[first]+IntToStr(aIPthHold.ICst)+aSPthHold.SPath[last];
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
  FreeAndNil(LPathPos);
  FreeAndNil(LOPS);
  FreeAndNil(TmpSL);
end;

constructor TACAtmPair.Create;
begin
  inherited Create;
end;

constructor TACAtmPair.Create(s: TStringList);
begin
  inherited Create;
  fFrgLst := s;
end;

destructor TACAtmPair.Destroy;
begin
  inherited Destroy;
end;

{type

  { TACAtmPair }

  TACAtmPair = class(TSequences)
  private
    RepAC: TAtomSymbol;
  protected
    fV: TNodeCost;
    procedure PathToString(s, t: Node; Mol: TMoleculeFrg; TSL: TStringList);
  public
    constructor Create;
    constructor Create(s: TStringList);
    destructor Destroy; override;
    procedure MolToFrgLst(Mol: TMoleculeFrg); override;
  end;

implementation

{ TACAtmPair }

procedure TACAtmPair.PathToString(s, t: Node; Mol: TMoleculeFrg; TSL: TStringList);
var
  u, v:  Node;
  FORWD: string;
begin
  RepAC.StereoParity := Mol.AtStereoParity(s);
  FORWD := RepAC.AtomString[s] + IntToStr(fV[t]);
  RepAC.StereoParity := Mol.AtStereoParity(t);
  FORWD := FORWD + RepAC.AtomString[t];

  TSL.Add('(' + FORWD + ')');
end;

constructor TACAtmPair.Create;
begin
  inherited Create;
  RepAC := TAtomSymbol.Create;
end;

constructor TACAtmPair.Create(s: TStringList);
begin
  inherited Create;
  RepAC   := TAtomSymbol.Create;
  fFrgLst := s;
end;

destructor TACAtmPair.Destroy;
begin
  FreeAndNil(RepAC);
  inherited Destroy;
end;

procedure TACAtmPair.MolToFrgLst(Mol: TMoleculeFrg);
var
  i, j:   AtomID;
  s, t, u, v: Node;
  P:      TNodeInfo;
  SFrg, strtmp:   string;
  ii, iw, jj, l, iimax:     integer;
  AW:     TArcCost;
  k:      ArcNum;
  TSLFrg, ms, mf, temp: TStringList;
begin
  //Initialization
  ms     := TStringList.Create;
  mf     := TStringList.Create;
  temp   := TStringList.Create;
  temp.StrictDelimiter := True;
  temp.Delimiter := ';';
  TSLFrg := TStringList.Create;
  TSLFrg.Sorted := True;
  TSLFrg.Duplicates := dupAccept;
  TSLFrg.Delimiter := ',';
  for l := 0 to ColorKeys.Count - 1 do
  begin
    ms.Clear;
    ms.Add('');
    if (ColorKeys[l] = 'Default') then
    begin
      strtmp := '';
      for i := 1 to Mol.nAtom do
        if (UseFormalCharge and Mol.IsChargedAt(i)) then
          //RepAC.AtomString[i] :=Mol.AtmSet[i]^.S + '&FC' + RepAC.ConvertFormalCharge(Mol.GetFormalCharge(i)) + '&'
          strtmp := strtmp + ';' + Mol.AtmSet[i]^.S + '&FC' +
            RepAC.ConvertFormalCharge(Mol.GetFormalCharge(i)) + '&'
        else
          //RepAC.AtomString[i] := Mol.AtmSet[i]^.S;
          strtmp := strtmp + ';' + Mol.AtmSet[i]^.S;
      ms[0] := strtmp;
      iw := 1;
    end
    else
    begin
      iw := RepAC.InitAtomStringWeight(ColorHash.Items[ColorKeys[l]]);
      for i := 1 to Mol.nAtom do
      begin
        InitMF(RepAC.AtomString[i], mf);
        if (UseFormalCharge and Mol.IsChargedAt(i)) then
        begin
          iimax := mf.Count - 1;
          for ii := 0 to iimax do
            mf.Add(mf[ii] + '&FC' + RepAC.ConvertFormalCharge(
              Mol.GetFormalCharge(i)) + '&');
          for ii := iimax downto 0 do
            mf.Delete(ii);
        end;
        iimax := ms.Count - 1;
        for ii := 0 to iimax do
          for jj := 0 to mf.Count - 1 do
            ms.Add(ms[ii] + ';' + mf[jj]);
        for ii := iimax downto 0 do
          ms.Delete(ii);
      end;
    end;

    for jj := 0 to ms.Count - 1 do
    begin
      temp.Clear;
      temp.DelimitedText := ms[jj];
      for ii := 1 to temp.Count - 1 do
      begin
        RepAC.AtomString[ii] := temp[ii];
      end;
      for i := 1 to Mol.nAtom do
        P[i] := 0;
      for k := 1 to Mol.p_M do
        AW[k] := 1;
      //Process each atom pair
      for s := 1 to Mol.nAtom do
      begin
        TSLFrg.Clear;
        if ((not UseMarkAtom) or (UseMarkAtom and Mol.IsMarkedAt(s))) then
        begin
          for t := s + 1 to Mol.nAtom do
            if ((not UseMarkAtom) or (UseMarkAtom and Mol.IsMarkedAt(t))) then
            begin
              Mol.DijHeap(AW, s, t, fV, P);
              //Path ending at t has a cost containing in fV[t]
              //Here dynamic bonds restrictions should be tested.
              if (fV[t] + 1 >= LenMin) and (fV[t] + 1 <= LenMax) then
                PathToString(s, t, Mol, TSLFrg);
            end;
          TSLFrg.Sort;
          if TSLFrg.Count > 0 then
          begin
            SFrg := TSLFrg.DelimitedText;
            RepAC.StereoParity := Mol.AtStereoParity(s);
            SFrg := SFrg + ',x' + RepAC.AtomString[s];
            AddFragment(SFrg);
          end;
        end;
      end;
    end;
  end;
  FreeAndNil(TSLFrg);
  FreeAndNil(ms);
  FreeAndNil(temp);
  FreeAndNil(mf);
end;}

end.

