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
unit UnitShortestPathAtmPair;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, strutils, UnitSequences, UnitShortestPath,
  unitatomsymbol, unitfragmentbase,UnitAtomBase,
  UnitMoleculeFrg, U_TYPE_GRAPHES, unitAtomAndBondType, UnitMoleculeBase;

type

  { TShrtPthAtmPair }

  //TShrtPthAtmPair = class(TSequences)
  TShrtPthAtmPair = class(TShortestPath)
  private
    fRepA: TAtomSymbol;
  protected
    //fV: TNodeCost;
    //function PathToString(Mol: TMoleculeFrg; s, t: AtomID;w: Cost): TStringList; virtual;
    function PathToString(Mol: TMoleculeFrg): TStringList; override;
    procedure PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList); override;
    procedure PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList; PathLst: TList); override;
  public
    constructor Create;
    constructor Create(s: TStringList);
    destructor Destroy; override;
    procedure InitAtomString(sdline: string); override;
    //procedure InitMF(mfstr: string; mf: TStringList);
    {procedure MolToFrgLst(Mol: TMoleculeFrg); override;}
  end;

implementation

{ TShrtPthAtmPair }
procedure TShrtPthAtmPair.PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList);
var
  u, v: Node;
  shead, stail, thead, ttail: AtomID;
  PA: PRAtom;
  PB: PRBond;
  FORWD, BCKWD: string;
  mf, fw, bw: TStringList;
  i, j, jj, ii, iimax: integer;
  piw: PInteger;
  piwused: boolean;
begin
  TSLFrg.Clear;
  for j := 0 to ColorKeys.Count - 1 do
  begin
    new(piw);
    piw^ := 1;
    piwused := False;
    if (ColorKeys[j] = 'Default') then
    begin
      //Itialisation of representation and weight
      for i := 1 to Mol.nAtom do
      begin
        fRepA.AtomString[i] := Mol.S_[i];
        if ((MarkAtom = 3) and (Mol.IsMarkedAt(i))) then
          fRepA.AtomString[i] := fRepA.AtomString[i] + '&MA&';
        if (UseFormalCharge and Mol.IsChargedAt(i)) then
          fRepA.AtomString[i] :=
            fRepA.AtomString[i] + '&FC' + fRepA.ConvertFormalCharge(
            Mol.GetFormalCharge(i)) + '&';
      end;
      //Generation of pair from path
      for u := 1 to fNP do
      begin
        if fPath[u,0]=1 then
        {if (fL[u] >= LenMin) then
        begin
          if (((DynBnd = 1) and (DynBondInFrg(fPath, fNPath + u, fL[u], Mol))) or
            ((DynBnd = 2) and (AllDynBondInFrg(fPath, fNPath + u, fL[u], Mol))) or
            (DynBnd = 0)) then}
          begin
            if ((MarkAtom = 0) or (MarkAtom = 1) or (MarkAtom = 3) or ((MarkAtom = 2) and
              (MrkAtmInFrg(fPath, u, fL[u], Mol)))) then
            begin
              piwused := True;
              stail := fPath[u, 1];
              ttail := fPath[u, fL[u]];

              PA := Mol.AtmSet[stail];
              fRepA.StereoParity := Mol.AtStereoParity(PA);
              FORWD := fRepA.AtomString[stail] + IntToStr(fV[u]);
              BCKWD := IntToStr(fV[u]) + fRepA.AtomString[stail];
              PA := Mol.AtmSet[ttail];
              fRepA.StereoParity := Mol.AtStereoParity(PA);
              FORWD := FORWD + fRepA.AtomString[ttail];
              BCKWD := fRepA.AtomString[ttail] + BCKWD;

          {PA    := Mol.AtmSet[ttail];
          fRepA.StereoParity := Mol.AtStereoParity(PA);
          BCKWD := fRepA.AtomString[ttail] + IntToStr(fV[ttail]);
          PA    := Mol.AtmSet[stail];
          fRepA.StereoParity := Mol.AtStereoParity(PA);
          BCKWD := BCKWD + fRepA.AtomString[stail];}

              if FORWD <= BCKWD then
                TSLFrg.AddObject(FORWD, TObject(piw))
              else
                TSLFrg.AddObject(BCKWD, TObject(piw));
            end;
          end;

        //end;
      end;
    end
    else
    begin
      mf := TStringList.Create;
      fw := TStringList.Create;
      bw := TStringList.Create;
      piw^ := fRepA.InitAtomStringWeight(ColorHash.Items[ColorKeys[j]]);
      for u := 1 to fNP do
      begin
        if fPath[u,0]=1 then
        {if (fL[u] >= LenMin) then
        begin
          if (((DynBnd = 1) and (DynBondInFrg(fPath, fNPath + u, fL[u], Mol))) or
            ((DynBnd = 2) and (AllDynBondInFrg(fPath, fNPath + u, fL[u], Mol))) or
            (DynBnd = 0)) then}
          begin
            if ((MarkAtom = 0) or (MarkAtom = 1) or (MarkAtom = 3) or((MarkAtom = 2) and
              (MrkAtmInFrg(fPath, u, fL[u], Mol)))) then
            begin
              piwused := True;
              fw.Clear;
              bw.Clear;

              stail := fPath[u, 1];
              ttail := fPath[u, fL[u]];

              PA := Mol.AtmSet[stail];
              fRepA.StereoParity := Mol.AtStereoParity(PA);
              InitMF(fRepA.AtomString[stail], mf);
              if (UseFormalCharge and Mol.IsChargedAt(stail)) then
              begin
                iimax := mf.Count - 1;
                for ii := 0 to iimax do
                  mf.Add(mf[ii] + '&FC' + fRepA.ConvertFormalCharge(
                    Mol.GetFormalCharge(stail)) + '&');
                for ii := iimax downto 0 do
                  mf.Delete(ii);
              end;
              if ((MarkAtom = 3) and Mol.IsMarkedAt(stail)) then
              begin
                iimax := mf.Count - 1;
                for ii := 0 to iimax do
                  mf.Add(mf[ii] + '&MA&');
                for ii := iimax downto 0 do
                  mf.Delete(ii);
              end;
              for jj := 0 to mf.Count - 1 do
              begin
                fw.Add(mf[jj] + IntToStr(fV[u]));
                bw.Add(mf[jj]);
              end;

              PA := Mol.AtmSet[ttail];
              fRepA.StereoParity := Mol.AtStereoParity(PA);
              InitMF(fRepA.AtomString[ttail], mf);
              if (UseFormalCharge and Mol.IsChargedAt(ttail)) then
              begin
                iimax := mf.Count - 1;
                for ii := 0 to iimax do
                  mf.Add(mf[ii] + '&FC' + fRepA.ConvertFormalCharge(
                    Mol.GetFormalCharge(ttail)) + '&');
                for ii := iimax downto 0 do
                  mf.Delete(ii);
              end;
              if ((MarkAtom = 3) and Mol.IsMarkedAt(ttail)) then
              begin
                iimax := mf.Count - 1;
                for ii := 0 to iimax do
                  mf.Add(mf[ii] + '&MA&');
                for ii := iimax downto 0 do
                  mf.Delete(ii);
              end;
              iimax := fw.Count - 1;
              for ii := 0 to iimax do
                for jj := 0 to mf.Count - 1 do
                begin
                  fw.Add(fw[ii] + mf[jj]);
                  bw.Add(mf[jj] + IntToStr(fV[u]) + bw[ii]);
                end;
              for ii := iimax downto 0 do
              begin
                fw.Delete(ii);
                bw.Delete(ii);
              end;

              for ii := 0 to fw.Count - 1 do
              begin
                if (fw[ii] <= bw[ii]) then
                  TSLFrg.AddObject(fw[ii], TObject(piw))
                else
                  TSLFrg.AddObject(bw[ii], TObject(piw));
              end;
            end;
          end;
        //end;
      end;
      FreeAndNil(fw);
      FreeAndNil(bw);
      FreeAndNil(mf);
    end;
    if (piwused = False) then
      dispose(piw);
  end;
end;

procedure TShrtPthAtmPair.PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList;
  PathLst: TList);
var
  piw: PInteger;
  j, k, iw, ii, jj, iimax, jjmax: integer;
  i, u, v, rv: Node;
  FORWD, BCKWD: string;
  THERE, BACK, Code: string;
  mfh, mft: TStringList;
  pu: PInteger;
  head, tail: AtomID;
begin
  for j := 0 to ColorKeys.Count - 1 do
  begin
    InitMF(ColorKeys[j],iw,TAtomBase(fRepA),Mol);
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
        tail:=fPath[u,1];
        head:=fPath[u,fL[u]];
        mft:=fColorPerAtom[tail] as TStringList;
        mfh:=fColorPerAtom[head] as TStringList;
        iimax := ffw.Count - 1;
        for ii := 0 to iimax do
          for jj := 0 to mft.Count - 1 do
          begin
            ffw.Add(ffw[ii] + mft[jj]+IntToStr(fV[u]));
            fbw.Add(fbw[ii] + IntToStr(fV[u])+mft[jj]);
          end;
        for ii := iimax downto 0 do
        begin
          ffw.Delete(ii);
          fbw.Delete(ii);
        end;
        iimax := ffw.Count - 1;
        for ii := 0 to iimax do
          for jj := 0 to mfh.Count - 1 do
          begin
            ffw.Add(ffw[ii] + mfh[jj]);
            fbw.Add(mfh[jj] + fbw[ii]);
          end;
        for ii := iimax downto 0 do
        begin
          ffw.Delete(ii);
          fbw.Delete(ii);
        end;
        for ii := 0 to ffw.Count - 1 do
        begin
          new(piw); piw^:=iw;
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

function TShrtPthAtmPair.PathToString(Mol: TMoleculeFrg): TStringList;
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

{function TShrtPthAtmPair.PathToString(Mol: TMoleculeFrg; s, t: AtomID;
  w: Cost): TStringList;
var
  i, u, v: AtomID;
  FORWD, BCKWD: string;
  mf, fw: TStringList;
  piw:    PInteger;
  j, jj, ii, iimax: integer;
begin
  Result := TStringList.Create;
  for j := 0 to ColorKeys.Count - 1 do
  begin
    if (ColorKeys[j] = 'Default') then
    begin //Valeur par defaut = symbol de l'atome
          //if (fRepA.UseSymbol) then
          //begin //Valeur par defaut = symbol de l'atome
      new(piw);
      piw^ := 1;
      for i := 1 to Mol.nAtom do
        fRepA.AtomString[i] := Mol.S_[i];//  AtmSet[i]^.S;
      //end;
      fRepA.StereoParity := Mol.AtStereoParity(s);
      FORWD := fRepA.AtomString[s];//Rep.GetSymbol(Mol.AtmSet[s]^.Z) + IntToStr(fV[t]);
      fRepA.StereoParity := Mol.AtStereoParity(t);
      FORWD := FORWD + IntToStr(w) + fRepA.AtomString[t];
      //Rep.GetSymbol(Mol.AtmSet[t]^.Z);

      //fRepA.StereoParity := Mol.AtStereoParity(t);
      //BCKWD := fRepA.AtomString[t];//Rep.GetSymbol(Mol.AtmSet[t]^.Z) + IntToStr(fV[t]);
      //fRepA.StereoParity := Mol.AtStereoParity(s);
      //BCKWD := BCKWD + IntToStr(w) + fRepA.AtomString[s];//Rep.GetSymbol(Mol.AtmSet[s]^.Z);

      BCKWD := ReverseString(FORWD);
      if FORWD <= BCKWD then
        Result.AddObject(FORWD, TObject(piw))
      else
        Result.AddObject(BCKWD, TObject(piw));
    end
    else
    begin
      mf := TStringList.Create;
      fw := TStringList.Create;
      new(piw);//assigne la mémoire! Ne pas oublier de la libérer!
      piw^ := fRepA.InitAtomStringWeight(ColorHash.Items[ColorKeys[j]]);
      fRepA.StereoParity := Mol.AtStereoParity(s);
      InitMF(fRepA.AtomString[s], mf);
      for jj := 0 to mf.Count - 1 do
        fw.Add(mf[jj]);
      fRepA.StereoParity := Mol.AtStereoParity(t);
      InitMF(fRepA.AtomString[t], mf);
      iimax := fw.Count - 1;
      for ii := 0 to iimax do
        for jj := 0 to mf.Count - 1 do
          fw.Add(fw[ii] + IntToStr(w) + mf[jj]);
      for ii := iimax downto 0 do
        fw.Delete(ii);
      for ii := 0 to fw.Count - 1 do
      begin
        BCKWD := ReverseString(fw[ii]);
        if (fw[ii] <= BCKWD) then
          Result.AddObject(fw[ii], TObject(piw))
        else
          Result.AddObject(BCKWD, TObject(piw));
      end;
      FreeAndNil(fw);
      FreeAndNil(mf);
    end;
  end;
end;}

constructor TShrtPthAtmPair.Create;
begin
  inherited Create;
  fRepA := TAtomSymbol.Create;
end;

constructor TShrtPthAtmPair.Create(s: TStringList);
begin
  inherited Create;
  fRepA := TAtomSymbol.Create;
  fFrgLst := s;
end;

destructor TShrtPthAtmPair.Destroy;
begin
  FreeAndNil(fRepA);
  inherited Destroy;
end;

procedure TShrtPthAtmPair.InitAtomString(sdline: string);
begin
  fRepA.InitAtomString(sdline);
end;

{procedure TShrtPthAtmPair.InitMF(mfstr: string; mf: TStringList);
begin
  mf.Clear;
  mf.StrictDelimiter := True;
  mf.Delimiter := '/';
  mf.DelimitedText := mfstr;
end;}

{procedure TShrtPthAtmPair.MolToFrgLst(Mol: TMoleculeFrg);
var
  i, j:   AtomID;
  s, t, u, v: Node;
  P:      TNodeInfo;
  SFrg:   TStringList;
  ii, iw: integer;
  AW:     TArcCost;
  k:      ArcNum;
  ptmpi:  Pinteger;
begin
  //Initialization
  for i := 1 to Mol.nAtom do
    P[i] := 0;  //array to contain a path
  for k := 1 to Mol.nBonds do
    AW[k] := 1; //array of atom weight - here all arcs have the same weight
  //Process each atom pair
  for s := 1 to Mol.nAtom do
  begin
    if ((not UseMarkAtom) or (UseMarkAtom and Mol.IsMarkedAt(s))) then
    begin
      for t := s + 1 to Mol.nAtom do
        if ((not UseMarkAtom) or (UseMarkAtom and Mol.IsMarkedAt(t))) then
        begin
          Mol.DijHeap(AW, s, t, fV, P);
          //Path ending at t has a cost containing in fV[t]
          //Here dynamic bonds restrictions should be tested.
          if (fV[t] + 1 >= LenMin) and (fV[t] + 1 <= LenMax) then
          begin
            SFrg := PathToString(Mol, s, t, fV[t]);
            for ii := 0 to SFrg.Count - 1 do
            begin
              iw := PInteger(SFrg.Objects[ii])^;
              AddFragment(SFrg[ii], iw);
            end;
            ptmpi := nil;
            for ii := 0 to SFrg.Count - 1 do
            begin
              if (ptmpi <> Pinteger(SFrg.Objects[ii])) then
              begin
                ptmpi := Pinteger(SFrg.Objects[ii]);
                dispose(ptmpi);
              end;
              SFrg.Objects[ii] := nil;
            end;
            FreeAndNil(SFrg);
          end;
        end;
    end;
  end;
end;}

end.

