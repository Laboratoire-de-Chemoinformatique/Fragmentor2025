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
unit UnitShortestPathAtmBndPair;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, strutils, UnitSequences, UnitShortestPath, unitatomsymbol,
  unitfragmentbase, UnitMoleculeFrg, U_TYPE_GRAPHES, UnitAtomBase,
  UnitMoleculeBase, unitAtomAndBondType; //UnitBondBase,

type

  { TShrtPthAtmPair }

  { TShrtPthAtmBndPair }

  TShrtPthAtmBndPair = class(TShortestPath)
  private
    fRepA: TAtomSymbol;
    //fRepB: TBondBase;
  protected
    function PathToString(Mol: TMoleculeFrg): TStringList; override;
    procedure PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList); override;
    procedure PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList; PathLst: TList); override;
  public
    constructor Create;
    constructor Create(s: TStringList);
    destructor Destroy; override;
    //procedure InitMF(mfstr: string; mf: TStringList);
  end;

implementation

{ TShrtPthAtmBndPair }
procedure TShrtPthAtmBndPair.PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList);
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
      //Initialisation
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
              shead := fPath[u, 2];
              ttail := fPath[u, fL[u]];
              thead := fPath[u, fL[u] - 1];

              PB := Mol.FindBond(stail, shead);
              PA := Mol.AtmSet[stail];
              //fRepA.StereoParity := Mol.AtStereoParity(PA);
              //fRepB.StereoBond := Mol.BdStereo(PB);
              FORWD := fRepA.AtomString[stail] + PB^.S + IntToStr(fV[u]);
              BCKWD := IntToStr(fV[u]) + PB^.S + fRepA.AtomString[stail];
              PB := Mol.FindBond(thead, ttail);
              PA := Mol.AtmSet[ttail];
              //fRepA.StereoParity := Mol.AtStereoParity(PA);
              //fRepB.StereoBond := Mol.BdStereo(PB);
              FORWD := FORWD + PB^.S + fRepA.AtomString[ttail];
              BCKWD := fRepA.AtomString[ttail] + PB^.S + BCKWD;

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
            if ((MarkAtom = 0) or (MarkAtom = 1) or (MarkAtom = 3) or ((MarkAtom = 2) and
              (MrkAtmInFrg(fPath, u, fL[u], Mol)))) then
            begin
              piwused := True;
              fw.Clear;
              bw.Clear;

              stail := fPath[u, 1];
              shead := fPath[u, 2];
              ttail := fPath[u, fL[u]];
              thead := fPath[u, fL[u] - 1];

              PB := Mol.FindBond(stail, shead);
              PA := Mol.AtmSet[stail];
              //fRepA.StereoParity := Mol.AtStereoParity(PA);
              //fRepB.StereoBond := Mol.BdStereo(PB);
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
                fw.Add(mf[jj] + PB^.S + IntToStr(fV[u]));
                bw.Add(IntToStr(fV[u]) + PB^.S + mf[jj]);
              end;
              PB := Mol.FindBond(thead, ttail);
              PA := Mol.AtmSet[ttail];
              //fRepA.StereoParity := Mol.AtStereoParity(PA);
              //fRepB.StereoBond := Mol.BdStereo(PB);
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
                  fw.Add(fw[ii] + PB^.S + mf[jj]);
                  bw.Add(mf[jj] + PB^.S + bw[ii]);
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

procedure TShrtPthAtmBndPair.PathToString(Mol: TMoleculeFrg;
  TSLFrg: TStringList; PathLst: TList);
var
  piw: PInteger;
  j, k, iw, ii, jj, iimax, jjmax: integer;
  i, u, v, rv: Node;
  FORWD, BCKWD: string;
  THERE, BACK, Code: string;
  mfh, mft: TStringList;
  pu: PInteger;
  shead,thead, stail,ttail: AtomID;
  PB1,PB2: PRBond;
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
        stail:=fPath[u,1];
        ttail:=fPath[u,2];
        PB1 := Mol.FindBond(stail, ttail);
        mft:=fColorPerAtom[stail] as TStringList;
        shead:=fPath[u,fL[u]];
        thead:=fPath[u,fL[u]-1];
        PB2 := Mol.FindBond(shead, thead);
        mfh:=fColorPerAtom[shead] as TStringList;
        iimax := ffw.Count - 1;
        for ii := 0 to iimax do
          for jj := 0 to mft.Count - 1 do
          begin
            ffw.Add(ffw[ii] + mft[jj]+PB1^.S+IntToStr(fV[u]));
            fbw.Add(fbw[ii] + IntToStr(fV[u])+PB1^.S+mft[jj]);
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
            ffw.Add(ffw[ii] + PB2^.S + mfh[jj]);
            fbw.Add(mfh[jj] + PB2^.S + fbw[ii]);
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


function TShrtPthAtmBndPair.PathToString(Mol: TMoleculeFrg): TStringList;
begin
  Result := TStringList.Create;
  PathToString(Mol, Result);
end;

constructor TShrtPthAtmBndPair.Create;
begin
  inherited Create;
  fRepA := TAtomSymbol.Create;
  //fRepB := TBondBase.Create;
end;

constructor TShrtPthAtmBndPair.Create(s: TStringList);
begin
  inherited Create;
  fRepA := TAtomSymbol.Create;
  //fRepB   := TBondBase.Create;
  fFrgLst := s;
end;

destructor TShrtPthAtmBndPair.Destroy;
begin
  FreeAndNil(fRepA);
  //FreeAndNil(fRepB);
  inherited Destroy;
end;

{procedure TShrtPthAtmBndPair.InitMF(mfstr: string; mf: TStringList);
begin
  mf.Clear;
  mf.StrictDelimiter := True;
  mf.Delimiter := '/';
  mf.DelimitedText := mfstr;
end;}


end.

