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
unit UnitAllPathAtmBond;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, strutils, unitatomsymbol, UnitAllPath,
  unitAtomAndBondType, UnitAtomBase, UnitMoleculeBase, UnitMoleculeFrg,
  U_TYPE_GRAPHES; //UnitBondBase,

type

  { TAllPathAtmBond }

  TAllPathAtmBond = class(TAllPath)
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
  end;

implementation

{ TAllPathAtmBond }
procedure TAllPathAtmBond.PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList);
var
  u, v: Node;
  head, tail: AtomID;
  PA: PRAtom;
  PB: PRBond;
  FORWD, BCKWD: string;
  i, j, ii, jj, iimax: integer;
  piw: PInteger;
  fw, bw: TStringList;
  mf: TStringList;
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
      for u := 1 to fNPath do
      begin
        if fPath[u,0]=1 then
        {if (fLen[u] >= LenMin) then
        begin
           if (((DynBnd = 1) and (DynBondInFrg(fPath, u, fLen[u], Mol))) or
            ((DynBnd = 2) and (AllDynBondInFrg(fPath, u, fLen[u], Mol))) or
            (DynBnd = 0)) then}
          begin
            if ((MarkAtom = 0) or (MarkAtom = 1) or (MarkAtom = 3) or((MarkAtom = 2) and
              (MrkAtmInFrg(fPath, u, fLen[u], Mol)))) then
            begin
              piwused := True;
              FORWD := '';
              BCKWD := '';
              for v := 1 to fLen[u] - 1 do
              begin
                tail := fPath[u, v];
                head := fPath[u, v + 1];
                PB := Mol.FindBond(tail, head);
                PA := Mol.AtmSet[tail];
                fRepA.StereoParity := Mol.AtStereoParity(PA);
                //fRepB.StereoBond := Mol.BdStereo(PB);
                FORWD := FORWD + fRepA.AtomString[tail] + PB^.S;
                BCKWD := PB^.S + fRepA.AtomString[tail] + BCKWD;
              end;
              PA := Mol.AtmSet[fPath[u, fLen[u]]];
              fRepA.StereoParity := Mol.AtStereoParity(PA);
              FORWD := FORWD + fRepA.AtomString[fPath[u, fLen[u]]];
              BCKWD := fRepA.AtomString[fPath[u, fLen[u]]] + BCKWD;
              if FORWD <= BCKWD then
              begin
                TSLFrg.AddObject(FORWD, TObject(piw));
              end
              else
              begin
                TSLFrg.AddObject(BCKWD, TObject(piw));
              end;
            end;
          end;
        //end;
      end;
    end
    else
    begin
      fw := TStringList.Create;
      bw := TStringList.Create;
      mf := TStringList.Create;
      piw^ := fRepA.InitAtomStringWeight(ColorHash.Items[ColorKeys[j]]);
      for u := 1 to fNPath do
      begin
        if fPath[u,0]=1 then
        {if (fLen[u] >= LenMin) then
        begin
          if (((DynBnd = 1) and (DynBondInFrg(fPath, u, fLen[u], Mol))) or
            ((DynBnd = 2) and (AllDynBondInFrg(fPath, u, fLen[u], Mol))) or
            (DynBnd = 0)) then}
          begin
            if ((MarkAtom = 0) or (MarkAtom = 1) or (MarkAtom = 3) or((MarkAtom = 2) and
              (MrkAtmInFrg(fPath, u, fLen[u], Mol)))) then
            begin
              piwused := True;
              fw.Clear;
              bw.Clear;
              FORWD := '';
              BCKWD := '';
              fw.Add(FORWD);
              bw.add(BCKWD);
              for v := 1 to fLen[u] - 1 do
              begin
                tail := fPath[u, v];
                head := fPath[u, v + 1];
                PB := Mol.FindBond(tail, head);
                PA := Mol.AtmSet[tail];
                fRepA.StereoParity := Mol.AtStereoParity(PA);
                //fRepB.StereoBond := Mol.BdStereo(PB);
                //FORWD := FORWD + fRepA.AtomString[tail] + fRepB.GetBondSymbol(PB^.B);
                InitMF(fRepA.AtomString[tail], mf);
                if (UseFormalCharge and Mol.IsChargedAt(fPath[u, v])) then
                begin
                  iimax := mf.Count - 1;
                  for ii := 0 to iimax do
                    mf.Add(mf[ii] + '&FC' + fRepA.ConvertFormalCharge(
                      Mol.GetFormalCharge(fPath[u, v])) + '&');
                  for ii := iimax downto 0 do
                    mf.Delete(ii);
                end;
                if ((MarkAtom = 3) and Mol.IsMarkedAt(fPath[u, v])) then
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
                    fw.Add(fw[ii] + mf[jj] + PB^.S);
                    bw.Add(PB^.S + mf[jj] + bw[ii]);
                  end;
                for ii := iimax downto 0 do
                begin
                  fw.Delete(ii);
                  bw.Delete(ii);
                end;
              end;
              PA := Mol.AtmSet[fPath[u, fLen[u]]];
              fRepA.StereoParity := Mol.AtStereoParity(PA);
              InitMF(fRepA.AtomString[fPath[u, fLen[u]]], mf);
              if (UseFormalCharge and Mol.IsChargedAt(fPath[u, fLen[u]])) then
              begin
                iimax := mf.Count - 1;
                for ii := 0 to iimax do
                  mf.Add(mf[ii] + '&FC' + fRepA.ConvertFormalCharge(
                    Mol.GetFormalCharge(fPath[u, fLen[u]])) + '&');
                for ii := iimax downto 0 do
                  mf.Delete(ii);
              end;
              if ((MarkAtom = 3) and Mol.IsMarkedAt(fPath[u, fLen[u]])) then
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
                  bw.Add(mf[jj] + bw[ii]);
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

procedure TAllPathAtmBond.PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList;
  PathLst: TList);
var
  piw: PInteger;
  j, k, iw, ii, jj, iimax, jjmax: integer;
  i, u, v, rv, tail, head: Node;
  FORWD, BCKWD: string;
  THERE, BACK, Code: string;
  PB: PRBond;
  mf: TStringList;
  pu: PInteger;
begin
  for j := 0 to ColorKeys.Count - 1 do
  begin
    InitMF(ColorKeys[j],iw,TAtomBase(fRepA),Mol);
    //initialization of the representation and extraction of weight associated to micro-species
    for u := 1 to fNPath do
      if fPath[u,0]=1 then
      begin
        ffw.Clear;
        fbw.Clear;
        FORWD := '';
        BCKWD := '';
        ffw.Add(FORWD);
        fbw.Add(BCKWD);
        for v := 1 to fLen[u]-1 do
        begin
          tail := fPath[u, v];
          head := fPath[u, v + 1];
          PB := Mol.FindBond(tail, head);
          mf:=fColorPerAtom[fPath[u, v]] as TStringList;
          iimax := ffw.Count - 1;
          for ii := 0 to iimax do
            for jj := 0 to mf.Count - 1 do
              ffw.Add(ffw[ii] + mf[jj]+PB^.S);
          for ii := iimax downto 0 do
            ffw.Delete(ii);
          //Generation of backward sequence - needed to have canonical sequences
          tail := fPath[u, fLen[u] - v + 1];
          head := fPath[u, fLen[u] - v];
          PB := Mol.FindBond(tail, head);
          mf:=fColorPerAtom[fPath[u, fLen[u] - v + 1]] as TStringList;
          iimax := fbw.Count - 1;
          for ii := 0 to iimax do
            for jj := 0 to mf.Count - 1 do
              fbw.Add(fbw[ii] + mf[jj]+PB^.S);
          for ii := iimax downto 0 do
            fbw.Delete(ii);
        end;
        //close forward and backward strings
        tail:=fPath[u, fLen[u]];
        mf:=fColorPerAtom[tail] as TStringList;
        iimax := ffw.Count - 1;
        for ii := 0 to iimax do
          for jj := 0 to mf.Count - 1 do
            ffw.Add(ffw[ii] + mf[jj]);
        for ii := iimax downto 0 do
          ffw.Delete(ii);
        tail:=fPath[u, 1];
        mf:=fColorPerAtom[tail] as TStringList;
        iimax := fbw.Count - 1;
        for ii := 0 to iimax do
          for jj := 0 to mf.Count - 1 do
            fbw.Add(fbw[ii] + mf[jj]);
        for ii := iimax downto 0 do
          fbw.Delete(ii);
        //Count fragment
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

function TAllPathAtmBond.PathToString(Mol: TMoleculeFrg): TStringList;
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

constructor TAllPathAtmBond.Create;
begin
  inherited Create;
  fRepA := TAtomSymbol.Create;
  //fRepB := TBondBase.Create;
end;

constructor TAllPathAtmBond.Create(s: TStringList);
begin
  inherited Create;
  fRepA := TAtomSymbol.Create;
  //fRepB   := TBondBase.Create;
  fFrgLst := s;
end;

destructor TAllPathAtmBond.Destroy;
begin
  FreeAndNil(fRepA);
  //FreeAndNil(fRepB);
  inherited Destroy;
end;

end.
