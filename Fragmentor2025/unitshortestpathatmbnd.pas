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
unit UnitShortestPathAtmBnd;

{$mode objfpc}{$H+}

interface

uses
//  Classes, SysUtils, strutils, unitfragmentbase, unitsequences, unitshortestpath,
//  UnitMoleculeFrg, unitAtomAndBondType, unitatomsymbol, UnitAtomBase,
//  U_TYPE_GRAPHES, UnitMoleculeBase; //UnitBondBase,
  Classes, SysUtils, strutils, unitfragmentbase, UnitMoleculeBase, UnitMoleculeFrg, UnitShortestPath,
  unitatombase, UnitAtomSymbol, unitAtomAndBondType, U_TYPE_GRAPHES;

type

  { TShrtPthAtmBnd }

  TShrtPthAtmBnd = class(TShortestPath)
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
//    procedure InitMF(mfstr: string; mf: TStringList);
  end;

implementation

{ TShrtPthAtmBnd }
procedure TShrtPthAtmBnd.PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList);
var
  u, v: Node;
  head, tail: AtomID;
  PA: PRAtom;
  PB: PRBond;
  FORWD, BCKWD: string;
  THERE, BACK, Code: string;
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
      //Initialization of weight and representation
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
      //Generating sequence from path
      for u := 1 to fNP do
      begin
        if fPath[u,0]=1 then
        {if (fL[u] >= LenMin) then
        begin
          if (((DynBnd = 1) and (DynBondInFrg(fPath, fNPath + u, fL[u], Mol))) or
            ((DynBnd = 2) and (AllDynBondInFrg(fPath, fNPath + u, fL[u], Mol))) or
            (DynBnd = 0)) then}
          begin
            if ((MarkAtom = 0) or (MarkAtom = 1) or (MarkAtom = 3) or((MarkAtom = 2) and (MrkAtmInFrg(fPath, u, fL[u], Mol)))) then
            begin
              piwused := True;
              FORWD := '';
              BCKWD := '';
              THERE := '';
              BACK := '';
              for v := 1 to fL[u] - 1 do
              begin
                tail := fPath[u, v];
                head := fPath[u, v + 1];
                PB := Mol.FindBond(tail, head);
                PA := Mol.AtmSet[tail];
                fRepA.StereoParity := Mol.AtStereoParity(PA);
                //fRepB.StereoBond := Mol.BdStereo(PB);
                //FORWD := FORWD + fRepA.AtomString[tail] + fRepB.GetBondSymbol(PB^.B);
                FORWD := FORWD + fRepA.AtomString[tail] + PB^.S;

                //Generating backward sequence for canonical sequences
                tail := fPath[u, fL[u] - v + 1];
                head := fPath[u, fL[u] - v];
                PB := Mol.FindBond(tail, head);
                PA := Mol.AtmSet[tail];
                fRepA.StereoParity := Mol.AtStereoParity(PA);
                //fRepB.StereoBond := Mol.BdStereo(PB);
                //BCKWD:=BCKWD+RepA.GetSymbol(PA^.Z)+RepB.GetBondSymbol(PB^.B);
                //BCKWD := BCKWD + fRepA.AtomString[tail] + fRepB.GetBondSymbol(PB^.B);
                BCKWD := BCKWD + fRepA.AtomString[tail] + PB^.S;

                Str(Mol.Z_[fPath[u, v]], Code);
                THERE := THERE + Code;
                Str(Mol.Z_[fPath[u, fL[u] - v + 1]], Code);
                BACK := BACK + Code;
              end;
              PA := Mol.AtmSet[fPath[u, fL[u]]];
              fRepA.StereoParity := Mol.AtStereoParity(PA);
              FORWD := FORWD + fRepA.AtomString[fPath[u, fL[u]]];

              PA := Mol.AtmSet[fPath[u, 1]];
              fRepA.StereoParity := Mol.AtStereoParity(PA);
              //BCKWD:=BCKWD+RepA.GetSymbol(PA^.Z);
              BCKWD := BCKWD + fRepA.AtomString[fPath[u, 1]];

              Str(Mol.Z_[fPath[u, fL[u]]], Code);
              THERE := THERE + Code;
              Str(Mol.Z_[fPath[u, 1]], Code);
              BACK := BACK + Code;
              if THERE <= BACK then
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
      fw := TStringList.Create;
      bw := TStringList.Create;
      mf := TStringList.Create;
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
            if ((MarkAtom = 0) or (MarkAtom = 1) or (MarkAtom = 3) or ((MarkAtom = 2) and (MrkAtmInFrg(fPath, u, fL[u], Mol)))) then
            begin
              piwused := True;
              fw.Clear;
              bw.Clear;
              FORWD := '';
              BCKWD := '';
              fw.Add(FORWD);
              bw.add(BCKWD);
              for v := 1 to fL[u] - 1 do
              begin
                tail := fPath[u, v];
                head := fPath[u, v + 1];
                PB := Mol.FindBond(tail, head);
                PA := Mol.AtmSet[tail];
                fRepA.StereoParity := Mol.AtStereoParity(PA);
                //fRepB.StereoBond := Mol.BdStereo(PB);
                //FORWD := FORWD + fRepA.AtomString[tail] + fRepB.GetBondSymbol(PB^.B);
                InitMF(fRepA.AtomString[tail], mf);
                if (UseFormalCharge and Mol.IsChargedAt(tail)) then
                begin
                  iimax := mf.Count - 1;
                  for ii := 0 to iimax do
                    mf.Add(mf[ii] + '&FC' + fRepA.ConvertFormalCharge(
                      Mol.GetFormalCharge(tail)) + '&');
                  for ii := iimax downto 0 do
                    mf.Delete(ii);
                end;
                if ((MarkAtom = 3) and Mol.IsMarkedAt(tail)) then
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
                    fw.Add(fw[ii] + mf[jj] + PB^.S);
                for ii := iimax downto 0 do
                  fw.Delete(ii);

                //Generating backward sequence
                tail := fPath[u, fL[u] - v + 1];
                head := fPath[u, fL[u] - v];
                PB := Mol.FindBond(tail, head);
                PA := Mol.AtmSet[tail];
                fRepA.StereoParity := Mol.AtStereoParity(PA);
                //fRepB.StereoBond := Mol.BdStereo(PB);
                InitMF(fRepA.AtomString[tail], mf);
                if (UseFormalCharge and Mol.IsChargedAt(tail)) then
                begin
                  iimax := mf.Count - 1;
                  for ii := 0 to iimax do
                    mf.Add(mf[ii] + '&FC' + fRepA.ConvertFormalCharge(
                      Mol.GetFormalCharge(tail)) + '&');
                  for ii := iimax downto 0 do
                    mf.Delete(ii);
                end;
                if ((MarkAtom = 3) and Mol.IsMarkedAt(tail)) then
                begin
                  iimax := mf.Count - 1;
                  for ii := 0 to iimax do
                    mf.Add(mf[ii] + '&MA&');
                  for ii := iimax downto 0 do
                    mf.Delete(ii);
                end;
                iimax := bw.Count - 1;
                for ii := 0 to iimax do
                  for jj := 0 to mf.Count - 1 do
                    bw.Add(bw[ii] + mf[jj] + PB^.S);
                for ii := iimax downto 0 do
                  bw.Delete(ii);
              end;
              PA := Mol.AtmSet[fPath[u, fL[u]]];
              fRepA.StereoParity := Mol.AtStereoParity(PA);
              InitMF(fRepA.AtomString[fPath[u, fL[u]]], mf);
              if (UseFormalCharge and Mol.IsChargedAt(tail)) then
              begin
                iimax := mf.Count - 1;
                for ii := 0 to iimax do
                  mf.Add(mf[ii] + '&FC' + fRepA.ConvertFormalCharge(
                    Mol.GetFormalCharge(tail)) + '&');
                for ii := iimax downto 0 do
                  mf.Delete(ii);
              end;
              if ((MarkAtom = 3) and Mol.IsMarkedAt(tail)) then
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
                  fw.Add(fw[ii] + mf[jj]);
              for ii := iimax downto 0 do
                fw.Delete(ii);
              //FORWD := FORWD + fRepA.AtomString[fPath[u + fNPath, fL[u]]];

              PA := Mol.AtmSet[fPath[u, 1]];
              fRepA.StereoParity := Mol.AtStereoParity(PA);
              //BCKWD:=BCKWD+fRepA.AtomString[fPath[u + fNPath, 1]];
              InitMF(fRepA.AtomString[fPath[u, 1]], mf);
              if (UseFormalCharge and Mol.IsChargedAt(PA)) then
              begin
                iimax := mf.Count - 1;
                for ii := 0 to iimax do
                  mf.Add(mf[ii] + '&FC' + fRepA.ConvertFormalCharge(
                    Mol.GetFormalCharge(PA)) + '&');
                for ii := iimax downto 0 do
                  mf.Delete(ii);
              end;
              if ((MarkAtom = 3) and Mol.IsMarkedAt(PA)) then
              begin
                iimax := mf.Count - 1;
                for ii := 0 to iimax do
                  mf.Add(mf[ii] + '&MA&');
                for ii := iimax downto 0 do
                  mf.Delete(ii);
              end;
              iimax := bw.Count - 1;
              for ii := 0 to iimax do
                for jj := 0 to mf.Count - 1 do
                  bw.Add(bw[ii] + mf[jj]);
              for ii := iimax downto 0 do
                bw.Delete(ii);

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

procedure TShrtPthAtmBnd.PathToString(Mol: TMoleculeFrg; TSLFrg: TStringList;
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
  s: string;
begin
  for j := 0 to ColorKeys.Count - 1 do
  begin
    s:=ColorKeys[j];
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
        for v := 1 to fL[u]-1 do
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
          tail := fPath[u, fL[u] - v + 1];
          head := fPath[u, fL[u] - v];
          PB := Mol.FindBond(tail, head);
          mf:=fColorPerAtom[fPath[u, fL[u] - v + 1]] as TStringList;
          iimax := fbw.Count - 1;
          for ii := 0 to iimax do
            for jj := 0 to mf.Count - 1 do
              fbw.Add(fbw[ii] + mf[jj]+PB^.S);
          for ii := iimax downto 0 do
            fbw.Delete(ii);
        end;
        //close forward and backward strings
        tail:=fPath[u, fL[u]];
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
  //for i:=0 to TSLFrg.Count-1 do writeln(TSLFrg.Text);
end;

function TShrtPthAtmBnd.PathToString(Mol: TMoleculeFrg): TStringList;
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

constructor TShrtPthAtmBnd.Create;
begin
  inherited Create;
  fRepA := TAtomSymbol.Create;
  //fRepB := TBondBase.Create;
end;

constructor TShrtPthAtmBnd.Create(s: TStringList);
begin
  inherited Create(s);
  fRepA := TAtomSymbol.Create;
  //fRepB   := TBondBase.Create;
  //fFrgLst := s;
end;

destructor TShrtPthAtmBnd.Destroy;
begin
  FreeAndNil(fRepA);
  //FreeAndNil(fRepB);
  inherited Destroy;
end;

{procedure TShrtPthAtmBnd.InitMF(mfstr: string; mf: TStringList);
begin
  mf.Clear;
  mf.StrictDelimiter := True;
  mf.Delimiter := '/';
  mf.DelimitedText := mfstr;
end;}

end.
