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
unit UnitTriplet;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Math, strutils, unitfragmentbase, UnitMoleculeFrg,
  UnitMoleculeBase, unitatomsymbol, unitAtomAndBondType, U_TYPE_GRAPHES,
  unitsequences, UnitAtomBase; //UnitBondBase;

type

  { TTriplet }
  //ligne de test

  TTriplet = class(TSequences)
  private
    fRepA: TAtomSymbol;
  protected
    fDistM: CostMatrix; // Matrix containing the distances btw each atom pairs
    fNPath, fNP: Node;  //Number of discovered paths, Number of t terminated path
    fLen, fL: TNodeInfo;
    fV: TNodeCost;  //Cost, Length of each path, Lenght of t terminated path
    fPath: NodeMatrix; //Array of each path
  public
    constructor Create;
    constructor Create(s: TStringList);
    destructor Destroy; override;
    //procedure InitMF(mfstr: string; mf: TStringList);
    procedure MolToFrgLst(Mol: TMoleculeFrg); override;
  end;

implementation

{ TTriplet }

constructor TTriplet.Create;
begin
  inherited Create;
  fRepA := TAtomSymbol.Create;
end;

constructor TTriplet.Create(s: TStringList);
begin
  inherited Create;
  fRepA := TAtomSymbol.Create;
  fFrgLst := s;
end;

destructor TTriplet.Destroy;
begin
  FreeAndNil(fRepA);
  inherited Destroy;
end;

{procedure TTriplet.InitMF(mfstr: string; mf: TStringList);
begin
  mf.Clear;
  mf.StrictDelimiter := True;
  mf.Delimiter := '/';
  mf.DelimitedText := mfstr;
end;}

procedure TTriplet.MolToFrgLst(Mol: TMoleculeFrg);
var
  s, t: Node;
  i, iw, ia, ib, ic, j, a, b, c: integer;
  k: ArcNum;
  AW: TArcCost;
  P: TNodeInfo;
  tmp: TStringList;
  mfa, mfb, mfc: TStringList;
  FORWD: string;
  idx: Integer;
  pidx: PRAtmFrg;
begin
  tmp := TStringList.Create;

  //Initialisation
  for k := 1 to Mol.p_M do
    AW[k] := 1;
  for i := 1 to Mol.nAtom do
    P[i] := 0;
  for i := 1 to Mol.nAtom do
    for j := 1 to Mol.nAtom do
      fDistM[i, j] := 0;

  fColorPerAtom.Add(TStringList.Create);
  for s := 1 to Mol.nAtom do
  begin
    fColorPerAtom.Add(TStringList.Create);
  end;

  //Filling the CostMatrix
  for s := 1 to Mol.nAtom do
  begin
    t := 0;
    Mol.DijHeap(AW, s, t, fV, P);
    for t := s + 1 to Mol.nAtom do
    begin
      if (fV[t] + 1 >= LenMin) and (fV[t] + 1 <= LenMax) then
      begin
        fDistM[s, t] := fV[t];
        fDistM[t, s] := fV[t];
      end;
    end;
  end;

  //Initializing the coloration
  for j := 0 to ColorKeys.Count - 1 do
  begin
    InitMF(ColorKeys[j],iw,TAtomBase(fRepA),Mol);
    for a := 1 to Mol.nAtom - 2 do
    begin
      for b := a + 1 to Mol.nAtom - 1 do
      begin
        for c := b + 1 to Mol.nAtom do
        begin
          if ((MarkAtom = 0) or (MarkAtom = 3) or
            ((MarkAtom = 1) and (Mol.IsMarkedAt(a) or Mol.IsMarkedAt(b) or
            Mol.IsMarkedAt(c)))) or ((MarkAtom = 2) and (Mol.IsMarkedAt(a)) or
            (Mol.IsMarkedAt(b)) or (Mol.IsMarkedAt(c))) then
          begin
            if ((fDistM[a, c] <> 0) and (fDistM[b, c] <> 0) and (fDistM[a, b] <> 0)) then
              //test "flat" triplets here with: if ((fDistM[a, c]= fDistM[b, c] + fDistM[a, b]) or (fDistM[b, c] = fDistM[a, c] + fDistM[a, b]) or (fDistM[a, b] = fDistM[b, c] + fDistM[a, c]))
            begin
              //canonicalization
              mfa:=fColorPerAtom[a] as TStringList;
              mfb:=fColorPerAtom[b] as TStringList;
              mfc:=fColorPerAtom[c] as TStringList;
              for ia := 0 to mfa.Count - 1 do
              begin
                for ib := 0 to mfb.Count - 1 do
                begin
                  for ic := 0 to mfc.Count - 1 do
                  begin
                    //color and marked atom
                    tmp.Add(mfa[ia] + Format('%.3D', [fDistM[b, c]]));
                    tmp.Add(mfb[ib] + Format('%.3D', [fDistM[a, c]]));
                    tmp.Add(mfc[ic] + Format('%.3D', [fDistM[a, b]]));
                    tmp.Sort;
                    //Adding Fragment
                    idx:=AddFragment(tmp[0] + '-' + tmp[1] + '-' + tmp[2],iw);
                    if GetFrgPerAtom then
                    begin
                      new(pidx); pidx^.idx:=idx; pidx^.len:=3;
                      (fFrgPerAtom[a] as TList).Add(TObject(pidx));//Attention à la destruction
                      new(pidx); pidx^.idx:=idx; pidx^.len:=3;
                      (fFrgPerAtom[b] as TList).Add(TObject(pidx));//Attention à la destruction
                      new(pidx); pidx^.idx:=idx; pidx^.len:=3;
                      (fFrgPerAtom[c] as TList).Add(TObject(pidx));//Attention à la destruction
                    end;
                    //Adding the fragments
                    tmp.Clear;
                  end;
                end;
              end;
            end;
          end;
        end;
      end;
    end;
  end;
  FreeAndNil(tmp);
end;

end.
