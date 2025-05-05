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
unit UnitAtmCnt;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, UnitMoleculeFrg, UnitFragmentBase, UnitAtomBase,
  UnitAtomSymbol, UnitMolecule;

type

  { TAtmCnt }

  TAtmCnt = class(TFrgBase)
  private
    fRepA: TAtomSymbol;
  public
    constructor Create;
    constructor Create(s: TStringList);
    constructor Create(s: TStringList; sdline: string);
    destructor Destroy; override;
    procedure InitAtomString(sdline: string); override;
    // and mf for multi-flagging as an atom may have several properties associated to it
    procedure MolToFrgLst(Mol: TMoleculeFrg); override;
  end;

implementation

{ TAtmCnt }

constructor TAtmCnt.Create;
begin
  inherited Create;
  fRepA := TAtomSymbol.Create;
end;

constructor TAtmCnt.Create(s: TStringList);
begin
  inherited Create(s);
  fRepA   := TAtomSymbol.Create;
end;

constructor TAtmCnt.Create(s: TStringList; sdline: string);
begin
  inherited Create(s);
  fRepA := TAtomSymbol.Create;
  InitAtomString(sdline);
end;

destructor TAtmCnt.Destroy;
begin
  inherited Destroy;
  FreeAndNil(fRepA);
end;

procedure TAtmCnt.InitAtomString(sdline: string);
begin
  fRepA.InitAtomString(sdline);
end;

procedure TAtmCnt.MolToFrgLst(Mol: TMoleculeFrg);
var
  i:   AtomID;
  iw, j, l: integer;
  mf: TStringList;
  pidx: PRAtmFrg;
begin
  iw:=1;
  //Prepare for repository of atom labels
  fColorPerAtom.Add(TStringList.Create);
  for i := 1 to Mol.nAtom do
    fColorPerAtom.Add(TStringList.Create);
  for j := 0 to ColorKeys.Count - 1 do
  begin
    InitMF(ColorKeys[j],iw,TAtomBase(fRepA),Mol);
    for i := 1 to Mol.nAtom do
    begin
      mf:=fColorPerAtom[i] as TStringList;
      for l := 0 to mf.Count - 1 do
      begin
        new(pidx); pidx^.idx:=AddFragment(mf[l], iw); pidx^.len:=1;
        if fGetFrgPerAtom then
          (fFrgPerAtom[i] as TList).Add(pidx)
        else
            dispose(pidx);
      end;
    end;
  end;
end;

end.

