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
unit UnitFragment;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, contnrs, UnitMoleculeFrg, UnitFragmentBase, UnitAtmCnt,
  UnitAllPathAtom, UnitAllPathBond, UnitAllPathAtmBond, UnitShortestPathAtom,
  UnitShortestPathBond, UnitShortestPathAtmbnd, UnitShortestPathAtmPair,
  UnitShortestPathAtmBndPair, UnitACPathAtom, UnitACPathBnd,
  UnitACPathAtmBnd, UnitACAtmPair, UnitACAtmBndPair, UnitACFXPathAtom,
  UnitACFXPathBnd, UnitACFXPathAtmBnd, UnitACFXAtmPair, UnitACFXAtmBndPair,
  UnitTriplet;

type
  PDefFrg = ^DefFrg;

  DefFrg = record
    FrgType, lmin, lmax, MarkAtom, DynBnd: byte;
    DoAllWays, EqFuzAt, UseFormalCharge, UseRadical, UseIsotope, AtomPairs,
    StrictFrg, GetAtomFrg, Cycle: boolean;
    ColorAFields, ColorBFields: string; //SDF fields for atom and bond coloring
  end;

  { LDefFrg - Store informations about fragmentors}

  LDefFrg = class(TList)
  private
    procedure SetDefFrg(Index: integer; Value: DefFrg);
    function GetDefFrg(Index: integer): DefFrg;
  public
    constructor Create;
    destructor Destroy; override;
    procedure Clear; override;
    property Objects[Index: integer]: DefFrg Read GetDefFrg Write SetDefFrg;
  end;

  TUFragmentException = class(Exception)
  end;

  { TFragment }

  TFragment = class(TObject)
  private
    fTotFrgList: TStringList;
    fFrgmnts: TObjectList;
    fMol: TMoleculeFrg;
    fFrgPerAtom: TObjectList;//List of pointers to fragment each atom belong to.
    fGetFrgPerAtom: Boolean;
    fbCycle: Boolean;
    procedure SetfGetFrgPerAtom(Value: boolean);
  public
    constructor Create; virtual;
    destructor Destroy; override;
    property TotFrgList: TStringList Read fTotFrgList Write fTotFrgList;
    property FrgPerAtom: TObjectList Read fFrgPerAtom;
    property GetFrgPerAtom: Boolean read fGetFrgPerAtom write SetfGetFrgPerAtom;
    property bCycle: Boolean read fbCycle write fbCycle;
    property Frgmnts: TObjectList read fFrgmnts;
    procedure Clear;
    procedure FrgReset;
    procedure SetMol(MolList: TStringList);
    procedure AddFragmentor(DFrg: DefFrg);
    procedure AddFragmentor(PDFrg: PDefFrg);
    procedure AddFragmentor(ftype, LenMinValue, LenMaxValue, MarkAtomValue, DynBnd: byte;
      EqFuzAt, UseFormalCharge, DoAllWays, AtomPairs, StrictFrg: boolean);
    procedure MolToFrgLst(MolList: TStringList); virtual;
    procedure SetFrgLst(fraglist: TStringList);
    function GetFrgLst: TStringList;
    function IsAnyNewFrg: boolean;
  end;

procedure mol_list_to_frag(MolList: TStringList; var fraglist: TStringList;
  ftype, LenMinValue, LenMaxValue, MarkAtomValue, DynBnd: byte;
  EqFuzAt, UseFormalCharge, DoAllWays, AtomPairs, StrictFrg: boolean);

procedure InitDefFrg(pf: PDefFrg);//Initialize a DefFrg with default values
procedure InitDefFrg(out af: DefFrg);//Initialize a DefFrg with default values

implementation

procedure mol_list_to_frag(MolList: TStringList; var fraglist: TStringList;
  ftype, LenMinValue, LenMaxValue, MarkAtomValue, DynBnd: byte;
  EqFuzAt, UseFormalCharge, DoAllWays, AtomPairs, StrictFrg: boolean);
var
  i:     integer;
  frag:  TFragment;
  stmp:  TStringList;
  pitem: PInteger;
begin
  frag := TFragment.Create;
  stmp := TStringList.Create;
  frag.AddFragmentor(ftype, LenMinValue, LenMaxValue, MarkAtomValue,
    DynBnd, EqFuzAt, UseFormalCharge, DoAllWays, AtomPairs, StrictFrg);
  frag.MolToFrgLst(MolList);
  stmp := frag.GetFrgLst;
  fraglist.Clear;
  for i := 0 to stmp.Count - 1 do
  begin
    new(pitem);
    pitem^ := PInteger(stmp.Objects[i])^;
    fraglist.AddObject(stmp.Strings[i], TObject(pitem));
  end;
  stmp.Free;
  frag.Free;
end;

procedure InitDefFrg(pf: PDefFrg);
begin
  InitDefFrg(pf^);
end;

procedure InitDefFrg(out af: DefFrg);
begin
  af.FrgType:=0;
  af.lmin:=0;
  af.lmax:=0;
  af.MarkAtom:=0;
  af.DynBnd:=0;
  af.DoAllWays:=False;
  af.EqFuzAt:=False;
  af.UseFormalCharge:=False;
  af.UseRadical:=False;
  af.UseIsotope:=False;
  af.AtomPairs:=False;
  af.StrictFrg:=False;
  af.GetAtomFrg:=False;
  af.Cycle:=False;
  af.ColorAFields:='Default';
  af.ColorBFields:='Default';
end;

{ LDefFrg }

procedure LDefFrg.SetDefFrg(Index: integer; Value: DefFrg);
var
  pitem: PDefFrg;
begin
  pitem := Items[Index];
  dispose(pitem);
  new(pitem);
  pitem^ := Value;
  Items[Index] := pitem;
end;

function LDefFrg.GetDefFrg(Index: integer): DefFrg;
begin
  Result := PDefFrg(Items[Index])^;
end;

constructor LDefFrg.Create;
begin
  inherited Create;
end;

destructor LDefFrg.Destroy;
begin
  Clear;
  inherited Destroy;
end;

procedure LDefFrg.Clear;
var
  i: integer;
begin
  for i := 0 to Count - 1 do
    dispose(PDefFrg(Items[i]));
  inherited Clear;
end;

{ TFragment }

procedure TFragment.SetfGetFrgPerAtom(Value: boolean);
//If the fragment ID to which belong an atom is always computed for all fragmentations
var
  i: integer;
begin
  if Value=True then
  begin
    for i:=0 to fFrgmnts.Count-1 do
        (fFrgmnts[i] as TFrgBase).GetFrgPerAtom:=True;
  end else
  begin
    i:=0;
    while (i<fFrgmnts.Count) and (Value=False) do
    begin
      Value:=(fFrgmnts[i] as TFrgBase).GetFrgPerAtom;
      inc(i);
    end;
  end;
  fGetFrgPerAtom:=Value;
end;

constructor TFragment.Create;
begin
  inherited Create;
  fTotFrgList := TStringList.Create;
  fTotFrgList.CaseSensitive:=True;
  fFrgmnts := TObjectList.Create;
  fFrgmnts.OwnsObjects := True;
  fMol := TMoleculeFrg.Create;
  fFrgPerAtom:=TObjectList.create;
  fFrgPerAtom.OwnsObjects:=True;
  fGetFrgPerAtom:=False;
  fbCycle:=False;
end;

destructor TFragment.Destroy;
var
  i,j: integer;
  pidx: PInteger;
  tli: TList;
begin

  Clear;
  FreeAndNil(fFrgmnts);
  FreeAndNil(fTotFrgList);
  FreeAndNil(fMol);
  FreeAndNil(fFrgPerAtom);
  inherited Destroy;
end;

procedure TFragment.Clear;
var
  i,j: integer;
  cnt: PInteger;
  pidx: PRAtmFrg;
  tli: TList;
begin
  for i := 0 to fTotFrgList.Count - 1 do
  begin
    if (PInteger(fTotFrgList.Objects[i]) <> nil) then
      dispose(PInteger(fTotFrgList.Objects[i]));
  end;
  fTotFrgList.Clear;
  fTotFrgList.CaseSensitive:=True;
  for i:=0 to fFrgmnts.Count-1 do
      (fFrgmnts.Items[i] as TFrgBase).NewFrg:=False;
  for i:=0 to fFrgPerAtom.Count-1 do
  begin
    tli:=fFrgPerAtom[i] as TList;
    for j:=0 to tli.Count-1 do
    begin
      if (tli.Items[j]<>nil) then
      begin
        pidx:=PRAtmFrg(tli.Items[j]);
        dispose(pidx);
      end;
      tli.Items[j]:=TObject(nil);
    end;
  end;
  fFrgPerAtom.Clear;
  fGetFrgPerAtom:=False;
  fbCycle:=False;
  fMol.Clear;
end;

procedure TFragment.FrgReset;
var
  i,j: integer;
  cnt: PInteger;
  pidx: PRAtmFrg;
  tli: TList;
begin
  //Reset fTotFrgList
  for i := 0 to fTotFrgList.Count - 1 do
  begin
    if (PInteger(fTotFrgList.Objects[i]) <> nil) then
      PInteger(fTotFrgList.Objects[i])^ := 0
    else
    begin
      new(cnt);
      cnt^ := 0;
      fTotFrgList.Objects[i] := TObject(cnt);
    end;
  end;
  for i:=0 to fFrgmnts.Count-1 do
      (fFrgmnts.Items[i] as TFrgBase).NewFrg:=False;

  //Reset fFrgPerAtom
  for i:=0 to fFrgPerAtom.Count-1 do
  begin
    tli:=fFrgPerAtom[i] as TList;
    for j:=0 to tli.Count-1 do
      if (tli.Items[j]<>nil) then
      begin
        pidx:=PRAtmFrg(tli.Items[j]);
        dispose(pidx);
        pidx:=nil;
        tli.Items[j]:=TObject(nil);
      end;
  end;
  fFrgPerAtom.Clear;
  fMol.Clear;
end;

procedure TFragment.SetMol(MolList: TStringList);
var
  i: integer;
begin
  fMol.Clear;
  fMol.LoadSDF(MolList);
  if fGetFrgPerAtom then
    for i:=0 to fMol.nAtom do
      fFrgPerAtom.Add(TList.Create);
end;

procedure TFragment.AddFragmentor(DFrg: DefFrg);
var
  i: integer;
begin
  with DFrg do
  begin
    //if (UseBenson) then
    //  raise TUFragmentException.Create('ERROR: Deprecated use of Benson option');

    //writeln('Adding a fragmentor **************************');

    if (FrgType = 0) then
    begin
      i := fFrgmnts.Add(TAtmCnt.Create(fTotFrgList));
      (fFrgmnts.Items[i] as TAtmCnt).MarkAtom := MarkAtom;
      (fFrgmnts.Items[i] as TAtmCnt).UseFormalCharge := UseFormalCharge;
      (fFrgmnts.Items[i] as TAtmCnt).UseRadical := UseRadical;
      (fFrgmnts.Items[i] as TAtmCnt).UseIsotope := UseIsotope;
    end
    else if (FrgType = 1) then
    begin
      if AtomPairs and DoAllWays then
        raise TUFragmentException.Create(
          'ERROR: Unsupported combination of arguments')
      else if AtomPairs then
      begin
        i := fFrgmnts.Add(TShrtPthAtmPair.Create(TotFrgList));
        (fFrgmnts.Items[i] as TShrtPthAtmPair).LenMin := lmin;
        (fFrgmnts.Items[i] as TShrtPthAtmPair).LenMax := lmax;
        (fFrgmnts.Items[i] as TShrtPthAtmPair).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TShrtPthAtmPair).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TShrtPthAtmPair).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TShrtPthAtmPair).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TShrtPthAtmPair).DynBnd := DynBnd;
        (fFrgmnts.Items[i] as TShrtPthAtmPair).IsPair := True;
      end
      else if DoAllWays then
      begin
        i := fFrgmnts.Add(TAllPathAtm.Create(TotFrgList));
        (fFrgmnts.Items[i] as TAllPathAtm).LenMin := lmin;
        (fFrgmnts.Items[i] as TAllPathAtm).LenMax := lmax;
        (fFrgmnts.Items[i] as TAllPathAtm).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TAllPathAtm).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TAllPathAtm).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TAllPathAtm).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TAllPathAtm).DynBnd := DynBnd;
      end
      else
      begin
        i := fFrgmnts.Add(TShrtPthAtm.Create(TotFrgList));
        (fFrgmnts.Items[i] as TShrtPthAtm).LenMin := lmin;
        (fFrgmnts.Items[i] as TShrtPthAtm).LenMax := lmax;
        (fFrgmnts.Items[i] as TShrtPthAtm).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TShrtPthAtm).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TShrtPthAtm).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TShrtPthAtm).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TShrtPthAtm).DynBnd := DynBnd;
      end;
    end
    else if (FrgType = 2) then
    begin
      if AtomPairs then
        raise TUFragmentException.Create(
          'ERROR: Unsupported combination of arguments')
      else if DoAllWays then
      begin
        i := fFrgmnts.Add(TAllPathBond.Create(TotFrgList));
        (fFrgmnts.Items[i] as TAllPathBond).LenMin := lmin;
        (fFrgmnts.Items[i] as TAllPathBond).LenMax := lmax;
        (fFrgmnts.Items[i] as TAllPathBond).MarkAtom := MarkAtom; //Atoms are not represented for this fragmentation. Is it usedful to set MarkAtom, etc?
        (fFrgmnts.Items[i] as TAllPathBond).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TAllPathBond).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TAllPathBond).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TAllPathBond).DynBnd := DynBnd;
      end
      else
      begin
        i := fFrgmnts.Add(TShrtPthBnd.Create(TotFrgList));
        (fFrgmnts.Items[i] as TShrtPthBnd).LenMin := lmin;
        (fFrgmnts.Items[i] as TShrtPthBnd).LenMax := lmax;
        (fFrgmnts.Items[i] as TShrtPthBnd).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TShrtPthBnd).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TShrtPthBnd).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TShrtPthBnd).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TShrtPthBnd).DynBnd := DynBnd;
      end;
    end
    else if (FrgType = 3) then
    begin
      if AtomPairs then
      begin
        i := fFrgmnts.Add(TShrtPthAtmBndPair.Create(TotFrgList));
        (fFrgmnts.Items[i] as TShrtPthAtmBndPair).LenMin := lmin;
        (fFrgmnts.Items[i] as TShrtPthAtmBndPair).LenMax := lmax;
        (fFrgmnts.Items[i] as TShrtPthAtmBndPair).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TShrtPthAtmBndPair).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TShrtPthAtmBndPair).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TShrtPthAtmBndPair).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TShrtPthAtmBndPair).DynBnd := DynBnd;
        (fFrgmnts.Items[i] as TShrtPthAtmBndPair).IsPair := True;
      end
      else if DoAllWays then
      begin
        i := fFrgmnts.Add(TAllPathAtmBond.Create(TotFrgList));
        (fFrgmnts.Items[i] as TAllPathAtmBond).LenMin := lmin;
        (fFrgmnts.Items[i] as TAllPathAtmBond).LenMax := lmax;
        (fFrgmnts.Items[i] as TAllPathAtmBond).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TAllPathAtmBond).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TAllPathAtmBond).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TAllPathAtmBond).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TAllPathAtmBond).DynBnd := DynBnd;
      end
      else
      begin
        i := fFrgmnts.Add(TShrtPthAtmBnd.Create(TotFrgList));
        (fFrgmnts.Items[i] as TShrtPthAtmBnd).LenMin := lmin;
        (fFrgmnts.Items[i] as TShrtPthAtmBnd).LenMax := lmax;
        (fFrgmnts.Items[i] as TShrtPthAtmBnd).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TShrtPthAtmBnd).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TShrtPthAtmBnd).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TShrtPthAtmBnd).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TShrtPthAtmBnd).DynBnd := DynBnd;
      end;
    end
    else if (FrgType = 4) then
    begin
      if (DoAllWays) then
        raise TUFragmentException.Create(
          'ERROR: Unsupported combination of arguments')
      else if AtomPairs then
      begin
        i := fFrgmnts.Add(TACAtmPair.Create(TotFrgList));
        (fFrgmnts.Items[i] as TACAtmPair).LenMin := lmin;
        (fFrgmnts.Items[i] as TACAtmPair).LenMax := lmax;
        (fFrgmnts.Items[i] as TACAtmPair).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TACAtmPair).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TACAtmPair).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TACAtmPair).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TACAtmPair).DynBnd := DynBnd;
        (fFrgmnts.Items[i] as TACAtmPair).IsPair := True;
      end
      else
      begin
        i := fFrgmnts.Add(TACPathAtom.Create(TotFrgList));
        (fFrgmnts.Items[i] as TACPathAtom).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TACPathAtom).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TACPathAtom).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TACPathAtom).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TACPathAtom).LenMin := lmin;
        (fFrgmnts.Items[i] as TACPathAtom).LenMax := lmax;
        (fFrgmnts.Items[i] as TACPathAtom).DynBnd := DynBnd;
      end;
    end
    else if (FrgType = 5) then
    begin
      if (AtomPairs or DoAllWays) then
        raise TUFragmentException.Create(
          'ERROR: Unsupported combination of arguments')
      else
      begin
        i := fFrgmnts.Add(TACPathBnd.Create(TotFrgList));
        (fFrgmnts.Items[i] as TACPathBnd).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TACPathBnd).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TACPathBnd).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TACPathBnd).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TACPathBnd).LenMin := lmin;
        (fFrgmnts.Items[i] as TACPathBnd).LenMax := lmax;
        (fFrgmnts.Items[i] as TACPathBnd).DynBnd := DynBnd;
      end;
    end
    else if (FrgType = 6) then
    begin
      if (DoAllWays) then
        raise TUFragmentException.Create(
          'ERROR: Unsupported combination of arguments')
      else if AtomPairs then
      begin
        i := fFrgmnts.Add(TACAtmBndPair.Create(TotFrgList));
        (fFrgmnts.Items[i] as TACAtmBndPair).LenMin := lmin;
        (fFrgmnts.Items[i] as TACAtmBndPair).LenMax := lmax;
        (fFrgmnts.Items[i] as TACAtmBndPair).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TACAtmBndPair).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TACAtmBndPair).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TACAtmBndPair).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TACAtmBndPair).DynBnd := DynBnd;
        (fFrgmnts.Items[i] as TACAtmBndPair).IsPair := True;
      end
      else
      begin
        i := fFrgmnts.Add(TACPathAtmBnd.Create(TotFrgList));
        (fFrgmnts.Items[i] as TACPathAtmBnd).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TACPathAtmBnd).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TACPathAtmBnd).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TACPathAtmBnd).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TACPathAtmBnd).LenMin := lmin;
        (fFrgmnts.Items[i] as TACPathAtmBnd).LenMax := lmax;
        (fFrgmnts.Items[i] as TACPathAtmBnd).DynBnd := DynBnd;
      end;
    end
    else if (FrgType = 7) then
    begin
      if (DoAllWays) then
        raise TUFragmentException.Create(
          'ERROR: Unsupported combination of arguments')
      else if (AtomPairs) then
      begin
        i := fFrgmnts.Add(TACFXAtmPair.Create(TotFrgList));
        (fFrgmnts.Items[i] as TACFXAtmPair).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TACFXAtmPair).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TACFXAtmPair).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TACFXAtmPair).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TACFXAtmPair).LenMin := lmin;
        (fFrgmnts.Items[i] as TACFXAtmPair).LenMax := lmax;
        (fFrgmnts.Items[i] as TACFXAtmPair).DynBnd := DynBnd;
        (fFrgmnts.Items[i] as TACFXAtmPair).IsPair := True;
      end
      else
      begin
        i := fFrgmnts.Add(TACFXPathAtom.Create(TotFrgList));
        (fFrgmnts.Items[i] as TACFXPathAtom).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TACFXPathAtom).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TACFXPathAtom).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TACFXPathAtom).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TACFXPathAtom).LenMin := lmin;
        (fFrgmnts.Items[i] as TACFXPathAtom).LenMax := lmax;
        (fFrgmnts.Items[i] as TACFXPathAtom).DynBnd := DynBnd;
      end;
    end
    else if (FrgType = 8) then
    begin
      if (AtomPairs or DoAllWays) then
        raise TUFragmentException.Create(
          'ERROR: Unsupported combination of arguments')
      else
      begin
        i := fFrgmnts.Add(TACFXPathBnd.Create(TotFrgList));
        (fFrgmnts.Items[i] as TACFXPathBnd).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TACFXPathBnd).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TACFXPathBnd).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TACFXPathBnd).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TACFXPathBnd).LenMin := lmin;
        (fFrgmnts.Items[i] as TACFXPathBnd).LenMax := lmax;
        (fFrgmnts.Items[i] as TACFXPathBnd).DynBnd := DynBnd;
      end;
    end
    else if (FrgType = 9) then
    begin
      if (DoAllWays) then
        raise TUFragmentException.Create(
          'ERROR: Unsupported combination of arguments')
      else if (AtomPairs) then
      begin
        i := fFrgmnts.Add(TACFXAtmBndPair.Create(TotFrgList));
        (fFrgmnts.Items[i] as TACFXAtmBndPair).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TACFXAtmBndPair).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TACFXAtmBndPair).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TACFXAtmBndPair).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TACFXAtmBndPair).LenMin := lmin;
        (fFrgmnts.Items[i] as TACFXAtmBndPair).LenMax := lmax;
        (fFrgmnts.Items[i] as TACFXAtmBndPair).DynBnd := DynBnd;
        (fFrgmnts.Items[i] as TACFXAtmBndPair).IsPair := True;
      end
      else
      begin
        i := fFrgmnts.Add(TACFXPathAtmBnd.Create(TotFrgList));
        (fFrgmnts.Items[i] as TACFXPathAtmBnd).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TACFXPathAtmBnd).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TACFXPathAtmBnd).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TACFXPathAtmBnd).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TACFXPathAtmBnd).LenMin := lmin;
        (fFrgmnts.Items[i] as TACFXPathAtmBnd).LenMax := lmax;
        (fFrgmnts.Items[i] as TACFXPathAtmBnd).DynBnd := DynBnd;
      end;
    end
    else if (FrgType = 10) then
    begin
      if (DoAllWays) then
        raise TUFragmentException.Create(
          'ERROR: Unsupported combination of arguments')
      else if (AtomPairs) then
        raise TUFragmentException.Create(
          'ERROR: Unsupported combination of arguments')
      else
      begin
        i := fFrgmnts.Add(TTriplet.Create(TotFrgList));
        (fFrgmnts.Items[i] as TTriplet).MarkAtom := MarkAtom;
        (fFrgmnts.Items[i] as TTriplet).UseFormalCharge := UseFormalCharge;
        (fFrgmnts.Items[i] as TTriplet).UseRadical := UseRadical;
        (fFrgmnts.Items[i] as TTriplet).UseIsotope := UseIsotope;
        (fFrgmnts.Items[i] as TTriplet).LenMin := lmin;
        (fFrgmnts.Items[i] as TTriplet).LenMax := lmax;
        (fFrgmnts.Items[i] as TTriplet).DynBnd := DynBnd;
      end;
    end
    else
      raise TUFragmentException.Create(
        'ERROR: Unsupported combination of arguments');
    (fFrgmnts.Last as TFrgBase).StrictFrg := StrictFrg;
    (fFrgmnts.Last as TFrgBase).ColorTerms.DelimitedText := DFrg.ColorAFields;
    (fFrgmnts.Last as TFrgBase).ColorBondSDField:=DFrg.ColorBFields;
    if (fFrgmnts.Last as TFrgBase).FrgPerAtom=nil then
      (fFrgmnts.Last as TFrgBase).FrgPerAtom:=fFrgPerAtom;
    (fFrgmnts.Last as TFrgBase).GetFrgPerAtom:=fGetFrgPerAtom;
    (fFrgmnts.Last as TFrgBase).bCycle:=Cycle;
  end;
end;

procedure TFragment.AddFragmentor(PDFrg: PDefFrg);
begin
  AddFragmentor(PDFrg^);
end;

procedure TFragment.AddFragmentor(ftype, LenMinValue, LenMaxValue,
  MarkAtomValue, DynBnd: byte;
  EqFuzAt, UseFormalCharge, DoAllWays, AtomPairs, StrictFrg: boolean);
var
  DFrg: DefFrg;
  i:    integer;
begin
  DFrg.FrgType   := ftype;
  DFrg.lmin      := LenMinValue;
  DFrg.lmax      := LenMaxValue;
  DFrg.MarkAtom  := MarkAtomValue;
  DFrg.EqFuzAt   := EqFuzAt;
  DFrg.UseFormalCharge := UseFormalCharge;
  DFrg.DoAllWays := DoAllWays;
  DFrg.AtomPairs := AtomPairs;
  DFrg.StrictFrg := StrictFrg;
  DFrg.DynBnd    := DynBnd;
  DFrg.ColorAFields := 'Default';
  DFrg.ColorBFields := 'Default';
  AddFragmentor(DFrg);
end;

procedure TFragment.MolToFrgLst(MolList: TStringList);
var
  i, j: integer;
  stmp: string;
begin
  SetMol(MolList);
  for i := 0 to fFrgmnts.Count - 1 do
  begin
    //(fFrgmnts.Items[i] as TFrgBase).ColorTerms.DelimitedText:=(fFrgmnts.Items[i] as TFrgBase).ColorFields;
    fMol.LoadSDFField(MolList, (fFrgmnts.Items[i] as TFrgBase).ColorTerms,
      (fFrgmnts.Items[i] as TFrgBase).ColorKeys, (fFrgmnts.Items[i] as TFrgBase).ColorHash);
    if (fFrgmnts.Items[i] as TFrgBase).bCycle then
      fMol.AnnotateCycleBnd(False)
    else
      fMol.AnnotateCycleBnd(True);
    fMol.LoadSDFFieldBdColor(MolList,(fFrgmnts.Items[i] as TFrgBase).ColorBondSDField);//Shall be included in LoadSDFField?
    if fFrgmnts.Items[i] is TAtmCnt then
      TAtmCnt(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TShrtPthAtm then
      TShrtPthAtm(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TAllPathAtm then
      TAllPathAtm(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TShrtPthAtmPair then
      TShrtPthAtmPair(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TShrtPthBnd then
      TShrtPthBnd(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TAllPathBond then
      TAllPathBond(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TShrtPthAtmBnd then
      TShrtPthAtmBnd(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TShrtPthAtmBndPair then
      TShrtPthAtmBndPair(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TAllPathAtmBond then
      TAllPathAtmBond(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TACPathAtom then
      TACPathAtom(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TACPathBnd then
      TACPathBnd(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TACPathAtmBnd then
      TACPathAtmBnd(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TACAtmPair then
      TACAtmPair(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TACAtmBndPair then
      TACAtmBndPair(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TACFXPathAtom then
      TACFXPathAtom(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TACFXPathBnd then
      TACFXPathBnd(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TACFXPathAtmBnd then
      TACFXPathAtmBnd(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TACFXAtmPair then
      TACFXAtmPair(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TACFXAtmBndPair then
      TACFXAtmBndPair(fFrgmnts.Items[i]).MolToFrgLst(fMol)
    else if fFrgmnts.Items[i] is TTriplet then
      TTriplet(fFrgmnts.Items[i]).MolToFrgLst(fMol);
  end;
end;

procedure TFragment.SetFrgLst(fraglist: TStringList);
var
  i: integer;
begin
  fTotFrgList.Assign(fraglist);
end;

function TFragment.GetFrgLst: TStringList;
begin
  Result := fTotFrgList;
end;

function TFragment.IsAnyNewFrg: boolean;
var
  i: integer;
begin
  Result := False;
  for i := 0 to fFrgmnts.Count - 1 do
    if (TFrgBase(fFrgmnts[i]).NewFrg) then
      Result := True;
end;

end.

