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
unit UnitAtomBase;
//unit to define atom based representations
{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, {UnitFragmentBase,} UnitRepresentationBase;

const
  AtSymMax = 4;                         // max letters for atom symbol
  strdflt  = '§';

type

  { TAtomBase }

  TAtomBase = class(TRepresentationBase)
  private
    fAtomString:   TStringList;
    fUseStereo:    boolean;
    fStereoParity: byte;
    function GetAtomString(indx: integer): string;
    procedure SetAtomString(indx: integer; const AValue: string);
  protected
    atomsTable_size: integer;           // periodic table size
    atomsTable:      array[1..MaxAtomSum] of string[AtSymMax]; // periodic table
  public
    constructor Create;
    destructor Destroy; override;
    procedure Clear; override;
    procedure init_atomsTable; virtual; abstract;
    procedure InitAtomString(sdline: string);
    procedure InitAtomString(sdline: TStringList);
    function InitWeight(msstr: string): integer;
    function InitAtomStringWeight(msstr: string): integer;
    // atom representation table
    function myAtomString(str: string): string;
    function GetSymbol(Id: integer): string; virtual;
    function ConvertFormalCharge(Id: integer):string; virtual;
    property AtomString[indx: integer]: string Read GetAtomString Write SetAtomString;
    property UseStereo: boolean Read fUseStereo Write fUseStereo;
    property StereoParity: byte Read fStereoParity Write fStereoParity;
  end;
  function ConvertFormalChargep(Id: integer):string;
  function UnConvertFormalCharge(Chg: Integer):byte;

implementation

{ TAtomBase }

function TAtomBase.GetAtomString(indx: integer): string;
begin
  if indx >= fAtomString.Count then
    Result := strdflt
  else
    Result := fAtomString[indx];
end;

procedure TAtomBase.SetAtomString(indx: integer; const AValue: string);
var
  j: integer;
begin
  if (fAtomString.Count - 1 < indx) then
    for j := fAtomString.Count to indx do
      fAtomString.Add(strdflt);
  fAtomString[indx] := AValue;
end;

constructor TAtomBase.Create;
var
  i: integer;
begin
  inherited Create;
  for i := Low(atomsTable) to High(atomsTable) do
    atomsTable[i] := '';
  fAtomString     := TStringList.Create;
  //fUseSymbol      := True;
  fUseStereo      := False;
  //fUseFormalCharge := False;
  //fFormalCharge   := 0;
  fStereoParity   := 0;
end;

destructor TAtomBase.Destroy;
begin
  inherited Destroy;
  FreeAndNil(fAtomString);
end;

procedure TAtomBase.Clear;
var
  i: integer;
begin
  inherited Clear;
  for i := Low(atomsTable) to High(atomsTable) do
    atomsTable[i] := '';
  fAtomString.Clear;
  fUseStereo    := False;
  //fUseFormalCharge := False;
  fStereoParity := 0;
end;

procedure TAtomBase.InitAtomString(sdline: string);
var
  tuples:  TStringList;
  tuple:   TStringList;
  i, j, n: integer;
begin
  fAtomString.Clear;
  tuples := TStringList.Create;
  tuple  := TStringList.Create;
  tuples.Delimiter := ' ';
  tuple.Delimiter := ':';
  tuples.DelimitedText := sdline;
  for i := 1 to tuples.Count - 1 do
  begin
    tuple.DelimitedText := tuples[i];
    n := StrToInt(tuple[0]);
    if (fAtomString.Count - 1 < n) then
      for j := fAtomString.Count to n do
        fAtomString.Add(strdflt);
    fAtomString[n] := tuple[1];
  end;
  //exemple: 1:A/D 2:A 4:H donne une liste fAtSym[0]="A/D" et fAtSym[1]="A" et fAtSym[3]=" H"
  FreeAndNil(tuples);
  FreeAndNil(tuple);
end;

procedure TAtomBase.InitAtomString(sdline: TStringList);
begin
  fAtomString.Assign(sdline);
end;


function TAtomBase.InitWeight(msstr: string): integer;
  // function in order to read micro-species line in the colored sdf field (ColorHash) and extract population level for incrementation of fragments
var
  iw:    integer;
  tuple: TStringList;
begin
  if (msstr = '') then
    iw := 1
  else
  begin
    tuple := TStringList.Create;
    tuple.StrictDelimiter := True;
    tuple.Delimiter := ' ';
    tuple.DelimitedText := msstr;
    try    //try to initialize population which is first given integer in line. If not given it will be set to 1
      iw := StrToInt(tuple[0]);
    except
      On E: EConvertError do begin
        iw := 1;
      end;
    end;
    FreeAndNil(tuple);
  end;
  Result := iw;
end;

function TAtomBase.InitAtomStringWeight(msstr: string): integer;
var
  tuples: TStringList;
  tuple:  TStringList;
  i, j, n, iw: integer;
begin
  fAtomString.Clear;
  tuples := TStringList.Create;
  tuples.StrictDelimiter := True;
  tuples.Delimiter := ' ';
  tuple  := TStringList.Create;
  tuple.StrictDelimiter := True;
  tuple.Delimiter := ':';
  tuples.DelimitedText := msstr; //tuples.DelimitedText := Trim(msstr);
  for i := 1 to tuples.Count - 1 do
  begin
    tuple.DelimitedText := tuples[i];
    n := StrToInt(tuple[0]);
    if (fAtomString.Count - 1 < n) then
      for j := fAtomString.Count to n do
        fAtomString.Add(strdflt);
    fAtomString[n] := tuple[1];
  end;
  //exemple: 1:A/D 2:A 4:H donne une liste fAtSym[0]="A/D" et fAtSym[1]="A" et fAtSym[3]=" H"
  //set iw=weight
  if (msstr = '') then
    iw := 1
  else if (tuples[0] = '?') then
    iw := 1
  else
  begin
    try    //try to initialize population which is first given integer in line. If not given it will be set to 1
      iw := StrToInt(tuples[0]);
    except
      On E: EConvertError do
      begin
        writeln('Error - could not recognise weight format of microspecies: should be either "?" or an integer');
        Halt;
      end;

    end;
  end;
  FreeAndNil(tuples);
  FreeAndNil(tuple);
  Result := iw;
end;

function TAtomBase.myAtomString(str: string): string;
var
  ip: integer;
begin
  ip     := Pos('_', str);
  Result := str;
  if (ip > 0) then
    Result := Copy(str, 1, ip - 1);
end;

function TAtomBase.GetSymbol(Id: integer): string;
begin
  Result := AtomString[Id];
  {if (fUseStereo and (fStereoParity <> 0)) then
    Result := Result + '_P' + IntToStr(fStereoParity);
  if (fUseFormalCharge) then
    Result := Result + '_C' + IntToStr();}
end;

function TAtomBase.ConvertFormalCharge(Id: integer):string;
begin
     Result := ConvertFormalChargep(Id);
end;

function ConvertFormalChargep(Id: integer):string;
begin
     if (Id = 1) then Result := '+3'
     else if (Id = 2) then Result := '+2'
     else if (Id = 3) then Result := '+1'
     else if (Id = 4) then Result := 'R' //doublet radical
     else if (Id = 5) then Result := '-1'
     else if (Id = 6) then Result := '-2'
     else if (Id = 7) then Result := '-3';
end;

function UnConvertFormalCharge(Chg: Integer): byte;
begin
  case Chg of
    -15 ..-4: Result := 9;
    4 ..15: Result := 8;
    -3: Result := 7;
    -2: Result := 6;
    -1: Result := 5;
    1:  Result := 3;
    2:  Result := 2;
    3:  Result := 1;
    else Result := 0;
  end;
end;
end.

