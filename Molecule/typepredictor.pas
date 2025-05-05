{ Fragmentor of the ISIDA Project

  Copyright (C) 2024 Laboratoire de Chemoinformatique, UMR 7140 CNRS (https://complex-matter.unistra.fr/equipes-de-recherche/laboratoire-de-chemoinformatique/home/)
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
unit TypePredictor;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils;

type
  MethodException = class(Exception);

  PEnumFileType = ^EnumFileType;
  EnumModelCategory=(REG,CLS,CLS5,CLU,UNKCategory);
  EnumModelType = (MLR, ASNN, MTL, SVM, CMD, UNKMethod);
  EnumADType = (zkNN, FrgCtrl, CmnFrgt, MinMax, UNKAD);
  EnumFileType = (CFR, ASNNModel, ASNNCfg, ASNNTrain, ASNNHead);
  EnumFrg = (AC, IA, IB, IAB, IIA, IIB, IIAB, IIHy, IIIA, IIIB, IIIAB, IVA, IVB, IVAB, UNKFRG);
  EnumTopFrg = (E,I,II,III, UNKTOPFRG);
  EnumExtTopFrg= (FX,NoFX,UNKFX);
  EnumColAtm = (A,Ph,Ep,Ba,NoCA,UNKCOLATM);
  EnumColBnd = (B,NoCB,UNKCOLBND);
  EnumCounting = (ms,NoCNT,UNKCNT);
  EnumOptFrg = (P, R, AP, FC, MA, MP, SF, AD, OD, NoOPT, UNKOPTFRG);

  PFrgtn = ^Frgtn;

  Frgtn = record
    FrgType: byte;
    Min:     byte;
    Max:     byte;
  end;

  { TModelFiles }

  TModelFiles = class(TStringList)
  private
  public
    constructor Create;
    destructor Destroy; override;
    procedure Clear; override;
  end;

  RASNNFiles = record
    ModelFile: string;
    CFGFile:   string;
    Moltest:   integer;
    TrainASNNFile: string;
    HdrTrainingfile: string;
    PropCol:   integer;
    NbProp:    integer;
  end;

  PRModel = ^RModel;

  RModel = record
    ModelType: EnumModelType;
    NbMFiles: integer;
    ModelFiles: TModelFiles;
    PropCol: integer;
    NbProp: integer;
    Frg: array of Frgtn;
  end;

  PRAD = ^RAD;

  RAD = record
    NbTFiles: integer;
    HdrFiles: TStringList;
    SmfFiles: TStringList;
    AD_typ:   array of EnumADType;
    AD_k:     array of integer;
    AD_z:     array of double;
    AD_min:   array of integer;
    AD_max:   array of integer;
  end;

  PRModelAD = ^RModelAD;

  RModelAD = record
    Model: RModel;
    AD:    RAD;
  end;

  { LModel }

  LModel = class(TList)
  private
    function GetRModel(Index: integer): RModel;
    procedure SetRModel(Index: integer; Value: RModel);
  public
    destructor Destroy; override;
    procedure Clear; override;
    procedure AddRModel(m: RModel);
    property Objects[Index: integer]: RModel Read GetRModel Write SetRModel;
  end;

  { LAD }

  LAD = class(TList)
  private
    function GetRAD(Index: integer): RAD;
    procedure SetRAD(Index: integer; Value: RAD);
  public
    destructor Destroy; override;
    procedure Clear; override;
    procedure AddRAD(m: RAD);
    property Objects[Index: integer]: RAD Read GetRAD Write SetRAD;
  end;

  { LModelAD }

  LModelAD = class(TList)
  private
    function GetRModelAD(Index: integer): RModelAD;
    procedure SetRModelAD(Index: integer; Value: RModelAD);
  public
    destructor Destroy; override;
    procedure Clear; override;
    procedure AddRModelAD(m: RModelAD);
    property Objects[Index: integer]: RModelAD Read GetRModelAD Write SetRModelAD;
  end;

function StrToEnumModelCategory(S: string): EnumModelCategory;
function StrToEnumModelType(S: string): EnumModelType;
function StrToEnumADType(S: string): EnumADType;
function StrToEnumFrg(S: string): EnumFrg;
function StrToEnumTopFrg(S: string): EnumTopFrg;
function StrToEnumColAtm(S: string): EnumColAtm;
function StrToEnumColBnd(S: string): EnumColBnd;
function StrToEnumCounting(S: string): EnumCounting;
function StrToEnumOptFrg(S: string): EnumOptFrg;
function StrToEnumExtTopFrg(S: string): EnumExtTopFrg;
function EnumModelCategoryToStr(EMCat: EnumModelCategory): string;
function EnumModelTypeToStr(EMType: EnumModelType): string;
function EnumADTypeToStr(EADType: EnumADType): string;
function EnumTopFrgToStr(En: EnumTopFrg): string;
function EnumColAtmToStr(E: EnumColAtm): string;
function EnumColBndToStr(E: EnumColBnd): string;
function EnumCountingToStr(E: EnumCounting): string;
function EnumOptFrgToStr(E: EnumOptFrg): string;
function EnumEnumExtTopFrg(E: EnumExtTopFrg): string;

implementation

function StrToEnumModelCategory(S: string): EnumModelCategory;
begin
  if (S = 'REG') then
    Result := REG
  else if (S = 'CLS') then
    Result := CLS
  else if (S = 'CLU') then
    Result := CLU
  else
    Result := UNKCategory;
end;

function StrToEnumModelType(S: string): EnumModelType;
begin
  if (S = 'MLR') then
    Result := MLR
  else if (S = 'ASNN') then
    Result := ASNN
  else if (S = 'MTL') then
    Result := MTL
  else if (S = 'SVM') then
    Result := SVM
  else if (S = 'CMD') then
    Result := CMD
  else
    Result := UNKMethod;
end;

function StrToEnumADType(S: string): EnumADType;
begin
  if (S = 'zkNN') then
    Result := zkNN
  else if (S = 'FrgCtrl') then
    Result := FrgCtrl
  else if (S = 'CmnFrgt') then
    Result := CmnFrgt
  else if (S = 'MinMax') then
    Result := MinMax
  else
    Result := UNKAD;
end;

function StrToEnumFrg(S: string): EnumFrg;
begin
  if (S = 'AC') then
    Result := AC
  else if (S = 'IA') then
    Result := IA
  else if (S = 'IB') then
    Result := IB
  else if (S = 'IAB') then
    Result := IAB
  else if (S = 'IIA') then
    Result := IIA
  else if (S = 'IIB') then
    Result := IIB
  else if (S = 'IIAB') then
    Result := IIAB
  else if (S = 'IIHy') then
    Result := IIHy
  else if (S = 'IIIA') then
    Result := IIIA
  else if (S = 'IIIB') then
    Result := IIIB
  else if (S = 'IIIAB') then
    Result := IIIAB
  else if (S = 'IVA') then
    Result := IVA
  else if (S = 'IVB') then
    Result := IVB
  else if (S = 'IVAB') then
    Result := IVAB
  else
    Result := UNKFRG;
end;

function StrToEnumTopFrg(S: string): EnumTopFrg;
begin
  if (S='E') then
    Result:=E
  else if (S='I') then
    Result:=I
  else if (S='II') then
    Result:=II
  else if (S='III') then
    Result:=III
  else
    Result:=UNKTOPFRG;
end;

function StrToEnumColAtm(S: string): EnumColAtm;
begin
  if (S='A') then
    Result:=A
  else if (S='Ph') then
    Result:=Ph
  else if (S='Ep') then
    Result:=Ep
  else if (S='Ba') then
    Result:=Ba
  else if (S='NoCA') then
    Result:=NoCA
  else
    Result:=UNKCOLATM;
end;

function StrToEnumColBnd(S: string): EnumColBnd;
begin
  if (S='B') then
    Result:=B
  else if (S='NoCB') then
    Result:=NoCB
  else
    Result:=UNKCOLBND;
end;

function StrToEnumCounting(S: string): EnumCounting;
begin
  if (S='ms') then
    Result:=ms
  else if (S='NoCNT') then
    Result:=NoCNT
  else
    Result:=UNKCNT;
end;

function StrToEnumOptFrg(S: string): EnumOptFrg;
begin
  if (S='P') then
    Result:=P
  else if (S='R') then
    Result:=R
  else if (S='AP') then
    Result:=AP
  else if (S='FC') then
    Result:=FC
  else if (S='MA') then
    Result:=MA
  else if (S='MP') then
    Result:=MP
  else if (S='SF') then
    Result:=SF
  else if (S='AD') then
    Result:=AD
  else if (S='OD') then
    Result:=OD
  else if (S='NoOPT') then
    Result:=NoOPT
  else
    Result:=UNKOPTFRG;
end;

function StrToEnumExtTopFrg(S: string): EnumExtTopFrg;
begin
  if (S='FX') then
    Result:=FX
  else if (S='NoFX') then
    Result:=NoFX
  else
    Result:=UNKFX;
end;

function EnumModelCategoryToStr(EMCat: EnumModelCategory): string;
begin
  case EMCat of
    REG: Result  := 'REG';
    CLS: Result := 'CLS';
    CLU: Result  := 'CLU';
    UNKCategory: Result := 'Unknown';
    else
      Result := 'Not supported';
  end;
end;

function EnumModelTypeToStr(EMType: EnumModelType): string;
begin
  case EMType of
    MLR: Result  := 'MLR';
    ASNN: Result := 'ASNN';
    SVM: Result  := 'SVM';
    MTL: Result  := 'MTL';
    CMD: Result  := 'CMD';
    UNKMethod: Result := 'Unknown';
    else
      Result := 'Not supported';
  end;
end;

function EnumADTypeToStr(EADType: EnumADType): string;
begin
  case EADType of
    zkNN: Result    := 'zkNN';
    FrgCtrl: Result := 'FrgCtrl';
    CmnFrgt: Result := 'CmnFrgt';
    MinMax: Result  := 'MinMax';
    UNKAD: Result   := 'Unknown';
    else
      Result := 'Not supported';
  end;
end;

function EnumTopFrgToStr(En: EnumTopFrg): string;
begin
  case En of
    E: Result:='E';
    I: Result:='I';
    II: Result:='II';
    III: Result:='III';
    else
      Result:='UNKTOPFRG';
  end;
end;

function EnumColAtmToStr(E: EnumColAtm): string;
begin
  case E of
    A: Result:='A';
    Ph: Result:='Ph';
    Ep: Result:='Ep';
    Ba: Result:='Ba';
    NoCA: Result:='NoCA';
    else
      Result:='UNKCOLATM';
  end;
end;

function EnumColBndToStr(E: EnumColBnd): string;
begin
  case E of
    B: Result:='B';
    NoCB: Result:='NoCB';
    else
      Result:='UNKCOLBND';
  end;
end;

function EnumCountingToStr(E: EnumCounting): string;
begin
  case E of
    ms: Result:='ms';
    NoCNT: Result:='NoCNT';
    else
      Result:='UNKCNT';
  end;
end;

function EnumOptFrgToStr(E: EnumOptFrg): string;
begin
  case E of
    P: Result:='P';
    R: Result:='R';
    AP: Result:='AP';
    FC: Result:='FC';
    MA: Result:='MA';
    MP: Result:='MP';
    SF: Result:='SF';
    AD: Result:='AD';
    OD: Result:='OD';
    NoOPT: Result:='NoOPT';
    else
      Result:='UNKOPTFRG';
  end;
end;

function EnumEnumExtTopFrg(E: EnumExtTopFrg): string;
begin
  case E of
    FX: Result:='FX';
    NoFX: Result:='NoFX';
    else
      Result:='UNKFX';
  end;
end;

{ LModel }

function LModel.GetRModel(Index: integer): RModel;
var
  pitem: PRModel;
begin
  pitem  := Items[Index];
  Result := pitem^;
end;

procedure LModel.SetRModel(Index: integer; Value: RModel);
var
  pitem: PRModel;
begin
  pitem := Items[Index];
  dispose(pitem);
  new(pitem);
  pitem^ := Value;
  Items[Index] := pitem;
end;

destructor LModel.Destroy;
begin
  Clear;
  inherited Destroy;
end;

procedure LModel.Clear;
var
  pitem: PRModel;
  i:     integer;
begin
  for i := 0 to Count - 1 do
  begin
    pitem := Items[i];
    dispose(pitem);
  end;
  inherited Clear;
end;

procedure LModel.AddRModel(m: RModel);
var
  pitem: PRModel;
begin
  new(pitem);
  pitem^ := m;
  Add(pitem);
end;

{ LAD }

function LAD.GetRAD(Index: integer): RAD;
var
  pitem: PRAD;
begin
  pitem  := Items[Index];
  Result := pitem^;
end;

procedure LAD.SetRAD(Index: integer; Value: RAD);
var
  pitem: PRAD;
begin
  pitem := Items[Index];
  dispose(pitem);
  new(pitem);
  pitem^ := Value;
  Items[Index] := pitem;
end;

destructor LAD.Destroy;
begin
  Clear;
  inherited Destroy;
end;

procedure LAD.Clear;
var
  pitem: PRAD;
  i:     integer;
begin
  for i := 0 to Count - 1 do
  begin
    pitem := Items[i];
    dispose(pitem);
  end;
  inherited Clear;
end;

procedure LAD.AddRAD(m: RAD);
var
  pitem: PRAD;
begin
  new(pitem);
  pitem^ := m;
  Add(pitem);
end;

{ LModelAD }

function LModelAD.GetRModelAD(Index: integer): RModelAD;
var
  pitem: PRModelAD;
begin
  pitem  := Items[Index];
  Result := pitem^;
end;

procedure LModelAD.SetRModelAD(Index: integer; Value: RModelAD);
var
  pitem: PRModelAD;
begin
  pitem := Items[Index];
  dispose(pitem);
  new(pitem);
  pitem^ := Value;
  Items[Index] := pitem;
end;

destructor LModelAD.Destroy;
begin
  Clear;
  inherited Destroy;
end;

procedure LModelAD.Clear;
var
  pitem: PRModelAD;
  i:     integer;
begin
  for i := 0 to Count - 1 do
  begin
    pitem := Items[i];
    dispose(pitem);
  end;
  inherited Clear;
end;

procedure LModelAD.AddRModelAD(m: RModelAD);
var
  pitem: PRModelAD;
begin
  new(pitem);
  pitem^ := m;
  Add(pitem);
end;

{ TModelFiles }

constructor TModelFiles.Create;
begin
  inherited Create;
end;

destructor TModelFiles.Destroy;
begin
  Clear;
  inherited Destroy;
end;

procedure TModelFiles.Clear;
var
  i: integer;
begin
  for i := 0 to Count - 1 do
    dispose(PEnumFileType(Objects[i]));
  inherited Clear;
end;

end.

