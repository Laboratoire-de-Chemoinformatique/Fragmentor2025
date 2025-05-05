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
unit UnitBondBase;
//Unit to describe bond representation
{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, UnitRepresentationBase;

const
     BTypeMax = 84;                        // max BondSymbol or BondCode
     BoSymMax = 1;                         // max letters for bond symbol
type

    { TBondBase }

    TBondBase = class(TRepresentationBase)
    private
           fUseStereo: boolean;
           fStereoBond: byte;
    protected
             BondSymbol: array [1..BTypeMax] of string[2];
    public
          constructor Create;
          destructor Destroy; override;
          procedure Clear; virtual;
          procedure init_BondTable; virtual;
          function myBondSymbol(str: string): integer; virtual;
          function GetBondSymbol(a :integer): string; virtual;
          property UseStereo: boolean read fUseStereo write fUseStereo;
          property StereoBond: byte read fStereoBond write fStereoBond;
    end;

implementation

{ TBondBase }

constructor TBondBase.Create;
var
   i: integer;
begin
     inherited Create;
     init_BondTable;
     fUseStereo:=False;
     fStereoBond:=0;
end;

destructor TBondBase.Destroy;
begin
  inherited Destroy;
end;

procedure TBondBase.Clear;
begin
     init_BondTable;
     fUseStereo:=False;
     fStereoBond:=0;
end;

procedure TBondBase.init_BondTable;
var
   i: integer;
begin
     for i:=1 to BTypeMax do BondSymbol[i]:='YY';
     BondSymbol[1]:='-';
     BondSymbol[2]:='=';
     BondSymbol[3]:='+';
     BondSymbol[4]:='*';
     BondSymbol[5]:='.';
     BondSymbol[6]:=':';
     BondSymbol[7]:='"';
     BondSymbol[8]:='~';
     BondSymbol[9]:='&';
     BondSymbol[12]:='2';
     BondSymbol[13]:='4';
     BondSymbol[18]:='5';
     BondSymbol[21]:='6';
     BondSymbol[23]:='3';
     BondSymbol[28]:='7';
     BondSymbol[31]:='9';
     BondSymbol[32]:='8';
     BondSymbol[81]:='1';
     BondSymbol[82]:='D';
     BondSymbol[83]:='T';
end;

function TBondBase.myBondSymbol(str: string): integer;
var
   ip: integer;
begin
     ip:=Pos('_',str);
     if (ip>0) then str:=Copy(str,1,ip-1);
     Result := Low(BondSymbol);
     repeat
           inc(Result);
     until ((Result>=High(BondSymbol)) or (BondSymbol[Result]=str));
end;

function TBondBase.GetBondSymbol(a: integer): string;
begin
     Result:=BondSymbol[a];
     if (fUseStereo and (fStereoBond<>0)) then Result:=Result+'_'+IntToStr(fStereoBond);
end;

end.

