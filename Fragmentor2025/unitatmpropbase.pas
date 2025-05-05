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
unit unitatmpropbase;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, UnitRepresentationBase;
  
type

{ TAtmProp }

TAtmPropBase = class(TRepresentationBase)
  private
  protected
           atmPrTable: array[1..MaxAtomSum] of double;
  public
           //atmPrTable: array[1..MaxAtomSum] of double;
        constructor Create;
        destructor Destroy; override;
        procedure Clear; override;
        procedure init_atmPrTable; virtual; abstract;        // atom representation table
        function GetAtmPrp(Id: integer): double; virtual;
        function AtmPrCompare(Id1, Id2: integer): double; virtual; abstract; //method to compare to atoms
  end;

implementation

{ TAtmProp }

constructor TAtmPropBase.Create;
var
   i: integer;
begin
     inherited Create;
     for i:=Low(atmPrTable) to High(atmPrTable) do
         atmPrTable[i]:=1;
end;

destructor TAtmPropBase.Destroy;
begin
     inherited Destroy;
end;

procedure TAtmPropBase.Clear;
var
   i: integer;
begin
     inherited Clear;
     for i:=Low(atmPrTable) to High(atmPrTable) do
         atmPrTable[i]:=1;
end;

function TAtmPropBase.GetAtmPrp(Id: integer): double;
begin
     Result:=atmPrTable[Id];
end;

end.

