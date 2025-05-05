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
unit UnitAtmPrpWeight;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, UnitRepresentationBase ,unitatmpropbase;

type

{ TAtmPropWeight }

TAtmPropWeight = class(TAtmPropBase)
    private
           prpTable_size: integer;
    public
          constructor Create;
          destructor Destroy; override;
          procedure init_atmPrTable; override;
          function AtmPrCompare(Id1, Id2: integer): double; override;
    end;

implementation

{ TAtmPropWeight }

constructor TAtmPropWeight.Create;
begin
     inherited Create;
end;

destructor TAtmPropWeight.Destroy;
begin
     inherited Destroy;
end;

procedure TAtmPropWeight.init_atmPrTable;
begin
     atmPrTable[1] := 1.00794;
     atmPrTable[2] := 4.00260;
     atmPrTable[3] := 6.941;
     atmPrTable[4] := 9.01218;
     atmPrTable[5] := 10.81;
     atmPrTable[6] := 12.011;
     atmPrTable[7] := 14.0067;
     atmPrTable[8] := 15.9994;
     atmPrTable[9] := 18.998403;
    atmPrTable[10] := 20.179;
    atmPrTable[11] := 22.98977;
    atmPrTable[12] := 24.30;
    atmPrTable[13] := 26.98154;
    atmPrTable[14] := 28.0855;
    atmPrTable[15] := 30.97376;
    atmPrTable[16] := 32.06;
    atmPrTable[17] := 35.453;
    atmPrTable[18] := 39.948;
    atmPrTable[19] := 39.0983;
    atmPrTable[20] := 40.08;
    atmPrTable[21] := 44.9559;
    atmPrTable[22] := 47.88;
    atmPrTable[23] := 50.9415;
    atmPrTable[24] := 51.996;
    atmPrTable[25] := 54.9380;
    atmPrTable[26] := 55.847;
    atmPrTable[27] := 58.9332;
    atmPrTable[28] := 55.69;
    atmPrTable[29] := 63.546;
    atmPrTable[30] := 65.38;
    atmPrTable[31] := 69.72;
    atmPrTable[32] := 72.59;
    atmPrTable[33] := 74.9216;
    atmPrTable[34] := 78.96;
    atmPrTable[35] := 79.904;
    atmPrTable[36] := 83.80;
    atmPrTable[37] := 85.4678;
    atmPrTable[38] := 87.62;
    atmPrTable[39] := 88.9069;
    atmPrTable[40] := 91.22;
    atmPrTable[41] := 92.9064;
    atmPrTable[42] := 95.94;
    atmPrTable[43] := 98.00;
    atmPrTable[44] := 101.07;
    atmPrTable[45] := 102.9055;
    atmPrTable[46] := 106.42;
    atmPrTable[47] := 104.8682;
    atmPrTable[48] := 112.41;
    atmPrTable[49] := 114.82;
    atmPrTable[50] := 118.69;
    atmPrTable[51] := 121.75;
    atmPrTable[52] := 127.60;
    atmPrTable[53] := 126.9045;
    atmPrTable[54] := 131.29;
    atmPrTable[55] := 132.9054;
    atmPrTable[56] := 137.33;
    atmPrTable[57] := 138.9055;
    atmPrTable[58] := 140.12;
    atmPrTable[59] := 140.9077;
    atmPrTable[60] := 144.24;
    atmPrTable[61] := 145.00;
    atmPrTable[62] := 150.36;
    atmPrTable[63] := 151.96;
    atmPrTable[64] := 157.25;
    atmPrTable[65] := 158.9254;
    atmPrTable[66] := 162.50;
    atmPrTable[67] := 164.9304;
    atmPrTable[68] := 167.26;
    atmPrTable[69] := 168.9342;
    atmPrTable[70] := 173.04;
    atmPrTable[71] := 174.967;
    atmPrTable[72] := 178.49;
    atmPrTable[73] := 180.9479;
    atmPrTable[74] := 183.85;
    atmPrTable[75] := 186.207;
    atmPrTable[76] := 190.2;
    atmPrTable[77] := 192.22;
    atmPrTable[78] := 195.08;
    atmPrTable[79] := 196.9665;
    atmPrTable[80] := 200.59;
    atmPrTable[81] := 204.383;
    atmPrTable[82] := 207.2;
    atmPrTable[83] := 208.9804;
    atmPrTable[84] := 209.00;
    atmPrTable[85] := 210.00;
    atmPrTable[86] := 222.00;
    atmPrTable[87] := 223.00;
    atmPrTable[88] := 226.0254;
    atmPrTable[89] := 227.0278;
    atmPrTable[90] := 232.0381;
    atmPrTable[91] := 231.0359;
    atmPrTable[92] := 238.0289;
    atmPrTable[93] := 237.0482;
    atmPrTable[94] := 244.00;
    atmPrTable[95] := 243.00;
    atmPrTable[96] := 247.00;
    atmPrTable[97] := 247.00;
    atmPrTable[98] := 251.00;
    atmPrTable[99] := 252.00;
    atmPrTable[100] := 257.00;
    atmPrTable[101] := 258.00;
    atmPrTable[102] := 259.00;
    atmPrTable[103] := 260.00;
    prpTable_size:=MaxAtomSum;
end;

function TAtmPropWeight.AtmPrCompare(Id1, Id2: integer): double;
begin
     Result:=atmPrTable[Id1]+atmPrTable[Id2];
end;

end.

