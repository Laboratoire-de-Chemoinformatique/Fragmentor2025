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
unit UnitAtomColumn;
//Unit describing atoms as the column of periodic table they belong to
{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, UnitAtomSymbol;

type

    { TAtomColumn }

    TAtomColumn = class(TAtomSymbol)
    private

    public
          constructor Create;
          destructor Destroy; override;
          procedure init_atomsTable; override;
    end;

implementation

{ TAtomColumn }

constructor TAtomColumn.Create;
begin
     inherited Create;
     init_atomsTable;
end;

destructor TAtomColumn.Destroy;
begin
     inherited Destroy;
end;

procedure TAtomColumn.init_atomsTable;
begin
     atomsTable[1] := '1';
     atomsTable[2] := '8';
     atomsTable[3] := '1';
     atomsTable[4] := '2';
     atomsTable[5] := '3';
     atomsTable[6] := '4';
     atomsTable[7] := '5';
     atomsTable[8] := '6';
     atomsTable[9] := '7';
    atomsTable[10] := '8';
    atomsTable[11] := '1';
    atomsTable[12] := '2';
    atomsTable[13] := '3';
    atomsTable[14] := '4';
    atomsTable[15] := '5';
    atomsTable[16] := '6';
    atomsTable[17] := '7';
    atomsTable[18] := '8';
    atomsTable[19] := '1';
    atomsTable[20] := '2';
    atomsTable[21] := 'M1';
    atomsTable[22] := 'M2';
    atomsTable[23] := 'M3';
    atomsTable[24] := 'M4';
    atomsTable[25] := 'M5';
    atomsTable[26] := 'M6';
    atomsTable[27] := 'M7';
    atomsTable[28] := 'M8';
    atomsTable[29] := 'M9';
    atomsTable[30] := 'M0';
    atomsTable[31] := '3';
    atomsTable[32] := '4';
    atomsTable[33] := '5';
    atomsTable[34] := '6';
    atomsTable[35] := '7';
    atomsTable[36] := '8';
    atomsTable[37] := '1';
    atomsTable[38] := '2';
    atomsTable[39] := 'M1';
    atomsTable[40] := 'M2';
    atomsTable[41] := 'M3';
    atomsTable[42] := 'M4';
    atomsTable[43] := 'M5';
    atomsTable[44] := 'M6';
    atomsTable[45] := 'M7';
    atomsTable[46] := 'M8';
    atomsTable[47] := 'M9';
    atomsTable[48] := 'M0';
    atomsTable[49] := '3';
    atomsTable[50] := '4';
    atomsTable[51] := '5';
    atomsTable[52] := '6';
    atomsTable[53] := '7';
    atomsTable[54] := '8';
    atomsTable[55] := '1';
    atomsTable[56] := '2';
    atomsTable[57] := 'L1';
    atomsTable[58] := 'L2';
    atomsTable[59] := 'L3';
    atomsTable[60] := 'L4';
    atomsTable[61] := 'L5';
    atomsTable[62] := 'L6';
    atomsTable[63] := 'L7';
    atomsTable[64] := 'L8';
    atomsTable[65] := 'L9';
    atomsTable[66] := 'L0';
    atomsTable[67] := 'K1';
    atomsTable[68] := 'K2';
    atomsTable[69] := 'K3';
    atomsTable[70] := 'K4';
    atomsTable[71] := 'M1';
    atomsTable[72] := 'M2';
    atomsTable[73] := 'M3';
    atomsTable[74] := 'M4';
    atomsTable[75] := 'M5';
    atomsTable[76] := 'M6';
    atomsTable[77] := 'M7';
    atomsTable[78] := 'M8';
    atomsTable[79] := 'M9';
    atomsTable[80] := 'M0';
    atomsTable[81] := '3';
    atomsTable[82] := '4';
    atomsTable[83] := '5';
    atomsTable[84] := '6';
    atomsTable[85] := '7';
    atomsTable[86] := '8';
    atomsTable[87] := '1';
    atomsTable[88] := '2';
    atomsTable[89] := 'L1';
    atomsTable[90] := 'L2';
    atomsTable[91] := 'L3';
    atomsTable[92] := 'L4';
    atomsTable[93] := 'L5';
    atomsTable[94] := 'L6';
    atomsTable[95] := 'L7';
    atomsTable[96] := 'L8';
    atomsTable[97] := 'L9';
    atomsTable[98] := 'L0';
    atomsTable[99] := 'K1';
    atomsTable[100] := 'K2';
    atomsTable[101] := 'K3';
    atomsTable[102] := 'K4';
    atomsTable[103] := 'M1';
    atomsTable[104] := 'CD';
    atomsTable[105] := 'CT';
    atomsTable[106] := 'CB';
    atomsTable[107] := 'CA';
    atomsTable[108] := 'CO';
    atomsTable[109] := 'CN';
    atomsTable[110] := 'NI';
    atomsTable[111] := 'NA';
end;

end.

