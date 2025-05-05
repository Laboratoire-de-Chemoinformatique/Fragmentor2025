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
unit UnitAtomSymbol;
//Unit to describe atoms as their two letters symbol
{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, UnitRepresentationBase,UnitAtomBase;

type

    { TAtomSymbol }

    TAtomSymbol = class(TAtomBase)
    private

    public
          constructor Create;
          destructor Destroy; override;
          procedure init_atomsTable; override;
    end;
    
implementation

{ TAtomSymbol }

constructor TAtomSymbol.Create;
begin
     inherited Create;
     init_atomsTable;
end;

destructor TAtomSymbol.Destroy;
begin
     inherited Destroy;
end;

procedure TAtomSymbol.init_atomsTable;
begin
     atomsTable[1] := 'H';
     atomsTable[2] := 'He';
     atomsTable[3] := 'Li';
     atomsTable[4] := 'Be';
     atomsTable[5] := 'B';
     atomsTable[6] := 'C';
     atomsTable[7] := 'N';
     atomsTable[8] := 'O';
     atomsTable[9] := 'F';
    atomsTable[10] := 'Ne';
    atomsTable[11] := 'Na';
    atomsTable[12] := 'Mg';
    atomsTable[13] := 'Al';
    atomsTable[14] := 'Si';
    atomsTable[15] := 'P';
    atomsTable[16] := 'S';
    atomsTable[17] := 'Cl';
    atomsTable[18] := 'Ar';
    atomsTable[19] := 'K';
    atomsTable[20] := 'Ca';
    atomsTable[21] := 'Sc';
    atomsTable[22] := 'Ti';
    atomsTable[23] := 'V';
    atomsTable[24] := 'Cr';
    atomsTable[25] := 'Mn';
    atomsTable[26] := 'Fe';
    atomsTable[27] := 'Co';
    atomsTable[28] := 'Ni';
    atomsTable[29] := 'Cu';
    atomsTable[30] := 'Zn';
    atomsTable[31] := 'Ga';
    atomsTable[32] := 'Ge';
    atomsTable[33] := 'As';
    atomsTable[34] := 'Se';
    atomsTable[35] := 'Br';
    atomsTable[36] := 'Kr';
    atomsTable[37] := 'Rb';
    atomsTable[38] := 'Sr';
    atomsTable[39] := 'Y';
    atomsTable[40] := 'Zr';
    atomsTable[41] := 'Nb';
    atomsTable[42] := 'Mo';
    atomsTable[43] := 'Tc';
    atomsTable[44] := 'Ru';
    atomsTable[45] := 'Rh';
    atomsTable[46] := 'Pd';
    atomsTable[47] := 'Ag';
    atomsTable[48] := 'Cd';
    atomsTable[49] := 'In';
    atomsTable[50] := 'Sn';
    atomsTable[51] := 'Sb';
    atomsTable[52] := 'Te';
    atomsTable[53] := 'I';
    atomsTable[54] := 'Xe';
    atomsTable[55] := 'Cs';
    atomsTable[56] := 'Ba';
    atomsTable[57] := 'La';
    atomsTable[58] := 'Ce';
    atomsTable[59] := 'Pr';
    atomsTable[60] := 'Nd';
    atomsTable[61] := 'Pm';
    atomsTable[62] := 'Sm';
    atomsTable[63] := 'Eu';
    atomsTable[64] := 'Gd';
    atomsTable[65] := 'Tb';
    atomsTable[66] := 'Dy';
    atomsTable[67] := 'Ho';
    atomsTable[68] := 'Er';
    atomsTable[69] := 'Tm';
    atomsTable[70] := 'Yb';
    atomsTable[71] := 'Lu';
    atomsTable[72] := 'Hf';
    atomsTable[73] := 'Ta';
    atomsTable[74] := 'W';
    atomsTable[75] := 'Re';
    atomsTable[76] := 'Os';
    atomsTable[77] := 'Ir';
    atomsTable[78] := 'Pt';
    atomsTable[79] := 'Au';
    atomsTable[80] := 'Hg';
    atomsTable[81] := 'Tl';
    atomsTable[82] := 'Pb';
    atomsTable[83] := 'Bi';
    atomsTable[84] := 'Po';
    atomsTable[85] := 'At';
    atomsTable[86] := 'Rn';
    atomsTable[87] := 'Fr';
    atomsTable[88] := 'Ra';
    atomsTable[89] := 'Ac';
    atomsTable[90] := 'Th';
    atomsTable[91] := 'Pa';
    atomsTable[92] := 'U';
    atomsTable[93] := 'Np';
    atomsTable[94] := 'Pu';
    atomsTable[95] := 'Am';
    atomsTable[96] := 'Cm';
    atomsTable[97] := 'Bk';
    atomsTable[98] := 'Cf';
    atomsTable[99] := 'Es';
    atomsTable[100] := 'Fm';
    atomsTable[101] := 'Md';
    atomsTable[102] := 'No';
    atomsTable[103] := 'Lr';
    atomsTable[104] := 'CD';
    atomsTable[105] := 'CT';
    atomsTable[106] := 'CB';
    atomsTable[107] := 'CA';
    atomsTable[108] := 'CO';
    atomsTable[109] := 'CN';
    atomsTable[110] := 'NI';
    atomsTable[111] := 'NA';
    atomsTable[112] := 'LP';
    atomsTable_size := MaxAtomSum;
end;

end.

