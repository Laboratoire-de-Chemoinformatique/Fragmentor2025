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
unit UnitBondComplexe;
//Unit to manage representation of non covalent bonds
{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, UnitBondOrder;

type

    { TBondComplexe }

    TBondComplexe = class(TBondOrder)
    private

    public
      constructor Create;
      destructor Destroy; override;
      procedure init_bondstable; override;
    end;

implementation

{ TBondComplexe }

constructor TBondComplexe.Create;
begin
     inherited Create;
     init_bondstable;
end;

destructor TBondComplexe.Destroy;
begin
  inherited Destroy;
end;

procedure TBondComplexe.init_bondstable;
begin
     BondSymbol[8]:='~';
     BondSymbol[9]:='&';
     inherited init_bondstable;
end;

end.

