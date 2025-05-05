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
unit UnitBondReactions;
//Unit to manage condensed graph of reaction (CGR) bond representation
{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, UnitBondOrder;

type

    { TBondReactions }

    TBondReactions = class(TBondOrder)
    private

    public
      constructor Create;
      destructor Destroy; override;
      procedure init_bondstable;
    end;
    
implementation

{ TBondReactions }

constructor TBondReactions.Create;
begin
     inherited Create;
     init_bondstable;
end;

destructor TBondReactions.Destroy;
begin
  inherited Destroy;
end;

procedure TBondReactions.init_bondstable;
begin
     BondSymbol[8]:='~';
     BondSymbol[9]:='&';
     //
     BondSymbol[81]:='1';  // made single bond;
     BondSymbol[12]:='2';  // made double bond from single bond
     BondSymbol[23]:='3';  // made triple bond from double bond
     BondSymbol[13]:='4';  // made triple bond from single bond
     BondSymbol[18]:='5';  // broken single bond
     BondSymbol[21]:='6';  // broken double bond to single bond
     BondSymbol[28]:='7';  // broken double bond
     BondSymbol[32]:='8';  // broken triple bond to double bond
     BondSymbol[31]:='9';  // broken triple bond to single bond
     BondSymbol[82]:='D';  // broken single bond
     BondSymbol[83]:='T';  // broken single bond
     inherited init_bondstable;
end;

end.

