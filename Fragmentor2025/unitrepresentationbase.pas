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
unit UnitRepresentationBase;
//Unit to define atom and bond representation
{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils; 

const
     MaxAtomSum = 118;

type

    { TRepresentationBase }

    TRepresentationBase=Class(TObject)
    private
           fUseMarkAtom: boolean;
           fUseFormalCharge: boolean;
    public
          property UseMarkAtom: boolean read fUseMarkAtom write fUseMarkAtom;
          property UseFormalCharge: boolean read fUseFormalCharge write fUseFormalCharge;
          constructor Create;
          destructor Destroy; override;
          procedure Clear; virtual;
    end;

implementation

{ TRepresentationBase }

constructor TRepresentationBase.Create;
begin
     inherited Create;
     fUseMarkAtom:=false;
     fUseFormalCharge:=false;
end;

destructor TRepresentationBase.Destroy;
begin
     inherited Destroy;
end;

procedure TRepresentationBase.Clear;
begin
     fUseMarkAtom:=false;
     fUseFormalCharge:=false;
end;

end.

