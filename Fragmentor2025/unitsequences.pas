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
unit UnitSequences;
//Unit to fragment molecules according to sequences of atoms and/or bonds
{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,UnitFragmentBase, UnitMoleculeFrg, unitAtomAndBondType, U_TYPE_GRAPHES;


const
     LAP     = 15;                         // length of atoms path
     MinFL   = 2;                          // min atoms path
     MaxFL   = 15;                         // max atoms path

type

    { TSequences }

    TSequences = class(TFrgBase)
    private
           fLenMin, fLenMax: integer;
    protected
             Start, Finish: integer;
             NAS: integer;
    public
          property LenMin: integer read fLenMin write fLenMin;
          property LenMax: integer read fLenMax write fLenMax;
          constructor Create;
          constructor Create(s: TStringList);
          destructor Destroy; override;
          procedure PathFilter(Mol:TMoleculeFrg); virtual; abstract;
          //procedure MrkAtFlter(Mol: TMoleculeFrg; s, t: Node); virtual; abstract;
          //PathFilter: In progress. Filtering of fragments should be done at a lower level.
          //MrkAtFlter: In progress. Filtering of fragments should be done at a lower level.
    end;

implementation

{ TSequences }

constructor TSequences.Create;
begin
     inherited Create;
end;

constructor TSequences.Create(s: TStringList);
begin
  inherited Create(s);
end;

destructor TSequences.Destroy;
begin
     inherited Destroy;
end;

end.

