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
unit UnitShortestPathPair;
//Unit to fragment based on shortest path sequences pairs
{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, math, unitfragmentbase, UnitMoleculeFrg, UnitMoleculeBase,
  unitAtomAndBondType, U_TYPE_GRAPHES, unitsequences, UnitBondBase;

type

    TShrtPthPairException = class(Exception)
    end;

    { TShortestPathPair }

    TShortestPathPair = class(TSequences)
    private
    protected
             fV: TNodeCost;
             function PathToString(Mol: TMoleculeFrg; s,t: Node): string; virtual; abstract;
    public
          constructor Create;
          destructor Destroy; override;
          procedure MolToFrgLst(Mol: TMoleculeFrg); override;
    end;

implementation

{ TShortestPathPair }

constructor TShortestPathPair.Create;
begin

end;

destructor TShortestPathPair.Destroy;
begin
  inherited Destroy;
end;

procedure TShortestPathPair.MolToFrgLst(Mol: TMoleculeFrg);
var
   i, j: AtomID;
   s, t, u, v: Node;
   P: TNodeInfo;
   SFrg: string;
   ii: integer;
   AW: TArcCost;
   k: ArcNum;
begin
     //Initialization
     fNPath:=0;
     for s:=1 to Mol.nAtom do begin
         fLen[s]:=0;
         for j:=1 to Mol.nAtom do fPath[s,t]:=0;
     end;
     for k:=1 to Mol.p_M do AW[k]:=1;
     //Process each atom pair
     for s:=1 to Mol.nAtom do begin
         if ((not UseMarkAtom) or (UseMarkAtom and Mol.IsMarkedAt(s))) then begin
            for t:=s+1 to Mol.nAtom do
                if ((not UseMarkAtom) or (UseMarkAtom and Mol.IsMarkedAt(t))) then begin
                   Mol.DijHeap(AW,s,t,fV,P);
                   //Path ending at t has a cost containing in fV[t]
                   //Here dynamic bonds restrictions should be tested.
                   SFrg:=PathToString(Mol,s,t);
                   AddFragment(SFrg);
                end;
         end;
         for u:=1 to fNPath do begin
             for v:=1 to fLen[u] do fPath[u,v]:=0;
                 fLen[u]:=0;
         end;
         fNPath:=0;
     end;
end;

end.

