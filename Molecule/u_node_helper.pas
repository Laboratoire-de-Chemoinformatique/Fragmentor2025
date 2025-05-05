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
unit U_NODE_HELPER;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, U_GRAPHES, U_TYPE_GRAPHES;

const
  margin=50;

type

PTNodePrp = ^TNodePrp;
TNodePrp = record
  Id: Node; //Index of the node in a graph
  Deg: Integer; //Degree of the node in a graph
  Col: array of Node; //Color of the node
end;

{ TListNode }

TListNode = class(TList)
private
  flvl,fcapacity: integer; //The current level to look in tables
  function GetCol(Index: integer): integer;
  function GetDeg(Index: integer): integer;
  function GetId(Index: integer): integer;
  function GetLvl: integer;
  function GetNode(Index: integer): PTNodePrp;
  procedure SetCol(Index: integer; AValue: integer);
  procedure SetDeg(Index: integer; AValue: integer);
  procedure SetId(Index: integer; AValue: integer);
  procedure SetLvl(AValue: integer);
  procedure SetNode(Index: integer; AValue: PTNodePrp);
public
  procedure SetListNode(G: T_GRAPHE_LISTE);
  procedure AddNode(aNode: Node; G: T_GRAPHE_LISTE);
  function FindById(Id: Node): PTNodePrp;
  procedure ClearP;
  property Items[Index: integer]:PTNodePrp read GetNode write SetNode;
  property lvl: integer read GetLvl write SetLvl;
  property Col[Index: integer]: integer read GetCol write SetCol;
  property Id[Index: integer]: integer read GetId write SetId;
  property Deg[Index: integer]: integer read GetDeg write SetDeg;
  constructor Create;
  destructor Destroy; override; //Use ClearP to free stored pointers before destruction of the object.
end;

function sortNodeDeg(pitem1,pitem2: Pointer): integer;
function sortNodeCol(pitem1,pitem2: Pointer): integer;
function sortNodeId(pitem1,pitem2: Pointer): integer;

implementation

function sortNodeDeg(pitem1, pitem2: Pointer): integer;
//Sort as decreasing node degree
//-1 head of list; +1 tail of list
//sorting by decreasing outer degree value
begin
  if PTNodePrp(pitem1)^.Deg>PTNodePrp(pitem2)^.Deg then Result:=-1
  else if PTNodePrp(pitem1)^.Deg<PTNodePrp(pitem2)^.Deg then Result:=1
  else Result:=0;
end;

function sortNodeCol(pitem1, pitem2: Pointer): integer;
//Sort as increasing color at level 0.
//-1 head of list; +1 tail of list
begin
  if PTNodePrp(pitem1)^.Col[0]>PTNodePrp(pitem2)^.Col[0] then Result:=1
  else if PTNodePrp(pitem1)^.Col[0]<PTNodePrp(pitem2)^.Col[0] then Result:=-1
  else Result:=0;
end;

function sortNodeId(pitem1, pitem2: Pointer): integer;
//Sort as increasing Id
//-1 head of list; +1 tail of list
begin
  if PTNodePrp(pitem1)^.Id>PTNodePrp(pitem2)^.Id then Result:=1
  else if PTNodePrp(pitem1)^.Id<PTNodePrp(pitem2)^.Id then Result:=-1
  else Result:=0;
end;

{ TListNode }

function TListNode.GetCol(Index: integer): integer;
begin
  Result:=Items[index]^.Col[flvl];
end;

function TListNode.GetDeg(Index: integer): integer;
begin
  Result:=Items[index]^.Deg;
end;

function TListNode.GetId(Index: integer): integer;
begin
  Result:=Items[index]^.Id;
end;

function TListNode.GetLvl: integer;
begin
  Result:=flvl;
end;

function TListNode.GetNode(Index: integer): PTNodePrp;
begin
  Result:=PTNodePrp(inherited Items[index]);
end;

procedure TListNode.SetCol(Index: integer; AValue: integer);
begin
  Items[index]^.Col[flvl]:=AValue;
end;

procedure TListNode.SetDeg(Index: integer; AValue: integer);
begin
  Items[index]^.Deg:=AValue;
end;

procedure TListNode.SetId(Index: integer; AValue: integer);
begin
  Items[index]^.Id:=AValue;
end;

procedure TListNode.SetLvl(AValue: integer);
var
  i: integer;
  pi: PTNodePrp;
begin
  if AValue>=fcapacity then //Resize Col if needed
  begin
    fcapacity:=AValue+margin;
    for i:=0 to Count-1 do
    begin
      pi:=Items[i];
      SetLength(pi^.Col,fcapacity);
    end;
  end;
  flvl:=AValue;
end;

procedure TListNode.SetNode(Index: integer; AValue: PTNodePrp);
begin
  Items[Index]:=Pointer(AValue);
end;

procedure TListNode.SetListNode(G: T_GRAPHE_LISTE);
var
  i: Node;
begin
  Capacity:=G.p_NX;//R does not necessitate more space than the number of nodes in G
  for i:=1 to G.p_NX do
    AddNode(i,G);
end;

procedure TListNode.AddNode(aNode: Node; G: T_GRAPHE_LISTE);
var
  pitem: PTNodePrp;
begin
  new(pitem);
  pitem^.Id:=aNode;
  pitem^.Deg:=G.p_OutDeg[aNode];
  SetLength(pitem^.Col,fcapacity);
  Add(pitem);
end;

function TListNode.FindById(Id: Node): PTNodePrp;
var
  i: integer;
  bFound: Boolean;
begin
  i:=0;
  bFound:=False;
  while (i<Count) and not bFound do
  begin
    If (PTNodePrp(Items[i])^.Id=Id) then bFound:=True;
    Inc(i);
  end;
  Dec(i);
  if bFound then Result:=PTNodePrp(Items[i]) else Result:=nil;
end;

procedure TListNode.ClearP;
var
  i: integer;
  pitem: PTNodePrp;
begin
  for i:=0 to Count-1 do
  begin
    pitem:=PTNodePrp(Items[i]);
    if pitem<>nil then dispose(pitem);
  end;
  Clear;
end;

constructor TListNode.Create;
begin
  flvl:=0;
  fcapacity:=20;
  inherited Create;
end;

destructor TListNode.Destroy;
begin
  inherited Destroy;
end;

end.

