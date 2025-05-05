{ Fragmentor of the ISIDA Project

  Copyright (C) 2024 Laboratoire de Chemoinformatique, UMR 7140 CNRS (https://complex-matter.unistra.fr/equipes-de-recherche/laboratoire-de-chemoinformatique/home/)
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
unit UnitMoleculeBase;

{$mode objfpc}{$H+}
//{$mode DELPHI}{$H+}

interface

uses
  Classes, SysUtils, U_BASE_GRAPHES, U_TYPE_GRAPHES, U_GRAPHES, U_NODE_HELPER,
  unitAtomAndBondType, contnrs;

const
  MaxAtom = MaxNode;
  MaxBond = MaxArcNum div 2;

type
  AtomID = Node;
  BondID = ArcNum;

  AAtmSet = array[AtomID] of PRAtom;
  ABndSet = array[BondID] of PRBond;
  //APrpSet = array[string] of string;

  AAoB = array of array of boolean;

  PRSAL=^RSAL;
  RSAL = record
    etype: EDynWrd;
    asize: byte;
    alist: array[1..2] of AtomID;
    aword: string;
  end;

  { TMoleculeBase }

  TMoleculeBase = class(T_GRAPHE_LISTE)
  private
    fMolName: string;  //molecule identifier or name
    fAtmSet: AAtmSet;
    // array containing the ID of atoms: AtomID constituting the molecule
    fBndSet: ABndSet;  // array containing the ID of bonds: BondID
    //fLBndSet: TList; //In progress, a bond pointer list to manage bonds
    //fPrpSet: APrpSet;
    fAPrpSze, fBPrpSze: integer;
    fABytSze, fBBytSze: integer;
    //freadpM: integer;
    //freadpNX: integer;
    //readpM: integer;
    //readpNX: integer;
    procedure DisposeAtmSet;
    procedure DisposeBndSet;
    function GetAtom(i: AtomID): PRAtom;
    function GetI(i: AtomID): AByt;
    function GetP(i: AtomID): APrp;
    function GetS(i: AtomID): TS;
    function GetW(i: AtomID): double;
    function GetZ(i: AtomID): TZ;
    procedure SetAtom(i: AtomID; PAt: PRAtom);
    function GetBond(i: BondID): PRBond;
    procedure SetBond(i: BondID; PBo: PRBond);
    function readpNX(): integer;
    function readpM(): integer;
    procedure SetS(i: AtomID; const AValue: TS);
    procedure SetZ(i: AtomID; const AValue: TZ);
    procedure SetW(i: AtomID; const AValue: double);
    procedure writepNX(i: integer);
    procedure writepM(i: integer);
  protected
    function int_readpos(str: string; bgn, lng: integer): integer;
    function dbl_readpos(str: string; bgn, lng: integer): double;
  public
    constructor Create;
    destructor Destroy; override;
    procedure Clear; virtual;
    procedure LoadSDF(sdfstr: TStringList); virtual;
    procedure LoadSDFString(sdfstr: string); virtual;
    function LoadSDFTStringList(sdfstr: TStringList;
      prpstr: string): TFPStringHashTable; virtual;
    function LoadSDFField(sdfstr: TStringList;
      prpstr: string): TFPStringHashTable; virtual;
    function LoadSDFField(sdfstr: TStringList;
      prpstr: TStringList): TFPStringHashTable; virtual;
    procedure LoadSDFField(sdfstr: TStringList; prpstr: string;
      var strhash: TFPStringHashTable); virtual;
    procedure LoadSDFField(sdfstr: TStringList; prpstr: TStringList;
      strhash: TFPStringHashTable); virtual;
    function LoadSDFSGroup(sdfstr:TStringList): TList; virtual;
    procedure LoadSDFSGroup(sdfstr:TStringList; var datalist: TList); virtual;
    function RemoveBondDir(PBo: PRBond): boolean;
    function RemoveBondDir(s, t: AtomID): boolean;
    function RemoveBond(s, t: AtomID): boolean;
    function RemoveBond(m: BondID): boolean;
    function convert_dict_hash(str: string): TFPStringHashTable;
    function convert_dict_list(str: string): TStringList;
    procedure RemoveAtom(x: AtomID);
    function AddAtomPT: PRAtom;
    function AddAtomID: AtomID;
    function AddBondDir(s, t: AtomID): PRBond;
    procedure AddBondPT(s, t: AtomID; var b1, b2: PRBond);
    procedure AddBondID(s, t: AtomID; var b1, b2: BondID);
    procedure SetAPrpSize(size: integer);
    procedure SetBPrpSize(size: integer);
    procedure SetABytSize(size: integer);
    procedure SetBBytSize(size: integer);
    function FindAtom(PAt: PRAtom): AtomID;
    function FindBond(PBo: PRBond): BondID;
    function FindBond(t, h: AtomID): PRBond;
    function FindBondID(t, h: AtomID): BondID;
    property nAtom: integer read readpNX write writepNX;
    property nBonds: integer read readpM write writepM;

    property MolName: string read fMolName write fMolName;
    property APrpSze: integer read fAPrpSze write fAPrpSze;
    property ABytSze: integer read fABytSze write fABytSze;
    property BPrpSze: integer read fBPrpSze write fBPrpSze;
    property BBytSze: integer read fBBytSze write fBBytSze;
    property AtmSet[i: AtomID]: PRAtom read GetAtom write SetAtom;
    property BndSet[i: BondID]: PRBond read GetBond write SetBond;
    property Z_[i: AtomID]: TZ read GetZ write SetZ;
    property S_[i: AtomID]: TS read GetS write SetS;
    property W_[i: AtomID]: double read GetW write SetW;
    property P_[i: AtomID]: APrp read GetP;
    property I_[i: AtomID]: AByt read GetI;
    //function IntToTB(idx: integer): TB;
    function IntToTB(col3, col7: integer): TB;
    function IntToS(col3, col7: integer): TS;
  end;

function CountOn(BV: TBits): integer;
function searchClique(graph: AAoB; sizeGraph: AtomID; var maxClique: TNodeInfo): AtomID;
procedure recursiveSearchClique(graph: AAoB; clique: TNodeInfo;
  sizeClique: AtomID; candidates: TNodeInfo; numCandidates: AtomID;
  var maxClique: TNodeInfo ; var maxSizeClique: AtomID);
procedure ResultInSdf(mol: TMoleculeBase ;First: TNodeInfo;szeMCS: AtomID;  var fileout: text);
procedure ResultInSdf(mol: TMoleculeBase ;First: TNodeInfo;szeMCS: AtomID;  var fileout: string);
//function MCS(mol1, mol2: TMoleculeBase; out First, Second: TNodeInfo): AtomID;
function MCSslow(mol1, mol2: TmoleculeBase; bUseDist: boolean; out First, Second: TNodeInfo): AtomID;
function MCS(mol1, mol2: TmoleculeBase; bUseDist: boolean; out First, Second: TNodeInfo): AtomID;
function mMCS(mol1, mol2: TmoleculeBase; bUseDist: boolean; out First, Second: TNodeInfo; QmaxDelta:ShortInt=0): AtomID;
procedure GetSSList(mol1, mol2: TmoleculeBase; minSze, maxSze: integer; bUseDist: boolean; out SSList: aiarray); // Search for list of common substructures in mol1 and mol2 of size between minSze nd maxSze

implementation

function CountOn(BV: TBits): integer;
var
  i: integer;
begin
  Result := 0;
  for i := 0 to BV.Size - 1 do
    if (BV.Bits[i]) then
      Inc(Result);
end;

function searchClique(graph: AAoB; sizeGraph: AtomID; var maxClique: TNodeInfo): AtomID;
//Search in the compatibility graph a maximal clique. Can be replaced by en EulerPath algorithm
//Input:
//- graph : a compatibility matrix; it is the setup of the search
//- sizeGraph: the size of the compatibility graph
//Output:
// - maxClique: the list of nodes in the clique, i.e. the list of atom pairs in the MCS
// - Result: the clique size, i.e. the MCS size
var
  i, j :integer;
  sizeClique, numCandidates, maxSizeClique: AtomID;
  clique, candidates: TNodeInfo;
  //dynDef:Boolean;
  //début calcul de la taille de la soustructure maximum commune
begin
  //initialization
  for i:=Low(clique) to High(clique) do clique[i]:=0;
  for i:=Low(candidates) to High(candidates) do candidates[i]:=0;
  maxSizeClique := 1;
  for i := 1 to sizeGraph do // on sélectionne un atome dans la molécule 1 ayant un atome image dans 2
  begin
    // on initialise les valeurs
    //dynDef:=False;
    clique[1] := i;
    sizeClique := 1;
    numCandidates := 0;
    if (i < sizeGraph) then
    begin
      //on liste tous les seconds atomes ayant une image compatible( true dans compatgraph)
      for j := i + 1 to sizeGraph do begin
        if (graph[i, j]) then
        begin
          //a un atomes i dans candidates, on a tous les atomes j compatible
          Inc(numCandidates);
          candidates[numCandidates] := j;
          // ici, on a regardé chaque atomes ayant une image et leur voisin ayant une image à un rang
          //on procède au recursiveSearchClique afin d'augmenter le rang jusqu'à atteindre le rang max
          if (sizeClique + numCandidates > maxSizeClique) then
            recursiveSearchClique(graph, clique, sizeClique,
              candidates, numCandidates,
              maxClique, maxSizeClique);
        end;
      end;
    end;
  end;
  //on renvoie la valeur de la taille de la soustructure max commune
  Result := maxSizeClique;
end;

procedure recursiveSearchClique(graph: AAoB; clique: TNodeInfo;
  sizeClique: AtomID; candidates: TNodeInfo; numCandidates: AtomID;
  var maxClique: TNodeInfo; var maxSizeClique: AtomID);
//From a clique and a set of candidate nodes, create a larger clique
//Input:
//- graph: the graph in which the clique is searched
//- clique: the current clique
//- sizeClique: the current clique size
//- candidates: the list of nodes that can increase the clique
//- numCandidates: the size of the list of candidate nodes
//Output:
//- maxClique: the largest clique
//- maxSizeClique: the size of the largest clique
var
  i, j, newNumCandidates: AtomID;
  newCandidates: TNodeInfo;
begin
  if (numCandidates=0)then
  begin // a clique found; if multiple clique are search -> this is the place to extract it
    if(sizeClique>maxSizeClique)then //it is the largest clique
    begin
      maxSizeClique:= sizeClique;
      for i:=1 to sizeClique do //Store the clique in maxClique
        maxClique[i]:= clique[i];
    end;
  end
  else
  begin
    Inc(sizeClique);
    for i := 1 to numCandidates do
    begin
      //The next node of the list is included to the clique
      clique[sizeClique] := candidates[i];
      newNumCandidates := 0;
      //From this node store the candidate nodes to increase the clique
      if (i < numCandidates) then
      begin
      for j := i + 1 to numCandidates do
        begin
          if (graph[candidates[i],candidates[j]]) then
          begin
            Inc(newNumCandidates);
            newCandidates[newNumCandidates] := candidates[j];
          end;
        end;
      end;
      if (sizeClique + newNumCandidates > maxSizeClique) then
        //Search using the new clique and set of candidates; it stops if the clique cannot increase.
        recursiveSearchClique(graph, clique, sizeClique,
          newCandidates, newNumCandidates,
          maxClique, maxSizeClique);
    end;
  end;
end;

procedure ResultInSdf(mol:TMoleculeBase;First:TNodeInfo; szeMCS: AtomID; var fileout: text);
var
   i,j:integer;
   mMCS:TMoleculeBase;
   Atmp:PRAtom;
   Btmp:PRBond;
   BdID,bd1,bd2:BondID;
   fmt,S:string;
begin
  mMCS:=TMoleculeBase.Create;
  bd1:=mol.nBonds; bd2:=mol.nBonds;
  //setAtom
  for i:=1 to szeMCS do
  begin
    Atmp:=mol.GetAtom(First[i]) ;
    mMCS.AddAtomID; //This adds an atom pointer to the molecule
    mMCS.SetAtom(i,Atmp); //Copy atom content
  end;
  //set Bond
  BdID:=0;
  for i:= 1 to szeMCS-1 do
    for j:= i+1 to szeMCS do
    begin
      Btmp:=mol.FindBond(First[i],First[j]);
      if(Btmp<>nil)then
      begin
        Inc(BdID);
        //changement des indices d'atomes1
        Btmp^.t:=i;
        Btmp^.h:=j;
        mMCS.AddBondID(i,j,bd1,bd2);
        mMCS.SetBond(BdID,Btmp);
      end;
    end;
  //write sdf
  //header
  writeln(fileout);
  writeln(fileout,'  MaxCommonSubstructure');
  writeln(fileout);
  //counts line
  writeln(fileout, IntToStr(szeMCS) :3,IntToStr(BdID):3,0:3,0:3,0:3,0:3, '':12,999:3, 'V2000':6);
  //Atom Block
  for i:=1 to szeMCS do begin
    //formatting Atom's name part
    Fmt:='%0:-3s';
    S:=Format(Fmt,[mMCS.S_[i]]);
    writeln(fileout,0.0:10:4,0.0:10:4,0.0:10:4,' ',S,0:2,0:3,0:3,0:3,0:3,0:3,0:3,0:3,0:3,IntToStr(First[i]):3,0:3,0:3);
  end;
  //Bond Block
  for i:= 1 to BdID do begin
    Btmp:=mMCS.GetBond(i);
    case Btmp^.B of
      14: writeln(fileout,Btmp^.t :3,Btmp^.h:3,81:3,0:3,0:3,0:3,4:3);
      15: writeln(fileout,Btmp^.t :3,Btmp^.h:3,82:3,0:3,0:3,0:3,4:3);
      16: writeln(fileout,Btmp^.t :3,Btmp^.h:3,83:3,0:3,0:3,0:3,4:3);
      17: writeln(fileout,Btmp^.t :3,Btmp^.h:3,84:3,0:3,0:3,0:3,4:3);
      18: writeln(fileout,Btmp^.t :3,Btmp^.h:3,18:3,0:3,0:3,0:3,4:3);
      19: writeln(fileout,Btmp^.t :3,Btmp^.h:3,28:3,0:3,0:3,0:3,4:3);
      20: writeln(fileout,Btmp^.t :3,Btmp^.h:3,38:3,0:3,0:3,0:3,4:3);
      21: writeln(fileout,Btmp^.t :3,Btmp^.h:3,48:3,0:3,0:3,0:3,4:3);
      22: writeln(fileout,Btmp^.t :3,Btmp^.h:3,12:3,0:3,0:3,0:3,8:3);
      23: writeln(fileout,Btmp^.t :3,Btmp^.h:3,13:3,0:3,0:3,0:3,8:3);
      24: writeln(fileout,Btmp^.t :3,Btmp^.h:3,14:3,0:3,0:3,0:3,8:3);
      25: writeln(fileout,Btmp^.t :3,Btmp^.h:3,21:3,0:3,0:3,0:3,8:3);
      26: writeln(fileout,Btmp^.t :3,Btmp^.h:3,23:3,0:3,0:3,0:3,8:3);
      27: writeln(fileout,Btmp^.t :3,Btmp^.h:3,24:3,0:3,0:3,0:3,8:3);
      28: writeln(fileout,Btmp^.t :3,Btmp^.h:3,31:3,0:3,0:3,0:3,8:3);
      29: writeln(fileout,Btmp^.t :3,Btmp^.h:3,32:3,0:3,0:3,0:3,8:3);
      30: writeln(fileout,Btmp^.t :3,Btmp^.h:3,34:3,0:3,0:3,0:3,8:3);
      31: writeln(fileout,Btmp^.t :3,Btmp^.h:3,41:3,0:3,0:3,0:3,8:3);
      32: writeln(fileout,Btmp^.t :3,Btmp^.h:3,42:3,0:3,0:3,0:3,8:3);
      33: writeln(fileout,Btmp^.t :3,Btmp^.h:3,43:3,0:3,0:3,0:3,8:3);
      else writeln(fileout,Btmp^.t :3,Btmp^.h:3,Btmp^.B:3,0:3,0:3,0:3,0:3);
    end;
  end;
  //End block
  writeln(fileout,'M  END');
  writeln(fileout);
  writeln(fileout,'$$$$');
  //
  FreeAndNil(mMCS);
end;

// méthode 1: sans matrice de distance
{%H-}{function MCS(mol1, mol2: TMoleculeBase; out First, Second: TNodeInfo): AtomID;
var
  i, j:integer;
  compatGraphSize,sizeMCS: AtomID;
  compatGraph: AAoB;
  index1, index2,maxClique: TNodeInfo;
  i1,j1: integer;

  function calcCompatibilityGraphSize(mol1, mol2: TMoleculeBase): AtomID;
  var
    compat_graph_size: AtomID;
    i, j: integer;
    //début du calcul du nombre de pairs d'atomes commun entre les deux molécules
  begin
    //Initialisation
    compat_graph_size := 0;
    for i := 1 to mol1.nAtom do
      for j := 1 to mol2.nAtom do
        //Penser à utiliser une fonction générique de comparaison de noeuds
        //Comparaison des atomes présents dans chaque molécules avec p_NX le nombre d'atomes dans chaques molécules
        if mol1.S_[i] = mol2.S_[j] then
          //Incrémente à chaque fois que 2 atomes similaires
          Inc(compat_graph_size);
    //Renvoie du résultat
    Result := compat_graph_size;
  end;

  procedure makeCompatibilityGraph(mol1, mol2: TMoleculeBase;
     compatGraph: AAoB; var index1, index2: TNodeInfo);
  //sans matrice de distance
  var
    i, j: integer;
    PRBnd1, PRBnd2: PRBond;
    compat_graph_size: AtomID;
    BondCompat: boolean;
  begin
    compat_graph_size := 0;
    writeln('size1: ' + IntToStr(mol1.nAtom) + '  size2: ' + IntToStr(mol2.nAtom));
    //Début génération du graph de compatibilité
    for i := 1 to mol1.nAtom do
      for j := 1 to mol2.nAtom do
        //Penser à généraliser par une fonction d'égalité portant sur les noeuds
        //Comparaison des atomes présents dans chaque molécules avec nAtom le nombre d'atomes dans chaques molécules
        if mol1.S_[i] = mol2.S_[j] then
          //Si deux atomes commun, on incrémente la taille du graph de compatibilité et on garde leur numero en index
        begin
          Inc(compat_graph_size);
          index1[compat_graph_size] := i;
          index2[compat_graph_size] := j;
        end;

    // On regarde les pair d'atomes : 1atome dans mol1 et son image dans mol2
    for i := 1 to compat_graph_size - 1 do
      //on prend une seconde pair
      for j := i + 1 to compat_graph_size do
      begin
        PRBnd1:=mol1.FindBond(index1[i],index1[j]);
        PRBnd2:=mol2.FindBond(index2[i],index2[j]);
        BondCompat:=False;
        if (index1[i] <> index1[j]) and (index2[i] <> index2[j]) then
        begin
             if ((PRBnd1=nil) and (PRBnd2=nil)) then
                BondCompat:=True
             else if ((PRBnd1<>nil) and (PRBnd2<>nil)) then
                  if (PRBnd1^.B=PRBnd2^.B) then
                     BondCompat:=True;
             if BondCompat then
             begin
                  compatGraph[i, j] := True;
                  compatGraph[j, i] := True;
             end
             else begin
                  compatGraph[i, j] := False;
                  compatGraph[j, i] := False;
             end;
        end;
      end;
  end;
  //Début programme baskin
begin
  //initialisation sizeMCS
  sizeMCS:= 0 ;
  //Calcule du nombre de pair d'atomes commun entre les deux molécules
  compatGraphSize := calcCompatibilityGraphSize(mol1, mol2);
  //verifie que au moins 2 atomes commun
  if compatGraphSize <> 0 then
  begin
    writeln('compatGraphSize= ' + IntToStr(compatGraphSize));
    //définie la taille de compat graph, indexé par nb de pair potentiel d'atomes
    SetLength(compatGraph, compatGraphSize + 1, compatGraphSize + 1);
    //génération du graph de compatibilité
    makeCompatibilityGraph(mol1, mol2, compatGraph, index1, index2);
    {for i1:=1 to compatGraphSize do
    begin
      for j1:=1 to compatGraphSize do
          write(BoolToStr(compatGraph[i1,j1],'T','.')+' ');
      writeln(' % '+IntToStr(index1[i1])+' '+IntToStr(index2[i1]));;
    end;}
    //calcul de la taille de la soustructure maximum commune
    sizeMCS := searchClique(compatGraph, compatGraphSize, maxClique);
    //writeln('SizeMCS: '+IntToStr(SizeMCS));
    // on verifie que les 2 molécules ont une chaine de plus d'un atome en commun
      if(sizeMCS>1)then begin
      for i := 1 to sizeMCS do
      begin
        //on donne les numero d'atomes communs
        First[i] := index1[maxclique[i]];
        Second[i] := index2[maxClique[i]];
      end
    end
    // on donne le numero de l'atome commun
    else
      for i := 1 to sizeMCS do
      begin
        First[i] := index1[1];
        Second[i] := index2[1];
      end;
  end;
  //on renvoie la taille de la sousstruc commune
  Result := sizeMCS;
end;
}

procedure ResultInSdf(mol: TMoleculeBase; First: TNodeInfo; szeMCS: AtomID;
  var fileout: string);
var
   SDFText: TextFile;
begin
  AssignFile(SDFText,fileout);
  Rewrite(SDFText);
  ResultInSdf(mol,First,szeMCS,SDFText);
  CloseFile(SDFText);
end;

//méthode 2: avec matrice de distance
function MCSslow(mol1, mol2: TmoleculeBase; bUseDist: boolean;
  out First, Second: TNodeInfo): AtomID;
var
  i, sizeMCS, compatGraphSize: integer;
  //compatGraph: AAoB;
  compatGraph,dynCompatGraph: AAoB;
  index1, index2, maxClique: TNodeInfo;
  {i1,j1: integer;}

  function DijHeapTopDist(G: TMoleculeBase; s, t: Node; Vmax: Cost = 20): Cost;
  //Compute topological distance between s and t from graph G; do not consider distances above Vmax bonds
  //Input:
  //G, a chemical structure
  //s: starting node
  //t: terminal node
  //Vmax: maximum distance allowed
  //Output: the distance between s and t if less than Vmax, Vmax otherwise
  Var x,y: Node;
      H  : T_Heap;
      k  : ArcNum;
      V  : TNodeCost;
  Begin
    H:=T_Heap.CREATE;
    with G do
    begin
      For x := 1 to p_NX+p_NY do begin
        V[x] := MaxCost
      End;
      H.ClearHeap;
      V[s] := 0;
      H.HeapInsert (V,s);
      Repeat
        H.HeapMin (V,x);
        For k := p_Head[x] to p_Head[x+1]-1 do begin
          y := p_Succ[k];
          If (V[x] + 1 < V[y]) and (V[x] + 1 <= Vmax) then begin
            V[y] := V[x] + 1;
            If H.InHeap (y) then H.MoveUp (V,y) else H.HeapInsert (V,y)
          End
        End;
      Until (H.HeapIsEmpty) or (x=t);
    end;
    FreeAndNil(H);
  End;
  function calcCompatibilityGraphSize(mol1, mol2: TMoleculeBase): AtomID;
  //Input: molecule1, molecule2; output: the size of the compatibility graph
  //Compute the size of a compatibility graph.
  //Atoms are compatible if they are of same element
  var
    compat_graph_size: integer;
    i, j: AtomID;
  begin
    compat_graph_size := 0;
    for i := 1 to mol1.p_NX do
      for j := 1 to mol2.p_NX do
        //Penser à utiliser une fonction générique de comparaison de noeuds
        if mol1.S_[i] = mol2.S_[j] then
          Inc(compat_graph_size);
    Result := compat_graph_size;
  end;
//avec matrice de distance
  procedure makeCompatibilityGraph(mol1, mol2: TMoleculeBase;
    bUseDist: boolean; compatGraph: AAoB; var index1, index2: TNodeInfo);
  //Input:
  //mol1, mol2: two molecules
  //bUseDist: use a distance matrix (True/False)
  //Output:
  //compatGraph: the compatibility graph as a boolean bidimensional array
  //index1, index2: list of atoms in mol1 and mol2
  //Compute the compatibility graph
  var
    i, j: AtomID;
    PRBnd1, PRBnd2: PRBond;
    compat_graph_size: integer;
    BondCompat: boolean;
    dist1, dist2: Cost;
    W1, W2: TArcCost;
    V: TNodeCost;//Distance between nodes
    P: TNodeInfo;//Shortest path between nodes (not used)
    bDistCompat: boolean;
  begin
    //Initialize graph as empty
    compat_graph_size := 0;
    //Compute the list of compatible atoms in mol1 and mol2
    for i := 1 to mol1.nAtom do
      for j := 1 to mol2.nAtom do
        //Penser à généraliser par une fonction d'égalité portant sur les noeuds
        if mol1.S_[i] = mol2.S_[j] then
        begin
          Inc(compat_graph_size);
          index1[compat_graph_size] := i;
          index2[compat_graph_size] := j;
        end;
    //Set the cost of edges to 1 in mol1 and mol2 to compute topological distances
    if bUseDist then
    begin
      for i:=1 to mol1.p_M do W1[i]:=1;
      for i:=1 to mol2.p_M do W2[i]:=1;
    end;
    // On regarde les pair d'atomes : 1atome dans mol1 et son image dans mol2
    for i := 1 to compat_graph_size - 1 do begin //an atom in mol1
      for j := i + 1 to compat_graph_size do //the corresponding compatible atom in mol2
      begin
        PRBnd1:=mol1.FindBond(index1[i],index1[j]);
        PRBnd2:=mol2.FindBond(index2[i],index2[j]);
        BondCompat:=False;
        if bUseDist then
        begin
          //If G1 is subgraph of G2, then the distances in G1 are preserved in G2
          //This is an efficient filter if the distance calculation overhead is
          //compensated by the decreased complexity of the problem
          dist1:=DijHeapTopDist(mol1,index1[i],index1[j]);
          dist2:=DijHeapTopDist(mol1,index1[i],index1[j]);
          bDistCompat:=(dist1=dist2);
        end else
          bDistCompat:=True;
        //
        if (index1[i] <> index1[j]) and (index2[i] <> index2[j]) and bDistCompat then
        begin
             if ((PRBnd1=nil) and (PRBnd2=nil)) then
                BondCompat:=True //atoms i and j are not bonded in both mol1 and mol2
             else if ((PRBnd1<>nil) and (PRBnd2<>nil)) then
                  if (PRBnd1^.B=PRBnd2^.B) then begin
                     BondCompat:=True; //atoms i and j are bonded by the same kind of bonds in mol1 and mol2
                  end;
             if BondCompat then begin //record the compatibility in the matrix
                  compatGraph[i, j] := True;
                  compatGraph[j, i] := True;
             end
             else begin
                  compatGraph[i, j] := False;
                  compatGraph[j, i] := False;
             end;
        end else
        begin //atoms i and j are not part of the same common substructure, they are excluded from the compatibility matrix
          compatGraph[i,j]:=False;
          compatGraph[j,i]:=False;
        end;
      end;
    end;
  end;

begin
  sizeMCS := 0;
  for i:=Low(maxClique) to High(maxClique) do
  begin
    maxClique[i]:=0;
    index1[i]:=0;
    index2[i]:=0;
  end;
  compatGraphSize := calcCompatibilityGraphSize(mol1, mol2);
  if compatGraphSize <> 0 then
  begin
    SetLength(compatGraph, compatGraphSize + 1, compatGraphSize + 1);
    makeCompatibilityGraph(mol1, mol2, bUseDist,
      compatGraph, index1, index2);
    sizeMCS := searchClique(compatGraph, compatGraphSize, maxClique);
    if (sizeMCS > 1) then
      for i := 1 to sizeMCS do
      begin
        First[i] := index1[maxclique[i]];
        Second[i] := index2[maxClique[i]];
      end
    else begin
      First[1] := index1[1];
      Second[1] := index2[1];
    end;
  end;
  Result := sizeMCS;
end;

function MCS(mol1, mol2: TmoleculeBase; bUseDist: boolean; out First,
  Second: TNodeInfo): AtomID;
var
  i,j,k, sizeMCS, compatGraphSize: integer;
  G: T_GRAPHE_LISTE; //compatibility graph
  Rref, R, Q, Qmax: TListNode;//For maximum clique search. Rref, R: list of nodes; Q: current clique; Qmax: maximum clique
  //Rngh: TListNode;
  index1, index2: TNodeInfo;
  SLog: TStringList;
  C: TNodeInfo;
  iq: integer;
  Sclq: TStringList;
  stmp: string;
  maxdeg: Node;
  Colors, LRngh: TObjectList;
  aPNode: PTNodePrp;
  ix: integer;
  px: PTNodePrp;
  dpth,niter: integer;

  function DijHeapTopDist(Gloc: TMoleculeBase; s, t: Node; Vmax: Cost = 20): Cost;
  //Compute topological distance between s and t from graph G; do not consider distances above Vmax bonds
  //Input:
  //G, a chemical structure
  //s: starting node
  //t: terminal node
  //Vmax: maximum distance allowed
  //Output: the distance between s and t if less than Vmax, Vmax otherwise
  Var x,y: Node;
      H  : T_Heap;
      k  : ArcNum;
      V  : TNodeCost;
  Begin
    H:=T_Heap.CREATE;
    with Gloc do
    begin
      For x := 1 to p_NX+p_NY do begin
        V[x] := MaxCost
      End;
      H.ClearHeap;
      V[s] := 0;
      H.HeapInsert (V,s);
      Repeat
        H.HeapMin (V,x);
        For k := p_Head[x] to p_Head[x+1]-1 do begin
          y := p_Succ[k];
          If (V[x] + 1 < V[y]) and (V[x] + 1 <= Vmax) then begin
            V[y] := V[x] + 1;
            If H.InHeap (y) then H.MoveUp (V,y) else H.HeapInsert (V,y)
          End
        End;
      Until (H.HeapIsEmpty) or (x=t);
    end;
    FreeAndNil(H);
    Result:=V[t];
  End;
  function calcCompatibilityGraphSize(mol1, mol2: TMoleculeBase; var index1, index2: TNodeInfo): AtomID;
  //Input: molecule1, molecule2; output: the size of the compatibility graph
  //Compute the size of a compatibility graph.
  //Atoms are compatible if they are of same element
  var
    compat_graph_size: integer;
    i, j: AtomID;
  begin
    compat_graph_size := 0;
    for i := 1 to mol1.nAtom do
      for j := 1 to mol2.nAtom do
        //Penser à utiliser une fonction générique de comparaison de noeuds
        if mol1.S_[i] = mol2.S_[j] then begin
          Inc(compat_graph_size);
          index1[compat_graph_size] := i;
          index2[compat_graph_size] := j;
        end;
    Result := compat_graph_size;
  end;
//avec matrice de distance
  procedure makeCompatibilityGraph(mol1, mol2: TMoleculeBase;
    bUseDist: boolean; var compatGraph: T_GRAPHE_LISTE; var index1, index2: TNodeInfo);
  //Input:
  //mol1, mol2: two molecules
  //bUseDist: use a distance matrix (True/False)
  //Output:
  //compatGraph: the compatibility graph as a boolean bidimensional array
  //index1, index2: list of atoms in mol1 and mol2
  //Compute the compatibility graph
  var
    i, j: AtomID;
    PRBnd1, PRBnd2: PRBond;
    compat_graph_size: integer;
    BondCompat: boolean;
    dist1, dist2: Cost;
    W, W1, W2: TArcCost;
    V: TNodeCost;//Distance between nodes
    P: TNodeInfo;//Shortest path between nodes (not used)
    bDistCompat: boolean;
    tmpGraphMatrix: T_GRAPHE_MATRICIEL;
  begin
    tmpGraphMatrix:=T_GRAPHE_MATRICIEL.CREATE;
    compat_graph_size:=calcCompatibilityGraphSize(mol1,mol2,index1,index2);
    tmpGraphMatrix.p_NX:=compat_graph_size; tmpGraphMatrix.p_NY:=0; tmpGraphMatrix.p_NoArc:=0;
    tmpGraphMatrix.p_Simple:=True;
    if compat_graph_size=0 then Exit;
    //Set the cost of edges to 1 in mol1 and mol2 to compute topological distances
    if bUseDist then
    begin
      for i:=1 to mol1.nAtom do W1[i]:=1;
      for i:=1 to mol2.nAtom do W2[i]:=1;
    end;
    // On regarde les pair d'atomes : 1atome dans mol1 et son image dans mol2
    for i := 1 to compat_graph_size - 1 do begin //an atom in mol1
      for j := i + 1 to compat_graph_size do //the corresponding compatible atom in mol2
      begin
        if (index1[i] <> index1[j]) and (index2[i] <> index2[j]) then
        begin
          PRBnd1:=mol1.FindBond(index1[i],index1[j]);
          PRBnd2:=mol2.FindBond(index2[i],index2[j]);
          BondCompat:=False;
          if bUseDist then
          begin
            //If G1 is subgraph of G2, then the distances in G1 are preserved in G2
            //This is an efficient filter if the distance calculation overhead is
            //not compensated by the decreased complexity of the problem
            dist1:=DijHeapTopDist(mol1,index1[i],index1[j]);
            dist2:=DijHeapTopDist(mol2,index2[i],index2[j]);
            //writeln('i=',index1[i],'/',index2[i],' j=',index1[j],'/',index2[j],' dist1=',dist1,' dist2=',dist2);
            bDistCompat:=(dist1=dist2);
          end else
            bDistCompat:=True;
          //
          if bDistCompat then
          begin
               if ((PRBnd1=nil) and (PRBnd2=nil)) then
                  BondCompat:=True //atoms i and j are not bonded in both mol1 and mol2; indispensable condition
               else if ((PRBnd1<>nil) and (PRBnd2<>nil)) then
                    if (PRBnd1^.B=PRBnd2^.B) then begin
                       BondCompat:=True; //atoms i and j are bonded by the same kind of bonds in mol1 and mol2
                    end;
               if BondCompat then begin //record the compatibility in the matrix
                    tmpGraphMatrix.p_A[i, j] := 1;
                    tmpGraphMatrix.p_A[j, i] := 1;
               end else begin
                 tmpGraphMatrix.p_A[i,j]:=tmpGraphMatrix.p_NoArc;
                 tmpGraphMatrix.p_A[j,i]:=tmpGraphMatrix.p_NoArc;
               end;
          end;
        end;
      end;
    end;

    compatGraph:=tmpGraphMatrix.PackC(nil);
    FreeAndNil(tmpGraphMatrix);
  end;

  function ApproximateColor(Rloc: TListNode): Node;
  var
    nx, nc: Node;
    Cnum, Cmax: Node;
    ir, ic, Qmin: integer;
    CList,Rngh: TListNode;
    bSameColor: boolean;
    px, py: PTNodePrp;
    ix, iy, xy: integer;
  begin
    Cmax:=0;  //maxno
    Qmin:=Qmax.Count-Q.Count+1; //min_k
    TListNode(Colors[1]).Clear;
    {for nc:=0 to Colors.Count-1 do
    begin
      CList:=Colors[nc] as TListNode;           // C
      if CList<>nil then CList.Clear;
    end;}
    ir:=0; //j
    Cnum:=1; //k
    for ix:=0 to Rloc.Count-1 do
    begin
      px:=Rloc[ix];
      Rngh:=LRngh[px^.Id] as TListNode;
      bSameColor:=True; //bSameColor=True if a neighbor y of x is in the color set nc
      nc:=1;
      while bSameColor and (nc<=Cmax) do
      begin
        CList:=Colors[nc] as TListNode;
        iy:=0;
        bSameColor:=False;
        while (iy<CList.Count) and (not bSameColor) do begin
          if Rngh.IndexOf(CList[iy])>=0 then bSameColor:=True;
          Inc(iy);
        end;
        if bSameColor then Inc(nc);
      end;
      if nc>Cmax then
      begin
        Cmax:=nc;
        CList:=Colors[Cmax] as TListNode;
        CList.Clear;
      end;
      CList.Add(px);
      if nc<Qmin then //this branch is not a good start for a maximum clique search
      begin
        Rloc[ir]:=px; //neutralize by moving it to the begining of the list
        Inc(ir);      //notice that the pointer located initially at location ir is already stored in Colors somewhere
      end;
    end;
    if (ir>0) then Rloc.Col[ir-1]:=0;
    if Qmin<=0 then Qmin:=1;
    for nc:=Qmin to Cmax do //These are colors indexing nodes usefull for clique search
    begin
      CList:=Colors[nc] as TListNode;
      for ic:=0 to CList.Count-1 do
      begin
        Rloc[ir]:=CList[ic];
        Rloc.Col[ir]:=nc;
        Inc(ir);
      end;
    end;
    Result:=Cmax;

    {//C++ code
    int j = 0;
    int maxno = 1;
    int min_k = QMAX.size() - Q.size() + 1;
    C[1].rewind();
    C[2].rewind();
    int k = 1;
    for (int i=0; i < R.size(); i++) {
      int pi = R.at(i).get_i();
      k = 1;
      while (cut1(pi, C[k]))
        k++;
      if (k > maxno) {
        maxno = k;
        C[maxno + 1].rewind();
      }
      C[k].push(pi);
      if (k < min_k) {
        R.at(j++).set_i(pi);
      }
    }
    if (j > 0) R.at(j-1).set_degree(0);
    if (min_k <= 0) min_k = 1;
    for (k = min_k; k <= maxno; k++)
      for (int i = 0; i < C[k].size(); i++) {
        R.at(j).set_i(C[k].at(i));
        R.at(j++).set_degree(k);
      }
    }
  end;

  procedure MaxClique(var Rloc: TListNode);
  //Recursive search for a maximum clique
  //R is a sorted list of nodes, in increasing color order
  //Q is the output maximum clique; Capacity of Q shall be preset to largest color
  //Rref, Q and Qmax have global scope.
  //However, the color depends on the depth of the recursion

  var
    iq: integer;
    ny: Node;
    Rnew, Rngh: TListNode;
    Cmax: Node;
    px, py: PTNodePrp;
    ix, iy, xy: integer;
  begin
    Inc(dpth);
    Inc(niter);
    Cmax:=0;
    while (Rloc.Count>0) do
    begin
      px:=Rloc.Last;//Choose a node with maximum color
      //writeln('Q.Count+px^.Col[Rloc.lvl]=',Q.Count,'+',px^.Col[Rloc.lvl]);
      if Q.Count+px^.Col[Rloc.lvl] > Qmax.Count then//This branch may include a larger clique than the current one
      begin
        Q.Add(px);
        //
        Rnew:=TListNode.Create;//The intersection of adjacent nodes with list nodes
        Rngh:=LRngh[px^.Id] as TListNode;//List of adjacent nodes to x
        Rnew.Assign(Rloc,laAnd,Rngh);//Rnew = R and (adajcente)
        Rnew.lvl:=Rloc.lvl+1;
        If Rnew.Count > 0 then //The intersection of adjacency list and R is not empty -> recursion
        begin
          Rnew.Sort(@sortNodeDeg);
          Cmax:=ApproximateColor(Rnew);
          MaxClique(Rnew);
          Dec(dpth);
        end else if Q.Count>Qmax.Count then //The current clique is the largest one so far
        begin
          Qmax.Assign(Q,laCopy);//Forget content of Qmax and store the content of Q
          writeln('Depth=',dpth,', Iteration=',niter,', Current clique size=',Q.Count,', Largest clique size=',Qmax.Count);
        end;
        FreeAndNil(Rnew);//Release the memory of Rnew but do not destroy the pointers it contain
        Q.Remove(px);//backtrack
      end else
      begin
        Exit;//Exit procedure here
      end;
      Rloc.Remove(px);
    end;
  end;

begin
  //Initialize
  dpth:=0;
  niter:=0;
  sizeMCS := 0;
  for i:=Low(index1) to High(index1) do
  begin
    index1[i]:=0;
    index2[i]:=0;
  end;
  SLog:=TStringList.Create;
  G:=nil;//T_GRAPHE_LISTE.CREATE;
  R:=TListNode.Create; //Set of nodes
  Rref:=TListNode.Create; //Nodes storage
  LRngh:=TObjectList.Create; //Set of adjacent nodes
  LRngh.OwnsObjects:=True;
  //Compute the compability graph
  makeCompatibilityGraph(mol1, mol2, bUseDist, G, index1, index2);
  writeln('G.p_NX=',G.p_NX,' G.p_M=',G.p_M);
  // DEBUG Create a CLQ for comparison
  {Sclq:=TStringList.Create;
  Sclq.Add('p edge '+IntToStr(G.p_M));
  for i:=1 to G.p_NX do
  begin
    for k:=G.p_HEAD[i] to G.p_HEAD[i+1]-1 do
    begin
      j:=G.p_SUCC[k];
      stmp:='e '+IntToStr(i)+' '+IntToStr(j);
      Sclq.Add(stmp);
    end;
  end;
  Sclq.SaveToFile('dimac.clq');
  FreeAndNil(Sclq);
//  halt;}
  // DEBUG
  //Compute the maximul clique
  Rref.SetListNode(G); //Rref store all pointers and is immutable
  Rref.Sort(@sortNodeId); //Rref is sorted as in G: Rref[Id in G]->PRNodePrp(Id=Id in G)
  // LRngh contains at position i the list of pointers j of adjacent nodes
  LRngh.Capacity:=G.p_NX;
  LRngh.Add(nil);//make the list 1-based as in G
  for i:=1 to G.p_NX do
  begin
    LRngh.Add(TListNode.Create as TObject);
    for k:=G.p_HEAD[i] to G.p_HEAD[i+1]-1 do
    begin
      j:=G.p_SUCC[k]-1;//Rref is 0-based
      (LRngh[i] as TListNode).Add(Rref[j]);
    end;
  end;
  //
  for i:=0 to Rref.Count-1 do R.Add(Rref[i]); //R is the working copy to be edited
  R.Sort(@sortNodeDeg);
  //Initialize colors
  maxdeg:=G.p_MaxOutDeg;
  for i:=0 to maxdeg-1 do
    R.Col[i]:=i+1;
  for i:=maxdeg to R.Count-1 do
    R.Col[i]:=maxdeg+1;
  //
  Colors:=TObjectList.Create;//Store list of nodes with the same colors
  Colors.OwnsObjects:=True;
  Colors.Capacity:=G.p_NX;
  for i:=0 to G.p_NX do
    Colors.Add(TListNode.Create);
  Q:=TListNode.Create; //a clique
  Qmax:=TListNode.Create; // max clique
  MaxClique(R);
  write('Clique:');
  for iq:=0 to Qmax.Count-1 do
    write(' '+IntToStr(Qmax.Items[iq]^.Id));
  writeln;
//
  if (Qmax.Count > 0) then
    for iq := 0 to Qmax.Count-1 do
    begin
      First[iq+1] := index1[Qmax.Items[iq]^.Id];
      Second[iq+1] := index2[Qmax.Items[iq]^.Id];
    end
  else begin
    First[1] := index1[1];
    Second[1] := index2[1];
  end;
  Result:=Qmax.Count;
//
  FreeAndNil(G);
  FreeAndNil(R);
  FreeAndNil(LRngh);
  FreeAndNil(Q);
  FreeAndNil(Qmax);
  FreeAndNil(Colors);
  Rref.ClearP;
  FreeAndNil(Rref);
  FreeAndNil(SLog);
end;


function mMCS(mol1, mol2: TmoleculeBase; bUseDist: boolean;
  out First, Second: TNodeInfo; QmaxDelta: ShortInt=0): AtomID;
var
  i,j,k, sizeMCS, compatGraphSize: integer;
  G: T_GRAPHE_LISTE; //compatibility graph
  Rref, R, Q, Qmax: TListNode;//For maximum clique search. Rref, R: list of nodes; Q: current clique; Qmax: maximum clique
  LQmax: TObjectList;
  QmaxSze: integer;
  //Rngh: TListNode;
  index1, index2: TNodeInfo;
  SLog: TStringList;
  C: TNodeInfo;
  iq,ilq: integer;
  Sclq: TStringList;
  stmp: string;
  maxdeg: Node;
  Colors, LRngh: TObjectList;
  aPNode: PTNodePrp;
  ix: integer;
  px: PTNodePrp;
  dpth,niter: integer;

  function DijHeapTopDist(Gloc: TMoleculeBase; s, t: Node; Vmax: Cost = 20): Cost;
  //Compute topological distance between s and t from graph G; do not consider distances above Vmax bonds
  //Input:
  //G, a chemical structure
  //s: starting node
  //t: terminal node
  //Vmax: maximum distance allowed
  //Output: the distance between s and t if less than Vmax, Vmax otherwise
  Var x,y: Node;
      H  : T_Heap;
      k  : ArcNum;
      V  : TNodeCost;
  Begin
    H:=T_Heap.CREATE;
    with Gloc do
    begin
      For x := 1 to p_NX+p_NY do begin
        V[x] := MaxCost
      End;
      H.ClearHeap;
      V[s] := 0;
      H.HeapInsert (V,s);
      Repeat
        H.HeapMin (V,x);
        For k := p_Head[x] to p_Head[x+1]-1 do begin
          y := p_Succ[k];
          If (V[x] + 1 < V[y]) and (V[x] + 1 <= Vmax) then begin
            V[y] := V[x] + 1;
            If H.InHeap (y) then H.MoveUp (V,y) else H.HeapInsert (V,y)
          End
        End;
      Until (H.HeapIsEmpty) or (x=t);
    end;
    FreeAndNil(H);
    Result:=V[t];
  End;
  function calcCompatibilityGraphSize(mol1, mol2: TMoleculeBase; var index1, index2: TNodeInfo): AtomID;
  //Input: molecule1, molecule2; output: the size of the compatibility graph
  //Compute the size of a compatibility graph.
  //Atoms are compatible if they are of same element
  var
    compat_graph_size: integer;
    i, j: AtomID;
  begin
    compat_graph_size := 0;
    for i := 1 to mol1.nAtom do
      for j := 1 to mol2.nAtom do
        //Penser à utiliser une fonction générique de comparaison de noeuds
        if mol1.S_[i] = mol2.S_[j] then begin
          Inc(compat_graph_size);
          index1[compat_graph_size] := i;
          index2[compat_graph_size] := j;
        end;
    Result := compat_graph_size;
  end;
//avec matrice de distance
  procedure makeCompatibilityGraph(mol1, mol2: TMoleculeBase;
    bUseDist: boolean; var compatGraph: T_GRAPHE_LISTE; var index1, index2: TNodeInfo);
  //Input:
  //mol1, mol2: two molecules
  //bUseDist: use a distance matrix (True/False)
  //Output:
  //compatGraph: the compatibility graph as a boolean bidimensional array
  //index1, index2: list of atoms in mol1 and mol2
  //Compute the compatibility graph
  var
    i, j: AtomID;
    PRBnd1, PRBnd2: PRBond;
    compat_graph_size: integer;
    BondCompat: boolean;
    dist1, dist2: Cost;
    W, W1, W2: TArcCost;
    V: TNodeCost;//Distance between nodes
    P: TNodeInfo;//Shortest path between nodes (not used)
    bDistCompat: boolean;
    tmpGraphMatrix: T_GRAPHE_MATRICIEL;
  begin
    tmpGraphMatrix:=T_GRAPHE_MATRICIEL.CREATE;
    compat_graph_size:=calcCompatibilityGraphSize(mol1,mol2,index1,index2);
    tmpGraphMatrix.p_NX:=compat_graph_size; tmpGraphMatrix.p_NY:=0; tmpGraphMatrix.p_NoArc:=0;
    tmpGraphMatrix.p_Simple:=True;
    if compat_graph_size=0 then Exit;
    //Set the cost of edges to 1 in mol1 and mol2 to compute topological distances
    if bUseDist then
    begin
      for i:=1 to mol1.nAtom do W1[i]:=1;
      for i:=1 to mol2.nAtom do W2[i]:=1;
    end;
    // On regarde les pair d'atomes : 1atome dans mol1 et son image dans mol2
    for i := 1 to compat_graph_size - 1 do begin //an atom in mol1
      for j := i + 1 to compat_graph_size do //the corresponding compatible atom in mol2
      begin
        if (index1[i] <> index1[j]) and (index2[i] <> index2[j]) then
        begin
          PRBnd1:=mol1.FindBond(index1[i],index1[j]);
          PRBnd2:=mol2.FindBond(index2[i],index2[j]);
          BondCompat:=False;
          if bUseDist then
          begin
            //If G1 is subgraph of G2, then the distances in G1 are preserved in G2
            //This is an efficient filter if the distance calculation overhead is
            //not compensated by the decreased complexity of the problem
            dist1:=DijHeapTopDist(mol1,index1[i],index1[j]);
            dist2:=DijHeapTopDist(mol2,index2[i],index2[j]);
            //writeln('i=',index1[i],'/',index2[i],' j=',index1[j],'/',index2[j],' dist1=',dist1,' dist2=',dist2);
            bDistCompat:=(dist1=dist2);
          end else
            bDistCompat:=True;
          //
          if bDistCompat then
          begin
               if ((PRBnd1=nil) and (PRBnd2=nil)) then
                  BondCompat:=True //atoms i and j are not bonded in both mol1 and mol2; indispensable condition
               else if ((PRBnd1<>nil) and (PRBnd2<>nil)) then
                    if (PRBnd1^.B=PRBnd2^.B) then begin
                       BondCompat:=True; //atoms i and j are bonded by the same kind of bonds in mol1 and mol2
                    end;
               if BondCompat then begin //record the compatibility in the matrix
                    tmpGraphMatrix.p_A[i, j] := 1;
                    tmpGraphMatrix.p_A[j, i] := 1;
               end else begin
                 tmpGraphMatrix.p_A[i,j]:=tmpGraphMatrix.p_NoArc;
                 tmpGraphMatrix.p_A[j,i]:=tmpGraphMatrix.p_NoArc;
               end;
          end;
        end;
      end;
    end;

    compatGraph:=tmpGraphMatrix.PackC(nil);
    FreeAndNil(tmpGraphMatrix);
  end;
  function ApproximateColor(Rloc: TListNode): Node;
  var
    nx, nc: Node;
    Cnum, Cmax: Node;
    ir, ic, Qmin: integer;
    CList,Rngh: TListNode;
    bSameColor: boolean;
    px, py: PTNodePrp;
    ix, iy, xy: integer;
  begin
    Cmax:=0;  //maxno
    Qmin:=Qmax.Count-Q.Count+1; //min_k
    TListNode(Colors[1]).Clear;
    {for nc:=0 to Colors.Count-1 do
    begin
      CList:=Colors[nc] as TListNode;           // C
      if CList<>nil then CList.Clear;
    end;}
    ir:=0; //j
    Cnum:=1; //k
    for ix:=0 to Rloc.Count-1 do
    begin
      px:=Rloc[ix];
      Rngh:=LRngh[px^.Id] as TListNode;
      bSameColor:=True; //bSameColor=True if a neighbor y of x is in the color set nc
      nc:=1;
      while bSameColor and (nc<=Cmax) do
      begin
        CList:=Colors[nc] as TListNode;
        iy:=0;
        bSameColor:=False;
        while (iy<CList.Count) and (not bSameColor) do begin
          if Rngh.IndexOf(CList[iy])>=0 then bSameColor:=True;
          Inc(iy);
        end;
        if bSameColor then Inc(nc);
      end;
      if nc>Cmax then
      begin
        Cmax:=nc;
        CList:=Colors[Cmax] as TListNode;
        CList.Clear;
      end;
      CList.Add(px);
      if nc<Qmin then //this branch is not a good start for a maximum clique search
      begin
        Rloc[ir]:=px; //neutralize by moving it to the begining of the list
        Inc(ir);      //notice that the pointer located initially at location ir is already stored in Colors somewhere
      end;
    end;
    if (ir>0) then Rloc.Col[ir-1]:=0;
    if Qmin<=0 then Qmin:=1;
    for nc:=Qmin to Cmax do //These are colors indexing nodes usefull for clique search
    begin
      CList:=Colors[nc] as TListNode;
      for ic:=0 to CList.Count-1 do
      begin
        Rloc[ir]:=CList[ic];
        Rloc.Col[ir]:=nc;
        Inc(ir);
      end;
    end;
    Result:=Cmax;

    {//C++ code
    int j = 0;
    int maxno = 1;
    int min_k = QMAX.size() - Q.size() + 1;
    C[1].rewind();
    C[2].rewind();
    int k = 1;
    for (int i=0; i < R.size(); i++) {
      int pi = R.at(i).get_i();
      k = 1;
      while (cut1(pi, C[k]))
        k++;
      if (k > maxno) {
        maxno = k;
        C[maxno + 1].rewind();
      }
      C[k].push(pi);
      if (k < min_k) {
        R.at(j++).set_i(pi);
      }
    }
    if (j > 0) R.at(j-1).set_degree(0);
    if (min_k <= 0) min_k = 1;
    for (k = min_k; k <= maxno; k++)
      for (int i = 0; i < C[k].size(); i++) {
        R.at(j).set_i(C[k].at(i));
        R.at(j++).set_degree(k);
      }
    }
  end;
  procedure MaxClique(var Rloc: TListNode);
  //Recursive search for a maximum clique
  //R is a sorted list of nodes, in increasing color order
  //Q is the output maximum clique; Capacity of Q shall be preset to largest color
  //Rref, Q and Qmax have global scope.
  //However, the color depends on the depth of the recursion

  var
    iq: integer;
    ny: Node;
    Rnew, Rngh: TListNode;
    Cmax: Node;
    px, py: PTNodePrp;
    ix, iy, xy: integer;
  begin
    Inc(dpth);
    Inc(niter);
    Cmax:=0;
    while (Rloc.Count>0) do
    begin
      px:=Rloc.Last;//Choose a node with maximum color
      //writeln('Q.Count+px^.Col[Rloc.lvl]=',Q.Count,'+',px^.Col[Rloc.lvl]);
      if Q.Count+px^.Col[Rloc.lvl] > Qmax.Count then//This branch may include a larger clique than the current one
      begin
        Q.Add(px);
        //
        Rnew:=TListNode.Create;//The intersection of adjacent nodes with list nodes
        Rngh:=LRngh[px^.Id] as TListNode;//List of adjacent nodes to x
        Rnew.Assign(Rloc,laAnd,Rngh);//Rnew = R and (adajcente)
        Rnew.lvl:=Rloc.lvl+1;
        If Rnew.Count > 0 then //The intersection of adjacency list and R is not empty -> recursion
        begin
          Rnew.Sort(@sortNodeDeg);
          Cmax:=ApproximateColor(Rnew);
          MaxClique(Rnew);
          Dec(dpth);
        end else if Q.Count>=QmaxSze-QmaxDelta then //The current clique is the largest one so far
        begin
          Qmax.Assign(Q,laCopy);//Forget content of Qmax and store the content of Q
          if Qmax.Count>QmaxSze then begin
            QmaxSze:=Qmax.Count;
            writeln('Depth=',dpth,', Iteration=',niter,', Current clique size=',Q.Count,', Largest clique size=',QmaxSze);
          end;
          Qmax:=TListNode.Create;
          LQmax.Add(Qmax);
        end;
        FreeAndNil(Rnew);//Release the memory of Rnew but do not destroy the pointers it contain
        Q.Remove(px);//backtrack
      end else
      begin
        Exit;//Exit procedure here
      end;
      Rloc.Remove(px);
    end;
  end;

begin
  //Initialize
  dpth:=0;
  niter:=0;
  sizeMCS := 0;
  QmaxSze:=0;
  for i:=Low(index1) to High(index1) do
  begin
    index1[i]:=0;
    index2[i]:=0;
  end;
  SLog:=TStringList.Create;
  G:=nil;//T_GRAPHE_LISTE.CREATE;
  R:=TListNode.Create; //Set of nodes
  Rref:=TListNode.Create; //Nodes storage
  LRngh:=TObjectList.Create; //Set of adjacent nodes
  LRngh.OwnsObjects:=True;
  //Compute the compability graph
  makeCompatibilityGraph(mol1, mol2, bUseDist, G, index1, index2);
  writeln('G.p_NX=',G.p_NX,' G.p_M=',G.p_M);
  // DEBUG Create a CLQ for comparison
  {Sclq:=TStringList.Create;
  Sclq.Add('p edge '+IntToStr(G.p_M));
  for i:=1 to G.p_NX do
  begin
    for k:=G.p_HEAD[i] to G.p_HEAD[i+1]-1 do
    begin
      j:=G.p_SUCC[k];
      stmp:='e '+IntToStr(i)+' '+IntToStr(j);
      Sclq.Add(stmp);
    end;
  end;
  Sclq.SaveToFile('dimac.clq');
  FreeAndNil(Sclq);
//  halt;}
  // DEBUG
  //Compute the maximul clique
  Rref.SetListNode(G); //Rref store all pointers and is immutable
  Rref.Sort(@sortNodeId); //Rref is sorted as in G: Rref[Id in G]->PRNodePrp(Id=Id in G)
  // LRngh contains at position i the list of pointers j of adjacent nodes
  LRngh.Capacity:=G.p_NX;
  LRngh.Add(nil);//make the list 1-based as in G
  for i:=1 to G.p_NX do
  begin
    LRngh.Add(TListNode.Create as TObject);
    for k:=G.p_HEAD[i] to G.p_HEAD[i+1]-1 do
    begin
      j:=G.p_SUCC[k]-1;//Rref is 0-based
      (LRngh[i] as TListNode).Add(Rref[j]);
    end;
  end;
  //
  for i:=0 to Rref.Count-1 do R.Add(Rref[i]); //R is the working copy to be edited
  R.Sort(@sortNodeDeg);
  //Initialize colors
  maxdeg:=G.p_MaxOutDeg;
  for i:=0 to maxdeg-1 do
    R.Col[i]:=i+1;
  for i:=maxdeg to R.Count-1 do
    R.Col[i]:=maxdeg+1;
  //
  Colors:=TObjectList.Create;//Store list of nodes with the same colors
  Colors.OwnsObjects:=True;
  Colors.Capacity:=G.p_NX;
  for i:=0 to G.p_NX do
    Colors.Add(TListNode.Create);
  Q:=TListNode.Create; //a clique
  Qmax:=TListNode.Create; // max clique
  LQmax:=TObjectList.Create;
  LQmax.OwnsObjects:=True;
  LQmax.Add(Qmax);
  MaxClique(R);
  //Cleanup LQmax
  for ilq:=LQmax.Count-1 downto 0 do
  begin
    Qmax:=LQmax[ilq] as TListNode;
    if Qmax.Count<QmaxSze-QmaxDelta then begin
      LQmax.Remove(Qmax);
    end;
  end;
  writeln('Clique:');
  for ilq:=0 to LQmax.Count-1 do
  begin
    Qmax:=LQmax[ilq] as TListNode;
    for iq:=0 to Qmax.Count-1 do
      write(' '+IntToStr(Qmax.Items[iq]^.Id));
    writeln;
  end;
//
  if (Qmax.Count > 0) then
    for iq := 0 to Qmax.Count-1 do
    begin
      First[iq+1] := index1[Qmax.Items[iq]^.Id];
      Second[iq+1] := index2[Qmax.Items[iq]^.Id];
    end
  else begin
    First[1] := index1[1];
    Second[1] := index2[1];
  end;
  Result:=Qmax.Count;
//
  FreeAndNil(G);
  FreeAndNil(R);
  FreeAndNil(LRngh);
  FreeAndNil(Q);
  FreeAndNil(LQmax);
  FreeAndNil(Colors);
  Rref.ClearP;
  FreeAndNil(Rref);
  FreeAndNil(SLog);
end;





procedure GetSSList(mol1, mol2: TmoleculeBase; minSze, maxSze: integer;
  bUseDist: boolean; out SSList: aiarray);
// Search for list of common substructures in mol1 and mol2 of size between minSze nd maxSze
var
  i, sizeMCS, compatGraphSize: integer;
  compatGraph,dynCompatGraph: AAoB;
  index1, index2, maxClique: TNodeInfo;
  {i1,j1: integer;}

  function calcCompatibilityGraphSize(mol1, mol2: TMoleculeBase): AtomID;
  var
    compat_graph_size: integer;
    i, j: AtomID;
  begin
    compat_graph_size := 0;
    for i := 1 to mol1.p_NX do
      for j := 1 to mol2.p_NX do
        //Penser à utiliser une fonction générique de comparaison de noeuds
        if mol1.S_[i] = mol2.S_[j] then
          Inc(compat_graph_size);
    Result := compat_graph_size;
  end;
  //avec matrice de distance
  procedure makeCompatibilityGraph(mol1, mol2: TMoleculeBase;
    bUseDist: boolean;dynCompatGraph:AAoB; compatGraph: AAoB; var index1, index2: TNodeInfo);
  //INPUTS:
  // - mol1, mol2: MoleculeBase instances to compare
  // - bUseDist: using the inter-atomic topological distances both in mol1 and mol2 to discard irrelevant potential clique
  // - dynCompatGraph: compatibility graphs for dynamic bonds -> to discard?
  // - compatGraph: compatibility graph. compatGraph[i,j]=True if mol1 atom i and mol2 atom j are members of a common substructure
  // - index1, index2: mol1 atom indices and mol2 atom indices of common substructures.
  var
    i, j: AtomID;
    PRBnd1, PRBnd2: PRBond;
    compat_graph_size: integer;
    BondCompat,DynBondCompat: boolean;
    dist1, dist2: Cost;
    W1, W2: TArcCost;
    V: TNodeCost;
    P: TNodeInfo;
    bDistCompat: boolean;
  begin
    for i:=Low(V) to High(V) do V[i]:=0;
    for i:=Low(P) to High(P) do P[i]:=0;
    compat_graph_size := 0;
    for i := 1 to mol1.nAtom do
      for j := 1 to mol2.nAtom do //Build a table that convert compatGraph indices to tuple: (Atom in mol1, Atom in mol2)
        //Penser à généraliser par une fonction d'égalité portant sur les noeuds
        if mol1.S_[i] = mol2.S_[j] then
        begin
          Inc(compat_graph_size);
          index1[compat_graph_size] := i;
          index2[compat_graph_size] := j;
        end;

    if bUseDist then
    begin
      for i:=1 to mol1.p_M do W1[i]:=1;
      for i:=1 to mol2.p_M do W2[i]:=1;
    end;

    // On regarde les pair d'atomes : 1atome dans mol1 et son image dans mol2
    for i := 1 to compat_graph_size - 1 do begin
      //on prend une seconde pair
      for j := i + 1 to compat_graph_size do
      begin
        PRBnd1:=mol1.FindBond(index1[i],index1[j]);
        PRBnd2:=mol2.FindBond(index2[i],index2[j]);
        BondCompat:=False;
        DynBondCompat:=False;
        if bUseDist then
        begin
          mol1.DijHeap(W1,index1[i],index1[j],v,P);
          dist1:=V[index1[j]];
          mol2.DijHeap(W2,index2[i],index2[j],v,P);
          dist2:=V[index2[j]];
          bDistCompat:=(dist1=dist2);
        end else
          bDistCompat:=True;
        if (index1[i] <> index1[j]) and (index2[i] <> index2[j]) and bDistCompat then
        begin
             if ((PRBnd1=nil) and (PRBnd2=nil)) then
                BondCompat:=True
             else if ((PRBnd1<>nil) and (PRBnd2<>nil)) then
                  if (PRBnd1^.B=PRBnd2^.B) then begin
                     BondCompat:=True;
                     if (PrBnd1^.B >= 14) then
                        DynBondCompat:=True;
                  end;
             if BondCompat then begin
                  compatGraph[i, j] := True;
                  compatGraph[j, i] := True;
             end
             else begin
                  compatGraph[i, j] := False;
                  compatGraph[j, i] := False;
             end;
             if DynBondCompat then begin
                  dynCompatGraph[i, j] := True;
                  dynCompatGraph[j, i] := True;
             end
             else begin
                  dynCompatGraph[i, j] := False;
                  dynCompatGraph[j, i] := False;
             end;

        end else
        begin
          compatGraph[i,j]:=False;
          compatGraph[j,i]:=False;
          dynCompatGraph[i, j] := False;
          dynCompatGraph[j, i] := False;
        end;
      end;
    end;
  end;

begin
  sizeMCS := 0;
  for i:=Low(maxClique) to High(maxClique) do
  begin
    maxClique[i]:=0;
    index1[i]:=0;
    index2[i]:=0;
  end;
  compatGraphSize := calcCompatibilityGraphSize(mol1, mol2);
  if compatGraphSize <> 0 then
  begin
    SetLength(compatGraph, compatGraphSize + 1, compatGraphSize + 1);
    SetLength(dynCompatGraph, compatGraphSize +1, compatGraphSize +1);
    makeCompatibilityGraph(mol1, mol2, bUseDist,dynCompatGraph,
      compatGraph, index1, index2);
    //sizeMCS := searchClique(compatGraph,dynCompatGraph, compatGraphSize, maxClique);
    sizeMCS := searchClique(compatGraph, compatGraphSize, maxClique);
    Writeln('Size MCS : '+IntToStr(SizeMCS));
    if (sizeMCS > 1) then
      for i := 1 to sizeMCS do
      //begin
        Write('('+IntToStr(index1[maxclique[i]])+' | '+IntToStr(index2[maxClique[i]])+') ');
        //First[i] := index1[maxclique[i]];
        //Second[i] := index2[maxClique[i]];
    {  end
    else
      for i := 1 to sizeMCS do
      begin
        First[i] := index1[1];
        Second[i] := index2[1];
      end;}
  end;
  //Result := sizeMCS;
end;

{ TMoleculeBase }

procedure TMoleculeBase.DisposeAtmSet;
var
  i: integer;
begin
  for i := Low(fAtmSet) to p_NX do
    if (fAtmSet[i] <> nil) then
    begin
      dispose(fAtmSet[i]);
      fAtmSet[i] := nil;
    end;
end;

procedure TMoleculeBase.DisposeBndSet;
var
  i, j: integer;
  pi,pj: PRBond;
begin
  {//work in progress fLBndSet as a bond pointer manager
  for i:=Low(fBndSet) to p_M do
    fBndSet[i]:=nil;
  for i:=0 to fLBndSet.Count-1 do
  begin
    pi:=fLBndSet[i];
    if pi<>nil then
    begin
      dipose(pi);
      fLBndSet[i]:=nil;
    end;
  end;
  fLBndSet.Clear;}
  for i := Low(fBndSet) to p_M do
  begin
    pi:=fBndSet[i];
    //search for multiple copies of a bond.
    for j := i + 1 to p_M do
    begin
      pj:=fBndSet[j];
      if (pi = pj) then
        fBndSet[j] := nil;
    end;
  end;
  for i := Low(fBndSet) to p_M do
  begin
    pi:=fBndSet[i];
    if (pi <> nil) then
    begin
      dispose(pi);
      fBndSet[i] := nil;
    end;
  end;
end;

function TMoleculeBase.GetAtom(i: AtomID): PRAtom;
  //Beware, the owner does not own the pointer to the RAtom record!
begin
  Result := fAtmSet[i];
end;

function TMoleculeBase.GetI(i: AtomID): AByt;
begin
  Result := fAtmSet[i]^.I;
end;

function TMoleculeBase.GetP(i: AtomID): APrp;
begin
  Result := fAtmSet[i]^.P;
end;

function TMoleculeBase.GetS(i: AtomID): TS;
begin
  Result := fAtmSet[i]^.S;
end;

function TMoleculeBase.GetW(i: AtomID): double;
begin
  Result := fAtmSet[i]^.W;
end;

function TMoleculeBase.GetZ(i: AtomID): TZ;
begin
  Result := fAtmSet[i]^.Z;
end;

procedure TMoleculeBase.SetW(i: AtomID; const AValue: double);
begin
  fAtmSet[i]^.W := AValue;
end;

procedure TMoleculeBase.SetAtom(i: AtomID; PAt: PRAtom);
//Beware, the caller does not anymore owns the pointer to the RAtom record!
begin
  {
  if (fAtmSet[i] <> nil) then
    dispose(fAtmSet[i]);
  fAtmSet[i] := PAt;
  }
  if (fAtmSet[i] <> nil) then
    PRAtom(fAtmSet[i])^:=PAt^
  else fAtmSet[i] := PAt;
end;

function TMoleculeBase.GetBond(i: BondID): PRBond;
  //Beware, the owner does not own the pointer to the RBond record!
begin
  Result := fBndSet[i];
end;

procedure TMoleculeBase.SetBond(i: BondID; PBo: PRBond);
//Beware, the caller does not anymore owns the pointer to the RBond record!
begin
  if (fBndSet[i] <> nil) then
    PRBond(fBndSet[i])^:=PBo^
  else fBndSet[i] := PBo;
  {
  if (fBndSet[i] <> nil) then
    dispose(fBndSet[i]);
  fBndSet[i] := PBo;
  }
end;

function TMoleculeBase.int_readpos(str: string; bgn, lng: integer): integer;
var
  CodeErr: integer;
  wrd: string;
begin
  Result := 0;
  CodeErr := 0;  // fix error in input SDF file
  wrd := Trim(Copy(str, bgn, lng));
  if Length(wrd) > 0 then
    Val(wrd, Result, CodeErr);                         // number of atoms
  if CodeErr <> 0 then
  begin
    Exception.Create('ERROR: reading sdf entry ' + fMolName);
    halt(1);
  end;
end;

function TMoleculeBase.dbl_readpos(str: string; bgn, lng: integer): double;
var
  CodeErr: integer;
  wrd: string;
begin
  Result := 0;
  CodeErr := 0;  // fix error in input SDF file
  wrd := Trim(Copy(str, bgn, lng));
  if Length(wrd) > 0 then
    Val(wrd, Result, CodeErr);                         // number of atoms
  if CodeErr <> 0 then
  begin
    Exception.Create('ERROR: reading sdf entry ' + fMolName);
    halt(1);
  end;
end;

function TMoleculeBase.convert_dict_hash(str: string): TFPStringHashTable;
var
  strsplit: TStringList;
  strhash: TFPStringHashTable;
  i: integer;
begin
  strsplit := TStringList.Create;
  strhash := TFPStringHashTable.Create;
  Assert(Assigned(strsplit));
  strsplit.Clear;

  //Splits the string with the ";" delimiter
  try
    strsplit.StrictDelimiter := True;
    strsplit.Delimiter := ';';
    strsplit.DelimitedText := str;
    for i := 0 to strsplit.Count - 1 do
    begin
      strhash.Add(strsplit[i], '');
    end;
  finally
    FreeAndNil(strsplit);

  end;
  Result := strhash;
end;

function TMoleculeBase.convert_dict_list(str: string): TStringList;
var
  strsplit: TStringList;
  {i: integer;}
begin
  strsplit := TStringList.Create;
  Assert(Assigned(strsplit));
  strsplit.Clear;

  //Splits the string with the ";" delimiter
  try
    strsplit.StrictDelimiter := True;
    strsplit.Delimiter := ';';
    strsplit.DelimitedText := str;
  finally

  end;
  Result := strsplit;
end;

function TMoleculeBase.readpNX(): integer;
begin
  Result := p_NX;
end;

function TMoleculeBase.readpM(): integer;
begin
  Result := p_M;
end;

procedure TMoleculeBase.SetS(i: AtomID; const AValue: TS);
begin
  fAtmSet[i]^.S := AValue;
end;

procedure TMoleculeBase.SetZ(i: AtomID; const AValue: TZ);
begin
  fAtmSet[i]^.Z := AValue;
end;

procedure TMoleculeBase.writepNX(i: integer);
begin
  p_NX := i;
end;

procedure TMoleculeBase.writepM(i: integer);
begin
  p_M := i;
end;

constructor TMoleculeBase.Create;
var
  i: AtomID;
  j: BondID;
  {n: integer;}
begin
  fMolName := 'NoName';
  fAPrpSze := 0;
  fABytSze := 0;
  fBPrpSze := 0;
  fBBytSze := 0;
  //fLBndSet:=TList.Create; //Work in progress -> a bond pointer manager
  for i := Low(fAtmSet) to High(fAtmSet) do
    fAtmSet[i] := nil;
  for j := Low(fBndSet) to High(fBndSet) do
    fBndSet[j] := nil;
  inherited Create;
end;

destructor TMoleculeBase.Destroy;
begin
  DisposeAtmSet;
  DisposeBndSet;
  //FreeAndNil(fLBndSet); //Work in progress -> a bond pointer manager
  inherited Destroy;
end;

procedure TMoleculeBase.Clear;
begin
  fMolName := 'NoName';
  fAPrpSze := 0;
  fABytSze := 0;
  fBPrpSze := 0;
  fBBytSze := 0;
  DisposeAtmSet;
  DisposeBndSet;
  p_M := 0;
  p_NX := 0;
  p_NY := 0;
end;

//Read an sdf String
procedure TMoleculeBase.LoadSDFString(sdfstr: string);
var
  stmp: TStringList;
begin
  stmp := TStringList.Create;
  stmp.Text := sdfstr;
  LoadSDF(stmp);
  FreeAndNil(stmp);
end;

//Read an sdf TStringList
{procedure TMoleculeBase.LoadSDFTStringList(sdfstr: TStringList);
const
     Np = 10;
var
     sdfstrsize: integer;
     strlist: TStringList;
     i, j, k, indexof: integer;
     LineNo: integer;
     AddPrp: integer;
     CodeErr : Integer;
     wrd     : string;
     PAt     : PRAtom;
     PBo     : PRBond;
     atmp    : CostMatrix;
     wtmp    : array of PRBond;
     s,t     : Node;
     M       : ArcNum;
begin
     sdfstrsize:=sdfstr.Count;
     //User defined
     fAPrpSze:=0;
     fABytSze:=0;
     fBPrpSze:=0;
     fBBytSze:=0;
     //Read each line and create the base molecule
     LineNo := 0;
     //---Line 1 : Molecule name---
     fMolName := sdfstr[LineNo];
     //---Line 2 of MOLfile---
     //  User's first and last initials (UserI), program name (ProgN), date/time:
     //  M/D/Y,H:m (DatTim), dimensional codes (DimCo), scaling factors (ScaI, ScaR),
     //  energy (Ene) if modeling program input, internal registry number (RegNo)
     //  if input through MDL format:
     Inc(LineNo);
     //---Line 3 of MOLfile---
     // A line for comments. If no comment is entered, a blank line must be present.
     Inc(LineNo);
     //---Line 4 of MOLfile : The Counts Line---
     Inc(LineNo);
     p_NX := 0;
     p_NY:= 0; //Molecular graphs are not bipartite a priori
     p_M := 0;
     p_NX:=int_readpos(sdfstr[LineNo],1,3);     // number of atoms
     p_M:=int_readpos(sdfstr[LineNo],4,3);      // number of bonds
     //Position 31 : Number of lines of additional properties
     AddPrp:=int_readpos(sdfstr[LineNo],31,3);
     //---Lines 5-... of MOLfile : The Atom Block---
     fAtmSet[0]:=nil;
     for i:=1 to p_NX do begin
         Inc(LineNo);
         new(PAt);
         SetLength(PAt^.P,fAPrpSze);
         for j:=Low(PAt^.P) to High(PAt^.P) do PAt^.P[j]:=0;
         PAt^.S:=Trim(Copy(sdfstr[LineNo],32,3));        // atom symbol
         PAt^.Z:=AtomSymbolToInt(PAt^.S);                // atomic number
         AtmSet[i]:=PAt;
         //Only these infos are stored. Add whichever you would like to add.
     end;  // 1..p_NX
     //---Lines of MOLfile : The Bond Block---
     //Express the molecular graph as a simple matrix graph
     for s:=0 to MaxAtom+1 do for t:=0 to MaxAtom+1 do atmp[s,t]:=0;
     SetLength(wtmp,p_M+1);
     fBndSet[0]:=nil;
     for i:=1 to p_M do begin
         Inc(LineNo);
         new(PBo);
         SetLength(PBo^.P,fBPrpSze); SetLength(PBo^.I,fBBytSze);
         PBo^.t:=int_readpos(sdfstr[LineNo],1,3);//Tail of bond
         PBo^.h:=int_readpos(sdfstr[LineNo],4,3);//Head of bond
         PBo^.B:=int_readpos(sdfstr[LineNo],7,3);//Type of bond
         PBo^.S:=BondSymbol[PBo^.B];
         atmp[PBo^.t,PBo^.h]:=i; atmp[PBo^.h,PBo^.t]:=i;
         wtmp[i]:=PBo;
         //Only these infos are stored. Add whichever you would like to add.
     end;
     //Store the graph in a packed format
     M:=0;
     for s:=1 to p_NX do begin
         p_HEAD[s]:=M+1;
         for t:=1 to p_NX do if (atmp[s,t]<>0) then begin
             Inc(M);
             p_SUCC[M]:=t;
             BndSet[M]:=wtmp[atmp[s,t]];
         end;
     end;
     p_HEAD[p_NX+p_NY+1]:=M+1;
     p_M:=M;
end;}

//Read an sdf TStringList containing only one molecule
//Add a dictionary to retrieve interesting properties
function TMoleculeBase.LoadSDFTStringList(sdfstr: TStringList;
  prpstr: string): TFPStringHashTable;
begin
  LoadSDF(sdfstr);
  Result := LoadSDFField(sdfstr, prpstr);
end;

function TMoleculeBase.LoadSDFField(sdfstr: TStringList;
  prpstr: string): TFPStringHashTable;
var
  strlist: TStringList;
begin
  //---Lines of Molfile : the properties---
  if (Length(prpstr) > 0) then
  begin
    strlist := convert_dict_list(prpstr);
    Result := LoadSDFField(sdfstr, strlist);
  end;
end;

function TMoleculeBase.LoadSDFField(sdfstr: TStringList;
  prpstr: TStringList): TFPStringHashTable;
var
  strhash: TFPStringHashTable;
begin
  strhash := TFPStringHashTable.Create;
  //For each line of the dictionary, retrieve the prop in the tstringlist of the pdf
  LoadSDFField(sdfstr, prpstr, strhash);
  Result := strhash;
end;

procedure TMoleculeBase.LoadSDFField(sdfstr: TStringList; prpstr: string;
  var strhash: TFPStringHashTable);
begin
  strhash := LoadSDFField(sdfstr, prpstr);
end;

procedure TMoleculeBase.LoadSDFField(sdfstr: TStringList; prpstr: TStringList;
  strhash: TFPStringHashTable);
var
  i, indexof: integer;
begin
  strhash.Clear;
  //For each line of the dictionary, retrieve the prop in the tstringlist of the pdf
  for i := 0 to prpstr.Count - 1 do
  begin
    //Test search
    indexof := sdfstr.IndexOf('>  <' + prpstr[i] + '>');
    if (indexof > 0) then
    begin
      strhash.Add(prpstr[i], sdfstr[indexof + 1]);
    end;
  end;
end;

function TMoleculeBase.LoadSDFSGroup(sdfstr: TStringList): TList;
var
  datalist: TList;
  i, j1{, n, indexof}: integer;
  lne,aword1,aword2: string;
  pitem: PRSAL;
  bStop: boolean;
  e: EDynWrd;

  function NextS(offset: integer): integer;
  var
    i: integer;
    bStop: boolean;
  begin
    i:=offset;
    bStop:=False;
    repeat
      if (Pos('M  SAL',sdfstr[i])>0) then bStop:=True;
      if (Pos('M  SED',sdfstr[i])>0) then bStop:=True;
      if (Pos('M  SDT',sdfstr[i])>0) then bStop:=True;
      Inc(i);
      if (i>=sdfstr.Count) then bStop:=True;
    until bStop;
    if (i<sdfstr.Count) then NextS:=i-1 else NextS:=-1;
  end;

begin
  datalist := TList.Create;
  i:=0; aword1:='';aword2:=''; bStop:=False; e:=UNK;
  {repeat
    if (Pos('M  SAL',sdfstr[i])>0) or (i>=(sdfstr.Count-1)) then bStop:=True;
    Inc(i);
  until bStop;}
  i:=NextS(0); j1:=i;
  if (i<sdfstr.Count) then
  begin
    //i:=i-1;
    repeat
      {j1:=i;} aword1:='';aword2:=''; bStop:=False;// aword1:=sdfstr[i];
      repeat
        j1:=NextS(j1+1);
        if j1>0 then
        begin
          if (Pos('M  SED',sdfstr[j1])>0) then
          begin
            lne:=sdfstr[j1];
            aword2:=Copy(lne,12,3);
          end;
          if (Pos('M  SDT',sdfstr[j1])>0) then
          begin
            lne:=sdfstr[j1];
            aword1:=Copy(lne,12,13);
            e:=StringToEDynWrd(Trim(aword1));
          end;
          if (Pos('M  SAL',sdfstr[j1])>0) then bStop:=True;
        end else bStop:=True;
      until bStop;
      if (e<UNK) then
      begin
        new(pitem);
        lne:=sdfstr[i];
        pitem^.etype:=e;
        pitem^.asize:=int_readpos(lne,11,3);
        pitem^.alist[1]:=int_readpos(lne,15,3);
        pitem^.alist[2]:=int_readpos(lne,19,3);
        pitem^.aword:=aword2;
        datalist.Add(pitem);
        //writeln(EDynWrdToString(e)+' '+IntToStr(pitem^.asize)+' '+IntToStr(pitem^.alist[1])+' '+IntToStr(pitem^.alist[2])+' '+pitem^.aword);
      end;
      i:=j1;
    until i<0;
  end;
  Result:=datalist;
end;

procedure TMoleculeBase.LoadSDFSGroup(sdfstr: TStringList;
  var datalist: TList);
//Add data to an existing datalist
var
  ToAdd: TList;
begin
  ToAdd := LoadSDFSGroup(sdfstr);
  datalist.AddList(ToAdd);
  FreeAndNil(ToAdd);
end;

//Copy and paste the parts of this procedure that are of interest
procedure TMoleculeBase.LoadSDF(sdfstr: TStringList);
{const
  Np = 10;}
var
  {%H-}sdfstrsize: integer;
  i, j{, k}: integer;
  LineNo: integer;
  {CodeErr: integer;
  atmL, cf: integer;
  stxt, ap: integer;
  wrd: string;}
  PAt: PRAtom;
  PBo: PRBond;
  atmp: CostMatrix;
  wtmp: array of PRBond;
  s, t: Node;
  M: ArcNum;

begin
  //------Line 1 of MOLfile------------------------------------------------------
  sdfstrsize := sdfstr.Count;
  LineNo := 0;
  MolName := sdfstr[LineNo];
  //------Line 2 of MOLfile------------------------------------------------------
  //  User's first and last initials (UserI), program name (ProgN), date/time:
  //  M/D/Y,H:m (DatTim), dimensional codes (DimCo), scaling factors (ScaI, ScaR),
  //  energy (Ene) if modeling program input, internal registry number (RegNo)
  //  if input through MDL format:
  Inc(LineNo);
  //------Line 3 of MOLfile------------------------------------------------------
  // A line for comments. If no comment is entered, a blank line must be present.
  Inc(LineNo);
  //------Line 4 of MOLfile : The Counts Line------------------------------------
  Inc(LineNo);
  p_NX := 0;
  p_NY := 0; //Molecular graphs are not bipartite a priori
  p_M := 0;
  p_NX := int_readpos(sdfstr[LineNo], 1, 3);      // number of atoms
  p_M := int_readpos(sdfstr[LineNo], 4, 3);      // number of bonds
  j := 0;

  //Position 7 : Number of atom lists
  //begin atmL:=int_readpos(sdfstr[LineNo],7,3); inc(j); end;   // number of atom lists
  //position 10 is obsolete
  //Position 13 : Chiral flag
  //begin cf:=int_readpos(sdfstr[LineNo],13,3); inc(j); end;  //chiral flag

  //Position 16 : Number of stext entries (ISIS)
  //begin stxt:=int_readpos(sdfstr[LineNo],16,3); inc(j); end;  //number of stext entries (ISIS)
  //position 19 is obsolete
  //position 22 is obsolete
  //position 25 is obsolete
  //position 28 is obsolete

  //Position 31 : Number of lines of additional properties
  //begin ap:=int_readpos(sdfstr[LineNo],31,3); inc(j); end;  // number of lines of additional properties

  //------Lines 5-... of MOLfile : The Atom Block--------------------------------
  AtmSet[0] := nil;
  for i := 1 to p_NX do
  begin
    Inc(LineNo);
    new(PAt);
    //SetAPrpSize(2); //set the atom properties array size (double)
    //SetABytSize(2); //set the atom properties array size (Byte)
    SetLength(PAt^.P, APrpSze);
    SetLength(PAt^.I, ABytSze);
    for j := Low(PAt^.P) to High(PAt^.P) do
      PAt^.P[j] := 0;
    for j := Low(PAt^.I) to High(PAt^.I) do
      PAt^.I[j] := 0;
    j := 0;
    PAt^.S := Trim(Copy(sdfstr[LineNo], 32, 3));        // atom symbol
         {begin
            PAt^.P[j]:=dbl_readpos(sdfstr[LineNo],1,10); //X
            inc(j);
         end;
         begin
            PAt^.P[j]:=dbl_readpos(sdfstr[LineNo],11,10);//Y
            inc(j);
         end;
         begin
            PAt^.P[j]:=dbl_readpos(sdfstr[LineNo],21,10);//Z
            inc(j);
         end;
         PAt^.S:=Trim(Copy(sdfstr[LineNo],32,3));        // atom symbol
         PAt^.Z:=AtomSymbolToInt(PAt^.S);                // atomic number
         j:=0;
         begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],35,2); //mass difference
            Inc(j);
         end;
         begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],37,3); //charge
            Inc(j);
         end;
         begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],40,3); //atom stereo parity
            Inc(j);
         end;
         begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],43,3); //hydrogen count + 1
            Inc(j);
         end;
         if fKeepAtmInfo[4] then begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],46,3); //stereo care box
            Inc(j);
         end;
         begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],49,3); //valence
            Inc(j);
         end;
         begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],52,3); //H0 designator
            Inc(j);
         end;
         begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],61,3); //atom-atom mapping number
            Inc(j);
         end;
          begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],64,3); //inversion/retention flag
            Inc(j);
         end;
          begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],67,3); //exact change flag
            Inc(j);
         end;}
    AtmSet[i] := PAt;
  end;  // 1..p_NX
  //------Lines of MOLfile : The Bond Block--------------------------------------
  //Express the molecular graph as a simple matrix graph
  for s := 0 to MaxAtom + 1 do
    for t := 0 to MaxAtom + 1 do
      atmp[s, t] := 0;
  SetLength(wtmp, p_M + 1);
  BndSet[0] := nil;
  for i := 1 to p_M do
  begin
    Inc(LineNo);
    new(PBo);
    //SetBPrpSize(2); //set the bond properties array size (double)
    //SetBBytSize(2); //set the bond properties array size (Byte)
    SetLength(PBo^.P, BPrpSze);
    SetLength(PBo^.I, BBytSze);
    PBo^.t := int_readpos(sdfstr[LineNo], 1, 3);//Tail of bond
    PBo^.h := int_readpos(sdfstr[LineNo], 4, 3);//Head of bond
    PBo^.B := IntToTB(int_readpos(sdfstr[LineNo], 7, 3),
      int_readpos(sdfstr[LineNo], 19, 3));//Type of bond
    PBo^.S := BondSymbol[PBo^.B];
    j := 0;
         {begin
            PBo^.I[j]:=int_readpos(sdfstr[LineNo],10,3);
            inc(j);
         end;
         begin
            PBo^.I[j]:=int_readpos(sdfstr[LineNo],16,3);
            inc(j);
         end;}
    atmp[PBo^.t, PBo^.h] := i;
    atmp[PBo^.h, PBo^.t] := i;
    wtmp[i] := PBo;
  end;
  //Store the graph in a packed format
  M := 0;
  for s := 1 to p_NX do
  begin
    p_HEAD[s] := M + 1;
    for t := 1 to p_NX do
      if (atmp[s, t] <> 0) then
      begin
        Inc(M);
        p_SUCC[M] := t;
        BndSet[M] := wtmp[atmp[s, t]];
      end;
  end;
  p_HEAD[p_NX + p_NY + 1] := M + 1;
  p_M := M;
end;

function TMoleculeBase.RemoveBondDir(PBo: PRBond): boolean;
begin
  Result := RemoveBondDir(PBo^.t, PBo^.h);
end;

function TMoleculeBase.RemoveBondDir(s, t: AtomID): boolean;
var
  i,j: integer;
  M, N: BondID;
  x: AtomID;
  PBd: PRBond;
begin
  //fBndSet is updated...
  Result := False;
  M := p_Head[s];
  i := 0;
  while ((i < p_OutDeg[s]) and (p_SUCC[M + i] <> t)) do
    Inc(i);
  M := M + i;
  if (p_SUCC[M] = t) then
  begin
    Result := True;
    PBd:=fBndSet[M];//This adress will be invalid
    if PBd<>nil then begin
      for j:=1 to nBonds do//Search for this adress to point to nil pointer
        if fBndSet[j]=PBd then //One record per edge = two arcs; if one arc has been destroyed the object may not be valid anymore
          fBndSet[j]:=nil;
      dispose(PBd); PBd:=nil;
    end;
    for N := M to nBonds do
    begin
      p_SUCC[N] := p_SUCC[N + 1];
      fBndSet[N] := fBndSet[N + 1];
    end;
    //p_SUCC[nBonds] := 0;
    //fBndSet[nBonds] := nil;
    for x := s + 1 to GraphOrder + 1 do
      p_HEAD[x] := p_HEAD[x] - 1;
    nBonds := nBonds - 1; //nBonds is in fact p_M
  end;
  p_HEAD[p_NX + p_NY + 1] := p_M + 1;
end;

function TMoleculeBase.RemoveBond(s, t: AtomID): boolean;
var
  b1, b2: boolean;
begin
  //fBndSet is not managed here
  Result := False;
  b1 := RemoveBondDir(s, t);
  b2 := RemoveBondDir(t, s);
  if (b1 and b2) then
    Result := True;
end;

function TMoleculeBase.RemoveBond(m: BondID): boolean;
var
  s, t: AtomID;
  mm: BondID;
  PBnd: PRBond;
begin
  s := BndSet[m]^.h;
  t := BndSet[m]^.t;
  Result := RemoveBond(s, t);
end;

procedure TMoleculeBase.RemoveAtom(x: AtomID);
var
  i,j: integer;
  M,MM: BondID;
  aa,bb: integer;
  y: Node;
  bmod: array of array of boolean;
  PBnd: PRBond;
  t,h:integer;
  b1,b2: boolean;
begin
  bmod:=nil;
  SetLength(bmod,nAtom+1,nAtom+1);
  for i:=1 to nAtom do
    for j:=1 to nAtom do
      bmod[i,j]:=False;//A given pointer to a bond may appear several times, but must be modified only once.
  //Remove bonds from x
  M := p_HEAD[x];
  for MM := 1 to p_OutDeg[x] do
    RemoveBond(M);//Removed bond, the remaining bond indexes are left shifted
  //The atom to remove is disconnected
  dispose(fAtmSet[x]);
  for i := x to nAtom do
  begin
    p_HEAD[i] := p_HEAD[i + 1];
    fAtmSet[i] := fAtmSet[i + 1];
  end;
  nAtom:=nAtom-1;
  //All references over x have to be shifted to the left
  for M := 1 to nBonds do
    if (p_SUCC[M] > x) then
      p_SUCC[M] := p_SUCC[M] - 1;
  for M := 1 to nBonds do begin
    PBnd:=fBndSet[M];
    if PBnd<>nil then begin
      t:=PBnd^.t; h:=PBnd^.h;
      if (t>x) or (h>x) then begin
  //      if ((bmod[PBnd^.t,PBnd^.h]=False) or (bmod[PBnd^.h,PBnd^.t]=False)) then begin
        if ((bmod[t,h]=False) or (bmod[h,t]=False)) then begin
          if t>x then PBnd^.t:=t-1;
          if h>x then PBnd^.h:=h-1;
          bmod[PBnd^.t,PBnd^.h]:=True;
          bmod[PBnd^.h,PBnd^.t]:=True;//Do not touch these bonds (with updated indices) anymore
        end;
      end;
    end;
  end;
  SetLength(bmod,0,0);
end;

function TMoleculeBase.AddAtomPT: PRAtom;
var
  PAt: PRAtom;
  i: integer;
begin
  p_NX := p_NX + 1;
  p_HEAD[p_NX] := p_M + 1;
  new(PAt);
  PAt^.S := 'XX';
  PAt^.Z := ZMax;
  PAt^.W := 0;
  SetLength(PAt^.P, fAPrpSze);
  SetLength(PAt^.I, fBBytSze);
  for i := Low(PAt^.P) to High(PAt^.P) do
    PAt^.P[i] := 0;
  for i := Low(PAt^.I) to High(PAt^.I) do
    PAt^.I[i] := 0;
  fAtmSet[p_NX] := PAt;
  p_HEAD[p_NX + p_NY + 1] := p_M + 1;
  Result := PAt;
end;

function TMoleculeBase.AddAtomID: AtomID;
var
  PAt: PRAtom;
begin
  PAt := AddAtomPT;
  Result := FindAtom(PAt);
end;

function TMoleculeBase.AddBondDir(s, t: AtomID): PRBond;
var
  M, B: BondID;
  x: AtomID;
  PBo: PRBond;
  i: integer;
begin
  new(PBo);
  PBo^.t := s;
  PBo^.h := t;
  PBo^.S := 'YY';
  PBo^.B := 21;
  SetLength(PBo^.P, fBPrpSze);
  SetLength(PBo^.I, fBBytSze);
  for i := Low(PBo^.P) to High(PBo^.P) do
    PBo^.P[i] := 0;
  for i := Low(PBo^.I) to High(PBo^.I) do
    PBo^.I[i] := 0;
  p_M := p_M + 1;
  M := p_HEAD[s + 1];
  for B := p_M downto M + 1 do
    p_SUCC[B] := p_SUCC[B - 1];
  for x := s + 1 to p_NX + 1 do
    p_HEAD[x] := p_HEAD[x] + 1;
  for B := p_M downto M + 1 do
    fBndSet[B] := fBndSet[B - 1];
  p_SUCC[M] := t;
  fBndSet[M] := PBo;
  Result := PBo;
end;

procedure TMoleculeBase.AddBondPT(s, t: AtomID; var b1, b2: PRBond);
begin
  b1 := AddBondDir(s, t);
  b2 := AddBondDir(t, s);
end;

procedure TMoleculeBase.AddBondID(s, t: AtomID; var b1, b2: BondID);
var
  PB1, PB2: PRBond;
begin
  PB1:=nil; PB2:=nil;
  AddBondPT(s, t, PB1, PB2);
  b1 := FindBond(PB1);
  b2 := FindBond(PB2);
end;

procedure TMoleculeBase.SetAPrpSize(size: integer);
begin
  fAPrpSze := size;
end;

procedure TMoleculeBase.SetBPrpSize(size: integer);
begin
  fBPrpSze := size;
end;

procedure TMoleculeBase.SetABytSize(size: integer);
begin
  fABytSze := size;
end;

procedure TMoleculeBase.SetBBytSize(size: integer);
begin
  fBBytSze := size;
end;

function TMoleculeBase.FindAtom(PAt: PRAtom): AtomID;
var
  x: AtomID;
begin
  x := 1;
  while ((PAt <> fAtmSet[x]) and (x <= p_NX)) do
    Inc(x);
  if (x > p_NX) then
    raise Exception.Create('TMoleculeBase - ERROR: Atom not found');
  Result := x;
end;

function TMoleculeBase.FindBond(PBo: PRBond): BondID;
var
  M: BondID;
begin
  M := 1;
  while ((PBo <> fBndSet[M]) and (M <= p_M)) do
    Inc(M);
  if (M > p_M) then
    raise Exception.Create('TMoleculeBase - ERROR: Bond not found');
  Result := M;
end;

function TMoleculeBase.FindBond(t, h: AtomID): PRBond;
var
  M: BondID;
begin
  Result := nil;
  M := p_Head[t];
  while ((M < p_Head[t+1]) and (p_SUCC[M] <> h)) do
    Inc(M);
  if (M<p_Head[t+1]) and (p_SUCC[M] = h) then
    Result := fBndSet[M];
end;

function TMoleculeBase.FindBondID(t, h: AtomID): BondID;
var
  M: BondID;
begin
  Result := MaxArcNum;
  M := p_Head[t];
  while ((M < p_Head[t+1]) and (p_SUCC[M] <> h)) do
    Inc(M);
  if (M<p_Head[t+1]) and (p_SUCC[M] = h) then
    Result := M;
end;

function TMoleculeBase.IntToTB(col3, col7: integer): TB;
begin
  if (col7 = 0) then
  begin
    Result := col3; //for case 1 to 9
    case col3 of
      50: Result := 10;  //bond .
      60: Result := 11;  //bond :
      70: Result := 12;  //bond #
      80: Result := 13;  //bond ~
    end;
  end
  else if (col3 = 1) then
  begin
    case col7 of
      8: Result := 25;   //bond 81
      -1: Result := 18;  //bond 18
      4: Result := 14;   //bond 21
      12: Result := 28;  //bond 31
      1: Result := 31;   //bond 41
      {BUG in CGR interpretation. Original version here.
      8: Result := 14;   //bond 81
      -1: Result := 18;  //bond 18
      4: Result := 25;   //bond 21
      12: Result := 28;  //bond 31
      1: Result := 31;   //bond 41}
    end;
  end
  else if (col3 = 2) then
  begin
    case col7 of
      4: Result := 15;   //bond 82
      -1: Result := 19;  //bond 28
      8: Result := 22;   //bond 12
      12: Result := 29;  //bond 32
      1: Result := 32;   //bond 42
    end;
  end
  else if (col3 = 3) then
  begin
    case col7 of
      12: Result := 16;  //bond 83
      -1: Result := 20;  //bond 38
      8: Result := 23;   //bond 13
      4: Result := 26;   //bond 23
      1: Result := 33;   //bond 43
    end;
  end
  else if (col3 = 4) then
  begin
    case col7 of
      1: Result := 17;   //bond 84
      -1: Result := 21;  //bond 48
      8: Result := 24;   //bond 14
      4: Result := 27;   //bond 24
      12: Result := 30;  //bond 34
    end;
  end
  else
  begin
    Result := BMax;
  end;
end;

function TMoleculeBase.IntToS(col3, col7: integer): TS;
begin
  Result := BondSymbol[IntToTB(col3, col7)];
end;

end.

