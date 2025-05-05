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
unit UnitMolecule;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,U_TYPE_GRAPHES,U_GRAPHES,unitAtomAndBondType;//, contnrs;

const
     MaxAtom = MaxNode;
     MaxBond = MaxArcNum/2;

type
    AtomID=Node;
    BondID=ArcNum;
    //
    AAtmSet = array[AtomID] of PRAtom;
    ABndSet = array[BondID] of PRBond;
    //

    { TMolecule }

    TMolecule = class(T_GRAPHE_LISTE)
    private
           fMolName: string;
           fAtmSet: AAtmSet;
           fBndSet: ABndSet;
           fMPrp: APrp;
           fMByt: AByt;
           fMPrpSze, fMBytSze: integer;
           fAPrpSze, fBPrpSze: integer;
           fABytSze, fBBytSze: integer;
           fKeepMolInfo, fKeepCrdInfo, fKeepAtmInfo, fKeepBndInfo: TBits;//Sze 4; 3; 6; 2
           procedure DisposeAtmSet;
           procedure DisposeBndSet;
           function GetAtom(i: AtomID): PRAtom;
           procedure SetAtom(i: AtomID; PAt: PRAtom);
           function GetBond(i: BondID): PRBond;
           procedure SetBond(i: BondID; PBo: PRBond);
           function GetMPrp(i: integer): double;
           procedure SetMPrp(i: integer; prp: double);
           function GetMByt(i: integer): byte;
           procedure SetMByt(i: integer; c: byte);
           procedure PermuteBnd(i, j: BondID);
           procedure PermuteAtm(i, j: AtomID);
    public
          constructor Create;
          destructor Destroy; override;
          procedure Clear;
          procedure LoadSDF(sdfstr: TStringList);
          function RemoveBondDir(PBo: PRBond): boolean;
          function RemoveBondDir(s,t:AtomID):boolean;
          function RemoveBond(s,t:AtomID):boolean;
          function RemoveBond(m:BondID):boolean;
          procedure RemoveAtom(x: AtomID);
          function AddAtomPT: PRAtom;
          function AddAtomID: AtomID;
          function AddBondDir(s,t: AtomID): PRBond;
          procedure AddBondPT(s,t: AtomID;var b1,b2: PRBond);
          procedure AddBondID(s,t: AtomID;var b1,b2: BondID);
          function FindAtom(PAt: PRAtom): AtomID;
          function FindBond(PBo: PRBond): BondID;
          function FindBond(t,h: AtomID): PRBond;
          function FindBondID(t,h: AtomID): BondID;
          property MolName: string read fMolName write fMolName;
          property MPrpSze: integer read fMPrpSze write fMPrpSze;
          property MBytSze: integer read fMBytSze write fMBytSze;
          property APrpSze: integer read fAPrpSze write fAPrpSze;
          property ABytSze: integer read fABytSze write fABytSze;
          property BPrpSze: integer read fBPrpSze write fBPrpSze;
          property BBytSze: integer read fBBytSze write fBBytSze;
          property KeepMolInfo: TBits read fKeepMolInfo write fKeepMolInfo;
          property KeepCrdInfo: TBits read fKeepCrdInfo write fKeepCrdInfo;
          property KeepAtmInfo: TBits read fKeepAtmInfo write fKeepAtmInfo;
          property KeepBndInfo: TBits read fKeepBndInfo write fKeepBndInfo;
          property AtmSet[i: AtomID]:PRAtom read GetAtom write SetAtom;
          property BndSet[i: BondID]:PRBond read GetBond write SetBond;
          property MPrp[i: integer]:double read GetMPrp write SetMPrp;
          property MByt[i: integer]:byte read GetMByt write SetMByt;
    end;

    function CountOn(BV: TBits): integer;

implementation

function CountOn(BV: TBits): integer;
var
   i:integer;
begin
     Result:=0;
     for i:=0 to BV.Size-1 do if (BV.Bits[i]) then inc(Result);
end;

{ TMolecule }

procedure TMolecule.DisposeAtmSet;
var
   i: integer;
begin
     for i:=Low(fAtmSet) to p_NX do if (fAtmSet[i]<>nil) then begin
         dispose(fAtmSet[i]);
         fAtmSet[i]:=nil;
     end;
end;

procedure TMolecule.DisposeBndSet;
var
   i,j: integer;
begin
     for i:=Low(fBndSet) to p_M do begin
         //writeln('*'+IntToStr(i)+'*');
     //end;
     //for i:=Low(fBndSet) to p_M do begin
         //search for multiple copies of a bond.
         for j:=i+1 to p_M do
             if (fBndSet[i]=fBndSet[j]) then fBndSet[j]:=nil;
         if (fBndSet[i]<>nil) then begin
            dispose(fBndSet[i]);
            fBndSet[i]:=nil;
         end;
     end;
end;

function TMolecule.GetAtom(i: AtomID): PRAtom;
//Beware, the owner does not own the pointer to the RAtom record!
begin
     Result:=fAtmSet[i];
end;

procedure TMolecule.SetAtom(i: AtomID; PAt: PRAtom);
//Beware, the caller does not anymore owns the pointer to the RAtom record!
begin
     if (fAtmSet[i]<>nil) then dispose(fAtmSet[i]);
     fAtmSet[i]:=PAt;
end;

function TMolecule.GetBond(i: BondID): PRBond;
//Beware, the owner does not own the pointer to the RBond record!
begin
     Result:=fBndSet[i];
end;

procedure TMolecule.SetBond(i: BondID; PBo: PRBond);
//Beware, the caller does not anymore owns the pointer to the RBond record!
begin
     if (fBndSet[i]<>nil) then dispose(fBndSet[i]);
     fBndSet[i]:=PBo;
end;

function TMolecule.GetMPrp(i: integer): double;
begin
     Result:=fMPrp[i];
end;

procedure TMolecule.SetMPrp(i: integer; prp: double);
begin
     fMPrp[i]:=prp;
end;

function TMolecule.GetMByt(i: integer): byte;
begin
     Result:=fMByt[i];
end;

procedure TMolecule.SetMByt(i: integer; c: byte);
begin
     fMByt[i]:=c;
end;

procedure TMolecule.PermuteBnd(i, j: BondID);
var
   PBo: PRBond;
begin
     PBo:=fBndSet[i];
     fBndSet[i]:=fBndSet[j];
     fBndSet[j]:=PBo;
end;

procedure TMolecule.PermuteAtm(i, j: AtomID);
var
   PAt: PRAtom;
begin
     PAt:=fAtmSet[i];
     fAtmSet[i]:=fAtmSet[j];
     fAtmSet[j]:=PAt;
end;

constructor TMolecule.Create;
var
   i: AtomID;
   j: BondID;
   n: integer;
begin
     fMolName:='NoName';
     fMPrpSze:=0;
     fMBytSze:=0;
     fAPrpSze:=0;
     fABytSze:=0;
     fBPrpSze:=0;
     fBBytSze:=0;
     fKeepMolInfo:=TBits.Create(4);
     fKeepCrdInfo:=TBits.Create(3);
     fKeepAtmInfo:=TBits.Create(6);
     fKeepBndInfo:=TBits.Create(2);
     for i:=Low(fAtmSet) to High(fAtmSet) do fAtmSet[i]:=nil;
     for j:=Low(fBndSet) to High(fBndSet) do fBndSet[j]:=nil;
     inherited Create;
end;

destructor TMolecule.Destroy;
begin
     DisposeAtmSet;
     DisposeBndSet;
     FreeAndNil(fKeepMolInfo);
     FreeAndNil(fKeepCrdInfo);
     FreeAndNil(fKeepAtmInfo);
     FreeAndNil(fKeepBndInfo);
     inherited Destroy;
end;

procedure TMolecule.Clear;
begin
     fMolName:='NoName';
     fMPrpSze:=0;
     fMBytSze:=0;
     fAPrpSze:=0;
     fABytSze:=0;
     fBPrpSze:=0;
     fBBytSze:=0;
     fKeepMolInfo.Clearall;
     fKeepCrdInfo.Clearall;
     fKeepAtmInfo.Clearall;
     fKeepBndInfo.Clearall;
     DisposeAtmSet;
     DisposeBndSet;
     p_M:=0;
     p_NX:=0;
     p_NY:=0;
end;

procedure TMolecule.LoadSDF(sdfstr: TStringList);
const
     Np = 10;
var
   i, j, k : Integer;
   LineNo  : Integer;
   CodeErr : Integer;
   wrd     : string;
   PAt     : PRAtom;
   PBo     : PRBond;
   atmp    : CostMatrix;
   wtmp    : array of PRBond;
   s,t     : Node;
   M       : ArcNum;
   function int_readpos(str: string; bgn, lng: integer): integer;
     var
        CodeErr: integer;
        wrd: string;
     begin
          Result:=0;
          CodeErr := 0;  // fix error in input SDF file
          wrd := Trim(Copy(str,bgn,lng));
          if Length(wrd)>0 then Val(wrd,Result,CodeErr);                         // number of atoms
          if CodeErr<>0 then begin Exception.Create('ERROR: reading sdf entry '+fMolName); halt(1); end;
     end;
   function dbl_readpos(str: string; bgn, lng: integer): double;
     var
        CodeErr: integer;
        wrd: string;
     begin
          Result:=0;
          CodeErr := 0;  // fix error in input SDF file
          wrd := Trim(Copy(str,bgn,lng));
          if Length(wrd)>0 then Val(wrd,Result,CodeErr);                         // number of atoms
          if CodeErr<>0 then begin Exception.Create('ERROR: reading sdf entry '+fMolName); halt(1); end;
     end;
begin
     // set tables sizes - increment size based on possible user's request
     fMBytSze:=fMBytSze+CountOn(fKeepMolInfo);
     fAPrpSze:=fAPrpSze+CountOn(fKeepCrdInfo);
     fABytSze:=fABytSze+CountOn(fKeepAtmInfo);
     fBBytSze:=fBBytSze+CountOn(fKeepBndInfo);
     //------Line 1 of MOLfile------------------------------------------------------
     LineNo:=0;
     fMolName := sdfstr[LineNo];
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
     p_NY:= 0; //Molecular graphs are not bipartite a priori
     p_M := 0;
     SetLength(fMByt,fMBytSze);
     p_NX:=int_readpos(sdfstr[LineNo],1,3);     // number of atoms
     p_M:=int_readpos(sdfstr[LineNo],4,3);      // number of bonds
     j:=0;
     if fKeepMolInfo[0] then begin fMByt[j]:=int_readpos(sdfstr[LineNo],7,3); inc(j); end;   // number of atom lists
                                                 //position 10 is obsolete
     if fKeepMolInfo[1] then begin fMByt[j]:=int_readpos(sdfstr[LineNo],13,3); inc(j); end;  //chiral flag
     if fKeepMolInfo[2] then begin fMByt[j]:=int_readpos(sdfstr[LineNo],16,3); inc(j); end;  //number of stext entries (ISIS)
                                                 //position 19 is obsolete
                                                 //position 22 is obsolete
                                                 //position 25 is obsolete
                                                 //position 28 is obsolete
     if fKeepMolInfo[3] then begin fMByt[j]:=int_readpos(sdfstr[LineNo],31,3); inc(j); end;  // number of lines of additional properties
     //------Lines 5-... of MOLfile : The Atom Block--------------------------------
     fAtmSet[0]:=nil;
     for i:=1 to p_NX do begin
         Inc(LineNo);
         new(PAt);
         SetLength(PAt^.P,fAPrpSze); SetLength(PAt^.I,fABytSze);
         for j:=Low(PAt^.P) to High(PAt^.P) do PAt^.P[j]:=0;
         for j:=Low(PAt^.I) to High(PAt^.I) do PAt^.I[j]:=0;
         j:=0;
         if fKeepCrdInfo[0] then begin
            PAt^.P[j]:=dbl_readpos(sdfstr[LineNo],1,10); //X
            inc(j);
         end;
         if fKeepCrdInfo[1] then begin
            PAt^.P[j]:=dbl_readpos(sdfstr[LineNo],11,10);//Y
            inc(j);
         end;
         if fKeepCrdInfo[2] then begin
            PAt^.P[j]:=dbl_readpos(sdfstr[LineNo],21,10);//Z
            inc(j);
         end;
         PAt^.S:=Trim(Copy(sdfstr[LineNo],32,3));        // atom symbol
         j:=0;
         if fKeepAtmInfo[0] then begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],35,2);
            Inc(j);
         end;
         if fKeepAtmInfo[1] then begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],37,3);
            Inc(j);
         end;
         if fKeepAtmInfo[2] then begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],40,3);
            Inc(j);
         end;
         if fKeepAtmInfo[3] then begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],43,3);
            Inc(j);
         end;
         if fKeepAtmInfo[4] then begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],46,3);
            Inc(j);
         end;
         if fKeepAtmInfo[5] then begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],49,3);
            Inc(j);
         end;
         AtmSet[i]:=PAt;
     end;  // 1..p_NX
     //------Lines of MOLfile : The Bond Block--------------------------------------
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
         j:=0;
         if fKeepBndInfo[0] then begin
            PBo^.I[j]:=int_readpos(sdfstr[LineNo],10,3);
            inc(j);
         end;
         if fKeepBndInfo[1] then begin
            PBo^.I[j]:=int_readpos(sdfstr[LineNo],16,3);
            inc(j);
         end;
         atmp[PBo^.t,PBo^.h]:=i; atmp[PBo^.h,PBo^.t]:=i;
         wtmp[i]:=PBo;
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
end;

function TMolecule.RemoveBondDir(PBo: PRBond): boolean;
begin
     Result:=RemoveBondDir(PBo^.t,PBo^.h);
end;

function TMolecule.RemoveBondDir(s, t: AtomID): boolean;
var
   i: integer;
   M,N: BondID;
   x: AtomID;
begin
     Result:=False;
     M:=p_Head[s];
     i:=0;
     while ((i<p_OutDeg[s]) and (p_SUCC[M+i]<>t)) do Inc(i);
     M:=M+i;
     if (p_SUCC[M]=t) then begin
        Result:=True;
        dispose(fBndSet[M]);
        for N:=M to p_M-1 do begin
            p_SUCC[N]:=p_SUCC[N+1];
            fBndSet[N]:=fBndSet[N+1];
        end;
        p_SUCC[p_M]:=0;
        fBndSet[p_M]:=nil;
        for x:=s+1 to GraphOrder+1 do p_HEAD[x]:=p_HEAD[x]-1;
        p_M:=p_M-1;
     end;
     p_HEAD[p_NX+p_NY+1]:=p_M+1;
end;

function TMolecule.RemoveBond(s, t: AtomID): boolean;
var
   b1,b2:boolean;
begin
     Result:=False;
     b1:=RemoveBondDir(s,t);
     b2:=RemoveBondDir(t,s);
     if (b1 and b2) then Result:=True;
end;

function TMolecule.RemoveBond(m: BondID): boolean;
var
   s,t:AtomID;
begin
     s:=BndSet[m]^.h; t:=BndSet[m]^.t;
     Result:=RemoveBond(s,t);
end;

procedure TMolecule.RemoveAtom(x: AtomID);
var
   i: integer;
   M: BondID;
begin
     M:=p_HEAD[x];
     for i:=1 to p_OutDeg[x] do RemoveBond(M);
     for i:=x to p_NX do p_HEAD[i]:=p_HEAD[i+1];
     for M:=1 to p_M do begin
         if (p_SUCC[M]>x) then p_SUCC[M]:=p_SUCC[M]-1;
         if (fBndSet[M]^.t>x) then fBndSet[M]^.t:=fBndSet[M]^.t-1;
         if (fBndSet[M]^.h>x) then fBndSet[M]^.h:=fBndSet[M]^.h-1;
     end;
     dispose(fAtmSet[x]);
     for i:=x to p_NX do fAtmSet[i]:=fAtmSet[i+1];
     fAtmSet[p_NX]:=nil;
     p_HEAD[p_NX+p_NY+1]:=p_M+1;
     p_NX:=p_NX-1;
end;

function TMolecule.AddAtomPT: PRAtom;
var
   PAt: PRAtom;
   i: integer;
begin
     p_NX:=p_NX+1;
     p_HEAD[p_NX]:=p_M+1;
     new(PAt);
     PAt^.S:='XX';
     PAt^.Z:=ZMax;
     PAt^.W:=0;
     SetLength(PAt^.P,fAPrpSze); SetLength(PAt^.I,fBBytSze);
     for i:=Low(PAt^.P) to High(PAt^.P) do PAt^.P[i]:=0;
     for i:=Low(PAt^.I) to High(PAt^.I) do PAt^.I[i]:=0;
     fAtmSet[p_NX]:=PAt;
     p_HEAD[p_NX+p_NY+1]:=p_M+1;
     Result:=PAt;
end;

function TMolecule.AddAtomID: AtomID;
var
   PAt: PRAtom;
begin
     PAt:=AddAtomPT;
     Result:=FindAtom(PAt);
end;

function TMolecule.AddBondDir(s, t: AtomID): PRBond;
var
   M,B: BondID;
   x: AtomID;
   PBo: PRBond;
   i: integer;
begin
     new(PBo);
     PBo^.t:=s;PBo^.h:=t;
     PBo^.S:='YY';
     PBo^.B:=21;
     SetLength(PBo^.P,fBPrpSze); SetLength(PBo^.I,fBBytSze);
     for i:=Low(PBo^.P) to High(PBo^.P) do PBo^.P[i]:=0;
     for i:=Low(PBo^.I) to High(PBo^.I) do PBo^.I[i]:=0;
     p_M:=p_M+1;
     M:=p_HEAD[s+1];
     for B:=p_M downto M+1 do p_SUCC[B]:=p_SUCC[B-1];
     for x:=s+1 to p_NX+1 do p_HEAD[x]:=p_HEAD[x]+1;
     for B:=p_M downto M+1 do fBndSet[B]:=fBndSet[B-1];
     p_SUCC[M]:=t; fBndSet[M]:=PBo;
     Result:=PBo;
end;

procedure TMolecule.AddBondPT(s, t: AtomID;var b1,b2: PRBond);
begin
     b1:=AddBondDir(s,t);
     b2:=AddBondDir(t,s);
end;

procedure TMolecule.AddBondID(s, t: AtomID; var b1, b2: BondID);
var
   PB1,PB2: PRBond;
begin
     AddBondPT(s,t,PB1,PB2);
     b1:=FindBond(PB1);
     b2:=FindBond(PB2);
end;

function TMolecule.FindAtom(PAt: PRAtom): AtomID;
var
   x: AtomID;
begin
     x:=1;
     while ((PAt<>fAtmSet[x]) and (x<=p_NX)) do Inc(x);
     if (x>p_NX) then raise Exception.Create('TMolecule - ERROR: Atom not found');
     Result:=x;
end;

function TMolecule.FindBond(PBo: PRBond): BondID;
var
   M: BondID;
begin
     M:=1;
     while ((PBo<>fBndSet[M]) and (M<=p_M)) do Inc(M);
     if (M>p_M) then raise Exception.Create('TMolecule - ERROR: Bond not found');
     Result:=M;
end;

function TMolecule.FindBond(t, h: AtomID): PRBond;
var
   i: Node;
   M: BondID;
begin
     Result:=nil;
     i:=0;
     M:=p_Head[t];
     while ((i<p_OutDeg[t]) and (p_SUCC[M+i]<>h)) do Inc(i);
     M:=M+i;
     if (p_SUCC[M]=h) then Result:=fBndSet[M];
end;

function TMolecule.FindBondID(t, h: AtomID): BondID;
var
   i: Node;
   M: BondID;
begin
     Result:=MaxArcNum;
     i:=0;
     M:=p_Head[t];
     while ((i<p_OutDeg[t]) and (p_SUCC[M+i]<>h)) do Inc(i);
     M:=M+i;
     if (p_SUCC[M]=h) then Result:=M;
end;

end.

