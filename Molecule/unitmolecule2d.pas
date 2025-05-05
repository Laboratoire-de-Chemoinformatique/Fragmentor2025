unit UnitMolecule2D;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, U_TYPE_GRAPHES, U_GRAPHES, UnitMoleculeBase,
  unitAtomAndBondType, contnrs;

type

  { TMolecule2D }

  TMolecule2D = class(TMoleculeBase)
    private

    public
      procedure Clear; override;
      procedure SaveSDF(sdfout: string);
      procedure LoadSDF(sdfstr: TStringList);
      function GetX(AIndex: AtomID): double;
      function GetY(AIndex: AtomID): double;
      function GetZ(AIndex: AtomID): double;
      procedure SetX(AIndex: AtomID; AValue: double);
      procedure SetY(AIndex: AtomID; AValue: double);
      procedure SetZ(AIndex: AtomID; AValue: double);
      constructor Create;
      destructor Destroy; override;
    end;

implementation

{ TMolecule2D }

procedure TMolecule2D.Clear;
begin
  inherited Clear;
  APrpSze := APrpSze + 3;//for X,Y,Z
  ABytSze := ABytSze + 10;//1 Byte for formal charge
  BBytSze := BBytSze + 3;//for bond direction (up/down) and orientation
end;

procedure TMolecule2D.SaveSDF(sdfout: string);
var
  SLtmp: TStringList;
  stmp: string;
  PAt: PRAtom;
  PBd: PRBond;
//  btmp: array of array of integer;
  btmp: array of array of boolean;
  i,j,LineNo: integer;
  s,t: Node;
  ii,jj:integer;
begin
  SLtmp := TStringList.Create;
  //------Line 1 of MOLfile------------------------------------------------------
  LineNo := 0;
  SLtmp.Add(MolName);
  //------Line 2 of MOLfile------------------------------------------------------
  //  User's first and last initials (UserI), program name (ProgN), date/time:
  //  M/D/Y,H:m (DatTim), dimensional codes (DimCo), scaling factors (ScaI, ScaR),
  //  energy (Ene) if modeling program input, internal registry number (RegNo)
  //  if input through MDL format:
  Inc(LineNo);
  SLtmp.Add('GMISIDA'+DateToStr(Date)+'111NA'+MolName);
  //------Line 3 of MOLfile------------------------------------------------------
  // A line for comments. If no comment is entered, a blank line must be present.
  Inc(LineNo);
  SLtmp.Add('');
  //------Line 4 of MOLfile : The Counts Line------------------------------------
  Inc(LineNo);
  stmp:=Format('%0:3d%1:3d',[nAtom,round(nBonds/2)])+'  0  0  0  0             999 V2000';
  SLtmp.Add(stmp);
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
  for i := 1 to p_NX do
  begin
    PAt:=AtmSet[i];
    Inc(LineNo);
    stmp:=Format('%0:10.4F%1:10.4F%2:10.4F',[GetX(i),GetY(i),GetZ(i)]);
    stmp:=stmp+' '+Format('%-3S',[PAt^.S]);        // atom symbol
    j:=0;
    stmp:=stmp+Format('%3D',[PAt^.I[j]]); Inc(j);     //mass difference
    stmp:=stmp+Format('%3D',[PAt^.I[j]]); Inc(j);     //charge
    stmp:=stmp+Format('%3D',[PAt^.I[j]]); Inc(j);     //atom stereo parity
    stmp:=stmp+Format('%3D',[PAt^.I[j]]); Inc(j);     //hydrogen count + 1
    stmp:=stmp+Format('%3D',[PAt^.I[j]]); Inc(j);     //stereo care box
    stmp:=stmp+Format('%3D',[PAt^.I[j]]); Inc(j);     //valence
    stmp:=stmp+Format('%3D',[PAt^.I[j]]); Inc(j);     //H0 designator
    stmp:=stmp+'  0  0';
    stmp:=stmp+Format('%3D',[PAt^.I[j]]); Inc(j);     //atom-atom mapping number
    stmp:=stmp+Format('%3D',[PAt^.I[j]]); Inc(j);     //inversion/retention flag
    stmp:=stmp+Format('%3D',[PAt^.I[j]]); Inc(j);     //exact change flag
    SLtmp.Add(stmp);
  end;  // 1..p_NX
  //------Lines of MOLfile : The Bond Block--------------------------------------
  //ii:=natom;
  SetLength(btmp,nAtom+1,nAtom+1);
  for i:=1 to nAtom do
    for j:=1 to nAtom do
      btmp[i,j]:=False;
  for i:=1 to nBonds do
  begin
    PBd:=BndSet[i];
    if PBd<>nil then
      if btmp[PBd^.h,PBd^.t]=False then
      begin
        stmp:=Format('%0:3D%1:3D%2:3D%3:3D',[PBd^.t,PBd^.h,PBd^.B,PBd^.I[0]])+
        '  0'+Format('%0:3D%1:3D',[PBd^.I[1],PBd^.I[2]]);
        SLtmp.Add(stmp);
        btmp[PBd^.h,PBd^.t]:=True;
        btmp[PBd^.t,PBd^.h]:=True;
      end;
  end;
  //End molecule
  SLtmp.Add('M  END');
  SLtmp.SaveToFile(sdfout);
  FreeAndNil(SLtmp);
end;

procedure TMolecule2D.LoadSDF(sdfstr: TStringList);
const
  Np = 10;
var
  sdfstrsize: integer;
  i, j, k: integer;
  LineNo: integer;
  CodeErr: integer;
  atmL, cf: integer;
  stxt, ap: integer;
  wrd: string;
  PAt: PRAtom;
  PBo: PRBond;
  atmp: CostMatrix;
  wtmp: array of PRBond;
  s, t: Node;
  M: ArcNum;
begin
  Clear;
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
    SetLength(PAt^.P, APrpSze); //set the atom properties array size (double)
    SetLength(PAt^.I, ABytSze); //set the atom properties array size (Byte)
    for j := Low(PAt^.P) to High(PAt^.P) do
      PAt^.P[j] := 0;
    for j := Low(PAt^.I) to High(PAt^.I) do
      PAt^.I[j] := 0;
    j := 0;
    begin
      PAt^.P[j] := dbl_readpos(sdfstr[LineNo], 1, 10); //X
      Inc(j);
    end;
    begin
      PAt^.P[j] := dbl_readpos(sdfstr[LineNo], 11, 10);//Y
      Inc(j);
    end;
    begin
      PAt^.P[j] := dbl_readpos(sdfstr[LineNo], 21, 10);//Z
      Inc(j);
    end;
    PAt^.S := Trim(Copy(sdfstr[LineNo], 32, 3));        // atom symbol
    PAt^.Z := AtomSymbolToInt(PAt^.S);                // atomic number
    j := 0;
    begin
      PAt^.I[j]:=int_readpos(sdfstr[LineNo],35,2); //mass difference
      Inc(j);
    end;
    begin
      PAt^.I[j] := int_readpos(sdfstr[LineNo], 37, 3); //charge;
      Inc(j);
    end;
    begin
      PAt^.I[j]:=int_readpos(sdfstr[LineNo],40,3); //atom stereo parity
      Inc(j);
    end;
    begin
      PAt^.I[j] := int_readpos(sdfstr[LineNo], 43, 3); //hydrogen count + 1
      Inc(j);
    end;
    begin
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
    end;
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
    SetLength(PBo^.P, BPrpSze); //set the bond properties array size (double)
    SetLength(PBo^.I, BBytSze); //set the bond properties array size (Byte)
    PBo^.t := int_readpos(sdfstr[LineNo], 1, 3);//Tail of bond
    PBo^.h := int_readpos(sdfstr[LineNo], 4, 3);//Head of bond
    PBo^.B := IntToTB(int_readpos(sdfstr[LineNo],7,3),int_readpos(sdfstr[LineNo],19,3));//Using col3 and col7
    PBo^.S := BondSymbol[PBo^.B];
    j := 0;
    begin
      PBo^.I[j] := int_readpos(sdfstr[LineNo], 10, 3);//bond orientation
      Inc(j);
    end;
    begin
      PBo^.I[j]:=int_readpos(sdfstr[LineNo],16,3); //xxx
      Inc(j);
    end;
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

function TMolecule2D.GetX(AIndex: AtomID): double;
begin
  Result := AtmSet[AIndex]^.P[0];
end;

function TMolecule2D.GetY(AIndex: AtomID): double;
begin
  Result := AtmSet[AIndex]^.P[1];
end;

function TMolecule2D.GetZ(AIndex: AtomID): double;
begin
  Result := AtmSet[AIndex]^.P[2];
end;

procedure TMolecule2D.SetX(AIndex: AtomID; AValue: double);
begin
  AtmSet[AIndex]^.P[0] := AValue;
end;

procedure TMolecule2D.SetY(AIndex: AtomID; AValue: double);
begin
  AtmSet[AIndex]^.P[1] := AValue;
end;

procedure TMolecule2D.SetZ(AIndex: AtomID; AValue: double);
begin
  AtmSet[AIndex]^.P[2] := AValue;
end;

constructor TMolecule2D.Create;
begin
  inherited Create;
  APrpSze := APrpSze + 3;//for X,Y,Z
  ABytSze := ABytSze + 10;//1 Byte for formal charge
  BBytSze := BBytSze + 3;//for bond direction (up/down) and orientation
end;

destructor TMolecule2D.Destroy;
begin
  inherited Destroy;
end;

end.

