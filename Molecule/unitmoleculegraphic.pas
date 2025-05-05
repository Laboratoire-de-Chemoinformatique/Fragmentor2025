unit UnitMoleculeGraphic;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, UnitMoleculeBase, unitAtomAndBondType, FPimage,
  Graphics, FPCanvas, U_TYPE_GRAPHES, U_GRAPHES;

const

  // MOLECULAR CONSTANTS

  MaxVale = 5;        //  Dimension of atom valence array
  SumElem = 103;      //  Sum of elements for periodic table
  MaxSumMol = 100000;   //  max sum of molecules in SDF file

type

  { ELEMENT TYPE }

  ELEMENT = record                            {  Chemical elements  }
    ES: string[2];                         {  Element symbol  }
    EN: string[12];                        {  Element name  }
    V: array[1..5] of byte;               {  Valence  }
    P: byte;                              {  Period in Mendeleev Table  }
    AW: real;                              {  Atom weit  }
  end;

  { TMoleculeGraphic }

  TMoleculeGraphic = class(TMoleculeBase)
  private
    fAFontOffset, fAStyleOffSet, fBStyleOffSet: integer;
    fDefaultFontSize: byte;
    fDefaultFontCol: TFPColor;
    procedure GetGeoCenter(var GX, GY, GZ: double);
  public
    constructor Create;
    destructor Destroy; override;
    procedure Reset;
    procedure Clear;
    procedure LoadSDF(sdfstr: TStringList);
    procedure SetAtomFont(AIndex: AtomID; AFontColor: TFPColor; AFontSize: byte);
    procedure SetAtomFont(PAt: PRAtom; AFontColor: TFPColor; AFontSize: byte);
    procedure SetAtomPenAndBrush(AIndex: AtomID; APenColor: TFPColor;
      APenStyle: TFPPenStyle; ABrushColor: TFPColor; ABrushStyle: TFPBrushStyle);
    procedure SetBondPenAndBrush(AIndex: BondID; APenColor: TFPColor;
      APenStyle: TFPPenStyle; ABrushColor: TFPColor; ABrushStyle: TFPBrushStyle);
    procedure GetAtomFont(AIndex: AtomID; out AFontColor: TFPColor; out AFontSize: byte);
    procedure GetAtomFont(PAt: PRAtom; out AFontColor: TFPColor; out AFontSize: byte);
    procedure GetAtomPenAndBrush(AIndex: AtomID; APenColor: TFPColor;
      APenStyle: TFPPenStyle; ABrushColor: TFPColor; ABrushStyle: TFPBrushStyle);
    procedure GetBondPenAndBrush(AIndex: BondID; APenColor: TFPColor;
      APenStyle: TFPPenStyle; ABrushColor: TFPColor; ABrushStyle: TFPBrushStyle);
    procedure SetDefaultStyle();
    function GetX(AIndex: AtomID): double;
    function GetY(AIndex: AtomID): double;
    function GetZ(AIndex: AtomID): double;
    procedure SetX(AIndex: AtomID; AValue: double);
    procedure SetY(AIndex: AtomID; AValue: double);
    procedure SetZ(AIndex: AtomID; AValue: double);
    procedure Center(CX, CY, CZ: double);
    procedure Scale(SX, SY, SZ: double);
    function DimX(): double;
    function DimY(): double;
    function DimZ(): double;
    function DimMaxXY(): double;
    function DimMaxXYZ(): double;
    function HCount(AIndex: AtomID): byte;
    property DefaultFontSize: byte read fDefaultFontSize write fDefaultFontSize;
    property DefaultFontCol: TFPColor read fDefaultFontCol write fDefaultFontCol;
  end;

function SDFChargeToInt(c: byte): integer;
{  Read the information about the elements  }
const
  Elem: array[1..SumElem] of ELEMENT = (
    (ES: 'H'; EN: 'Hydrogen'; V: (0, 1, 0, 0, 0); P: 1; AW: 1.0079),
    (ES: 'He'; EN: 'Helium'; V: (0, 0, 0, 0, 0); P: 0; AW: 4.0026),
    (ES: 'Li'; EN: 'Lithium'; V: (0, 1, 0, 0, 0); P: 1; AW: 6.941),
    (ES: 'Be'; EN: 'Beryllium'; V: (0, 2, 0, 0, 0); P: 2; AW: 9.0122),
    (ES: 'B'; EN: 'Boron'; V: (0, 3, 0, 0, 0); P: 3; AW: 10.811),
    (ES: 'C'; EN: 'Carbon'; V: (0, 4, 0, 0, 0); P: 4; AW: 12.011),
    (ES: 'N'; EN: 'Nitrogen'; V: (0, 3, 5, 0, 0); P: 5; AW: 14.0067),
    (ES: 'O'; EN: 'Oxygen'; V: (0, 2, 0, 0, 0); P: 6; AW: 15.9994),
    (ES: 'F'; EN: 'Fluorine'; V: (0, 1, 0, 0, 0); P: 7; AW: 18.9984),
    (ES: 'Ne'; EN: 'Neon'; V: (0, 0, 0, 0, 0); P: 0; AW: 20.179),
    (ES: 'Na'; EN: 'Sodium'; V: (0, 1, 0, 0, 0); P: 1; AW: 22.9898),
    (ES: 'Mg'; EN: 'Magnesium'; V: (0, 2, 0, 0, 0); P: 2; AW: 24.305),
    (ES: 'Al'; EN: 'Aluminum'; V: (0, 3, 0, 0, 0); P: 3; AW: 26.9815),
    (ES: 'Si'; EN: 'Silicon'; V: (0, 4, 0, 0, 0); P: 4; AW: 28.0855),
    (ES: 'P'; EN: 'Phosphorus'; V: (0, 3, 5, 0, 0); P: 5; AW: 30.9738),
    (ES: 'S'; EN: 'Sulfur'; V: (0, 2, 4, 6, 0); P: 6; AW: 32.066),
    (ES: 'Cl'; EN: 'Chlorine'; V: (0, 1, 3, 5, 7); P: 7; AW: 35.453),
    (ES: 'Ar'; EN: 'Argon'; V: (0, 0, 0, 0, 0); P: 0; AW: 39.948),
    (ES: 'K'; EN: 'Potassium'; V: (0, 1, 0, 0, 0); P: 1; AW: 39.0983),
    (ES: 'Ca'; EN: 'Calcium'; V: (0, 2, 0, 0, 0); P: 2; AW: 40.078),
    (ES: 'Sc'; EN: 'Scandium'; V: (0, 3, 0, 0, 0); P: 3; AW: 44.9559),
    (ES: 'Ti'; EN: 'Titanium'; V: (0, 3, 4, 0, 0); P: 4; AW: 47.88),
    (ES: 'V'; EN: 'Vanadium'; V: (0, 2, 3, 4, 5); P: 5; AW: 50.9415),
    (ES: 'Cr'; EN: 'Chromium'; V: (0, 2, 3, 6, 0); P: 6; AW: 51.996),
    (ES: 'Mn'; EN: 'Manganese'; V: (0, 2, 3, 4, 6); P: 7; AW: 54.938),
    (ES: 'Fe'; EN: 'Iron'; V: (0, 2, 3, 4, 6); P: 8; AW: 55.847),
    (ES: 'Co'; EN: 'Cobalt'; V: (0, 2, 3, 0, 0); P: 8; AW: 58.9332),
    (ES: 'Ni'; EN: 'Nickel'; V: (0, 2, 3, 0, 0); P: 8; AW: 58.69),
    (ES: 'Cu'; EN: 'Copper'; V: (0, 1, 2, 0, 0); P: 1; AW: 63.546),
    (ES: 'Zn'; EN: 'Zinc'; V: (0, 2, 0, 0, 0); P: 2; AW: 65.39),
    (ES: 'Ga'; EN: 'Gallium'; V: (0, 3, 0, 0, 0); P: 3; AW: 69.723),
    (ES: 'Ge'; EN: 'Germanium'; V: (0, 2, 4, 0, 0); P: 4; AW: 72.59),
    (ES: 'As'; EN: 'Arsenic'; V: (0, 3, 5, 0, 0); P: 5; AW: 74.9216),
    (ES: 'Se'; EN: 'Selenium'; V: (0, 2, 4, 6, 0); P: 6; AW: 78.96),
    (ES: 'Br'; EN: 'Bromine'; V: (0, 1, 3, 5, 7); P: 7; AW: 79.904),
    (ES: 'Kr'; EN: 'Krypton'; V: (0, 0, 0, 0, 0); P: 0; AW: 83.80),
    (ES: 'Rb'; EN: 'Rubidium'; V: (0, 1, 0, 0, 0); P: 1; AW: 85.4678),
    (ES: 'Sr'; EN: 'Strontium'; V: (0, 2, 0, 0, 0); P: 2; AW: 87.62),
    (ES: 'Y'; EN: 'Yttrium'; V: (0, 3, 0, 0, 0); P: 3; AW: 88.9059),
    (ES: 'Zr'; EN: 'Zirconium'; V: (0, 4, 0, 0, 0); P: 4; AW: 91.224),
    (ES: 'Nb'; EN: 'Niobium'; V: (0, 3, 5, 0, 0); P: 5; AW: 92.9064),
    (ES: 'Mo'; EN: 'Molybdenum'; V: (0, 3, 4, 5, 6); P: 6; AW: 95.94),
    (ES: 'Tc'; EN: 'Technetium'; V: (0, 7, 0, 0, 0); P: 7; AW: 98.9062),
    (ES: 'Ru'; EN: 'Ruthenium'; V: (0, 2, 3, 4, 6); P: 8; AW: 101.07),
    (ES: 'Rh'; EN: 'Rhodium'; V: (0, 2, 3, 4, 0); P: 8; AW: 102.9055),
    (ES: 'Pd'; EN: 'Palladium'; V: (0, 2, 4, 0, 0); P: 8; AW: 106.42),
    (ES: 'Ag'; EN: 'Silver'; V: (0, 1, 0, 0, 0); P: 1; AW: 107.868),
    (ES: 'Cd'; EN: 'Cadmium'; V: (0, 2, 0, 0, 0); P: 2; AW: 112.41),
    (ES: 'In'; EN: 'Indium'; V: (0, 3, 0, 0, 0); P: 3; AW: 114.82),
    (ES: 'Sn'; EN: 'Tin'; V: (0, 2, 4, 0, 0); P: 4; AW: 118.710),
    (ES: 'Sb'; EN: 'Antimony'; V: (0, 3, 5, 0, 0); P: 5; AW: 121.75),
    (ES: 'Te'; EN: 'Tellurium'; V: (0, 2, 4, 6, 0); P: 6; AW: 127.60),
    (ES: 'I'; EN: 'Iodine'; V: (0, 1, 3, 5, 7); P: 7; AW: 126.9045),
    (ES: 'Xe'; EN: 'Xenon'; V: (0, 0, 0, 0, 0); P: 0; AW: 131.29),
    (ES: 'Cs'; EN: 'Cesium'; V: (0, 1, 0, 0, 0); P: 1; AW: 132.905),
    (ES: 'Ba'; EN: 'Barium'; V: (0, 2, 0, 0, 0); P: 2; AW: 137.33),
    (ES: 'La'; EN: 'Lanthanum'; V: (0, 3, 0, 0, 0); P: 3; AW: 138.9055),
    (ES: 'Ce'; EN: 'Cerium'; V: (0, 3, 4, 0, 0); P: 3; AW: 140.12),
    (ES: 'Pr'; EN: 'Praseodymium'; V: (0, 3, 4, 0, 0); P: 3; AW: 140.9077),
    (ES: 'Nd'; EN: 'Neodymium'; V: (0, 3, 0, 0, 0); P: 3; AW: 144.24),
    (ES: 'Pm'; EN: 'Promethium'; V: (0, 3, 0, 0, 0); P: 3; AW: 147.0),
    (ES: 'Sm'; EN: 'Samarium'; V: (0, 2, 3, 0, 0); P: 3; AW: 150.36),
    (ES: 'Eu'; EN: 'Europium'; V: (0, 2, 3, 0, 0); P: 3; AW: 151.96),
    (ES: 'Gd'; EN: 'Gadolinium'; V: (0, 3, 0, 0, 0); P: 3; AW: 157.25),
    (ES: 'Tb'; EN: 'Terbium'; V: (0, 3, 4, 0, 0); P: 3; AW: 158.9254),
    (ES: 'Dy'; EN: 'Dysprosium'; V: (0, 3, 0, 0, 0); P: 3; AW: 162.50),
    (ES: 'Ho'; EN: 'Holmium'; V: (0, 3, 0, 0, 0); P: 3; AW: 164.9304),
    (ES: 'Er'; EN: 'Erbium'; V: (0, 3, 0, 0, 0); P: 3; AW: 167.26),
    (ES: 'Tm'; EN: 'Thulium'; V: (0, 2, 3, 0, 0); P: 3; AW: 168.9342),
    (ES: 'Yb'; EN: 'Ytterbium'; V: (0, 2, 3, 0, 0); P: 3; AW: 173.04),
    (ES: 'Lu'; EN: 'Lutetium'; V: (0, 3, 0, 0, 0); P: 3; AW: 174.967),
    (ES: 'Hf'; EN: 'Hafnium'; V: (0, 4, 0, 0, 0); P: 4; AW: 178.49),
    (ES: 'Ta'; EN: 'Tantalum'; V: (0, 5, 0, 0, 0); P: 5; AW: 180.9479),
    (ES: 'W'; EN: 'Tungsten'; V: (0, 3, 4, 5, 6); P: 6; AW: 183.85),
    (ES: 'Re'; EN: 'Rhenium'; V: (0, 2, 4, 6, 7); P: 7; AW: 186.207),
    (ES: 'Os'; EN: 'Osmium'; V: (0, 2, 3, 4, 6); P: 8; AW: 190.2),
    (ES: 'Ir'; EN: 'Iridium'; V: (0, 2, 3, 4, 6); P: 8; AW: 192.22),
    (ES: 'Pt'; EN: 'Platinum'; V: (0, 2, 4, 0, 0); P: 8; AW: 195.08),
    (ES: 'Au'; EN: 'Gold'; V: (0, 1, 3, 0, 0); P: 1; AW: 196.9665),
    (ES: 'Hg'; EN: 'Mercury'; V: (0, 1, 2, 0, 0); P: 2; AW: 200.59),
    (ES: 'Tl'; EN: 'Thallium'; V: (0, 1, 3, 0, 0); P: 3; AW: 204.383),
    (ES: 'Pb'; EN: 'Lead'; V: (0, 2, 4, 0, 0); P: 4; AW: 207.2),
    (ES: 'Bi'; EN: 'Bismuth'; V: (0, 3, 5, 0, 0); P: 5; AW: 208.9804),
    (ES: 'Po'; EN: 'Polonium'; V: (0, 2, 4, 0, 0); P: 6; AW: 209.0),
    (ES: 'At'; EN: 'Astatine'; V: (0, 1, 3, 5, 7); P: 7; AW: 210.0),
    (ES: 'Rn'; EN: 'Radon'; V: (0, 0, 0, 0, 0); P: 0; AW: 222.0),
    (ES: 'Fr'; EN: 'Francium'; V: (0, 1, 0, 0, 0); P: 1; AW: 223.0),
    (ES: 'Ra'; EN: 'Radium'; V: (0, 2, 0, 0, 0); P: 2; AW: 226.0254),
    (ES: 'Ac'; EN: 'Actinium'; V: (0, 3, 0, 0, 0); P: 3; AW: 227.0),
    (ES: 'Th'; EN: 'Thorium'; V: (0, 3, 4, 0, 0); P: 3; AW: 232.0381),
    (ES: 'Pa'; EN: 'Protactinium'; V: (0, 3, 4, 5, 0); P: 3; AW: 231.03559),
    (ES: 'U'; EN: 'Uranium'; V: (0, 3, 4, 5, 6); P: 3; AW: 238.0289),
    (ES: 'Np'; EN: 'Neptunium'; V: (0, 3, 4, 5, 6); P: 3; AW: 237.0482),
    (ES: 'Pu'; EN: 'Plutonium'; V: (0, 3, 4, 5, 6); P: 3; AW: 242.0),
    (ES: 'Am'; EN: 'Americium'; V: (0, 3, 4, 5, 6); P: 3; AW: 243.0),
    (ES: 'Cm'; EN: 'Curium'; V: (0, 3, 0, 0, 0); P: 3; AW: 247.0),
    (ES: 'Bk'; EN: 'Berkelium'; V: (0, 3, 4, 0, 0); P: 3; AW: 247.0),
    (ES: 'Cf'; EN: 'Californium'; V: (0, 3, 0, 0, 0); P: 3; AW: 249.0),
    (ES: 'Es'; EN: 'Einsteinium'; V: (0, 3, 0, 0, 0); P: 3; AW: 254.0),
    (ES: 'Fm'; EN: 'Fermium'; V: (0, 3, 0, 0, 0); P: 3; AW: 253.0),
    (ES: 'Md'; EN: 'Mendelevium'; V: (0, 3, 0, 0, 0); P: 3; AW: 256.0),
    (ES: 'No'; EN: 'Nobelium'; V: (0, 1, 0, 0, 0); P: 3; AW: 254.0),
    (ES: 'Lr'; EN: 'Lawrencium'; V: (0, 1, 0, 0, 0); P: 3; AW: 257.0));

implementation

function SDFChargeToInt(c: byte): integer;
begin
  case c of  //According to CTFiles standard
    1: Result := 3;
    2: Result := 2;
    3: Result := 1;
    4: Result := 0; //doublet radical
    5: Result := -1;
    6: Result := -2;
    7: Result := -3;
    else
      Result := 0;
  end;
end;

{ TMoleculeGraphic }

procedure TMoleculeGraphic.GetGeoCenter(var GX, GY, GZ: double);
var
  i: AtomID;
begin
  GX := 0;
  GY := 0;
  GZ := 0;
  for i := 1 to p_NX do
  begin
    GX := GX + GetX(i);
    GY := GY + GetY(i);
    GZ := GZ + GetZ(i);
  end;
  GX := GX / p_NX;
  GY := GY / p_NX;
  GZ := GZ / p_NX;
end;

constructor TMoleculeGraphic.Create;
begin
  inherited Create;
  Reset;
end;

destructor TMoleculeGraphic.Destroy;
begin
  inherited Destroy;
end;

procedure TMoleculeGraphic.Reset;
begin
  fDefaultFontSize := 10;
  fDefaultFontCol := colBlack;

  APrpSze := APrpSze + 3;//for X,Y,Z

  ABytSze := ABytSze + 1;//1 Byte for formal charge

  fAFontOffset := ABytSze - 1;
  ABytSze := ABytSze + 5;//5 Byte for Font color and size

  fAStyleOffSet := ABytSze - 1;
  ABytSze := ABytSze + 10;//5 Byte for Pen and 5 Byte for Brush

  BBytSze := BBytSze + 1;//for bond direction (up/down) and orientation

  fBStyleOffSet := BBytSze - 1;
  BBytSze := BBytSze + 10;//5 Byte for Pen and 5 Byte for Brush
end;

procedure TMoleculeGraphic.Clear;
begin
  inherited Clear;
  Reset;
end;

procedure TMoleculeGraphic.LoadSDF(sdfstr: TStringList);
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
  FntCol: TFPColor;
begin
  Clear;
  FntCol := fDefaultFontCol;
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
    if PAt^.S = 'C' then
      SetAtomFont(PAt, colBlack, fDefaultFontSize)
    else if PAt^.S = 'N' then
      SetAtomFont(PAt, colBlue, fDefaultFontSize)
    else if PAt^.S = 'O' then
      SetAtomFont(PAt, colRed, fDefaultFontSize)
    else if PAt^.S = 'S' then
      SetAtomFont(PAt, colDkYellow, fDefaultFontSize)
    else if PAt^.S = 'P' then
      SetAtomFont(PAt, colFuchsia, fDefaultFontSize)
    else if PAt^.S = 'F' then
      SetAtomFont(PAt, colGreen, fDefaultFontSize)
    else if PAt^.S = 'Cl' then
      SetAtomFont(PAt, colDkGreen, fDefaultFontSize)
    else if PAt^.S = 'Br' then
      SetAtomFont(PAt, colMaroon, fDefaultFontSize)
    else if PAt^.S = 'I' then
      SetAtomFont(PAt, colPurple, fDefaultFontSize)
    else if PAt^.S = 'H' then
      SetAtomFont(PAt, colTeal, fDefaultFontSize)
    else
      SetAtomFont(PAt, fDefaultFontCol, fDefaultFontSize);
    j := 0;
         {begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],35,2); //mass difference
            Inc(j);
         end;}
    begin
      PAt^.I[j] := int_readpos(sdfstr[LineNo], 37, 3); //charge;
      Inc(j);
    end;
         {begin
            PAt^.I[j]:=int_readpos(sdfstr[LineNo],40,3); //atom stereo parity
            Inc(j);
         end;}
    {begin
      PAt^.I[j] := int_readpos(sdfstr[LineNo], 43, 3); //hydrogen count + 1
      Inc(j);
    end;}
         {if fKeepAtmInfo[4] then begin
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
    SetLength(PBo^.P, BPrpSze); //set the bond properties array size (double)
    SetLength(PBo^.I, BBytSze); //set the bond properties array size (Byte)
    PBo^.t := int_readpos(sdfstr[LineNo], 1, 3);//Tail of bond
    PBo^.h := int_readpos(sdfstr[LineNo], 4, 3);//Head of bond
    PBo^.B := IntToTB(int_readpos(sdfstr[LineNo], 7, 3),
      int_readpos(sdfstr[LineNo], 19, 3));//Type of bond
    PBo^.S := BondSymbol[PBo^.B];
    j := 0;
    begin
      PBo^.I[j] := int_readpos(sdfstr[LineNo], 10, 3);//bond orientation
      Inc(j);
    end;
         {begin
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

procedure TMoleculeGraphic.SetAtomFont(AIndex: AtomID; AFontColor: TFPColor;
  AFontSize: byte);
var
  PAt: PRAtom;
begin
  PAt := AtmSet[AIndex];
  SetAtomFont(PAt, AFontColor, AFontSize);
end;

procedure TMoleculeGraphic.SetAtomFont(PAt: PRAtom; AFontColor: TFPColor;
  AFontSize: byte);
const
  WtoB = 257;
begin
  PAt^.I[fAFontOffset + 1] := AFontColor.red div WtoB;
  PAt^.I[fAFontOffSet + 2] := AFontColor.green div WtoB;
  PAt^.I[fAFontOffSet + 3] := AFontColor.blue div WtoB;
  PAt^.I[fAFontOffSet + 4] := AFontColor.alpha div WtoB;
  PAt^.I[fAFontOffset + 5] := AFontSize;
end;

procedure TMoleculeGraphic.SetAtomPenAndBrush(AIndex: AtomID;
  APenColor: TFPColor; APenStyle: TFPPenStyle; ABrushColor: TFPColor;
  ABrushStyle: TFPBrushStyle);
var
  PAt: PRAtom;
begin
  PAt := AtmSet[AIndex];
  PAt^.I[fAStyleOffSet + 1] := APenColor.red;
  PAt^.I[fAStyleOffSet + 2] := APenColor.green;
  PAt^.I[fAStyleOffSet + 3] := APenColor.blue;
  PAt^.I[fAStyleOffSet + 4] := APenColor.alpha;
  if APenStyle = psSolid then
    PAt^.I[fAStyleOffSet + 5] := 1
  else if APenStyle = psDash then
    PAt^.I[fAStyleOffSet + 5] := 2
  else if APenStyle = psDot then
    PAt^.I[fAStyleOffSet + 5] := 3
  else if APenStyle = psDashDot then
    PAt^.I[fAStyleOffSet + 5] := 4
  else if APenStyle = psDashDotDot then
    PAt^.I[fAStyleOffSet + 5] := 5
  else if APenStyle = psClear then
    PAt^.I[fAStyleOffSet + 5] := 6
  else if APenStyle = psInsideframe then
    PAt^.I[fAStyleOffSet + 5] := 7
  else if APenStyle = psPattern then
    PAt^.I[fAStyleOffSet + 5] := 8;

  PAt^.I[fAStyleOffSet + 6] := ABrushColor.red;
  PAt^.I[fAStyleOffSet + 7] := ABrushColor.green;
  PAt^.I[fAStyleOffSet + 8] := ABrushColor.blue;
  PAt^.I[fAStyleOffSet + 9] := ABrushColor.alpha;
  if ABrushStyle = bsSolid then
    PAt^.I[fAStyleOffSet + 10] := 1
  else if ABrushStyle = bsClear then
    PAt^.I[fAStyleOffSet + 10] := 2
  else if ABrushStyle = bsHorizontal then
    PAt^.I[fAStyleOffSet + 10] := 3
  else if ABrushStyle = bsVertical then
    PAt^.I[fAStyleOffSet + 10] := 4
  else if ABrushStyle = bsFDiagonal then
    PAt^.I[fAStyleOffSet + 10] := 5
  else if ABrushStyle = bsBDiagonal then
    PAt^.I[fAStyleOffSet + 10] := 6
  else if ABrushStyle = bsCross then
    PAt^.I[fAStyleOffSet + 10] := 7
  else if ABrushStyle = bsDiagCross then
    PAt^.I[fAStyleOffSet + 10] := 8;
end;

procedure TMoleculeGraphic.SetBondPenAndBrush(AIndex: BondID;
  APenColor: TFPColor; APenStyle: TFPPenStyle; ABrushColor: TFPColor;
  ABrushStyle: TFPBrushStyle);
var
  PBd: PRBond;
begin
  PBd := BndSet[AIndex];
  PBd^.I[fBStyleOffSet + 1] := APenColor.red;
  PBd^.I[fBStyleOffSet + 2] := APenColor.green;
  PBd^.I[fBStyleOffSet + 3] := APenColor.blue;
  PBd^.I[fBStyleOffSet + 4] := APenColor.alpha;
  if APenStyle = psSolid then
    PBd^.I[fBStyleOffSet + 5] := 1
  else if APenStyle = psDash then
    PBd^.I[fBStyleOffSet + 5] := 2
  else if APenStyle = psDot then
    PBd^.I[fBStyleOffSet + 5] := 3
  else if APenStyle = psDashDot then
    PBd^.I[fBStyleOffSet + 5] := 4
  else if APenStyle = psDashDotDot then
    PBd^.I[fBStyleOffSet + 5] := 5
  else if APenStyle = psClear then
    PBd^.I[fBStyleOffSet + 5] := 6
  else if APenStyle = psInsideframe then
    PBd^.I[fBStyleOffSet + 5] := 7
  else if APenStyle = psPattern then
    PBd^.I[fBStyleOffSet + 5] := 8;

  PBd^.I[fBStyleOffSet + 6] := ABrushColor.red;
  PBd^.I[fBStyleOffSet + 7] := ABrushColor.green;
  PBd^.I[fBStyleOffSet + 8] := ABrushColor.blue;
  PBd^.I[fBStyleOffSet + 9] := ABrushColor.alpha;
  if ABrushStyle = bsSolid then
    PBd^.I[fBStyleOffSet + 10] := 1
  else if ABrushStyle = bsClear then
    PBd^.I[fBStyleOffSet + 10] := 2
  else if ABrushStyle = bsHorizontal then
    PBd^.I[fBStyleOffSet + 10] := 3
  else if ABrushStyle = bsVertical then
    PBd^.I[fBStyleOffSet + 10] := 4
  else if ABrushStyle = bsFDiagonal then
    PBd^.I[fBStyleOffSet + 10] := 5
  else if ABrushStyle = bsBDiagonal then
    PBd^.I[fBStyleOffSet + 10] := 6
  else if ABrushStyle = bsCross then
    PBd^.I[fBStyleOffSet + 10] := 7
  else if ABrushStyle = bsDiagCross then
    PBd^.I[fBStyleOffSet + 10] := 8;
end;

procedure TMoleculeGraphic.GetAtomFont(AIndex: AtomID; out AFontColor: TFPColor;
  out AFontSize: byte);
var
  PAt: PRAtom;
begin
  PAt := AtmSet[AIndex];
  GetAtomFont(PAt, AFontColor, AFontSize);
end;

procedure TMoleculeGraphic.GetAtomFont(PAt: PRAtom; out AFontColor: TFPColor;
  out AFontSize: byte);
const
  WtoB = 257;
begin
  ;
  AFontColor.red := PAt^.I[fAFontOffset + 1] * WtoB;
  AFontColor.green := PAt^.I[fAFontOffSet + 2] * WtoB;
  AFontColor.blue := PAt^.I[fAFontOffSet + 3] * WtoB;
  AFontColor.alpha := PAt^.I[fAFontOffSet + 4] * WtoB;
  AFontSize := PAt^.I[fAFontOffset + 5];
end;

procedure TMoleculeGraphic.GetAtomPenAndBrush(AIndex: AtomID;
  APenColor: TFPColor; APenStyle: TFPPenStyle; ABrushColor: TFPColor;
  ABrushStyle: TFPBrushStyle);
var
  PAt: PRAtom;
begin
  PAt := AtmSet[AIndex];
  APenColor := FPColor(PAt^.I[fAStyleOffSet + 1], PAt^.I[fAStyleOffSet + 2],
    PAt^.I[fAStyleOffSet + 3], PAt^.I[fAStyleOffSet + 4]);
  case PAt^.I[fAStyleOffSet + 5] of
    1: APenStyle := psSolid;
    2: APenStyle := psDash;
    3: APenStyle := psDot;
    4: APenStyle := psDashDot;
    5: APenStyle := psDashDotDot;
    6: APenStyle := psClear;
    7: APenStyle := psInsideframe;
    8: APenStyle := psPattern;
    else
      APenStyle := psSolid;
  end;

  ABrushColor := FPColor(PAt^.I[fAStyleOffSet + 6], PAt^.I[fAStyleOffSet + 7],
    PAt^.I[fAStyleOffSet + 8], PAt^.I[fAStyleOffSet + 9]);
  case PAt^.I[fAStyleOffSet + 10] of
    1: ABrushStyle := bsSolid;
    2: ABrushStyle := bsClear;
    3: ABrushStyle := bsHorizontal;
    4: ABrushStyle := bsVertical;
    5: ABrushStyle := bsFDiagonal;
    6: ABrushStyle := bsBDiagonal;
    7: ABrushStyle := bsCross;
    8: ABrushStyle := bsDiagCross;
    else
      ABrushStyle := bsSolid;
  end;
end;

procedure TMoleculeGraphic.GetBondPenAndBrush(AIndex: BondID;
  APenColor: TFPColor; APenStyle: TFPPenStyle; ABrushColor: TFPColor;
  ABrushStyle: TFPBrushStyle);
var
  PBd: PRBond;
begin
  PBd := BndSet[AIndex];
  APenColor := FPColor(PBd^.I[fBStyleOffSet + 1], PBd^.I[fBStyleOffSet + 2],
    PBd^.I[fBStyleOffSet + 3], PBd^.I[fBStyleOffSet + 4]);
  case PBd^.I[fBStyleOffSet + 5] of
    1: APenStyle := psSolid;
    2: APenStyle := psDash;
    3: APenStyle := psDot;
    4: APenStyle := psDashDot;
    5: APenStyle := psDashDotDot;
    6: APenStyle := psClear;
    7: APenStyle := psInsideframe;
    8: APenStyle := psPattern;
    else
      APenStyle := psSolid;
  end;

  ABrushColor := FPColor(PBd^.I[fBStyleOffSet + 6], PBd^.I[fBStyleOffSet + 7],
    PBd^.I[fBStyleOffSet + 8], PBd^.I[fBStyleOffSet + 9]);
  case PBd^.I[fBStyleOffSet + 10] of
    1: ABrushStyle := bsSolid;
    2: ABrushStyle := bsClear;
    3: ABrushStyle := bsHorizontal;
    4: ABrushStyle := bsVertical;
    5: ABrushStyle := bsFDiagonal;
    6: ABrushStyle := bsBDiagonal;
    7: ABrushStyle := bsCross;
    8: ABrushStyle := bsDiagCross;
    else
      ABrushStyle := bsSolid;
  end;
end;

procedure TMoleculeGraphic.SetDefaultStyle();
var
  i: AtomID;
  k: BondID;
  AColor: TFPColor;
  APStyle: TPenStyle;
  ABStyle: TBrushStyle;
begin
  AColor := TColorToFPColor(clBlack);
  APStyle := psSolid;
  ABStyle := bsSolid;
  for i := 1 to p_NX do
  begin
    SetAtomPenAndBrush(i, AColor, APStyle, AColor, ABStyle);
  end;
  for i := 1 to p_M do
  begin
    SetBondPenAndBrush(i, AColor, APStyle, AColor, ABStyle);
  end;
end;

function TMoleculeGraphic.GetX(AIndex: AtomID): double;
begin
  Result := AtmSet[AIndex]^.P[0];
end;

function TMoleculeGraphic.GetY(AIndex: AtomID): double;
begin
  Result := AtmSet[AIndex]^.P[1];
end;

function TMoleculeGraphic.GetZ(AIndex: AtomID): double;
begin
  Result := AtmSet[AIndex]^.P[2];
end;

procedure TMoleculeGraphic.SetX(AIndex: AtomID; AValue: double);
begin
  AtmSet[AIndex]^.P[0] := AValue;
end;

procedure TMoleculeGraphic.SetY(AIndex: AtomID; AValue: double);
begin
  AtmSet[AIndex]^.P[1] := AValue;
end;

procedure TMoleculeGraphic.SetZ(AIndex: AtomID; AValue: double);
begin
  AtmSet[AIndex]^.P[2] := AValue;
end;

procedure TMoleculeGraphic.Center(CX, CY, CZ: double);
var
  GX, GY, GZ: double;
  i: AtomID;
  PAt: PRAtom;
begin
  GetGeoCenter(GX, GY, GZ);
  CX := CX - GX;
  CY := CY - GY;
  CZ := CZ - GZ;
  for i := 1 to p_NX do
  begin
    SetX(i, GetX(i) + CX);
    SetY(i, GetY(i) + CY);
    SetZ(i, GetZ(i) + CZ);
  end;
end;

procedure TMoleculeGraphic.Scale(SX, SY, SZ: double);
var
  GX, GY, GZ: double;
  i: AtomID;
  PAt: PRAtom;
begin
  GetGeoCenter(GX, GY, GZ);
  Center(0, 0, 0);
  for i := 1 to p_NX do
  begin
    SetX(i, GetX(i) * SX);
    SetY(i, GetY(i) * SY);
    SetZ(i, GetZ(i) * SZ);
  end;
  Center(GX,GY,GZ);
end;

function TMoleculeGraphic.DimX(): double;
var
  i, j: AtomID;
  PAt: PRAtom;
  d: double;
begin
  Result := 0;
  for i := 1 to p_NX do
    for j := i + 1 to p_NX do
    begin
      d := abs(GetX(i) - GetX(j));
      if Result < d then
        Result := d;
    end;
end;

function TMoleculeGraphic.DimY(): double;
var
  i, j: AtomID;
  PAt: PRAtom;
  d: double;
begin
  Result := 0;
  for i := 1 to p_NX do
    for j := i + 1 to p_NX do
    begin
      d := abs(GetY(i) - GetY(j));
      if Result < d then
        Result := d;
    end;
end;

function TMoleculeGraphic.DimZ(): double;
var
  i, j: AtomID;
  PAt: PRAtom;
  d: double;
begin
  Result := 0;
  for i := 1 to p_NX do
    for j := i + 1 to p_NX do
    begin
      d := abs(GetZ(i) - GetZ(j));
      if Result < d then
        Result := d;
    end;
end;

function TMoleculeGraphic.DimMaxXY(): double;
var
  dX, dY: double;
begin
  dX := DimX();
  dY := DimY();
  writeln('dX: '+FloatToStr(dX));
  writeln('dY: '+FloatToStr(dY));
  if dX > dY then
    Result := dX
  else
    Result := dY;
end;

function TMoleculeGraphic.DimMaxXYZ(): double;
var
  dX, dY, dZ: double;
begin
  dX := DimX();
  dY := DimY();
  dZ := DimZ();
  if (dX > dY) and (dX > dZ) then
    Result := dX
  else if (dY > dX) and (dY > dZ) then
    Result := dY
  else if (dZ > dX) and (dZ > dY) then
    Result := dZ;
end;

function TMoleculeGraphic.HCount(AIndex: AtomID): byte;
const
  MaxVal = 5;
var
  PAt: PRAtom;
  i, last: AtomID;
  j: BondID;
  Z: TZ;
  e: ELEMENT;
  succat: THead;
  nval: byte;
  B: TB;
  S: TS;
  ia, k, nh: byte;
  chg: integer;
begin
  PAt := AtmSet[AIndex];
  S := PAt^.S;
  Z := PAt^.Z;
  chg := 0;
  if PAt^.I[0] <> 0 then
    chg := SDFChargeToInt(PAt^.I[0]);
  e := elem[Z];
  for i:=1 to MaxVal do write(IntToStr(e.V[i])+' ');
  //If compiler stops here, check if you do not have multiple copies of U_GRAPH.pas in your unit paths.
  GetOutArc(AIndex, succat, last);
  nval := 0;
  nh := 0;
  ia := 0;
  for j := 1 to last do
  begin
    B := BndSet[succat[j]]^.B;
    case B of
      1: nval := nval + 1;
      2: nval := nval + 2;
      3: nval := nval + 3;
      4:
      begin
        nval := nval + 2;
        ia := ia + 1;
      end;
    end;
  end;
  if ia = 2 then
  begin
    nval := nval - 1;
    if (S = 'N') and (nval = 4) and (chg = 0) then
      nval := nval - 1;
    if ((S = 'S') or (S = 'O')) and (nval = 3) and (chg = 0) then
      nval := nval - 1;
  end;
  nval := nval + abs(chg);
  for k := MaxVale downto 1 do
    if e.V[k] >= nval then
      nh := e.V[k] - nval;
  Result := nh;
end;


end.

