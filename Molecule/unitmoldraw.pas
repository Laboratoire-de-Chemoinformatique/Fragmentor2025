unit UnitMolDraw;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Types, UnitMoleculeGraphic, UnitMoleculeBase, unitAtomAndBondType,
  fpvectorial, svgvectorialwriter, Graphics, math,
  {$IfNDef CGImode}
   fpvutils, fpvtocanvas, FPCanvas, FPimage, LCLProc, ExtCtrls;
  {$EndIf}
  {$IfDef CGImode}
   fpvutils, FPimage;
  {$EndIf}

type

  RCircles=record
    Cin,Cout: TvCircle;
  end;

  { TMol2DVec }

  TMol2DVec = class(TObject)
  private
    fMolG: TMoleculeGraphic;
    fVMolDoc: TvVectorialDocument;
    fVMolPage: TvVectorialPage;
    fVMolFmt: TvVectorialFormat;
    fVMolWidth, fVMolHeight: double;
    fDoubleDist: double;
    fTripleDist: double;
    fAromaticDist: double;
    fCircleRadius: double;
    fEpsilonRadius: double;
    fFontSize: integer;
    fnbin: integer;
    fRedMap: TFPPalette;
    fBlueMap: TFPPalette;
    fAtmCircles: array of RCircles;
    //    ftmpCanva: TCanvas;
    procedure SetWidth(d: double);
    procedure SetHeight(d: double);
    procedure MolCenter();
    procedure MolScale();
    procedure DrawBond(PBd: PRBond);
    procedure DrawBondD(PBd: PRBond);
    procedure DrawBond(X1, Y1, X2, Y2: double);
    procedure DrawBondD(X1, Y1, X2, Y2:double; T: TB);
    procedure DrawDoubleBond(PBd: PRBond);
    procedure DrawDoubleBond(X1, Y1, X2, Y2: double);
    procedure DrawDoubleBondD(PBd: PRBond);
    procedure DrawDoubleBondD(X1, Y1, X2, Y2: double; T: TB);
    procedure DrawTripleBond(PBd: PRBond);
    procedure DrawTripleBond(X1, Y1, X2, Y2: double);
    procedure DrawTripleBondD(PBd: PRBond);
    procedure DrawTripleBondD(X1, Y1, X2, Y2: double; T: TB);
    procedure DrawAromaticBond(PBd: PRBond);
    procedure DrawAromaticBond(X1, Y1, X2, Y2: double);
    procedure DrawAromaticBondD(PBd: PRBond);
    procedure DrawAromaticBondD(X1, Y1, X2, Y2: double; T: TB);
    procedure DrawBondWedge(PBd:PRBond; AColor: TFPColor);
    procedure DrawBondWedge(X1, Y1, X2, Y2: double; AColor: TFPColor);
    function FloatToFPColor(d, dmin, dmax: double): TFPColor;
  public
    constructor Create;
    destructor Destroy; override;
    procedure Clear;
    procedure SetMolGSDF(SDF: TStringList);
    procedure DrawMol();
    procedure DrawMol(SDF: TStringList);
    procedure DrawMolBckg(W: TDoubleDynArray);
    procedure DrawMolBckg(W: TDoubleDynArray; n: integer);
    {$IfNDef CGImode}
    procedure WriteToImage(AImage: TImage);
    {$EndIf}
    procedure WriteToFile(fnme: string);
    function WriteToCGI(): TStrings;
    procedure SetPalette(n: integer);
    procedure SetPalette();
    property nbin: integer read fnbin write fnbin;
    property MolG: TMoleculeGraphic read fMolG;
    property VMolDoc: TvVectorialDocument read fVMolDoc write fVMolDoc;
    property VMolPage: TvVectorialPage read fVMolPage write fVMolPage;
    property VMolFmt: TvVectorialFormat read fVMolFmt write fVMolFmt;
    property VMolWidth: double read fVMolWidth write SetWidth;
    property VMolHeight: double read fVMolHeight write SetHeight;
    property DoubleDist: double read fDoubleDist write fDoubleDist;
    property TripleDist: double read fTripleDist write fTripleDist;
    property AromaticDist: double read fAromaticDist write fAromaticDist;
    property CircleRadius: double read fCircleRadius write fCircleRadius;
    property EpsilonRadius: double read fEpsilonRadius write fEpsilonRadius;
    property FontSize: integer read fFontSize write fFontSize;
  end;

implementation



{ TMol2DVec }

procedure TMol2DVec.SetWidth(d: double);
begin
  fVMolWidth := d;
  fVMolPage.Width := fVMolWidth;
end;

procedure TMol2DVec.SetHeight(d: double);
begin
  fVMolHeight := d;
  fVMolPage.Height := fVMolHeight;
end;

procedure TMol2DVec.MolCenter();
var
   i: integer;
begin
  MolScale();
  fMolG.Center(fVMolWidth / 2.0, fVMolHeight / 2.0, 0);
{  fVMolPage.AddText(fVMolWidth / 2.0, fVMolHeight / 2.0,'X');
  for i:=1 to fMolG.nAtom do
      writeln('x: '+FloatToStr(fMolG.GetX(i))+' Y: '+FloatToStr(fMolG.GetX(i)));
  writeln('---------');}
end;

procedure TMol2DVec.MolScale();
var
  dX,dY,dMax, sf, MinDimOut: double;
begin
  dX:=fMolG.DimX();
  dY:=fMolG.DimY();
  if dX>dY then
     sf:=0.7 * fVMolWidth/dX
  else
     sf:=0.7 * fVMolHeight/dY;
  fMolG.Scale(sf, sf, sf);
end;

procedure TMol2DVec.DrawBond(PBd: PRBond);
var
  X1, Y1, X2, Y2: double;
begin
  X1 := fMolG.GetX(PBd^.t);
  Y1 := fMolG.GetY(PBd^.t);
  X2 := fMolG.GetX(PBd^.h);
  Y2 := fMolG.GetY(PBd^.h);
  DrawBond(X1, Y1, X2, Y2);
end;

procedure TMol2DVec.DrawBondD(PBd: PRBond);
var
  X1, Y1, X2, Y2: double;
begin
  X1 := fMolG.GetX(PBd^.t);
  Y1 := fMolG.GetY(PBd^.t);
  X2 := fMolG.GetX(PBd^.h);
  Y2 := fMolG.GetY(PBd^.h);
  DrawBondD(X1, Y1, X2, Y2, PBd^.B);
end;

procedure TMol2DVec.DrawBond(X1, Y1, X2, Y2: double);
begin
  with fVMolPage do
  begin
    StartPath(X1, Y1);
    AddLineToPath(X2, Y2);
    EndPath();
  end;
end;

procedure TMol2DVec.DrawBondD(X1, Y1, X2, Y2: double; T: TB);
var
  dX,dY,dX1,dX2,dY1,dY2: double;
begin
  DrawBond(X1,Y1,X2,Y2);
  dX:=abs(X2-X1)*0.25;
  dY:=abs(Y2-Y1)*0.25;
  if (X1<X2) then
  begin
    dX1:=X1+dX;
    dX2:=X2-dX;
  end else begin
    dX1:=X1-dX;
    dX2:=X2+dX;
  end;
  if (Y1<Y2) then
  begin
    dY1:=Y1+dY;
    dY2:=Y2-dY;
  end else begin
    dY1:=Y1-dY;
    dY2:=Y2+dY;
  end;
  with fVMolPage do
  begin
    StartPath(dX1, dY1);
    case T of
      14:SetPenColor(colGreen);
      18:SetPenColor(colRed);
      end;
    SetPenWidth(3);
    AddLineToPath(dX2, dY2);
    EndPath();
  end;
end;

procedure TMol2DVec.DrawDoubleBond(PBd: PRBond);
var
  X1, Y1, X2, Y2: double;
begin
  X1 := fMolG.GetX(PBd^.t);
  Y1 := fMolG.GetY(PBd^.t);
  X2 := fMolG.GetX(PBd^.h);
  Y2 := fMolG.GetY(PBd^.h);
  DrawDoubleBond(X1, Y1, X2, Y2);
end;

procedure TMol2DVec.DrawDoubleBond(X1, Y1, X2, Y2: double);
var
  ux, uy, vx, vy, nv: double;
  X1_P, Y1_P: double;
  X2_P, Y2_P: double;
begin
  ux := X2 - X1;
  uy := Y2 - Y1;
  vx:=0;
  vy:=0;
  if abs(ux)<abs(uy) then begin
    vx:=1;
    vy:=-ux/uy;
  end else begin
    vy:=1;
    vx:=-uy/ux;
  end;
  nv := sqrt(vx * vx + vy * vy);
  vx := vx / nv;
  vy := vy / nv;
  X1_P := X1 + DoubleDist * vx;
  Y1_P := Y1 + DoubleDist * vy;
  X2_P := X2 + DoubleDist * vx;
  Y2_P := Y2 + DoubleDist * vy;
  DrawBond(X1_P, Y1_P, X2_P, Y2_P);
  X1_P := X1 - DoubleDist * vx;
  Y1_P := Y1 - DoubleDist * vy;
  X2_P := X2 - DoubleDist * vx;
  Y2_P := Y2 - DoubleDist * vy;
  DrawBond(X1_P, Y1_P, X2_P, Y2_P);
end;

procedure TMol2DVec.DrawDoubleBondD(PBd: PRBond);
var
  X1, Y1, X2, Y2: double;
begin
  X1 := fMolG.GetX(PBd^.t);
  Y1 := fMolG.GetY(PBd^.t);
  X2 := fMolG.GetX(PBd^.h);
  Y2 := fMolG.GetY(PBd^.h);
  DrawDoubleBondD(X1, Y1, X2, Y2,PBd^.B);
end;

procedure TMol2DVec.DrawDoubleBondD(X1, Y1, X2, Y2: double; T: TB);
var
  ux, uy, vx, vy, nv: double;
  X1_P, Y1_P: double;
  X2_P, Y2_P: double;
begin
  ux := X2 - X1;
  uy := Y2 - Y1;
  vx:=0;
  vy:=0;
  if abs(ux)<abs(uy) then begin
    vx:=1;
    vy:=-ux/uy;
  end else begin
    vy:=1;
    vx:=-uy/ux;
  end;
  nv := sqrt(vx * vx + vy * vy);
  vx := vx / nv;
  vy := vy / nv;
  X1_P := X1 + DoubleDist * vx;
  Y1_P := Y1 + DoubleDist * vy;
  X2_P := X2 + DoubleDist * vx;
  Y2_P := Y2 + DoubleDist * vy;
  case T of
    25,19: DrawBondD(X1_P, Y1_P, X2_P, Y2_P, 18);
    22,15: DrawBondD(X1_P, Y1_P, X2_P, Y2_P, 14);
    else DrawBond(X1_P, Y1_P, X2_P, Y2_P);
  end;
  X1_P := X1 - DoubleDist * vx;
  Y1_P := Y1 - DoubleDist * vy;
  X2_P := X2 - DoubleDist * vx;
  Y2_P := Y2 - DoubleDist * vy;
  case T of
    19: DrawBondD(X1_P, Y1_P, X2_P, Y2_P, 18);
    15: DrawBondD(X1_P, Y1_P, X2_P, Y2_P, 14);
    else DrawBond(X1_P, Y1_P, X2_P, Y2_P);
  end;
end;

procedure TMol2DVec.DrawTripleBond(PBd: PRBond);
var
  X1, Y1, X2, Y2: double;
begin
  X1 := fMolG.GetX(PBd^.t);
  Y1 := fMolG.GetY(PBd^.t);
  X2 := fMolG.GetX(PBd^.h);
  Y2 := fMolG.GetY(PBd^.h);
  DrawTripleBond(X1, Y1, X2, Y2);
end;

procedure TMol2DVec.DrawTripleBond(X1, Y1, X2, Y2: double);
var
  ux, uy, vx, vy, nv: double;
  X1_P, Y1_P: double;
  X2_P, Y2_P: double;
begin
  DrawBond(X1, Y1, X2, Y2);
  ux := X2 - X1;
  uy := Y2 - Y1;
  vx:=0;
  vy:=0;
  if abs(ux)<abs(uy) then begin
    vx:=1;
    vy:=-ux/uy;
  end else begin
    vy:=1;
    vx:=-uy/ux;
  end;
  nv := sqrt(vx * vx + vy * vy);
  vx := vx / nv;
  vy := vy / nv;
  X1_P := X1 + TripleDist * vx;
  Y1_P := Y1 + TripleDist * vy;
  X2_P := X2 + TripleDist * vx;
  Y2_P := Y2 + TripleDist * vy;
  DrawBond(X1_P, Y1_P, X2_P, Y2_P);
  X1_P := X1 - TripleDist * vx;
  Y1_P := Y1 - TripleDist * vy;
  X2_P := X2 - TripleDist * vx;
  Y2_P := Y2 - TripleDist * vy;
  DrawBond(X1_P, Y1_P, X2_P, Y2_P);
end;

procedure TMol2DVec.DrawTripleBondD(PBd: PRBond);
var
  X1, Y1, X2, Y2: double;
begin
  X1 := fMolG.GetX(PBd^.t);
  Y1 := fMolG.GetY(PBd^.t);
  X2 := fMolG.GetX(PBd^.h);
  Y2 := fMolG.GetY(PBd^.h);
  DrawTripleBondD(X1, Y1, X2, Y2, PBd^.B);
end;

procedure TMol2DVec.DrawTripleBondD(X1, Y1, X2, Y2: double; T: TB);
var
  ux, uy, vx, vy, nv: double;
  X1_P, Y1_P: double;
  X2_P, Y2_P: double;
begin
  case T of
    20: DrawBondD(X1_P, Y1_P, X2_P, Y2_P,18);
    16: DrawBondD(X1_P, Y1_P, X2_P, Y2_P,14);
    else DrawBond(X1_P, Y1_P, X2_P, Y2_P);
  end;
  ux := X2 - X1;
  uy := Y2 - Y1;
  vx:=0;
  vy:=0;
  if abs(ux)<abs(uy) then begin
    vx:=1;
    vy:=-ux/uy;
  end else begin
    vy:=1;
    vx:=-uy/ux;
  end;
  nv := sqrt(vx * vx + vy * vy);
  vx := vx / nv;
  vy := vy / nv;
  X1_P := X1 + TripleDist * vx;
  Y1_P := Y1 + TripleDist * vy;
  X2_P := X2 + TripleDist * vx;
  Y2_P := Y2 + TripleDist * vy;
  case T of
    20,28,29: DrawBondD(X1_P, Y1_P, X2_P, Y2_P,18);
    26,23,16: DrawBondD(X1_P, Y1_P, X2_P, Y2_P,14);
    else DrawBond(X1_P, Y1_P, X2_P, Y2_P);
  end;
  X1_P := X1 - TripleDist * vx;
  Y1_P := Y1 - TripleDist * vy;
  X2_P := X2 - TripleDist * vx;
  Y2_P := Y2 - TripleDist * vy;
  case T of
    20,28: DrawBondD(X1_P, Y1_P, X2_P, Y2_P,18);
    23,16: DrawBondD(X1_P, Y1_P, X2_P, Y2_P,14);
    else DrawBond(X1_P, Y1_P, X2_P, Y2_P);
  end;
end;

procedure TMol2DVec.DrawAromaticBond(PBd: PRBond);
var
  At1, At2, At3: AtomID;
  M, s: integer;
  iBd: BondID;
  iPBd: PRBond;
  X1, Y1, X2, Y2: double;
  ux, uy, vx, vy, nv, scal1, scal2: double;
  X1_P, Y1_P: double;
  X2_P, Y2_P: double;
  bSet: boolean;
begin
  At1 := PBd^.t;
  At2 := PBd^.h;
  //orient vector v interior to aromaticity
  with fMolG do
  begin
    X1 := GetX(At1);
    Y1 := GetY(At1);
    X2 := GetX(At2);
    Y2 := GetY(At2);
    ux := X2 - X1;
    uy := Y2 - Y1;
    vx:=0;
    vy:=0;
    if abs(ux)<abs(uy) then begin
      vx:=1;
      vy:=-ux/uy;
    end else begin
      vy:=1;
      vx:=-uy/ux;
    end;
    nv := sqrt(vx * vx + vy * vy);
    vx := vx / nv;
    vy := vy / nv;
    M := p_HEAD[At1];
    bset := False;
    while (M < p_HEAD[At1 + 1]) and (not bset) do
    begin
      At3 := p_SUCC[M];
      iPBd := FindBond(At1, At3);
      if (iPBd <> PBd) and (iPBd <> nil) then
      begin
        if (iPBd^.B = 4) then
        begin
          scal1 := (GetX(At3) - X1) * (X2 - X1) + (GetY(At3) - Y1) * (Y2 - Y1);
          scal2 := (GetX(At3) - X1) * vx + (GetY(At3) - Y1) * vy;
          if scal1 * scal2 > 0 then
          begin
            vx := -vx;
            vy := -vy;
          end;
          bSet := True;
        end;
      end;
      Inc(M);
    end;
  end;
  X1_P := X1 + AromaticDist * vx;
  Y1_P := Y1 + AromaticDist * vy;
  X2_P := X2 + AromaticDist * vx;
  Y2_P := Y2 + AromaticDist * vy;
  DrawAromaticBond(X1_P, Y1_P, X2_P, Y2_P);
  DrawBond(X1, Y1, X2, Y2);
end;

procedure TMol2DVec.DrawAromaticBond(X1, Y1, X2, Y2: double);
begin
  with fVMolPage do
  begin
    StartPath(X1, Y1);
    SetPenStyle(psDot);
    AddLineToPath(X2, Y2);
    EndPath();
  end;
end;

procedure TMol2DVec.DrawAromaticBondD(PBd: PRBond);
var
  At1, At2, At3: AtomID;
  M, s: integer;
  iBd: BondID;
  iPBd: PRBond;
  X1, Y1, X2, Y2: double;
  ux, uy, vx, vy, nv, scal1, scal2: double;
  X1_P, Y1_P: double;
  X2_P, Y2_P: double;
  bSet: boolean;
  T: TB;
begin
  At1 := PBd^.t;
  At2 := PBd^.h;
  T:=PBd^.B;
  //orient vector v interior to aromaticity
  with fMolG do
  begin
    X1 := GetX(At1);
    Y1 := GetY(At1);
    X2 := GetX(At2);
    Y2 := GetY(At2);
    ux := X2 - X1;
    uy := Y2 - Y1;
    vx:=0;
    vy:=0;
    if abs(ux)<abs(uy) then begin
      vx:=1;
      vy:=-ux/uy;
    end else begin
      vy:=1;
      vx:=-uy/ux;
    end;
    nv := sqrt(vx * vx + vy * vy);
    vx := vx / nv;
    vy := vy / nv;
    M := p_HEAD[At1];
    bset := False;
    while (M < p_HEAD[At1 + 1]) and (not bset) do
    begin
      At3 := p_SUCC[M];
      iPBd := FindBond(At1, At3);
      if (iPBd <> PBd) and (iPBd <> nil) then
      begin
        if (iPBd^.B = 4) then
        begin
          scal1 := (GetX(At3) - X1) * (X2 - X1) + (GetY(At3) - Y1) * (Y2 - Y1);
          scal2 := (GetX(At3) - X1) * vx + (GetY(At3) - Y1) * vy;
          if scal1 * scal2 > 0 then
          begin
            vx := -vx;
            vy := -vy;
          end;
          bSet := True;
        end;
      end;
      Inc(M);
    end;
  end;
  X1_P := X1 + AromaticDist * vx;
  Y1_P := Y1 + AromaticDist * vy;
  X2_P := X2 + AromaticDist * vx;
  Y2_P := Y2 + AromaticDist * vy;
  case T of
    17,24,30: DrawAromaticBondD(X1_P, Y1_P, X2_P, Y2_P,14);
    27: begin
          DrawBondD(X1_P, Y1_P, X2_P, Y2_P,18);
          DrawAromaticBondD(X1_P, Y1_P, X2_P, Y2_P,14);
        end;
    21,31,33: DrawAromaticBondD(X1_P, Y1_P, X2_P, Y2_P,18);
    32: begin
          DrawBondD(X1_P, Y1_P, X2_P, Y2_P,14);
          DrawAromaticBondD(X1_P, Y1_P, X2_P, Y2_P,18);
        end;
    else DrawAromaticBond(X1_P, Y1_P, X2_P, Y2_P);
  end;
  case T of
    17: DrawBondD(X1, Y1, X2, Y2,14);
    21: DrawBondD(X1, Y1, X2, Y2,18);
    else DrawBond(X1, Y1, X2, Y2);
  end;
  X1_P := X1 - AromaticDist * vx;
  Y1_P := Y1 - AromaticDist * vy;
  X2_P := X2 - AromaticDist * vx;
  Y2_P := Y2 - AromaticDist * vy;
  case T of
    30: DrawBondD(X1_P, Y1_P, X2_P, Y2_P,18);
    33: DrawBondD(X1_P, Y1_P, X2_P, Y2_P,14);
  end;
end;

procedure TMol2DVec.DrawAromaticBondD(X1, Y1, X2, Y2: double; T: TB);
begin
  with fVMolPage do
  begin
    StartPath(X1, Y1);
    SetPenStyle(psDot);
    case T of
      14: SetPenColor(colGreen);
      18: SetPenColor(colRed);
    end;
    AddLineToPath(X2, Y2);
    EndPath();
  end;
end;

procedure TMol2DVec.DrawBondWedge(PBd: PRBond; AColor: TFPColor);
var
  X1, Y1, X2, Y2: double;
begin
  X1 := fMolG.GetX(PBd^.t);
  Y1 := fMolG.GetY(PBd^.t);
  X2 := fMolG.GetX(PBd^.h);
  Y2 := fMolG.GetY(PBd^.h);
  DrawBondWedge(X1, Y1, X2, Y2,AColor);
end;

procedure TMol2DVec.DrawBondWedge(X1, Y1, X2, Y2: double; AColor: TFPColor);
var
  ux, uy, vx, vy, nv: double;
  X1_P, Y1_P: double;
  X2_P, Y2_P: double;
begin
  ux := X2 - X1;
  uy := Y2 - Y1;
  vx:=0;
  vy:=0;
  if abs(ux)<abs(uy) then begin
    vx:=1;
    vy:=-ux/uy;
  end else begin
    vy:=1;
    vx:=-uy/ux;
  end;
  nv := sqrt(vx * vx + vy * vy);
  vx := vx / nv;
  vy := vy / nv;
  X1_P := X2 + DoubleDist * vx;
  Y1_P := Y2 + DoubleDist * vy;
  X2_P := X2 - DoubleDist * vx;
  Y2_P := Y2 - DoubleDist * vy;
  with fVMolPage do
  begin
    StartPath(X1, Y1);
    SetPenStyle(psSolid);
    SetPenColor(colBlack);
    SetBrushColor(AColor);
    SetBrushStyle(bsSolid);
    AddLineToPath(X1_P, Y1_P);
    AddLineToPath(X2_P, Y2_P);
    AddLineToPath(X1, Y1);
    EndPath();
  end;
end;

function TMol2DVec.FloatToFPColor(d, dmin, dmax: double): TFPColor;
var
  i,j: integer;
  upb: double;
begin
  upb:=0;
  if d=0 then
  begin
     j:=0;
     Result.alpha:=High(Word);
     Result.red:=High(Word);
     Result.green:=High(Word);
     Result.blue:=High(Word);
  end
  else if d>0 then
  begin
    j:=0;
    repeat
      Inc(j);
      upb:=j*dmax/nbin;
    until (d<=upb);
    Result:=fBlueMap.Color[j-1];
  end
  else if d<0 then
  begin
    j:=0;
    repeat
      Inc(j);
      upb:=j*dmin/nbin;
    until (d>=upb);
    Result:=fRedMap.Color[j-1];
  end;
end;

constructor TMol2DVec.Create;
begin
  fMolG := TMoleculeGraphic.Create;
  fVMolFmt := vfSVG;
  fVMolDoc := TvVectorialDocument.Create;
  fVMolPage := fVMolDoc.AddPage();
  fVMolWidth := 100;
  fVMolHeight := 100;
  //fVMolDoc.Width:=fVMolWidth;
  //fVMolDoc.Height:=fVMolHeight;
  //fVMolDoc.GuessDocumentSize();
  fVMolPage.Width := fVMolWidth;
  fVMolPage.Height := fVMolHeight;
  fDoubleDist := 3;
  fTripleDist := 2;
  fAromaticDist := 3;
  fCircleRadius := 7;
  fEpsilonRadius:=Ceil(0.4*fCircleRadius);
  fFontSize := 10;
  fRedMap:=TFPPalette.Create(0);
  fBlueMap:=TFPPalette.Create(0);
  fnbin:=0;
  //ftmpCanva:=TCanvas.Create;
  //fVMolPage.Clear;
end;

destructor TMol2DVec.Destroy;
begin
  Clear;
  FreeAndNil(fVMolDoc);
  FreeAndNil(fMolG);
  FreeAndNil(fRedMap);
  FreeAndNil(fBlueMap);
  //FreeAndNil(ftmpCanva);
  inherited Destroy;
end;

procedure TMol2DVec.Clear;
begin
  fMolG.Clear;
  fVMolDoc.Clear;
end;

procedure TMol2DVec.SetMolGSDF(SDF: TStringList);
begin
  fMolG.LoadSDF(SDF);
  //fMolG.SetDefaultStyle();
end;

procedure TMol2DVec.DrawMol();
var
  i, j, k: integer;
  iAt1, iAt2: AtomID;
  PAt, PAt1, PAt2: PRAtom;
  PBd: PRBond;
  dMax, sf, MinDimOut: double;
  Circle: TvCircle;
  Text: TvText;
  FontWidth, FontHeight: double;
  FntCol: TFPColor;
  FntSze: byte;
  ChgStr: string;
  chg: integer;
  dX, dY: double;

begin
  MolCenter();
  //SetWidth(1000);
  //SetHeight(1000);
  //fVMolPage.Width := 10;
  //fVMolPage.Height := 10;
  for k := 1 to fMolG.p_M do
  begin
    PBd := fMolG.BndSet[k];
    iAt1 := PBd^.h;
    iAt2 := PBd^.t;
    case PBd^.B of
      1: DrawBond(PBd);
      2: DrawDoubleBond(PBd);
      3: DrawTripleBond(PBd);
      4: DrawAromaticBond(PBd);
      14: DrawBondD(PBd);//single bond create
      18: DrawBondD(PBd);//single bond break
      25: DrawDoubleBondD(PBd);//double bond to single bond
      22 : DrawDoubleBondD(PBd);//single bond to double bond
      else;
    end;
    if (PBd^.B = 1) and (PBd^.I[0] = 1) then
      DrawBondWedge(PBd,colBlack);
    if (PBd^.B = 1) and (PBd^.I[0] = 6) then
      DrawBondWedge(PBd,colWhite);
    if (PBd^.B = 1) and (PBd^.I[0] = 4) then
      DrawBondWedge(PBd,colGray);
  end;
  for i := 1 to fMolG.p_NX do
  begin
    PAt := fMolG.AtmSet[i];
    chg:=0;
    ChgStr := '';
    if PAt^.I[0]<>0 then
     chg:=SDFChargeToInt(PAt^.I[0]);
    if (chg = 1) and (PAt^.I[0]>0) then
       ChgStr:='+'
    else if (chg = 1) and (PAt^.I[0]<0) then
     ChgStr:='-'
    else if (chg > 1) and (PAt^.I[0]>0) then
      ChgStr := IntToStr(chg) + '+'
    else if (chg < 1) and (PAt^.I[0]>0) then
      ChgStr := IntToStr(-chg) + '-';
    if fAtmCircles<>nil then
    begin
      fAtmCircles[i].Cin.X:=fMolG.GetX(i);
      fAtmCircles[i].Cin.Y:=fMolG.GetY(i);
      fAtmCircles[i].Cout.X:=fMolG.GetX(i);
      fAtmCircles[i].Cout.Y:=fMolG.GetY(i);
    end;
    if (PAt^.S <> 'C') or (ChgStr <> '') then
    begin
      fMolG.GetAtomFont(i, FntCol, FntSze);
      //
      fVMolPage.AddCircle(fMolG.GetX(i), fMolG.GetY(i), fCircleRadius);
      //writeln('dX: '+FloatToStr(fMolG.DimX())+' dY: '+FloatToStr(fMolG.DimY()));
      Circle := TvCircle(fVMolPage.GetLastEntity());
      Circle.Brush.Color := TColorToFPColor(clWhite);
      Circle.Brush.Style := bsSolid;
      Circle.Pen.Style := psClear;
      //We have to move the atoms in case of Linux machines
      //fVMolPage.AddText(fMolG.GetX(i) - 0.5 * FntSze, fMolG.GetY(i) -
      //  0.5 * FntSze, 0, 'Sans', FntSze, fMolG.S_[i]); //0.5->2.0 for version befor lazarus 1.6 and lazarus 2.0
      {$IfDef UNIX}{$IfnDef darwin}
      fVMolPage.AddText(fMolG.GetX(i) - 0.5 * FntSze, fMolG.GetY(i) -
        0.5 * FntSze, 0, 'Sans', FntSze, fMolG.S_[i]); //0.5->2.0 for version befor lazarus 1.6 and lazarus 2.0 and later
      {$EndIf}{$Endif}
      {$IfDef WINDOWS}
      fVMolPage.AddText(fMolG.GetX(i) - 0.5 * FntSze, fMolG.GetY(i) -
        0.8 * FntSze, 0, 'Sans', FntSze, fMolG.S_[i]); //0.8 for windows 10 (Lazarus 2.0.4)
      {$EndIf}
      {$IfDef darwin}
      fVMolPage.AddText(fMolG.GetX(i) - 0.5 * FntSze, fMolG.GetY(i) -
        2.0 * FntSze, 0, 'Sans', FntSze, fMolG.S_[i]); //0.5->2.0 for version befor lazarus 1.6
      {$ENDIF}
      Text := TvText(fVMolPage.GetLastEntity());
      Text.Font.Color := FntCol;//  FntCol;
      if ChgStr <> '' then
      begin
        //Do we have to move the atoms in case of Linux machines here ?
        fVMolPage.AddText(fMolG.GetX(i) + (0.5 + 0.5 * Length(fMolG.S_[i])) * FntSze -
          0.7 * FntSze, fMolG.GetY(i) -0.85 * FntSze, 0, 'Sans', (2 * FntSze div 3), ChgStr);
        Text := TvText(fVMolPage.GetLastEntity());
        Text.Font.Color := FntCol;//  FntCol;
      end;
    end;
  end;
  //Change the fVMolDoc Width and Height to fit the molecule
  MolCenter();
  //fVMolDoc.GuessDocumentSize(); Note: this function produces a SIGSEGV error on Ubuntu
  fVMolDoc.Width:=fVMolPage.Width;
  fVMolDoc.Height:=fVMolPage.Height;
  {WriteLn('fVMolPage.MaxX='+FloatToStr(fVMolPage.MaxX));
  WriteLn('fVMolPage.MaxY='+FloatToStr(fVMolPage.MaxY));
  WriteLn('fVMolPage.MaxZ='+FloatToStr(fVMolPage.MaxZ));

  WriteLn('fVMolPage.MinX='+FloatToStr(fVMolPage.MinX));
  WriteLn('fVMolPage.MinY='+FloatToStr(fVMolPage.MinY));
  WriteLn('fVMolPage.MinZ='+FloatToStr(fVMolPage.MinZ));

  WriteLn('fVMolPage.GetEntitiesCount='+IntToStr(fVMolPage.GetEntitiesCount));

  WriteLn('fVMolPage.GetEntity(1).X='+FloatToStr(fVMolPage.GetEntity(1).X));
  WriteLn('fVMolPage.GetEntity(1).Y='+FloatToStr(fVMolPage.GetEntity(1).Y));
  WriteLn('fVMolPage.GetEntity(1).Z='+FloatToStr(fVMolPage.GetEntity(1).Z));

  WriteLn('fVMolDoc.Width='+FloatToStr(fVMolDoc.Width)+' fVMolDoc.Height='+FloatToStr(fVMolDoc.Height));}


end;

procedure TMol2DVec.DrawMol(SDF: TStringList);
begin
  SetMolGSDF(SDF);
  DrawMol();
end;

procedure TMol2DVec.DrawMolBckg(W: TDoubleDynArray);
var
  i: integer;
  wmin, wmax: double;
  c: TFPColor;
  Circle: TvCircle;
begin
  wmin:=W[1];
  wmax:=W[1];
  for i:=1 to High(W) do
    if W[i]<wmin then
       wmin:=W[i]
    else if W[i]>wmax then
       wmax:=W[i];
  if wmin>0 then wmin:=0;
  if wmax<0 then wmax:=0;
  SetLength(fAtmCircles,fMolG.p_NX+1);
  for i := 1 to fMolG.p_NX do
  begin
    c:=FloatToFPColor(W[i],wmin,wmax);
    fVMolPage.AddCircle(0, 0, fCircleRadius+fEpsilonRadius);
    Circle := TvCircle(fVMolPage.GetLastEntity());
    Circle.Brush.Color := c;
    Circle.Brush.Style := bsSolid;
    Circle.Pen.Style := psSolid;
    Circle.Pen.Width:=1;
    fAtmCircles[i].Cout:=Circle;
    fVMolPage.AddCircle(0, 0, fCircleRadius);
    Circle := TvCircle(fVMolPage.GetLastEntity());
    Circle.Brush.Color := TColorToFPColor(clWhite);
    Circle.Brush.Style := bsSolid;
    Circle.Pen.Style := psSolid;
    Circle.Pen.Width:=1;
    fAtmCircles[i].Cin:=Circle;
  end;
end;

procedure TMol2DVec.DrawMolBckg(W: TDoubleDynArray; n: integer);
var
  i: integer;
begin
  SetPalette(n);
  DrawMolBckg(W);
end;

{$IfNDef CGImode}
procedure TMol2DVec.WriteToImage(AImage: TImage);
begin
  AImage.Canvas.Brush.Color := clWhite;
  AImage.Canvas.Brush.Style := bsSolid;
  AImage.Canvas.FillRect(0, 0, AImage.Width, AImage.Height);
  DrawFPVectorialToCanvas(fVMolPage, AImage.Canvas);//Creates a memory leak with fpvectorial
  //Bug already reported in November 2019. Target is Lazarus 2.0.8
  AImage.Invalidate;
end;
{$EndIf}

procedure TMol2DVec.WriteToFile(fnme: string);
begin
  fVMolDoc.WriteToFile(fnme, fVMolFmt);
end;

function TMol2DVec.WriteToCGI(): TStrings;
var
  AStrings: TStrings;
begin
  //First get the strings
  AStrings:=TStrings.Create;
  fVMolDoc.WriteToStrings(AStrings, fVMolFmt);
  //Then modify 3 things in the file:
  //1-Change the dimensions in <svg to w*0.25 and h*0.25
  //2-Add <g transform="scale(0.25)"> right after <svg>
  //3-Add </g> right before </svg>

  Result:=AStrings;
end;

procedure TMol2DVec.SetPalette(n: integer);
begin
  fnbin:=n;
  SetPalette();
end;

procedure TMol2DVec.SetPalette();
var
  i: integer;
  w: Word;
  c: TFPColor;
  ur,ub: array [1..3] of double;
  t: double;
begin
  ur[1]:=0;
  ur[2]:=0-High(Word);
  ur[3]:=0-High(Word);
  ub[1]:=0-High(Word);
  ub[2]:=0-High(Word);
  ub[3]:=0;
  c.alpha:=High(Word);
  for i:=1 to fnbin+1 do //Since a change in FPimage, need to add a bin...
  begin
    t:=i/fnbin;
    c.red:=Floor(ub[1]*t)+High(Word);
    c.green:=Floor(ub[2]*t)+High(Word);
    c.blue:=Floor(ub[3]*t)+High(Word);
    fBlueMap.Add(c);
    c.red:=Floor(ur[1]*t)+High(Word);
    c.green:=Floor(ur[2]*t)+High(Word);
    c.blue:=Floor(ur[3]*t)+High(Word);
    fRedMap.Add(c);
  end;
end;

end.
//procedure test;
{procedure TMol2DVec.test;
var
  VecDoc: TvVectorialDocument;
  Vec: TvVectorialPage;
begin
  VecDoc := TvVectorialDocument.Create;
  try
    Vec := VecDoc.AddPage();
    // All documents are 10cm x 10cm
    Vec.Width := 100;
    Vec.Height := 100;

    // ...

    // multi_test_1     Combines various elements
    //Vec.Clear;
    Vec.StartPath(0, 20);
    Vec.AddLineToPath(30, 30);
    Vec.EndPath();
    Vec.StartPath(0, 0);
    Vec.AddLineToPath(100, 0);
    Vec.AddLineToPath(100, 100);
    Vec.AddLineToPath(0, 100);
    Vec.AddLineToPath(0, 0);
    Vec.EndPath();
    Vec.StartPath(0, 0);
    Vec.AddLineToPath(10, 10);
    Vec.AddBezierToPath(10, 20, 20, 20, 20, 10);
    Vec.AddLineToPath(30, 0);
    Vec.EndPath();
    Vec.AddText(10, 10, 0, '10,10 Some text in english.');
    Vec.AddText(20, 20, 0, '20, 20 Mówić, cześć, Włosku, Parabéns.');
    Vec.AddText(30, 30, 0, '30, 30 森林，是一个高密');
    VecDoc.WriteToFile('multi_test_1.svg', vfSVG);
  finally
    FreeAndNil(VecDoc);
  end;
end;}

