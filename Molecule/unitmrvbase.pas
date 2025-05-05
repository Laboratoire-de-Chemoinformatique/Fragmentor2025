unit UnitMRVBase;

{Read and Write Chemaxon MRV format}

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, contnrs, DOM, XMLRead, XMLWrite, UnitMoleculeBase,
  unitAtomAndBondType;

type

  { TMoleculeMrv }

  TMoleculeMrv = class(TMoleculeBase)
  private
    Global: TFPStringHashTable;
    AtmID:  TStringList;
    procedure SetX(Index: integer; d: double);
    function GetX(Index: integer): double;
    procedure SetY(Index: integer; d: double);
    function GetY(Index: integer): double;
    procedure SetZ(Index: integer; d: double);
    function GetZ(Index: integer): double;
    procedure SetMrvAtmSet(Index: integer; i: byte);
    function GetMrvAtmSet(Index: integer): byte;
    procedure SetMrvBndSet(Index: integer; i: byte);
    function GetMrvBndSet(Index: integer): byte;
    procedure SetMrvStereo(Index: integer; i: byte);
    function GetMrvStereo(Index: integer): byte;
  public
    constructor Create;
    destructor Destroy; override;
    procedure ReadMRV(mrfle: string);
    procedure WriteMRV(mrfle: string);
    property X[Index: integer]: double Read GetX Write SetX;
    property Y[Index: integer]: double Read GetY Write SetY;
    property Z[Index: integer]: double Read GetZ Write SetZ;
    property MrvAtmSet[Index: integer]: byte Read GetMrvAtmSet Write SetMrvAtmSet;
    property MrvBndSet[Index: integer]: byte Read GetMrvBndSet Write SetMrvBndSet;
    property MrvStereo[Index: integer]: byte Read GetMrvStereo Write SetMrvStereo;
  end;

implementation

{ TMoleculeMrv }

procedure TMoleculeMrv.SetX(Index: integer; d: double);
begin
  AtmSet[Index]^.P[0] := d;
end;

function TMoleculeMrv.GetX(Index: integer): double;
begin
  Result := AtmSet[Index]^.P[0];
end;

procedure TMoleculeMrv.SetY(Index: integer; d: double);
begin
  AtmSet[Index]^.P[1] := d;
end;

function TMoleculeMrv.GetY(Index: integer): double;
begin
  Result := AtmSet[Index]^.P[1];
end;

procedure TMoleculeMrv.SetZ(Index: integer; d: double);
begin
  AtmSet[Index]^.P[2] := d;
end;

function TMoleculeMrv.GetZ(Index: integer): double;
begin
  Result := AtmSet[Index]^.P[2];
end;

procedure TMoleculeMrv.SetMrvAtmSet(Index: integer; i: byte);
begin
  AtmSet[Index]^.I[0] := i;
end;

function TMoleculeMrv.GetMrvAtmSet(Index: integer): byte;
begin
  Result := AtmSet[Index]^.I[0];
end;

procedure TMoleculeMrv.SetMrvBndSet(Index: integer; i: byte);
begin
  BndSet[Index]^.I[1] := i;
end;

function TMoleculeMrv.GetMrvBndSet(Index: integer): byte;
begin
  Result := BndSet[Index]^.I[1];
end;

procedure TMoleculeMrv.SetMrvStereo(Index: integer; i: byte);
begin
  BndSet[Index]^.I[0] := i;
end;

function TMoleculeMrv.GetMrvStereo(Index: integer): byte;
begin
  Result := BndSet[Index]^.I[0];
end;

constructor TMoleculeMrv.Create;
begin
  AtmID  := TStringList.Create;
  Global := TFPStringHashTable.Create;
  inherited Create;
end;

destructor TMoleculeMrv.Destroy;
begin
  FreeAndNil(AtmID);
  FreeAndNil(Global);
  inherited Destroy;
end;

procedure TMoleculeMrv.ReadMRV(mrfle: string);
var
  doc:    TXMLDocument;
  root, entry, MChemStruc, molecule, molitem, atom, bond, bondchld: TDOMNode;
  i, j, j1, j2, n, m, Mbnd, s, t: integer;
  PAt:    PRAtom;
  PBd:    PRBond;
  tmpSL:  TStringList;
  atmp:   array [1..MaxAtom] of array [1..MaxAtom] of integer;
  wtmp:   array [1..MaxBond] of PRBond;
  bfound: boolean;

  procedure InitLocalArrays(dima: integer);
  var
    i, j, dimb: integer;
  begin
    //Reset
    dimb := 4 * dima;
    //Init
    for i := 1 to dima do
      for j := 1 to dima do
        atmp[i, j] := 0;
    for i := 1 to dimb do
      wtmp[i] := nil;
  end;

begin
  tmpSL := TStringList.Create;
  tmpSL.Delimiter := ' ';
  ReadXMLFile(doc, mrfle);
  root  := doc.FindNode('cml');
  entry := root.FirstChild;

  APrpSze := 3;
  ABytSze := 1;
  BPrpSze := 0;
  BBytSze := 2;

  while entry <> nil do
  begin
    with entry.Attributes do
    begin
      for i := 0 to Length - 1 do
      begin
        writeln(Item[i].NodeName + ' = ' + Item[i].NodeValue);
        Global.Add(Item[i].NodeName, Item[i].NodeValue);
      end;
    end;
    MChemStruc := entry.FindNode('MChemicalStruct');
    molecule   := MChemStruc.FindNode('molecule');
    molitem    := molecule.FirstChild;
    while molitem <> nil do
    begin
      writeln('molitem: ' + molitem.NodeName + ' / ' + molitem.NodeValue);
      if molitem.NodeName = 'atomArray' then //Start of atom block
      begin
        n    := 0;
        atom := molitem.FirstChild;
        while (atom <> nil) do
        begin
          new(PAt);
          Inc(n);
          while (n >= AtmID.Count) do
            AtmID.Add('');
          SetLength(PAt^.P, APrpSze);
          SetLength(PAt^.I, ABytSze);
          PAt^.S    := '';
          PAt^.Z    := ZMax;
          PAt^.W    := 0;
          AtmSet[n] := PAt;
          X[n]      := 0;
          Y[n]      := 0;
          Z[n]      := 0;
          MrvAtmSet[n] := 0;
          with atom.Attributes do
          begin
            for i := 0 to Length - 1 do
            begin
              if Item[i].NodeName = 'id' then
                AtmID.Strings[n] := Item[i].NodeValue;
              if (Item[i].NodeName = 'x2') or (Item[i].NodeName = 'x3') then
                X[n] := StrToFloat(Item[i].NodeValue);
              if (Item[i].NodeName = 'y2') or (Item[i].NodeName = 'y3') then
                Y[n] := StrToFloat(Item[i].NodeValue);
              if (Item[i].NodeName = 'z3') then
                Z[n] := StrToFloat(Item[i].NodeValue);
              if Item[i].NodeName = 'elementType' then
              begin
                S_[n] := Item[i].NodeValue;
                Z_[n] := AtomSymbolToInt(PAt^.S);
              end;
              if Item[i].NodeName = 'mrvSetSeq' then
                MrvAtmSet[n] := StrToInt(Item[i].NodeValue);
              writeln('atom ' + Item[i].NodeName + ' = ' + Item[i].NodeValue);
            end;
          end;
          atom := atom.NextSibling;
        end;
        InitLocalArrays(n);
      end; //End of atom block
      if molitem.NodeName = 'bondArray' then //Start of bond block
      begin
        m    := 0;
        bond := molitem.FirstChild;
        while (bond <> nil) do
        begin
          new(PBd);
          Inc(m);
          SetLength(PBd^.P, BPrpSze);
          SetLength(PBd^.I, BBytSze);
          PBd^.S    := '';
          PBd^.B    := BMax;
          PBd^.t    := 0;
          PBd^.h    := 0;
          PBd^.I[0] := 0;
          PBd^.I[1] := 0;
          with bond.Attributes do
          begin
            for i := 0 to bond.Attributes.Length - 1 do
            begin
              if (Item[i].NodeName = 'order') then
              begin
                PBd^.B := StrToInt(Item[i].NodeValue);
                PBd^.S := BondSymbol[PBd^.B];
              end;
              if (Item[i].NodeName = 'atomRefs2') then
              begin
                tmpSL.Clear;
                tmpSL.DelimitedText := Item[i].NodeValue;
                j1 := AtmID.IndexOf(tmpSL[0]);
                j2 := AtmID.IndexOf(tmpSL[1]);
                if (j1 >= 0) then
                  PBd^.t := j1;
                if (j2 >= 0) then
                  PBd^.h := j2;
                //writeln(tmpSL[0]+' '+tmpSL[1]+' '+IntToStr(j1)+' '+IntToStr(j2)+' '+AtmID[j1]+' '+AtmID[j2]);
              end;
              writeln('bond attributes ' + Item[i].NodeName +
                ' = ' + Item[i].NodeValue);
            end;
            if Item[i].NodeName = 'mrvSetSeq' then
              PBd^.I[1] := StrToInt(Item[i].NodeValue);
          end;
          bondchld := bond.FirstChild;
          while (bondchld <> nil) do
          begin
            if bondchld.NodeName = 'bondStereo' then
              if bondchld.FirstChild.NodeValue = 'H' then
                PBd^.I[0] := 1
              else if bondchld.FirstChild.NodeValue = 'W' then
                PBd^.I[0] := 2
              else if bondchld.FirstChild.NodeValue = 'T' then
                PBd^.I[0] := 3
              else if bondchld.FirstChild.NodeValue = 'C' then
                PBd^.I[0] := 4;
            writeln('bond childs: ' + bondchld.NodeName + ' / ' +
              bondchld.FirstChild.NodeValue);
            bondchld := bondchld.NextSibling;
          end;

          atmp[PBd^.t, PBd^.h] := m;
          atmp[PBd^.h, PBd^.t] := m;
          wtmp[m] := PBd;
          bond    := bond.NextSibling;
        end;
      end;

      p_NX := n;
      p_NY := 0;
      p_M  := m;
      //Store the graph in a packed format
      M    := 0;
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

      molitem := molitem.NextSibling;
    end;
    entry := entry.NextSibling;
  end;
  FreeAndNil(doc);
  FreeAndNil(tmpSL);
end;

procedure TMoleculeMrv.WriteMRV(mrfle: string);
var
  doc:      TXMLDocument;
  root, entry, MChemStruc, molecule, atomarray, bondarray, atom, bond,
  bondchld: TDOMNode;
  i, j, j1, j2, n, m, Mbnd, s, t: integer;
  PAt:      PRAtom;
  PBd:      PRBond;
  tmpSL:    TStringList;
  atmp:     array [1..MaxAtom] of array [1..MaxAtom] of integer;

  procedure InitLocalArrays(dima: integer);
  var
    i, j, dimb: integer;
  begin
    //Init
    for i := 1 to dima do
      for j := 1 to dima do
        atmp[i, j] := 0;
  end;

begin
  doc   := TXMLDocument.Create;
  root  := doc.CreateElement('cml');
  //cml->Marvin Document
  entry := doc.CreateElement('MDocument');
  TDOMElement(entry).SetAttribute('atomSetRGB', Global.Items['atomSetRGB']);
  TDOMElement(entry).SetAttribute('atomSetFont', Global.Items['atomSetFont']);
  TDOMElement(entry).SetAttribute('bondSetRGB', Global.Items['bondSetRGB']);
  TDOMElement(entry).SetAttribute('bondSetLineThickness',
    Global.Items['bondSetLineThickness']);
  //Marvin Document -> Marvin Chemical Structure
  MChemStruc := doc.CreateElement('MChemicalStruct');
  molecule   := doc.CreateElement('molecule');
  TDOMElement(molecule).SetAttribute('molID', 'm1');
  //Marvin Chemical Structure -> atomArray -> atom
  atomarray := doc.CreateElement('atomArray');
  for i := 1 to p_NX do
  begin
    atom := doc.CreateElement('atom');
    TDOMElement(atom).SetAttribute('id', 'a' + IntToStr(i));
    TDOMElement(atom).SetAttribute('elementType', string(AtmSet[i]^.S));
    if (Z[i] = 0) then
    begin
      TDOMElement(atom).SetAttribute('x2', FloatToStr(X[i]));
      TDOMElement(atom).SetAttribute('y2', FloatToStr(Y[i]));
    end
    else
    begin
      TDOMElement(atom).SetAttribute('x3', FloatToStr(X[i]));
      TDOMElement(atom).SetAttribute('y3', FloatToStr(Y[i]));
      TDOMElement(atom).SetAttribute('z3', FloatToStr(Z[i]));
    end;
    if MrvAtmSet[i] <> 0 then
      TDOMElement(atom).SetAttribute('mrvSetSeq', IntToStr(MrvAtmSet[i]));
    atomarray.AppendChild(atom);
  end;
  //Marvin Chemical Structure -> bondArray -> bond
  InitLocalArrays(p_M);
  bondarray := doc.CreateElement('bondArray');
  for i := 1 to p_M do
  begin
    if (atmp[BndSet[i]^.t, BndSet[i]^.h] = 0) then
    begin
      bond := doc.CreateElement('bond');
      TDOMElement(bond).SetAttribute('atomRefs2', 'a' + IntToStr(
        BndSet[i]^.t) + ' a' + IntToStr(BndSet[i]^.h));
      TDOMElement(bond).SetAttribute('order', IntToStr(BndSet[i]^.B));
      if MrvStereo[i] <> 0 then
      begin
        bondchld := doc.CreateElement('bondStereo');
        case MrvStereo[i] of
          1: bondchld.AppendChild(doc.CreateTextNode('H'));
          2: bondchld.AppendChild(doc.CreateTextNode('W'));
          3: bondchld.AppendChild(doc.CreateTextNode('T'));
          4: bondchld.AppendChild(doc.CreateTextNode('C'));
        end;
        bond.AppendChild(bondchld);
      end;
      if MrvBndSet[i] <> 0 then
      begin
        TDOMElement(bond).SetAttribute('mrvSetSeq', IntToStr(MrvBndSet[i]));
      end;
      bondarray.AppendChild(bond);
      atmp[BndSet[i]^.t, BndSet[i]^.h] := i;
      atmp[BndSet[i]^.h, BndSet[i]^.t] := i;
    end;
  end;
  //Append childs
  molecule.AppendChild(atomarray);
  molecule.AppendChild(bondarray);
  MChemStruc.AppendChild(molecule);
  entry.AppendChild(MChemStruc);
  root.AppendChild(entry);
  doc.AppendChild(root);
  //Write XML and close the document
  WriteXMLFile(doc, mrfle);
  FreeAndNil(doc);
end;

end.

