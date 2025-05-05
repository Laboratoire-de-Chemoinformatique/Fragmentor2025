unit ReadSmiles;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, StrUtils, intlist, checktools, Math;

type

  Bondtype = record
    First: integer; { ID of first atom}
    second: integer; { ID of second atom}
    iorder: byte;    { Integer bond order:1, 2, 3 or 4 for aromatic}
    order: real;    { Real bond order value: 1.0, 1.5 for aromatic, 2.0 or 3.0}
  end;
  pBondtype = ^Bondtype;

  TAtom = class(TObject)
  public
    Symbol: string;
    fcharge: integer;
    isotope: integer;
    chiralflag, map: byte; {chirality info and atom mapping nr if provided}
    implicitH: integer;
    aromatic: boolean;
    ringConnect: TFPList;
    ringborder: Tintlist;
    parent: Tatom;
    parentbond: byte;
    this: integer; {index in atom list}
    constructor Create(const s: string);
    destructor Destroy;
  end;

  TMolecule = class(TObject)
  public
    atom, bond: TFPList;      { list of atoms and bonds of the molecule}
    constructor Create(const smi: string; dict: TAtomDict);
    destructor Destroy;
    function natoms: integer; {returns the number of atoms in the molecule}
    function nbonds: integer; {returns the number of atoms in the molecule}
    procedure save(var outfile: Text);
  end;

procedure ProcSmiles(const filename, outname: string; const smicol: integer);
function FirstUpper(const s: string): string;

implementation

function FirstUpper(const s: string): string;
var
  head: string;

begin
  head := uppercase(LeftStr(s, 1));
  if length(s) > 1 then
    Result := head + RightStr(s, length(s) - 1)
  else
    Result := head;
end;

constructor TAtom.Create(const s: string);
begin
  inherited Create();
  Symbol := s;
  fcharge := 0;
  isotope := 0;
  aromatic := (Symbol = lowercase(Symbol));
  ringConnect := TFPList.Create();
  ringborder := Tintlist.Create();
  parent := nil;
  parentbond := 0;
  this := 0;
  chiralflag := 0;
  map := 0;
  if aromatic then
    Symbol := FirstUpper(Symbol);
end;

destructor TAtom.Destroy;
begin
  FreeAndNil(ringConnect);
  FreeAndNil(ringborder);
  FreeAndNil(self);
end;

constructor TMolecule.Create(const smi: string; dict: TAtomDict);
var
  pos, isotope, currbord, r, rr: integer;
  symb: string;
  curratom, currparent: TAtom;
  inbra: byte;
  lastref, ringhead: TFPlist;
  ringlist: Tintlist;
  bp: pBondType;
begin
  inherited Create();
  atom := TFPList.Create();
  bond := TFPList.Create();
  lastref := TFPList.Create();
  ringhead := TFPList.Create();
  ringlist := TIntList.Create();
  pos := 1;
  currparent := nil;
  currbord := 0;
  isotope := 0; {temporary storage of isotope label if given}
  symb := '';
  inbra := 0; {true when currently reading within [atom] square brackets}
  repeat
    if inbra > 0 then
      Inc(inbra);
    if isAtom(smi, pos, symb, dict) and ((symb <> 'H') or (inbra = 2)) then
    begin
      {accept a H symbol as atom only if it follows immediately after the [ sign}
      curratom := TAtom.Create(symb);
      atom.add(curratom);
      curratom.parent := currparent;
      curratom.this := natoms;
      curratom.isotope := isotope;
      isotope := 0; {reset isotope label after assignment}
      currparent := curratom;
      curratom.parentbond := currbord;
      currbord := 0; {reset currbord after assigning to an atom}
    end

    else if not isBond(smi, pos, currbord, dict) then

      if inbra > 0 then {if we are within an [atom] specification block}
        if closebra(smi, pos) then
        begin
          inbra := 0;  {[atom] block terminated}
          if currbord > 2 then
            {must have captured stereobond in bracketed atom specification [X@H] - asign it to the last atom in the list, not to the following atom}
            TAtom(atom[natoms - 1]).chiralflag := currbord
          else if abs(currbord) < 3 then
            {we are inbra but found "-" or "+" signs, that is "bond orders" 1 or -1 - interpret as formal charge}
            TAtom(atom[natoms - 1]).fcharge :=
              TAtom(atom[natoms - 1]).fcharge - currbord;
          currbord := 0;
        end
        else
        begin
          if isMark(smi, pos) then
          begin
            {found a ":" - atom mapping applies}
            TAtom(atom[natoms - 1]).map := isNr(smi, pos, 2);
            if TAtom(atom[natoms - 1]).map = 0 then
            begin
              writeln('FATAL - expecting numerical atom map number after the : at character ', pos, ' in ', smi);
              Halt;
            end;
          end
          else if (not isH(smi, pos)) and (isNr(smi, pos, 1) = 0) then
          begin
            writeln('FATAL - illegal character at position ', pos, ' in ', smi);
            Halt;
          end;
        end  {end inbra scenario}

      else if openbra(smi, pos) then
      begin
        inbra := 1;  {[atom] block just opened - check if isotope data is given first}
        isotope := isNr(smi, pos, 2);
      end
      else
      begin
        {we are not in any [atom] block}
        {check parentheses}
        if isDot(smi, pos) then
          currparent := nil {disjoined pieces}
        else if openpar(smi, pos) then
          lastref.add(atom[natoms - 1])
        {by default, the last seen atom becomes the future parent after closepar}
        else if closepar(smi, pos) then
        begin
          {reset the reference atom to the latest in the queue lastref}
          if lastref.Count = 0 then
          begin
            writeln('FATAL - closing paranthesis at position ', pos -
              1, ' has no match');
            Halt;
          end;
          currparent := TAtom(lastref[lastref.Count - 1]);
          lastref.Delete(lastref.Count - 1);
          if openpar(smi, pos) then
            lastref.add(currparent);
          {if openpar follows immediately after closedpar, then the parent is not the last encountered atom but the last encountered parent}
        end
        else
        begin
          {the only other option to encounter here is a number: ring marker}
          r := isRingClosure(smi, pos);
          if r = 0 then
          begin
            writeln('FATAL - expecting a ring closure label at ', pos, ' in ', smi);
            Halt;
          end;
          rr := ringlist.find(r);
          if rr < 0 then
          begin
            {r is not an active ring label: add it to ringlist, and the latest encountered atom to the ringhead list}
            ringlist.add(r);
            ringhead.add(atom[natoms - 1]);
            currbord := 0; {ring label also resets latest bond order}
          end
          else
          begin
            {r is an open ring label: link closure to data in positions rr of the lists}
            TAtom(atom[natoms - 1]).ringconnect.add(ringhead[rr]);
            TAtom(atom[natoms - 1]).ringborder.add(currbord);
            ringhead.Delete(rr);
            ringlist.Delete(rr);
            currbord := 0;
          end;
        end;
      end;

  until pos > length(smi);
  if lastref.Count > 0 then
  begin
    writeln('FATAL - unclosed parantheses seem to exist in ', smi);
    Halt;
  end;
  FreeAndNil(lastref);
  if (ringhead.Count > 0) or (ringlist.IntCount > 0) then
  begin
    writeln('FATAL - something went wrong with ring closure marks in ', smi);
    Halt;
  end;
  FreeAndNil(ringhead);
  FreeAndNil(ringlist);
  if inbra > 0 then
  begin
    writeln('FATAL - something went wrong with [ ] delimitors in ', smi);
    Halt;
  end;

  {make bond list}
  for r := 0 to natoms - 1 do
  begin
    curratom := TAtom(atom[r]);
    {add direct bond to parent to the list}
    if curratom.parent <> nil then
    begin
      new(bp);
      bond.add(bp);
      bp^.First := min(curratom.this, curratom.parent.this);
      bp^.second := max(curratom.this, curratom.parent.this);
      if curratom.parentbond = 0 then
        if curratom.aromatic and curratom.parent.aromatic then
        begin
          bp^.iorder := 4;
          bp^.order := 1.5;
        end
        else
        begin
          bp^.iorder := 1;
          bp^.order := 1.0;
        end
      else
      begin
        bp^.iorder := curratom.parentbond mod 10;
        bp^.order := curratom.parentbond mod 10;
      end;
    end;
    {Add any ring closure bonds}
    for rr := 0 to curratom.ringconnect.Count - 1 do
    begin
      new(bp);
      bond.add(bp);
      bp^.First := min(curratom.this, TAtom(curratom.ringconnect[rr]).this);
      bp^.second := max(curratom.this, TAtom(curratom.ringconnect[rr]).this);
      if curratom.ringborder[rr] = 0 then
        if curratom.aromatic and TAtom(curratom.ringconnect[rr]).aromatic then
        begin
          bp^.iorder := 4;
          bp^.order := 1.5;
        end
        else
        begin
          bp^.iorder := 1;
          bp^.order := 1.0;
        end
      else
      begin
        bp^.iorder := curratom.ringborder[rr] mod 10;
        if bp^.iorder < 4 then
          bp^.order := bp^.iorder
        else
          bp^.order := 1.5;
      end;
    end;
  end;
end;



function TMolecule.natoms: integer;
begin
  exit(atom.Count);
end;

function TMolecule.nbonds: integer;
begin
  exit(bond.Count);
end;

procedure TMolecule.save(var outfile: Text);
var
  i, nfc, niso: integer;
  fcstring, isostring, acstring: string; {to resume formal charges & isotopes}
  a: TAtom;
begin
  writeln(outfile, '');
  writeln(outfile, '  SMI2SDF 00000000002D');
  writeln(outfile, '');
  writeln(outfile, PadLeft(IntToStr(natoms), 3), PadLeft(IntToStr(nbonds), 3),
    '  0  0  0  0            999 V2000');
  fcstring := '';
  nfc := 0;
  isostring := '';
  niso := 0;
  for i := 0 to natoms - 1 do
  begin
    acstring := '  0';
    {This adds atom mark labels in columns required by chemaxon and fragmentor}
    a := TAtom(atom[i]);
    if a.fcharge <> 0 then
    begin
      Inc(nfc);
      fcstring := fcstring + PadLeft(IntToStr(a.this), 4) +
        PadLeft(IntToStr(a.fcharge), 4);
      acstring := PadLeft(IntToStr(4 - a.fcharge), 3);
    end;
    writeln(outfile, FormatFloat('0.0000', 0): 10,
      FormatFloat('0.0000', 0): 10, FormatFloat('0.0000', 0): 10, ' ',
      PadRight(a.symbol, 3): 3, ' 0', acstring, '  0  0  0  0  0',
      PadLeft(IntToStr(min(1, a.map)), 3), '  0', PadLeft(
      IntToStr(a.map), 3), '  0  0');
    if a.isotope <> 0 then
    begin
      Inc(niso);
      isostring := isostring + PadLeft(IntToStr(a.this), 4) +
        PadLeft(IntToStr(a.isotope), 4);
    end;
  end;
  for i := 0 to nbonds - 1 do
    writeln(outfile, PadLeft(IntToStr(pBondtype(bond[i])^.First), 3): 3,
      PadLeft(IntToStr(pBondtype(bond[i])^.second), 3): 3,
      PadLeft(IntToStr(pBondtype(bond[i])^.iorder), 3): 3,
      '  0  0  0  0');
  if niso > 0 then
    writeln(outfile, 'M  ISO', PadLeft(IntToStr(niso), 3), isostring);
  if nfc > 0 then
    writeln(outfile, 'M  CHG', PadLeft(IntToStr(nfc), 3), fcstring);
  writeln(outfile, 'M  END');
  writeln(outfile, '$$$$');
end; { TMolecule.writemdl}



destructor TMolecule.Destroy;
var
  i: integer;
  bp: ^BondType;

begin
  for i := 0 to natoms - 1 do
    TAtom(atom[i]).Destroy;
  FreeAndNil(atom);
  for i := 0 to bond.Count - 1 do
  begin
    bp := bond[i];
    Dispose(bp);
  end;
  FreeAndNil(bond);

  FreeAndNil(self);
end;

procedure ProcSmiles(const filename, outname: string; const smicol: integer);
var
  inp, outp: Text;
  line, del: string;
  mol: TMolecule;
  dict: TAtomDict;
  n, i: integer;
begin
  if not FileExists(filename) then
  begin
    writeln('FATAL - input file ', filename, ' does not exist');
    Halt;
  end;
  dict := TAtomDict.Create();
  Assign(inp, filename);
  reset(inp);
  Assign(outp, outname);
  rewrite(outp);
  del := '';
  for i := 1 to 50 do
    del := del + chr(8);
  n := 0;
  while not EOF(inp) do
  begin
    Inc(n);
    Write('Processing item ', n);
    readln(inp, line);
    mol := TMolecule.Create(ExtractWord(smicol, line, [#32, #9]), dict);
    mol.save(outp);
    mol.Destroy;
    Write(del);
  end;
  writeln(' ');
  Close(outp);
  Close(inp);
end;


end.
