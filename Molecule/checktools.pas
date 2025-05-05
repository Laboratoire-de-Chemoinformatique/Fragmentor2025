unit checktools;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Strutils, intlist;

type
  TAtomDict = class(TObject)
  public
    symbol: TStringList;
    minValence: TIntlist;
    maxValence: TIntList;
    bsymbol: TStringList;
    bord: Tintlist;
    constructor Create();
  end;

function isAtom(const smi: string; var pos: integer; var symb: string;
  dict: TAtomDict): boolean;
function isBond(const smi: string; var pos, currbord: integer; dict: TAtomDict): boolean;
function isH(const smi: string; var pos: integer): boolean;
function isdot(const smi: string; var pos: integer): boolean;
function isMark(const smi: string; var pos: integer): boolean;
function openpar(const smi: string; var pos: integer): boolean;
function closepar(const smi: string; var pos: integer): boolean;
function openbra(const smi: string; var pos: integer): boolean;
function closebra(const smi: string; var pos: integer): boolean;
function isNr(const smi: string; var pos: integer; maxlen: byte): byte;
function isRingClosure(const smi: string; var pos: integer): byte;


implementation

function isRingClosure(const smi: string; var pos: integer): byte;
begin
  if MidStr(smi, pos, 1)='%' then begin
    inc(pos);
    Result := isNr(smi,pos,2);
  end
  else
    Result := isNr(smi,pos,1);
end;


function isAtom(const smi: string; var pos: integer; var symb: string;
  dict: TAtomDict): boolean;
begin
  Result := False;
  symb := MidStr(smi, pos, 2);
  if dict.symbol.indexOf(symb) > -1 then
  begin
    pos := pos + 2;
    Result := True;
  end
  else
  begin
    symb := MidStr(smi, pos, 1);
    if dict.symbol.indexOf(symb) > -1 then
    begin
      pos := pos + 1;
      Result := True;
    end;
  end;
end;


function isBond(const smi: string; var pos, currbord: integer; dict: TAtomDict): boolean;
var
  i: integer;
begin
  Result := false;
  i := dict.bsymbol.indexOf(MidStr(smi, pos, 2));
  if i > -1 then
  begin
    pos := pos + 2;
    currbord := dict.bord[i];
    Result := true;
  end
  else
  begin
    i := dict.bsymbol.indexOf(MidStr(smi, pos, 1));
    if i > -1 then
    begin
      pos := pos + 1;
      currbord := dict.bord[i];
      Result := true;
    end;
  end;
end;



function isH(const smi: string; var pos: integer): boolean;
begin
  Result := (MidStr(smi, pos, 1) = 'H');
  if Result then
    Inc(pos);
end;


function isDot(const smi: string; var pos: integer): boolean;
begin
  Result := (MidStr(smi, pos, 1) = '.');
  if Result then
    Inc(pos);
end;


function isMark(const smi: string; var pos: integer): boolean;
begin
  Result := (MidStr(smi, pos, 1) = ':');
  if Result then
    Inc(pos);
end;

function openpar(const smi: string; var pos: integer): boolean;
begin
  Result := (MidStr(smi, pos, 1) = '(');
  if Result then
    Inc(pos);
end;

function closepar(const smi: string; var pos: integer): boolean;
begin
  Result := (MidStr(smi, pos, 1) = ')');
  if Result then
    Inc(pos);
end;

function openbra(const smi: string; var pos: integer): boolean;
begin
  Result := (MidStr(smi, pos, 1) = '[');
  if Result then
    Inc(pos);
end;

function closebra(const smi: string; var pos: integer): boolean;
begin
  Result := (MidStr(smi, pos, 1) = ']');
  if Result then
    Inc(pos);
end;


function isNr(const smi: string; var pos: integer; maxlen: byte): byte;
var
  err: integer;
begin
  while maxlen > 0 do
  begin
    val(MidStr(smi, pos, maxlen), Result, err);
    if err <> 0 then
    begin
      Result := 0;
      Dec(maxlen);
    end
    else
    begin
      pos := pos + maxlen;
      exit;
    end;
  end;
end;


constructor TAtomDict.Create;
begin
  inherited Create();
  symbol := TStringList.Create;
  bsymbol := TStringList.Create;
  symbol.casesensitive := True;
  minValence := TIntList.Create;
  maxValence := TIntList.Create;
  bord := TIntList.Create;
  symbol.add('H');
  minValence.add(1);
  maxValence.add(1);
  symbol.add('C');
  minValence.add(4);
  maxValence.add(4);
  symbol.add('c');
  minValence.add(4);
  maxValence.add(4);
  symbol.add('N');
  minValence.add(3);
  maxValence.add(4);
  symbol.add('n');
  minValence.add(3);
  maxValence.add(4);
  symbol.add('O');
  minValence.add(2);
  maxValence.add(2);
  symbol.add('o');
  minValence.add(3);
  maxValence.add(3);
  symbol.add('P');
  minValence.add(3);
  maxValence.add(5);
  symbol.add('p');
  minValence.add(3);
  maxValence.add(5);
  symbol.add('F');
  minValence.add(1);
  maxValence.add(1);
  symbol.add('Cl');
  minValence.add(1);
  maxValence.add(1);
  symbol.add('Br');
  minValence.add(1);
  maxValence.add(1);
  symbol.add('I');
  minValence.add(1);
  maxValence.add(1);
  symbol.add('S');
  minValence.add(2);
  maxValence.add(6);
  symbol.add('s');
  minValence.add(2);
  maxValence.add(6);
  symbol.add('Se');
  minValence.add(2);
  maxValence.add(6);
  symbol.add('se');
  minValence.add(2);
  maxValence.add(6);
  symbol.add('Zn');
  minValence.add(0);
  maxValence.add(8);
  symbol.add('Mg');
  minValence.add(0);
  maxValence.add(8);
  symbol.add('Fe');
  minValence.add(0);
  maxValence.add(8);
  symbol.add('Co');
  minValence.add(0);
  maxValence.add(8);
  symbol.add('Hg');
  minValence.add(0);
  maxValence.add(8);
  symbol.add('Si');
  minValence.add(4);
  maxValence.add(4);
  symbol.add('Ca');
  minValence.add(2);
  maxValence.add(2);
  symbol.add('B');
  minValence.add(1);
  maxValence.add(4);
  symbol.add('As');
  minValence.add(3);
  maxValence.add(5);

  bsymbol.add('-');
  bord.add(1);
  bsymbol.add('+');
  {dummy "bond" - since "-" is context-dependently single bond or negative formal charge, let "+" also figure among possible bond characters}
  bord.add(-1);
  bsymbol.add('=');
  bord.add(2);
  bsymbol.add('#');
  bord.add(3);
  bsymbol.add('@');
  bord.add(11);
  bsymbol.add('@@');
  bord.add(21);
  bsymbol.add('\');
  bord.add(31);
  bsymbol.add('/');
  bord.add(41);
  bsymbol.add('--');
  bord.add(2);
  bsymbol.add('++');
  bord.add(-2);
end;

end.
