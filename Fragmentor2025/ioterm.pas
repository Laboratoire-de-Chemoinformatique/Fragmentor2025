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
unit IOterm;

{$mode objfpc}{$H+}

interface

uses
  getopts, strutils, Classes, SysUtils, UnitFragmentBase, UnitFragment, contnrs;

const
  ArgLst = 'i:o:t:l:u:f:s:m:h:b:c:d:x:z:';
  UseMsg = 'Usage: Fragmentor -i <input> -o <output> -x <int>:<xml> -t <type> ' +
    '-l <lower> -u <upper> -f <oformat> -s <sfield> -c <colorfields> -h <header> -m <markedatom> -d <dynamicbonds>' +
    ' --DoAllWays --UseFormalCharge --UseRadical --UseIsotope --AtomPairs --Cycle --StrictFrg --GetAtomFragment --Pipe';

type

  { TIOterm }

  TIOterm = class(TObject)
  private
    IFileName, OFileName, HFileName, XFileName: string;
    oformat, ifield: string;
    FrgType:     LDefFrg;
    RWSFrgXML: byte; //Manage the XML option: read, write or simulate XML file
    fTbLongOpt:  array[1..10] of TOption;
    fPTbLongOpt: POption;
    fbPipe: boolean; //Output on standard output
    fbStrict: boolean; //Strict Fragment
    fbGetA: boolean; //Get Atom Fragment
    fiRS: integer; //Save fragments with this periodicity
    procedure Interpret(a: char; opt: string; LOptIndex: integer);
    procedure Integrity();
  public
    constructor Create();
    destructor Destroy(); override;
    property SGIFileName: string Read IFileName Write IFileName;
    property SGOFileName: string Read OFileName Write OFileName;
    property SGHFileName: string Read HFileName Write HFileName;
    property SGXFileName: string Read XFileName Write XFileName;
    property SGRWSFrgXML: byte Read RWSFrgXML write RWSFrgXML;
    property SGFrgType: LDefFrg Read FrgType Write FrgType;
    property SGoformat: string Read oformat Write oformat;
    property SGifield: string Read ifield Write ifield;
    property SGbPipe: boolean Read fbPipe Write fbPipe;
    property SGbStrict: boolean Read fbStrict Write fbStrict;
    property SGbGetA: Boolean Read fbGetA Write fbGetA;
    property SGiRS: integer Read fiRS Write fiRS;
    procedure Display();
    procedure WriteHeader(frglst: TStringList);
    procedure WriteSVM(prop: string; frglst: TStringList);
    procedure WriteSParse(prop: string; frglst: TStringList);
    procedure CleanSParse(NFrg: integer);
    procedure WriteSMF(frglst: TStringList);
    procedure WriteProp(prop: string);
    procedure WriteAtmFrg(AtomFrg: TObjectList);
  end;

  TIOtermException = class(Exception)
  end;


implementation

constructor TIOterm.Create();
var
  LInd: longint;
  Arg:  char;
begin
  writeln('*********************************************');
  writeln('                                             ');
  writeln('                ISIDA Fragmentor2012         ');
  writeln('       This is a terminal interface for      ');
  writeln('       generating molecular fragments        ');
  writeln('                                             ');
  writeln('       G. Marcou, D. Horvath, A. Varnek      ');
  writeln('     F. Ruggiu, V. Solov''ev, E. Moyemont    ');
  writeln('                                             ');
  writeln('                      2022                   ');
  writeln('                                             ');
  writeln(' Universite de Strasbourg                    ');
  writeln(' Faculte de Chimie                           ');
  writeln(' Laboratoire d''infochimie                   ');
  writeln('                                             ');
  writeln('*********************************************');
  IFileName := '';
  OFileName := '';
  HFileName := '';
  XFileName := 'default.xml';
  RWSFrgXML := 0;
  oformat   := 'SVM';
  ifield    := '';
  fbPipe    := False;
  fbStrict  := False;
  fbGetA    := False;
  fiRS      := 0;

  with fTbLongOpt[1] do
  begin
    Name    := 'DoAllWays';
    Has_arg := 0;
    Flag    := nil;
    Value   := #0;
  end;
  with fTbLongOpt[2] do
  begin
    Name    := 'EqFuzAt';
    Has_arg := 0;
    Flag    := nil;
    Value   := #0;
  end;
  with fTbLongOpt[3] do
  begin
    Name    := 'AtomPairs';
    Has_arg := 0;
    Flag    := nil;
    Value   := #0;
  end;
  with fTbLongOpt[4] do
  begin
    Name    := 'StrictFrg';
    Has_arg := 0;
    Flag    := nil;
    Value   := #0;
  end;
  with fTbLongOpt[5] do
  begin
    Name    := 'UseFormalCharge';
    Has_arg := 0;
    Flag    := nil;
    Value   := #0;
  end;
  with fTbLongOpt[6] do
  begin
    Name    := 'UseRadical';
    Has_arg := 0;
    Flag    := nil;
    Value   := #0;
  end;
  with fTbLongOpt[7] do
  begin
    Name    := 'UseIsotope';
    Has_arg := 0;
    Flag    := nil;
    Value   := #0;
  end;
  with fTbLongOpt[8] do
  begin
    Name    := 'GetAtomFragment';
    Has_arg := 0;
    Flag    := nil;
    Value   := #0;
  end;
  with fTbLongOpt[9] do
  begin
    Name    := 'Pipe';
    Has_arg := 0;
    Flag    := nil;
    Value   := #0;
  end;
  with fTbLongOpt[10] do
  begin
    Name    := 'Cycle';
    Has_arg := 0;
    Flag    := nil;
    Value   := #0;
  end;
  fPTbLongOpt := @fTbLongOpt[1];

  FrgType := LDefFrg.Create;
  Arg     := GetLongOpts(ArgLst, fPTbLongOpt, LInd);
  //Arg:=GetOpt(ArgLst);
  if (Arg = EndOfOptions) then
  begin
    writeln(UseMsg);
    writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    raise TIOtermException.Create(UseMsg);
  end;
  while (Arg <> EndOfOptions) do
  begin
    Interpret(Arg, OptArg, LInd);
    Arg := GetLongOpts(ArgLst, fPTbLongOpt, LInd);
  end;
  Integrity;
end;

destructor TIOterm.Destroy();
begin
  FrgType.Destroy;
  inherited Destroy;
end;

procedure TIOterm.Interpret(a: char; opt: string; LOptIndex: integer);
var
  Frg: PDefFrg;

begin
  if (a = '?') then
  begin
    writeln(UseMsg);
    writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    raise TIOtermException.Create(UseMsg);
  end;
  if (a = 'i') then
  begin
    SGIFileName := opt;
    if SGOFileName = '' then
      SGOFileName:=ChangeFileExt(SGIFileName,'');
  end
  else if (a = 'o') then
  begin
    SGOFileName := opt;
  end
  else if (a = 't') then
  begin
    new(Frg);
    InitDefFrg(Frg);
    Frg^.FrgType   := Numb2Dec(opt, 10);
    {Frg^.lmin      := 0;
    Frg^.lmax      := 0;
    Frg^.MarkAtom := 0;
    Frg^.DynBnd    := 0;
    Frg^.ColorAFields := 'Default';
    Frg^.ColorBFields := 'Default';
    Frg^.DoAllWays := False;
    Frg^.UseFormalCharge := False;
    Frg^.UseRadical:=False;
    Frg^.UseIsotope:=False;
    Frg^.EqFuzAt := False;
    Frg^.AtomPairs := False;
    Frg^.StrictFrg := False;
    Frg^.GetAtomFrg:=False;
    Frg^.Cycle:=False;}
    FrgType.Add(Frg);
  end
  else if (a = 'l') then
  begin
    Frg := FrgType[FrgType.Count - 1];
    if (Frg <> nil) then
      Frg^.lmin := Numb2Dec(opt, 10);
  end
  else if (a = 'u') then
  begin
    Frg := FrgType[FrgType.Count - 1];
    if (Frg <> nil) then
      Frg^.lmax := Numb2Dec(opt, 10);
  end
  else if (a = #0) then
  begin
    try
      Frg := FrgType[FrgType.Count - 1];
      if (Frg <> nil) then
      begin
        if (fTbLongOpt[LOptIndex].Name = 'DoAllWays') then
          Frg^.DoAllWays := True
        else if (fTbLongOpt[LOptIndex].Name = 'EqFuzAt') then
          Frg^.EqFuzAt := True
        else if (fTbLongOpt[LOptIndex].Name = 'AtomPairs') then
          Frg^.AtomPairs := True
        else if (fTbLongOpt[LOptIndex].Name = 'StrictFrg') then
          fbStrict := True //Global fragment option
        else if (fTbLongOpt[LOptIndex].Name = 'UseFormalCharge') then
          Frg^.UseFormalCharge := True
        else if (fTbLongOpt[LOptIndex].Name = 'UseRadical') then
          Frg^.UseRadical := True
        else if (fTbLongOpt[LOptIndex].Name = 'UseIsotope') then
          Frg^.UseIsotope := True
        else if (fTbLongOpt[LOptIndex].Name = 'GetAtomFragment') then
          fbGetA := True //Global fragment option
        else if (fTbLongOpt[LOptIndex].Name = 'Pipe') then
          fbPipe:=True
        else if (fTbLongOpt[LOptIndex].Name = 'Cycle') then
          Frg^.Cycle :=True;
      end;
    except
      Writeln('ERROR: Long option '+fTbLongOpt[LOptIndex].Name+' could not be associated to a given fragment definition.');
      halt;
    end;
  end
  else if (a = 'f') then
  begin
    SGoformat := opt;
  end
  else if (a = 's') then
  begin
    SGifield := opt;
  end
  else if (a = 'm') then
  begin
    Frg := FrgType[FrgType.Count - 1];
    if (Frg <> nil) then
    begin
      Frg^.MarkAtom := Numb2Dec(opt, 10);
    end;
  end
  else if (a = 'd') then
  begin
    Frg := FrgType[FrgType.Count - 1];
    if (Frg <> nil) then
    begin
      Frg^.DynBnd := Numb2Dec(opt, 10);
    end;
  end
  else if (a = 'h') then
  begin
    SGHFileName := opt;
  end
  else if (a='b') then
  begin
    Frg:=FrgType[FrgType.Count-1];
    if (Frg <> nil) then
       Frg^.ColorBFields:=opt;
  end
  else if (a = 'c') then
  begin
    Frg := FrgType[FrgType.Count - 1];
    if (Frg <> nil) then
    begin
      Frg^.ColorAFields := opt;
    end;
  end else if (a='x') then
  begin //XML fragmentation file
    //default value "0"
    //"1:file.xml" -> output the fragmentation xml code in file.xml
    //"2:file.xml" -> read the fragmentation xml code in file.xml, takes priority
    //"3:file.xml" -> output the fragmentation xml code in file.xml. Fragmentation does not take place.
    opt:=StringReplace(opt,':',' ',[]);
    ReadStr(opt,RWSFrgXML,XFileName);
    XFileName:=Trim(XFileName);
  end
  else if (a='z') then fiRS:=StrToInt(opt);
end;

procedure TIOterm.Display();
var
  i:     integer;
  pitem: PDefFrg;
begin
  pitem := nil;
  writeln('Input file: ' + SGIFileName);
  writeln('Output file: ' + SGOFileName);
  writeln('Output format: ' + SGoformat);
  if fbPipe then
    writeln('REMINDER: the output will be concatenated to the output file');
  if (HFileName <> '') then
    writeln('Header file: ' + SGHFileName);
  case RWSFrgXML of
    0: writeln('No XML');
    1: writeln('Fragmentation setup stored in: '+XFileName);
    2: writeln('Fragmentation setup read from: '+XFileName);
    3: writeln('Fragmentation setup stored in: '+XFileName+
       ' (Fragmentation simulated)');
    else writeln('WARNING: XML file option unkwown');
  end;
  if fiRS>0 then writeln('Header file saved each '+IntToStr(fiRS)+' molecule.');
  writeln('/////////////////////////////////');
  for i := 0 to FrgType.Count - 1 do
  begin
    pitem := FrgType[i];
    if (pitem <> nil) then
    begin
      writeln('Fragmentation type: ' + IntToStr(pitem^.FrgType));
      if (pitem^.FrgType >= 1) then
      begin
        writeln('Lower boundary: ' + IntToStr(pitem^.lmin));
        writeln('Upper boundary: ' + IntToStr(pitem^.lmax));
      end;
      writeln('MarkAtom: ' + IntToStr(pitem^.MarkAtom));
      writeln('Atom colors: ' + pitem^.ColorAFields);
      writeln('Bond colors: ' + pitem^.ColorBFields);
      if (pitem^.Cycle) then writeln('Annotate bonds in cycle');
      writeln('Equivalent Fuzzy Atom count (EqFuzAt): ' + BoolToStr(pitem^.EqFuzAt));
      writeln('UseFormalCharge: ' + BoolToStr(pitem^.UseFormalCharge));
      writeln('UseRadical: ' + BoolToStr(pitem^.UseRadical));
      writeln('UseIsotope : ' + BoolToStr(pitem^.UseIsotope));
      writeln('DoAllWays: ' + BoolToStr(pitem^.DoAllWays));
      writeln('AtomPairs: ' + BoolToStr(pitem^.AtomPairs));
      writeln('StrictFragmentation: ' + BoolToStr(pitem^.StrictFrg));
      writeln('DynamicBond: ' + IntToStr(pitem^.DynBnd));
    end;
    writeln('/////////////////////////////////');
  end;
  if (Length(SGifield) > 0) then
  begin
    writeln('Searching for field: ' + SGifield);
  end;
end;

procedure TIOterm.WriteHeader(frglst: TStringList);
var
  Outfile: Text;
  i: integer;

begin
  AssignFile(OutFile, SGOFileName + '.hdr');
  ReWrite(OutFile);
  for i := 0 to frglst.Count - 1 do
  begin
    writeln(Outfile, format('%6s', [IntToStr(i + 1)]) + '. ' +
      format('%120s', [frglst.Strings[i]]));
          {writeln(format('%6s', [IntToStr(i + 1)]) + '. ' +
                           format('%120s', [frglst.Strings[i]]));}
  end;
  Close(Outfile);
end;

procedure TIOterm.WriteSVM(prop: string; frglst: TStringList);
var
  Outfile: TextFile;
  i, j:    integer;

begin
  AssignFile(OutFile, SGOFileName + '.svm');
  if (FileExists(SGOFileName + '.svm')) then
    Append(OutFile)
  else
    ReWrite(OutFile);
  if(prop='') then prop:='?';
  Write(OutFile, prop);

  for i := 0 to frglst.Count - 1 do
  begin
    j := PInteger(frglst.Objects[i])^;   //pointeur n'a pas de type donc on doit lui donner(PInteger), ^ dit qu'on est pas intéressé par l'adresse mais par la valeur
    if (j <> 0) then
    begin
      Write(OutFile, ' ' + IntToStr(i + 1) + ':' + IntToStr(j));
    end;
  end;

  //WriteLn('Max nb of frg : '+IntToStr(frglst.Count - 1));

  //Write the last element
  //if(i=frglst.Count - 1) then Write(OutFile, ' ' + IntToStr(i + 1) + ':' + IntToStr(j));

  writeln(OutFile);
  Close(Outfile);
end;

procedure TIOterm.WriteSParse(prop: string; frglst: TStringList);
var
  Outfile: TextFile;
  i, j:    integer;

begin
  AssignFile(OutFile, SGOFileName + '.csv.tmp');
  if (FileExists(SGOFileName + '.csv.tmp')) then
    Append(OutFile)
  else
    Rewrite(OutFile);
  Write(OutFile, prop);
  for i := 0 to frglst.Count - 1 do
  begin
    j := PInteger(frglst.Objects[i])^;
    if (j <> 0) then
    begin
      Write(OutFile,';'+IntToStr(j));
    end else Write(OutFile,';0');
  end;
  WriteLn(OutFile);
  Close(Outfile);
end;

procedure TIOterm.CleanSParse(NFrg: integer);
//Add trailing 0 to CSV output.
var
  Tmpfile, Outfile: TextFile;
  SLline: TStringList;
  sline: string;
  i:    integer;

begin
  SLline:=TStringList.Create;
  SLline.Delimiter:=';';
  AssignFile(TmpFile, SGOFileName + '.csv.tmp');
  AssignFile(Outfile, SGOFileName + '.csv');
  Reset(Tmpfile);
  Rewrite(Outfile);
  repeat
    SLline.Clear;
    ReadLn(Tmpfile,sline);
    SLline.DelimitedText:=sline;
    for i:=SLline.Count to NFrg do
      sline:=sline+';0';
    Writeln(Outfile,sline);
  until(EOF(TmpFile));
  Close(Outfile);
  Close(Tmpfile);
  FreeAndNil(SLline);
  DeleteFile(SGOFileName + '.csv.tmp');
end;

procedure TIOterm.WriteSMF(frglst: TStringList);
var
  Outfile: TextFile;
  i, j:    integer;

begin
  AssignFile(OutFile, SGOFileName + '.smf');
  if (FileExists(SGOFileName + '.smf')) then
    Append(OutFile)
  else
    ReWrite(OutFile);
  for i := 0 to frglst.Count - 1 do
  begin
    j := PInteger(frglst.Objects[i])^;
    if (j <> 0) then
    begin
      Write(OutFile, IntToStr(i + 1) + ' ' + IntToStr(j) + ' ');
    end;
  end;
  writeln(OutFile);
  Close(Outfile);
end;

procedure TIOterm.WriteProp(prop: string);
var
  Outfile: TextFile;
begin
  AssignFile(OutFile, SGOFileName + '.prp');
  if (FileExists(SGOFileName + '.prp')) then
    Append(OutFile)
  else
    ReWrite(OutFile);
  writeln(OutFile, prop);
  Close(OutFile);
end;

procedure TIOterm.WriteAtmFrg(AtomFrg: TObjectList);
var
  Outfile: TextFile;
  tmplst: TList;
  i,j,kk: integer;
begin
  AssignFile(OutFile, SGOFileName + '.atm');
  if (FileExists(SGOFileName + '.atm')) then
    Append(OutFile)
  else
    ReWrite(OutFile);
  for i:=1 to AtomFrg.Count-1 do
  begin
    tmplst:=AtomFrg.Items[i] as TList;
    if (tmplst.Count>0) or (i=(AtomFrg.Count-1)) then write(OutFile,IntToStr(i)+':');
    for j:=0 to tmplst.Count-1 do begin
      kk:=PRAtmFrg(tmplst.Items[j])^.idx;
      write(OutFile,IntToStr(kk+1));
      if j<tmplst.Count-1 then
        write(OutFile,',')
      else
        write(OutFile,' ');
    end;
  end;
  writeln(OutFile);
  Close(OutFile);
end;

procedure TIOterm.Integrity();
var
  pitem: PDefFrg;
  i:     integer;
begin
  if ((SGIFileName = '') or (FrgType.Count < 0)) then
  begin
    writeln('Missing input file name or fragmentation type');
    writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    raise TIOtermException.Create(UseMsg);
  end;
  for i := 0 to FrgType.Count - 1 do
  begin
    pitem := FrgType[i];
    //Guarantee that global options are transfered to each Fragment machine
    if fbGetA then pitem^.GetAtomFrg:=True;
    if fbStrict then pitem^.StrictFrg:=True;
    //Check for a Fragment machine that all requirements are valid
    if ((pitem^.FrgType < 0) or (pitem^.FrgType > 10)) then
    begin
      writeln('Wrong fragmentation type');
      writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      raise TIOtermException.Create(UseMsg);
    end;
    if ((pitem^.StrictFrg) and (SGHFileName = '')) then
    begin
      writeln('Missing header file for strict fragmentation');
      writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      raise TIOtermException.Create(UseMsg);
    end;
  end;
  if ((SGHFileName <> '') and not (FileExists(SGHFileName))) then
  begin
    writeln('Impossible to open header file: ' + SGHFileName);
    writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    raise TIOtermException.Create(UseMsg);
  end;
end;

end.

