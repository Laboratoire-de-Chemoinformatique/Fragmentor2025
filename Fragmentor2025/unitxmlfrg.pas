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
unit UnitXMLFrg;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, DOM, XMLRead, XMLWrite, TypePredictor, UnitFragment;

type

{EnumTopFrg = (E,I,II,III, UNKTOPFRG);
EnumColAtm = (A,Ph,Ep,Ba,NoCA,UNKCOLATM);
EnumColBnd = (B,NoCB,UNKCOLBND);
EnumExtTopFrg= (FX,NoFX,UNKFX);
EnumOptFrg = (P, R, AP, FC, MA, MP, SF, AD, OD, NoOPT, UNKOPTFRG)}


{ TXMLFrg }

TXMLFrg = class(TObject)
  private
    fdoc: TXMLDocument;
    froot: TDOMElement;
    ffrgnode: TDOMElement;
    fRW: boolean; //True Write, False Read
    fFileName: string;
    flof: string;
    fbStrictFrg: boolean;
    //foutstr: TStringList;
    fowndoc: boolean;
    fstop: boolean;
    fUseHDR: string;
  public
    constructor Create(s, sp: string; b: Boolean);
    constructor Create(s: string; b: Boolean);
    constructor Create;
    destructor Destroy; override;
    property RW: Boolean read fRW;
    property FileName: string read fFileName;
    property lof: string read flof;
    property bStrictFrg: boolean read fbStrictFrg;
    property UseHDR: string read fUseHDR;
    //property outstr: TStringList read foutstr;
    property stop: boolean read fstop write fstop;
    Procedure SetDocument(rn: TDOMNode; b:boolean);
    Procedure SetDocument(rn: TDOMNode; b:boolean; mixtureElement:integer);
    procedure Finalize;
    procedure WriteFragmentations();
    procedure ReadFragmentations();
    procedure ReadHeader();
    procedure WriteHeader(aHdr: string; abStrictFrg: boolean);
    procedure ReadFragmentation(pf: PDefFrg);
    procedure WriteFragmentation(pf: PDefFrg);
    procedure NextFragmentation();
    procedure XmlReadError(w1,w2: string);
    procedure WriteSDFnames(SDFnames: TStringList; const sprp: string);
    procedure WriteSDFnames(SDFnames: string; const sprp: string);
    procedure ReadSDFnames(SDFnames:TStringList; out sprp: string);
    procedure SetUseHDR(s: string);
  end;

implementation

{ TXMLFrg }

constructor TXMLFrg.Create(s, sp: string; b: Boolean);
begin
  Create;
  fowndoc:=True;
  fRW:=b;
  fFileName:=s; //XML name
  flof:=ChangeFileExt(sp,'.hdr');  //base name for the hdr file
  fstop:=False;
  if fRW then
  begin
     ReadXMLFile(fdoc, fFileName);
  end else
    begin
      fdoc:=TXMLDocument.Create;
    end;
end;

{constructor TXMLFrg.Create(s, sp: string; b: Boolean);
begin
  Create;
  fowndoc:=True;
  fRW:=b;
  fFileName:=s;
  flof:=sp;
//  flof:=ChangeFileExt(fFileName,'.hdr');
  fstop:=False;
  if fRW then
  begin
     ReadXMLFile(fdoc, fFileName);
     froot:=fdoc.FindNode('Fragmentations') as TDOMElement;
     if (froot=nil) then XmlReadError('Fragmentations','TXMLFrg.Create');
     ffrgnode:=froot.FindNode('Fragmentation') as TDOMElement;
     if (ffrgnode=nil) then XmlReadError('Fragmentation','TXMLFrg.Create');
     if froot.Attributes.GetNamedItem('UseHDR')<>nil then//<--TO CORRECT
        fUseHDR:=froot.Attributes.GetNamedItem('UseHDR').NodeValue;
  end else
    begin
      fdoc:=TXMLDocument.Create;
      froot:=fdoc.CreateElement('Fragmentations');
      fdoc.AppendChild(froot);
      WriteHeader(flof,False);
    end;
end;}

constructor TXMLFrg.Create(s: string; b: Boolean);
var
  sp: string;
begin
  sp:=ChangeFileExt(s,'.hdr');//Guess a name for the header file. Updated later if needed
  Create(s,sp,b);
end;

constructor TXMLFrg.Create;
begin
  inherited Create;
  fdoc:=nil;
  froot:=nil;
  ffrgnode:=nil;
  fRW:=True;
  fFileName:='';
  flof:='';
  fbStrictFrg:=False;
  fowndoc:=True;
  fstop:=True;
  fUseHDR:='';
end;

procedure TXMLFrg.Finalize;
begin
  if fRW then
  else
    WriteXMLFile(fdoc,fFileName);
  if fowndoc=True then FreeAndNil(fdoc);
end;

procedure TXMLFrg.WriteFragmentations();
//Write the node Fragmentations and subnode Header
begin
  froot:=fdoc.CreateElement('Fragmentations');
  fdoc.AppendChild(froot);
  WriteHeader(flof,False);
end;

procedure TXMLFrg.ReadFragmentations();
//Read the node Fragmentations and subnode Header
begin
  froot:=fdoc.FindNode('Fragmentations') as TDOMElement;
  if (froot=nil) then XmlReadError('Fragmentations','TXMLFrg.Create');
  ReadHeader();
  ffrgnode:=froot.FindNode('Fragmentation') as TDOMElement;
  if (ffrgnode=nil) then XmlReadError('Fragmentation','TXMLFrg.Create');
  if froot.Attributes.GetNamedItem('UseHDR')<>nil then//<--TO CORRECT
     fUseHDR:=froot.Attributes.GetNamedItem('UseHDR').NodeValue;
end;

procedure TXMLFrg.ReadHeader();
var
  hdrnode: TDOMElement;
  sStrictFrg: string;
begin
  hdrnode:=froot.FindNode('Header') as TDOMElement;
  if (hdrnode.FirstChild<>nil) then begin
     flof:=string(hdrnode.FirstChild.NodeValue);
     sStrictFrg:=string(hdrnode.GetAttribute('StrictFrg'));
     if (sStrictFrg='True') then fbStrictFrg:=True else fbStrictFrg:=False;
  end
  else XmlReadError('Header','TXMLFrg.ReadHeader');
end;

procedure TXMLFrg.WriteHeader(aHdr: string; abStrictFrg: boolean);
var
  hdrnode: TDOMElement;
  hdrtxtnode: TDOMText;
begin
  hdrnode:=fdoc.CreateElement('Header');
  //hdrtxtnode:=fdoc.CreateTextNode(aHdr);
  hdrtxtnode:=fdoc.CreateTextNode(ExtractFileName(aHdr)); //Only write the relative path - modified 05-08-2022
  if abStrictFrg then hdrnode.SetAttribute('StrictFrg','True')
  else hdrnode.SetAttribute('StrictFrg','False');
  hdrnode.AppendChild(hdrtxtnode);
  froot.AppendChild(hdrnode);
end;

destructor TXMLFrg.Destroy;
begin
  if (fdoc<>nil) then Finalize;
  inherited Destroy;
end;

procedure TXMLFrg.SetDocument(rn: TDOMNode; b: boolean);
//Set the Fragmentation XML parser from a <Fragmentations> node
//b: True, read mode; False, write mode
//lof: the header file.
//bStrict: strict fragmentation True/False
begin
  fowndoc:=False;
  fRW:=b;
  fdoc:=rn.OwnerDocument as TXMLDocument;
  froot:=rn as TDOMElement;
  //write('$$$ ');write(froot.Attributes.Item[0].NodeValue);writeln(' $$$');
  fFileName:=fdoc.NamespaceURI;
  flof:=ChangeFileExt(fFileName,'hdr');
  fstop:=False;
  if fRW then
  begin
     if (froot.NodeName <> 'Fragmentations') then XmlReadError('Fragmentations','TXMLFrg.Create');
     ffrgnode:=froot.FindNode('Fragmentation') as TDOMElement;
     if (ffrgnode=nil) then XmlReadError('Fragmentation','TXMLFrg.Create');
  end else
    begin
      froot:=fdoc.CreateElement('Fragmentations');
      fdoc.AppendChild(froot);
    end;
    //write('$$$$$$ ');write(ffrgnode.Attributes.GetNamedItem('ID').NodeValue);writeln(' $$$$$$');
end;

procedure TXMLFrg.SetDocument(rn: TDOMNode; b: boolean; mixtureElement: integer);
//Set the Fragmentation XML parser from a <Fragmentations> node
//b: True, read mode; False, write mode
//lof: the header file.
begin
  fowndoc:=False;
  fRW:=b;
  fdoc:=rn.OwnerDocument as TXMLDocument;
  froot:=rn as TDOMElement;
  //write('$$$ ');write(froot.Attributes.Item[0].NodeValue);writeln(' $$$');
  fFileName:=fdoc.NamespaceURI;
  flof:=ChangeFileExt(fFileName,'hdr');
  fstop:=False;
  if fRW then
  begin
     if (froot.NodeName <> 'Fragmentations_'+IntToStr(mixtureElement)) then XmlReadError('Fragmentations_'+IntToStr(mixtureElement),'TXMLFrg.Create');
     ffrgnode:=froot.FindNode('Fragmentation') as TDOMElement;
     if (ffrgnode=nil) then XmlReadError('Fragmentation','TXMLFrg.Create');
  end else
    begin
      froot:=fdoc.CreateElement('Fragmentations_'+IntToStr(mixtureElement));
      fdoc.AppendChild(froot);
    end;
    //write('$$$$$$ ');write(ffrgnode.Attributes.GetNamedItem('ID').NodeValue);writeln(' $$$$$$');
end;

procedure TXMLFrg.ReadFragmentation(pf: PDefFrg);
var
  t: byte;
  Etop: EnumTopFrg;
  EColA: EnumColAtm;
  EColB: EnumColBnd;
  EFX: EnumExtTopFrg;
  //ECnt: EnumCounting;
  EOpt: EnumOptFrg;
  LoFNode: TDOMText;
  stmp: string;
begin
  pf^.FrgType:=255;
  pf^.lmin:=0;
  pf^.lmax:=0;
  pf^.MarkAtom:=0;
  pf^.DynBnd:=0;
  pf^.ColorAFields:='Default';
  pf^.ColorBFields:='Default';
  pf^.UseFormalCharge:=False;
  pf^.AtomPairs:=False;
  pf^.DoAllWays:=False;
  pf^.EqFuzAt:=False;
  pf^.UseRadical:=False; // ?
  pf^.UseIsotope:=False; // ?
  if bStrictFrg then pf^.StrictFrg:=True else pf^.StrictFrg:=False;
  pf^.GetAtomFrg:=False;
  pf^.Cycle:=False;
  //Read byte nodes
  if (ffrgnode.Attributes.GetNamedItem('Min') = nil) or
     (ffrgnode.Attributes.GetNamedItem('Max') = nil) or
     (ffrgnode.Attributes.GetNamedItem('Topology') = nil) or
     (ffrgnode.Attributes.GetNamedItem('Atom') = nil) or
     (ffrgnode.Attributes.GetNamedItem('Bond') = nil) then begin
    XmlReadError(ffrgnode.TagName+' Missing keyword','TXMLFrg.ReadFragmentation');
  end;
  pf^.lmin := byte(StrToInt(ffrgnode.Attributes.GetNamedItem('Min').NodeValue));
  pf^.lmax := byte(StrToInt(ffrgnode.Attributes.GetNamedItem('Max').NodeValue));
  Etop     := StrToEnumTopFrg(ffrgnode.Attributes.GetNamedItem('Topology').NodeValue);
  EColA    := StrToEnumColAtm(ffrgnode.Attributes.GetNamedItem('Atom').NodeValue);
  EColB    := StrToEnumColBnd(ffrgnode.Attributes.GetNamedItem('Bond').NodeValue);
  if (ffrgnode.Attributes.GetNamedItem('Extended')<>nil) then
    EFX      := StrToEnumExtTopFrg(ffrgnode.Attributes.GetNamedItem('Extended').NodeValue);
  //Interpret FrgType
  if (Etop=E)then pf^.FrgType:=0
  else if (Etop=I) and (EColA=A) and (EColB=NoCB) then pf^.FrgType:=1
  else if (Etop=I) and (EColA=NoCA) and (EColB=B) then pf^.FrgType:=2
  else if (Etop=I) and (EColA=A) and (EColB=B) then pf^.FrgType:=3
  else if (Etop=II) and (EColA=A) and (EColB=NoCB) and (EFX=NoFX) then pf^.FrgType:=4
  else if (Etop=II) and (EColA=NoCA) and (EColB=B) and (EFX=NoFX) then pf^.FrgType:=5
  else if (Etop=II) and (EColA=A) and (EColB=B) and (EFX=NoFX) then pf^.FrgType:=6
  else if (Etop=II) and (EColA=A) and (EColB=NoCB) and (EFX=FX) then pf^.FrgType:=7
  else if (Etop=II) and (EColA=NoCA) and (EColB=B) and (EFX=FX) then pf^.FrgType:=8
  else if (Etop=II) and (EColA=A) and (EColB=B) and (EFX=FX) then pf^.FrgType:=9
  else if (Etop=III) then pf^.FrgType:=10;
  //
  if (ffrgnode.Attributes.GetNamedItem('MarkAtom')<>nil) then
    pf^.MarkAtom:=byte(StrToInt(ffrgnode.Attributes.GetNamedItem('MarkAtom').NodeValue));
  if (ffrgnode.Attributes.GetNamedItem('DynBnd')<>nil) then
    pf^.DynBnd:=byte(StrToInt(ffrgnode.Attributes.GetNamedItem('DynBnd').NodeValue));
  //Read string nodes
  if (ffrgnode.Attributes.GetNamedItem('ColorAField')<>nil) then
    pf^.ColorAFields:=ffrgnode.Attributes.GetNamedItem('ColorAField').NodeValue;
  if (ffrgnode.Attributes.GetNamedItem('ColorBField')<>nil) then
    pf^.ColorBFields:=ffrgnode.Attributes.GetNamedItem('ColorBField').NodeValue;
  if pf^.ColorAFields='' then pf^.ColorAFields:='Default';
  if pf^.ColorBFields='' then pf^.ColorBFields:='Default';
  //Read boolean nodes
  if (ffrgnode.Attributes.GetNamedItem('DoAllWays')<>nil) then
    pf^.DoAllWays:=StrToBool(ffrgnode.Attributes.GetNamedItem('DoAllWays').NodeValue);
  if (ffrgnode.Attributes.GetNamedItem('AtomPairs')<>nil) then
    pf^.AtomPairs:=StrToBool(ffrgnode.Attributes.GetNamedItem('AtomPairs').NodeValue);
  if (ffrgnode.Attributes.GetNamedItem('FuzzyAtom')<>nil) then
    pf^.EqFuzAt:=StrToBool(ffrgnode.Attributes.GetNamedItem('FuzzyAtom').NodeValue);
  if (ffrgnode.Attributes.GetNamedItem('FormalCharge')<>nil) then
    pf^.UseFormalCharge:=StrToBool(ffrgnode.Attributes.GetNamedItem('FormalCharge').NodeValue);
  if (ffrgnode.Attributes.GetNamedItem('Isotope')<>nil) then
    pf^.UseIsotope:=StrToBool(ffrgnode.Attributes.GetNamedItem('Isotope').NodeValue);
  if (ffrgnode.Attributes.GetNamedItem('Radical')<>nil) then
    pf^.UseRadical:=StrToBool(ffrgnode.Attributes.GetNamedItem('Radical').NodeValue);
  if (ffrgnode.Attributes.GetNamedItem('Cycle')<>nil) then
    pf^.Cycle:=StrToBool(ffrgnode.Attributes.GetNamedItem('Cycle').NodeValue);
  if (ffrgnode.Attributes.GetNamedItem('StrictFrg')<>nil) then
    pf^.StrictFrg:=StrToBool(ffrgnode.Attributes.GetNamedItem('StrictFrg').NodeValue);
  if (ffrgnode.Attributes.GetNamedItem('AtomFrg')<>nil) then
    pf^.GetAtomFrg:=StrToBool(ffrgnode.Attributes.GetNamedItem('AtomFrg').NodeValue);
  //Read header file name
  LoFNode:=ffrgnode.FirstChild as TDOMText;
  if LoFNode<>nil then flof:=ffrgnode.FirstChild.NodeValue;
end;

procedure TXMLFrg.WriteFragmentation(pf: PDefFrg);
var
  t: byte;
  Etop: EnumTopFrg;
  EColA: EnumColAtm;
  EColB: EnumColBnd;
  EFX: EnumExtTopFrg;
  //ECnt: EnumCounting;
  EOpt: EnumOptFrg;
  LoFNode: TDOMText;
begin
  Etop:=UNKTOPFRG;
  EColA:=NoCA;
  EColB:=NoCB;
  //ECnt:=NoCNT;
  EFX:=NoFX;
  EOpt:=UNKOPTFRG;
  //
  ffrgnode:=fdoc.CreateElement('Fragmentation');
  froot.AppendChild(ffrgnode);
  ffrgnode.SetAttribute('Min',IntToStr(pf^.lmin));
  ffrgnode.SetAttribute('Max',IntToStr(pf^.lmax));
  t:=pf^.FrgType;
  case t of
    0: Etop:=E;
    1: begin
         Etop:=I;
         EColA:=A;
         EColB:=NoCB;
       end;
    2: begin
         Etop:=I;
         EColA:=NoCA;
         EColB:=B;
       end;
    3: begin
         Etop:=I;
         EColA:=A;
         EColB:=B;
       end;
    4: begin
         Etop:=II;
         EColA:=A;
         EColB:=NoCB;
         EFX:=NoFX;
       end;
    5: begin
         Etop:=II;
         EColA:=NoCA;
         EColB:=B;
         EFX:=NoFX;
       end;
    6: begin
         Etop:=II;
         EColA:=A;
         EColB:=B;
         EFX:=NoFX;
       end;
    7: begin
         Etop:=II;
         EColA:=A;
         EColB:=NoCB;
         EFX:=FX;
       end;
    8: begin
         Etop:=II;
         EColA:=NoCA;
         EColB:=B;
         EFX:=FX;
       end;
    9: begin
         Etop:=II;
         EColA:=A;
         EColB:=B;
         EFX:=FX;
       end;
    10: begin
          Etop:=III;
        end;
  end;
  ffrgnode.SetAttribute('Topology',EnumTopFrgToStr(Etop));
  ffrgnode.SetAttribute('Atom',EnumColAtmToStr(EColA));
  ffrgnode.SetAttribute('Bond',EnumColBndToStr(EColB));
  ffrgnode.SetAttribute('Extended',EnumEnumExtTopFrg(EFX));
  ffrgnode.SetAttribute('MarkAtom',IntToStr(pf^.MarkAtom));
  ffrgnode.SetAttribute('DynBnd',IntToStr(pf^.DynBnd));
  ffrgnode.SetAttribute('ColorAField',pf^.ColorAFields);
  ffrgnode.SetAttribute('ColorBField',pf^.ColorBFields);
  if (pf^.DoAllWays) then ffrgnode.SetAttribute('DoAllWays','True')
  else ffrgnode.SetAttribute('DoAllWays','False');
  if (pf^.AtomPairs) then ffrgnode.SetAttribute('AtomPairs','True')
  else ffrgnode.SetAttribute('AtomPairs','False');
  if (pf^.EqFuzAt) then ffrgnode.SetAttribute('FuzzyAtom','True')
  else ffrgnode.SetAttribute('FuzzyAtom','False');
  if (pf^.UseFormalCharge)then ffrgnode.SetAttribute('FormalCharge','True')
  else ffrgnode.SetAttribute('FormalCharge','False');
  if (pf^.UseIsotope)then ffrgnode.SetAttribute('Isotope','True')
  else ffrgnode.SetAttribute('Isotope','False');
  if (pf^.UseRadical)then ffrgnode.SetAttribute('Radical','True')
  else ffrgnode.SetAttribute('Radical','False');
  if (pf^.Cycle) then ffrgnode.SetAttribute('Cycle','True')
  else ffrgnode.SetAttribute('Cycle','False');
  //if (pf^.StrictFrg)then ffrgnode.SetAttribute('StrictFrg','True')
  //else ffrgnode.SetAttribute('StrictFrg','False');
  if (pf^.GetAtomFrg)then ffrgnode.SetAttribute('AtomFrg','True')
  else ffrgnode.SetAttribute('AtomFrg','False');
  //Write header file name
  //LoFNode:=fdoc.CreateTextNode(flof);
  //ffrgnode.AppendChild(LoFNode);
  //lof (list of fragment) and StrictFrg are now global and managed by header
end;

procedure TXMLFrg.NextFragmentation();
var
  tmpNode: TDOMElement;
begin
  tmpNode:=ffrgnode.NextSibling as TDOMElement;
  if (tmpNode<>nil) and (tmpNode.NodeName='Fragmentation') then
    ffrgnode:=tmpNode
  else
    fstop:=True;
end;

procedure TXMLFrg.XmlReadError(w1, w2: string);
begin
  Writeln('XML read error at node: '+w1);
  raise Exception.Create('*** ERROR: '+w2+' ***');
end;

procedure TXMLFrg.WriteSDFnames(SDFnames: TStringList; const sprp: string);
var
  //tmpnode: TDOMElement;
  LoSDF: TDOMText;
  NdSDF: TDOMElement;
  i: integer;
begin
  //Save a list of SDF files. Useful for GUI.
  ffrgnode:=fdoc.CreateElement('SDFiles');
  if sprp<>'' then ffrgnode.SetAttribute('TargetField',sprp);
  froot.AppendChild(ffrgnode);
  for i:=0 to SDFnames.Count-1 do
  begin
    NdSDF:=fdoc.CreateElement('SDFile');
    ffrgnode.AppendChild(NdSDF);
    LoSDF:=fdoc.CreateTextNode(SDFnames[i]);
    NdSDF.AppendChild(LoSDF);
  end;
end;

procedure TXMLFrg.WriteSDFnames(SDFnames: string; const sprp: string);
var
  tmpSL: TStringList;
begin
  tmpSL:=TStringList.Create;
  tmpSL.Add(SDFnames);
  WriteSDFnames(tmpSL,sprp);
  FreeAndNil(tmpSL);
end;

procedure TXMLFrg.ReadSDFnames(SDFnames: TStringList; out sprp: string);
var
  s: string;
  tmpnode: TDOMNode;
begin
  //Read a list of SDF files. Useful for GUI
  SDFnames.Clear;
  ffrgnode:=froot.FindNode('SDFiles') as TDOMElement;
  if (ffrgnode.Attributes.GetNamedItem('TargetField') <> nil) then
     sprp:=ffrgnode.Attributes.GetNamedItem('TargetField').NodeValue
  else sprp:='';
  if (ffrgnode<>nil) then
  begin
    tmpnode:=ffrgnode.FirstChild;
    while (tmpnode<>nil) and (tmpnode.NodeName='SDFile') do
    begin
      s:=tmpnode.FirstChild.NodeValue;
      SDFnames.Add(s);
      tmpnode:=tmpnode.NextSibling;
    end;
  end;
end;

procedure TXMLFrg.SetUseHDR(s: string);
//Set the attribute UseHDR (meaningfull mainly for the xFragmentor)
begin
  froot.SetAttribute('UseHDR',s);
end;

end.

