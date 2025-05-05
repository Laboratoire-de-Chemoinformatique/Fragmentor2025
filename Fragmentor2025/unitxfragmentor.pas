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
unit unitxfragmentor;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, strutils, contnrs, FileUtil, LazFileUtils, Forms, Controls,
  Graphics, Dialogs, StdCtrls, Grids, FileCtrl, CheckLst, ExtCtrls, MaskEdit,
  Spin, Menus, UnitFragment, UnitFragmentBase, BasicSDF, UnitXMLFrg;

type

  { TFrmxFragmentor }

  TFrmxFragmentor = class(TForm)
    BtCheckHeader: TButton;
    BtCheckXML: TButton;
    BtSDFAdd: TButton;
    BtSDFDel: TButton;
    BtAddFrg: TButton;
    BtDeleteFrg: TButton;
    BtSaveXML: TButton;
    BtReadXML: TButton;
    BtReset: TButton;
    BtQuit: TButton;
    BtRun: TButton;
    CBStrictFrg: TCheckBox;
    CBTopo: TComboBox;
    CBAColor: TComboBox;
    CBBColor: TComboBox;
    CBActivateGBHeader: TCheckBox;
    CBMarkedAtom: TComboBox;
    EdHeader: TEdit;
    EdXML: TEdit;
    EdSDFCBfield: TEdit;
    EdSDFfield: TEdit;
    EdSDFCAfield: TEdit;
    GBxHeader: TGroupBox;
    GBxXML: TGroupBox;
    LbSDFCBField: TLabel;
    LbSDFfield: TLabel;
    LbSDFCAField: TLabel;
    Lbmax: TLabel;
    Lbmin: TLabel;
    LbBColor: TLabel;
    LbAColor: TLabel;
    LbTopo: TLabel;
    LbSDFList: TLabel;
    LBSDF: TListBox;
    LBFrgLst: TListBox;
    MemoLog: TMemo;
    ODSDFOpen: TOpenDialog;
    CBAllWays: TCheckBox;
    CBAtomPairs: TCheckBox;
    CBCharge: TCheckBox;
    CBGetFrag: TCheckBox;
    ODXMLRead: TOpenDialog;
    ODXMLSave: TSaveDialog;
    SBSLSDF: TScrollBox;
    SEmin: TSpinEdit;
    SEmax: TSpinEdit;
    procedure BtCheckHeaderClick(Sender: TObject);
    procedure BtDeleteFrgClick(Sender: TObject);
    procedure BtQuitClick(Sender: TObject);
    procedure BtReadXMLClick(Sender: TObject);
    procedure BtResetClick(Sender: TObject);
    procedure BtRunClick(Sender: TObject);
    procedure BtSaveXMLClick(Sender: TObject);
    procedure BtSDFAddClick(Sender: TObject);
    procedure BtSDFDelClick(Sender: TObject);
    procedure BtAddFrgClick(Sender: TObject);
    procedure CBAColorChange(Sender: TObject);
    procedure CBActivateGBHeaderChange(Sender: TObject);
    procedure CBBColorChange(Sender: TObject);
    procedure CBStrictFrgChange(Sender: TObject);
    procedure CBTopoChange(Sender: TObject);
  private
    { private declarations }
    IFileName, OFileName, HFileName, XFileName, WDFileName: string;
    //Input, Output, Header, XML and Workdir
    procedure RefreshStrictFrg(bStrictFrg: boolean);
    procedure LBSDFAddAndResize(flenme: string);
  public
    { public declarations }
  end;

procedure WriteAtmFrg(AtomFrg: TObjectList; atmfile: string);
procedure WriteSVM(prop: string; frglst: TStringList; svmfile: string);
procedure WriteARFFData(prop: string; frglst: TStringList; arfffile: string);
procedure WriteHeader(frglst: TStringList; hdrfile: string);
function CheckSDFField(filenme, field: string): boolean;
function DefFrgToString(aFrg: DefFrg): string;

var
  FrmxFragmentor: TFrmxFragmentor;

implementation

procedure WriteAtmFrg(AtomFrg: TObjectList; atmfile: string);
var
  Outfile: TextFile;
  tmplst: TList;
  i, j, kk: integer;
begin
  AssignFile(OutFile, atmfile);
  if FileExists(atmfile) then
    Append(OutFile)
  else
    ReWrite(OutFile);
  for i := 1 to AtomFrg.Count - 1 do
  begin
    tmplst := AtomFrg.Items[i] as TList;
    Write(OutFile, IntToStr(i) + ':');
    for j := 0 to tmplst.Count - 1 do
    begin
      kk := PRAtmFrg(tmplst.Items[j])^.idx;
      Write(OutFile, IntToStr(kk + 1));
      if j < tmplst.Count - 1 then
        Write(OutFile, ',')
      else
        Write(OutFile, ' ');
    end;
  end;
  writeln(OutFile);
  Close(OutFile);
end;

procedure WriteSVM(prop: string; frglst: TStringList; svmfile: string);
var
  Outfile: TextFile;
  i, j: integer;
begin
  AssignFile(OutFile, svmfile);
  if (FileExists(svmfile)) then
    Append(OutFile)
  else
    ReWrite(OutFile);
  if(prop='') then prop:='?';
  Write(OutFile, prop);
  for i := 0 to frglst.Count - 1 do
  begin
    j := PInteger(frglst.Objects[i])^;
    if (j <> 0) then
    begin
      Write(OutFile, ' ' + IntToStr(i + 1) + ':' + IntToStr(j));
    end;
  end;
  writeln(OutFile);
  Close(Outfile);
end;

procedure WriteARFFData(prop: string; frglst: TStringList; arfffile: string);
var
  Outfile: TextFile;
  i, j: integer;
begin
  AssignFile(OutFile, arfffile);
  if (FileExists(arfffile)) then
    Append(OutFile)
  else
    ReWrite(OutFile);
  Write(OutFile,'{');
  for i := 0 to frglst.Count - 1 do
  begin
    j := PInteger(frglst.Objects[i])^;
    if (j <> 0) then
    begin
      Write(OutFile, IntToStr(i) + ' ' + IntToStr(j)+', ');//attributes are 0 based
    end;
  end;
  if(prop='') then prop:='?';
  Write(OutFile,'last '+prop+'}');
  writeln(OutFile);
  Close(Outfile);
end;

procedure WriteHeader(frglst: TStringList; hdrfile: string);
var
  Outfile: Text;
  i: integer;
begin
  AssignFile(OutFile, hdrfile);
  ReWrite(OutFile);
  for i := 0 to frglst.Count - 1 do
  begin
    writeln(Outfile, format('%6s', [IntToStr(i + 1)]) + '. ' +
      format('%120s', [frglst.Strings[i]]));
  end;
  Close(Outfile);
end;

function CheckSDFField(filenme, field: string): boolean;
var
  sdf: SDFmini;
  smol: TStringList;
begin
  smol := nil;
  sdf := SDFmini.Create(filenme);
  sdf.prptag := field;
  sdf.NextTMol;
  if sdf.prpval = '?' then
    Result := False
  else
    Result := True;
  FreeAndNil(smol);
  FreeAndNil(sdf);
end;

function DefFrgToString(aFrg: DefFrg): string;
var
  sFrg: string;
begin
  sFrg := '';
  case aFrg.FrgType of
    0: sFrg := 'E';
    1: sFrg := 'IA';
    2: sFrg := 'IB';
    3: sFrg := 'IAB';
    4: sFrg := 'IIA';
    5: sFrg := 'IIB';
    6: sFrg := 'IIAB';
    7: sFrg := 'IIA_R';
    8: sFrg := 'IIB_R';
    9: sFrg := 'IIAB_R';
    10: sFrg := 'III';
    else
      sFrg := 'UNK';
  end;
  if aFrg.FrgType <> 0 then
    sFrg := sFrg + '(' + IntToStr(aFrg.lmin) + '-' + IntToStr(aFrg.lmax) + ')';
  if aFrg.DoAllWays then
    sFrg := sFrg + '_AP';
  if aFrg.AtomPairs then
    sFrg := sFrg + '_P';
  if aFrg.UseFormalCharge then
    sFrg := sFrg + '_FC';
  if aFrg.MarkAtom = 1 then
    sFrg := sFrg + '_MA1';
  if aFrg.MarkAtom = 2 then
    sFrg := sFrg + '_MA2';
  if aFrg.MarkAtom = 3 then
    sFrg := sFrg + '_MA3';
  if aFrg.ColorAFields <> 'Default' then
    sFrg := sFrg + '_CuA(' + aFrg.ColorAFields + ')';
  if aFrg.ColorAFields = '' then
     ShowMessage('BUG: atom coloring not defined. How could that be possible?');
  if aFrg.ColorBFields <> 'Default' then
    sFrg := sFrg + '_CuB(' + aFrg.ColorBFields + ')';
  if aFrg.ColorBFields = '' then
     ShowMessage('BUG: bond coloring not defined. How could that be possible?');
  Result := sFrg;
end;

{$R *.lfm}

{ TFrmxFragmentor }

procedure TFrmxFragmentor.BtSDFAddClick(Sender: TObject);
var
  isze: integer;
  counter: integer;
begin
  //Modified 10-04-2018
  if ODSDFOpen.Execute then
    IFileName := ODSDFOpen.FileName;
    for counter:=0 to ODSDFOpen.Files.Count-1 do
       begin
          IFileName:= ODSDFOpen.Files.Strings[counter];
          LBSDFAddAndResize(IFileName);
       end;
  //LBSDFAddAndResize(IFileName);
end;

procedure TFrmxFragmentor.BtDeleteFrgClick(Sender: TObject);
var
  pFrg: PDefFrg;
  i: integer;
begin
  i := LBFrgLst.ItemIndex;
  if i >= 0 then
  begin
    pFrg := PDefFrg(LBFrgLst.Items.Objects[i]);
    Dispose(pFrg);
    pFrg := nil;
    LBFrgLst.Items.Delete(i);
  end;
end;

procedure TFrmxFragmentor.BtQuitClick(Sender: TObject);
begin
  Close;
end;

procedure TFrmxFragmentor.BtReadXMLClick(Sender: TObject);
var
  xml: TXMLFrg;
  pFrg: PDefFrg;
  sxml, sFrg, sPrp,hBase: string;
  SDFList: TStringList;
  i: integer;
begin
  //Set the LBFrgList using data stored into an XML file
  for i:=LBFrgLst.Count-1 downto 0 do
  begin
    pFrg:=PDefFrg(LBFrgLst.Items.Objects[i]);
    dispose(pFrg);
    LBFrgLst.Items.Delete(i);
  end;
  for i:=LBSDF.Count-1 downto 0 do LBSDF.Items.Delete(i);
  //
  sxml := '';
  if ODXMLRead.Execute then
  begin
    sxml := ODXMLRead.FileName;
  end;
  if FileExists(sxml) then
  begin
    xml := TXMLFrg.Create(sxml, True);
    if xml.UseHDR<>'' then CBActivateGBHeader.Checked:=True;
    while (not xml.stop) do
    begin
      pFrg := new(PDefFrg);
      xml.ReadFragmentation(pFrg);
      sFrg := DefFrgToString(pFrg^);
      LBFrgLst.AddItem(sFrg, TObject(pFrg));
      if pFrg^.GetAtomFrg then CBGetFrag.Checked:=True;
      if pFrg^.StrictFrg then CBStrictFrg.Checked:=True;
      if xml.UseHDR<>'' then EdHeader.Text:=xml.UseHDR;
      xml.NextFragmentation;
    end;
    SDFList:=TStringList.Create;
    xml.ReadSDFnames(SDFList,sPrp);
    EdSDFfield.Text:=sPrp;
    for i:=0 to SDFList.Count-1 do
    begin
      LBSDFAddAndResize(SDFList[i]);
    end;
    FreeAndNil(SDFList);
    FreeAndNil(xml);
  end;
end;

procedure TFrmxFragmentor.BtResetClick(Sender: TObject);
var
  i: integer;
  pFrg: PDefFrg;
begin
  for i := LBSDF.Count - 1 downto 0 do
  begin
    LBSDF.Items.Delete(i);
  end;
  for i := LBFrgLst.Count - 1 downto 0 do
  begin
    pFrg := PDefFrg(LBFrgLst.Items.Objects[i]);
    dispose(pFrg);
    pFrg := nil;
    LBFrgLst.Items.Delete(i);
  end;
  CBAColor.ItemIndex := 1;
  CBBColor.ItemIndex := 1;
  CBBColor.ItemIndex := 0;
  SEmin.Value := 2;
  SEmax.Value := 2;
  EdSDFCAfield.Text := 'SDF field name...';
  EdSDFfield.Text := '';
  CBAtomPairs.Checked := False;
  CBAllWays.Checked := False;
  CBGetFrag.Checked := False;
  CBCharge.Checked := False;
  CBActivateGBHeader.Checked := False;
  CBStrictFrg.Checked := False;
  EdHeader.Text := 'Base name of header files...';
  MemoLog.Clear;
  MemoLog.Append('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
  MemoLog.Append('ISIDA/xFragmentor');
  MemoLog.Append('');
  MemoLog.Append('Université de Strasbourg, 2016');
  MemoLog.Append('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
end;

procedure TFrmxFragmentor.BtRunClick(Sender: TObject);
var
  i, ii, j, jj, k, NFrg: integer;
  testPos: integer;
  dtmp: double;
  stmp, stmp1, stmp2, sdname, hdrname, svmname, arffname, atmname, xmlname, sID, CWD, svalue,slog: string;
  frag: TFragment;
  pFrag: PDefFrg;
  mol: SDFmini;
  hdrlst, tmol, totfraglst, SVMSL, ARFFSL, SData, STuple: TStringList;
  HFile: TextFile;
  xml: TXMLFrg;
begin
  frag := TFragment.Create;
  sID := '';
  hdrname := '';
  svmname := '';
  atmname := '';
  for i := 0 to LBFrgLst.Count - 1 do
  begin
    sId := sID + '_' + LBFrgLst.Items[i];
    pFrag:=PDefFrg(LBFrgLst.Items.Objects[i]);
    if CBGetFrag.Checked then
    begin
      pFrag^.GetAtomFrg:=True;
      frag.GetFrgPerAtom:=True;
    end;
    pFrag^.StrictFrg:=CBStrictFrg.Checked;
    frag.AddFragmentor(pFrag);
  end;
  for i := 0 to LBSDF.Count - 1 do
  begin
    Application.ProcessMessages;
    sdname := LBSDF.Items[i];
    CWD := ExtractFilePath(sdname);
    if CBActivateGBHeader.Checked then
    begin
      if ((EdHeader.Text <> 'Base name of header files...') and (EdHeader.Text <> '')) then //Search a meaningful header name
        hdrname := ChangeFileExt(CWD + EdHeader.Text,'.hdr')
      else begin
        //default behavior:
        //if SDF files are *test*.sdf, search for the corresponding *train*.hdr
        testPos:=Pos('test',ExtractFileName(sdname));
        If (testPos>0) then
        begin
          hdrname:=StringReplace(ExtractFileNameOnly(sdname),'test','train',[]);
          hdrname := CWD + hdrname+sId;
          hdrname:=ChangeFileExt(hdrname,'.hdr');
          If not FileExists(hdrname) then
          begin
            Raise Exception.Create('ERROR: Header file '+hdrname+' does not exists.');
            Halt;
          end;
        end else
        begin
          Raise Exception.Create('ERROR: No valid header file.');
          Halt;
        end;
      end;
    end else hdrname := ChangeFileExt(sdname, sID + '.hdr');
    svmname := ChangeFileExt(sdname, sID + '.svm');
    arffname := ChangeFileExt(sdname, sID + '.arff');
    atmname := ChangeFileExt(sdname, sID + '.atm');
    //Save an XML file describing current fragmentation
    xmlname:=ChangeFileExt(sdname, sID + '.xml');
    if FileExists(xmlname) then
    begin
      MemoLog.Lines.add('WARNING: file ' + xmlname + ' overwritten');
      DeleteFile(xmlname);
    end;
    xml:=TXMLFrg.Create(xmlname,False);
    xml.WriteFragmentations();
    for ii := 0 to LBFrgLst.Items.Count - 1 do
    begin
      pFrag := PDefFrg(LBFrgLst.Items.Objects[ii]);
      xml.WriteFragmentation(pFrag);
    end;
    xml.WriteSDFnames(sdname,EdSDFfield.Text);
    xml.Finalize;
    FreeAndNil(xml);
    //
    mol := SDFmini.Create(sdname);
    mol.prptag := EdSDFfield.Text;
    if EdSDFfield.Text='' then MemoLog.Lines.add('WARNING: No SDF Field of interest has been given.');
    hdrlst := TStringList.Create;
    if FileExists(hdrname) and CBActivateGBHeader.Checked then
    begin
      AssignFile(HFile, hdrname);
      Reset(HFile);
      repeat
        readln(HFile, dtmp, stmp);
        stmp := DelSpace(stmp);
        hdrlst.Add(stmp);
        MemoLog.Lines.Add(stmp);
      until (EOF(HFile));
      frag.SetFrgLst(hdrlst);
      hdrname := ChangeFileExt(sdname, sID + '.hdr'); //Save a new HDR file
    end;
    if FileExists(svmname) or FileExists(arffname) then
    begin
      MemoLog.Lines.add('WARNING: file ' + svmname + ' overwritten');
      MemoLog.Lines.add('WARNING: file ' + arffname + ' overwritten');
      MemoLog.Lines.add('WARNING: file ' + hdrname + ' overwritten');
      DeleteFile(svmname);
      DeleteFile(arffname);
      DeleteFile(hdrname);
    end;
    if FileExists(atmname) then
    begin
      MemoLog.Lines.add('WARNING: file ' + atmname + ' overwritten');
      DeleteFile(atmname);
    end;
    slog:='/';
    k := 0;
    Application.ProcessMessages;
    repeat
      if (k mod 25 = 0) then
      begin
        MemoLog.Lines.Add(slog);
        slog:='';
        Application.ProcessMessages;
      end;
      tmol := mol.NextTMol;
      Inc(k);
      if tmol.Count < 4 then
        continue;
      slog:=slog+'.' + IntToStr(k);
      frag.FrgReset;
      frag.MolToFrgLst(tmol);
      totfraglst := frag.GetFrgLst;
      if (frag.IsAnyNewFrg) then
        slog:=slog+'!';
      if frag.GetFrgPerAtom then
      begin
        WriteAtmFrg(frag.FrgPerAtom, atmname);
      end;
      svalue := mol.prpval;
      svalue := Trim(svalue);
      WriteSVM(svalue, totfraglst, svmname);
    until (mol.fend);
    slog:=slog+'/';
    MemoLog.Lines.add(slog);
    slog:='';
    Application.ProcessMessages;
    //Write the header file
    WriteHeader(totfraglst,hdrname);
    // Finalize SVM file and ARFF conversion
    Nfrg := totfraglst.Count;
    SVMSL := TStringList.Create;
    ARFFSL := TStringList.Create;
    //Head of ARFF
    ARFFSL.Add('@RELATION "'+sdname+'"');
    ARFFSL.Add('');
    for ii:=0 to totfraglst.Count-1 do ARFFSL.Add('@ATTRIBUTE "'+totfraglst[ii]+'" NUMERIC');
    ARFFSL.Add('@ATTRIBUTE class NUMERIC');
    ARFFSL.Add('');
    ARFFSL.Add('@data');
    //
    SData := TStringList.Create;
    SData.Delimiter := ' ';
    Stuple := TStringList.Create;
    Stuple.Delimiter := ':';
    SVMSL.LoadFromFile(svmname);
    for ii := 0 to SVMSL.Count - 1 do
    begin
      stmp := SVMSL[ii];
      SData.DelimitedText := stmp;
      if SData.Count > 1 then
      begin
        Stuple.DelimitedText := SData[SData.Count - 1];
        j := StrToInt(Stuple[0]);
        if (j < Nfrg) then
          stmp := stmp + ' ' + IntToStr(Nfrg) + ':0';
      end;
      SVMSL[ii] := stmp;
      //ARFF output
      stmp:='{';
      for jj:=1 to SData.Count-1 do
      begin
          STuple.DelimitedText:=SData[jj];
          j:=StrToInt(STuple[0])-1;
          stmp:=stmp+IntToStr(j)+' '+STuple[1]+', ';
      end;
      //In ARFF no need to add a 0 for the last descriptor
      stmp:=stmp+IntToStr(Nfrg)+' '+SData[0]+'}';
      ARFFSL.Add(stmp);
    end;
    SVMSL.SaveToFile(svmname);
    ARFFSL.SaveToFile(arffname);
    FreeAndNil(Stuple);
    FreeAndNil(SData);
    FreeAndNil(SVMSL);
    FreeAndNil(ARFFSL);
    //
    FreeAndNil(hdrlst);
    FreeAndNil(mol);
  end;
  MemoLog.Lines.Add('!!!!!!!!! Fragmentation successfull !!!!!!!!!');
  FreeAndNil(frag);
end;

procedure TFrmxFragmentor.BtSaveXMLClick(Sender: TObject);
var
  xml: TXMLFrg;
  i: integer;
  pFrg: PDefFrg;
  sxml: string;
begin
  //Save the data stored into LBFrgLst into an XML file
  xml := nil;
  sxml := '';
  if ODXMLSave.Execute then
    sxml := ODXMLSave.FileName;
  if sxml <> '' then
  begin
    XFileName := sxml;
    xml := TXMLFrg.Create(XFileName, False);//Write mode
    if (CBActivateGBHeader.Checked) and (EdHeader.Text<>'Base name of header files...') then xml.SetUseHDR(EdHeader.Text);
    for i := 0 to LBFrgLst.Items.Count - 1 do
    begin
      pFrg := PDefFrg(LBFrgLst.Items.Objects[i]);
      pFrg^.GetAtomFrg:=CBGetFrag.Checked;
      pFrg^.StrictFrg:=CBStrictFrg.Checked;
      xml.WriteFragmentation(pFrg);
    end;
    xml.WriteSDFnames(TStringList(LBSDF.Items),EdSDFfield.Text);
    xml.Finalize;
    FreeAndNil(xml);
  end;
end;

procedure TFrmxFragmentor.BtCheckHeaderClick(Sender: TObject);
var
  i: integer;
  testPos: integer;
  bFound: boolean;
  basenme, sWD, hnme, sID: string;
begin
  basenme := EdHeader.Text;
  if (basenme <> 'Base name of header files...') and (basenme <> '') then
  begin
    basenme := EdHeader.Text;
    sID:='';
    for i := 0 to LBFrgLst.Count - 1 do sId := sID + '_' + LBFrgLst.Items[i];
    bFound := True;
    i := 0;
    while (i <= LBSDF.Items.Count - 1) and bFound do
    begin
      sWD := ExtractFilePath(LBSDF.Items.Strings[i]);
      hnme := sWD + basenme;
      hnme:=ChangeFileExt(hnme,'.hdr');
      bFound := FileExists(hnme);
      Inc(i);
    end;
    if (not bFound) then
    begin
      ShowMessage('ERROR: Header file not found. Was searching for : ' + hnme);
      EdHeader.Text := 'Base name of header files...';
    end
    else
      ShowMessage('All header files were found.');
  end else
  begin
    //default behavior:
    //if SDF files are *test*.sdf, search for the corresponding *train*.hdr
    sID:='';
    for i := 0 to LBFrgLst.Count - 1 do sId := sID + '_' + LBFrgLst.Items[i];
    bFound := True;
    i := 0;
    while (i <= LBSDF.Items.Count - 1) and bFound do
    begin
      sWD := ExtractFilePath(LBSDF.Items.Strings[i]);
      testPos:=Pos('test',ExtractFileName(LBSDF.Items.Strings[i]));
      If (testPos>0) then
      begin
        basenme:= StringReplace(ExtractFileNameOnly(LBSDF.Items.Strings[i]),'test','train',[]);
        hnme := sWD + basenme+sID;
        hnme:=ChangeFileExt(hnme,'.hdr');
        bFound := FileExists(hnme);
      end;
      Inc(i);
    end;
    if (not bFound) then
    begin
      ShowMessage('ERROR: Header file not found. Was searching for : ' + hnme);
      EdHeader.Text := 'Base name of header files...';
    end
    else
      ShowMessage('All header files were found.');
  end;
end;

procedure TFrmxFragmentor.BtSDFDelClick(Sender: TObject);
var
  i: integer;
begin
  //Modified 10/04/2018
  //i := LBSDF.ItemIndex;
  i := LBSDF.Items.Count-1;
  while i>=0 do begin
    if LBSDF.Selected[i] then
       LBSDF.Items.Delete(i);
    i:=i-1;
  end;
end;

procedure TFrmxFragmentor.BtAddFrgClick(Sender: TObject);
var
  aFrg: DefFrg;
  init_Frg: DefFrg = (FrgType: 0; lmin: 0; lmax: 0; MarkAtom: 0; DynBnd: 0;
  DoAllWays: False; EqFuzAt: False; UseFormalCharge: False; UseRadical: False; UseIsotope : False; AtomPairs: False;
  StrictFrg: False; GetAtomFrg: False; Cycle: False; ColorAFields: 'Default';
  ColorBFields: 'Default';);
  pFrg: PDefFrg;
  i: integer;
  SDFOK, NoError: boolean;
  sFrg: string;
begin
  //initialize values
  aFrg := init_Frg;
  SDFOK := True;
  NoError := True;
  aFrg.AtomPairs := CBAtomPairs.Checked;
  if CBAColor.ItemIndex = 2 then
  begin
    //Check validity of atom color field
    i := 0;
    if LBSDF.Count < 1 then
      SDFOK := False
    else
      SDFOK := True;
    while (i < LBSDF.Count) and (SDFOK = True) do
    begin
      SDFOK := CheckSDFField(LBSDF.Items[i], EdSDFCAfield.Text);
      Inc(i);
    end;
    if SDFOK = False then
    begin
      EdSDFCAfield.Text := 'SDF field name...';
      ShowMessage('ERROR: SDF does not contained the atom coloration sheme.');
      NoError := False;
    end
    else
      aFrg.ColorAFields := EdSDFCAfield.Text;
  end;
  if CBBColor.ItemIndex = 2 then
  begin
    //Check validity of bond color field
    i := 0;
    if LBSDF.Count < 1 then
      SDFOK := False
    else
      SDFOK := True;
    while (i < LBSDF.Count) and (SDFOK = True) do
    begin
      SDFOK := CheckSDFField(LBSDF.Items[i], EdSDFCBfield.Text);
      Inc(i);
    end;
    if SDFOK = False then
    begin
      EdSDFCAfield.Text := 'SDF field name...';
      ShowMessage('ERROR: SDF does not contained the bond coloration sheme.');
      NoError := False;
    end
    else
      aFrg.ColorBFields := EdSDFCBfield.Text;
  end;
  aFrg.DoAllWays := CBAllWays.Checked;
  if CBTopo.ItemIndex = 0 then
    aFrg.FrgType := 0 //t 0, E
  else if (CBTopo.ItemIndex = 1) and ((CBAColor.ItemIndex = 1) or (CBAColor.ItemIndex = 2)) and
    (CBBColor.ItemIndex = 0) then
    aFrg.FrgType := 1 //t 1, IA;
  else if (CBTopo.ItemIndex = 1) and (CBAColor.ItemIndex = 0) and ((CBBColor.ItemIndex = 1) or
    (CBBColor.ItemIndex = 2)) then
    aFrg.FrgType := 2 //t 2, IB;
  else if (CBTopo.ItemIndex = 1) and ((CBAColor.ItemIndex = 1) or (CBAColor.ItemIndex = 2)) and
    ((CBBColor.ItemIndex = 1) or (CBBColor.ItemIndex = 2)) then
    aFrg.FrgType := 3 //t 3, IAB;
  else if (CBTopo.ItemIndex = 2) and ((CBAColor.ItemIndex = 1) or (CBAColor.ItemIndex = 2)) and
    (CBBColor.ItemIndex = 0) then
    aFrg.FrgType := 4 //t 4, IIA;
  else if (CBTopo.ItemIndex = 2) and (CBAColor.ItemIndex = 0) and ((CBBColor.ItemIndex = 1) or
    (CBBColor.ItemIndex = 2)) then
    aFrg.FrgType := 5 //t 5, IIB;
  else if (CBTopo.ItemIndex = 2) and ((CBAColor.ItemIndex = 1) or (CBAColor.ItemIndex = 2)) and
    ((CBBColor.ItemIndex = 1) or (CBBColor.ItemIndex = 2)) then
    aFrg.FrgType := 6 //t 6, IIAB;
  else if (CBTopo.ItemIndex = 3) and ((CBAColor.ItemIndex = 1) or (CBAColor.ItemIndex = 2)) and
    (CBBColor.ItemIndex = 0) then
    aFrg.FrgType := 7 //t 7, IIRA;
  else if (CBTopo.ItemIndex = 3) and (CBAColor.ItemIndex = 0) and ((CBBColor.ItemIndex = 1) or
    (CBBColor.ItemIndex = 2)) then
    aFrg.FrgType := 8 //t 8, IIRB;
  else if (CBTopo.ItemIndex = 3) and ((CBAColor.ItemIndex = 1) or (CBAColor.ItemIndex = 2)) and
    ((CBBColor.ItemIndex = 1) or (CBBColor.ItemIndex = 2)) then
    aFrg.FrgType := 9 //t 9, IIRAB;
  else if (CBTopo.ItemIndex = 4) then
    aFrg.FrgType := 10 //t 10, IIIAB;
  else
  begin
    CBTopo.ItemIndex := 0;
    CBAColor.ItemIndex := 1;
    CBBColor.ItemIndex := 1;
    ShowMessage('ERROR: the combination of options selected is not currently supported.');
    NoError := False;
  end;
  aFrg.GetAtomFrg := CBGetFrag.Checked;
  if SEmax.Value < SEmin.Value then //Check validity of min/max values
  begin
    SEmax.Value := 2;
    SEmin.Value := 2;
    ShowMessage('ERROR: max fragment size cannot be smaller than min fragment size.');
    NoError := False;
  end
  else
  begin
    aFrg.lmin := SEmin.Value;
    aFrg.lmax := SEmax.Value;
  end;
  aFrg.MarkAtom := CBMarkedAtom.ItemIndex;
  aFrg.UseFormalCharge := CBCharge.Checked;
  aFrg.StrictFrg := False;
  //Build a string represention of the Fragmentation scheme
  sFrg := DefFrgToString(aFrg);
  // Add to List
  if NoError then
  begin
    New(pFrg);
    pFrg^ := aFrg;
    LBFrgLst.AddItem(sFrg, TObject(pFrg));
  end;
end;

procedure TFrmxFragmentor.CBAColorChange(Sender: TObject);
begin
  if CBAColor.ItemIndex = 2 then
  begin
    //Display interface dedicated to custom color definition
    LbSDFCAField.Show;
    EdSDFCAfield.Show;
  end
  else
  begin
    //Hide interface dedicated to custom color definition
    LbSDFCAField.Hide;
    EdSDFCAfield.Hide;
  end;
end;

procedure TFrmxFragmentor.CBActivateGBHeaderChange(Sender: TObject);
begin
  if CBActivateGBHeader.Checked then
    GBxHeader.Enabled := True
  else
    GBxHeader.Enabled := False;
    // Disable the StrictFrg button if this is not checked
    CBStrictFrg.Checked:=False;
end;

procedure TFrmxFragmentor.CBBColorChange(Sender: TObject);
begin
  if CBBColor.ItemIndex = 2 then
  begin
    //Display interface dedicated to custom color definition
    LbSDFCBField.Show;
    EdSDFCBfield.Show;
  end
  else
  begin
    //Hide interface dedicated to custom color definition
    LbSDFCBField.Hide;
    EdSDFCBfield.Hide;
  end;
end;

procedure TFrmxFragmentor.CBStrictFrgChange(Sender: TObject);
begin
  RefreshStrictFrg(CBStrictFrg.Checked);
end;

procedure TFrmxFragmentor.CBTopoChange(Sender: TObject);
begin
  if CBTopo.ItemIndex <> 0 then
  begin
    //Display size selection menu
    Lbmin.Show;
    Lbmax.Show;
    SEmin.Show;
    SEmax.Show;
    LbBColor.Show;
    CBBColor.Show;
  end
  else
  begin
    //Hide size selection menu
    Lbmin.Hide;
    Lbmax.Hide;
    SEmin.Hide;
    SEmax.Hide;
    LbBColor.Hide;
    CBBColor.Hide;
  end;
end;

procedure TFrmxFragmentor.RefreshStrictFrg(bStrictFrg: boolean);
var
  i: integer;
  pFrg: PDefFrg;
begin
  for i := 0 to LBFrgLst.Items.Count - 1 do
  begin
    pFrg := PDefFrg(LBFrgLst.Items.Objects[i]);
    pFrg^.StrictFrg := bStrictFrg;
  end;
end;

procedure TFrmxFragmentor.LBSDFAddAndResize(flenme: string);
var
  i,isze: integer;
begin
  if LBSDF.Items.IndexOf(flenme) < 0 then
  begin
    i:=LBSDF.Items.Add(flenme);
    isze := LBSDF.Canvas.TextWidth(LBSDF.Items[i] + '     ');
    if isze > LBSDF.ScrollWidth then
      LBSDF.Width := isze;
  end;
end;

end.
