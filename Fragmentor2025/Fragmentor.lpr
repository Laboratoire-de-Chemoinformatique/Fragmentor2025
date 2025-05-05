{ Fragmentor of the ISIDA Project

  Copyright (C) 2022 Laboratoire de Chemoinformatique, UMR 7140 CNRS (http://complex-matter.unistra.fr/equipes-de-recherche/laboratoire-de-chemoinformatique/home/)
  contact: Gilles Marcou g.marcou@unistra.fr

  This library is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
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
program Fragmentor;

{$mode objfpc}{$H+}
//{$M 2000000,2000000}

uses {$IFDEF UNIX} {$IFDEF UseCThreads}
  cthreads, {$ENDIF} {$ENDIF}
  Classes, getopts, strutils, SysUtils, contnrs, IOterm, BasicSDF,
  typepredictor, U_BASE_GRAPHES, U_GRAPHES, unitatomandbondtype,
  unitmoleculebase, U_TYPE_GRAPHES, unitmolecule, UnitFragment,
  UnitFragmentBase, UnitMoleculeFrg, UnitXMLFrg;

var
  h:      SDFmini;  // Objet qui gère le fichier sdf
  nme:    string;//PChar;
  stmp,sNFrg:   string;//temp string, last fragments on SVM files
  tmol:   TStringList; // Molécule au format sdfs
  pitem:  PDefFrg;   // Pointeur vers un objet qui décrit une méthode de fragmentation.
  sfield, svalue: string;
  totfraglst: TStringList;
  hdrlst,test: TStringList;
  SLStdo: TStringList;//used for standard output
  frag:   TFragment;   // S'occupe de fragmenter les fichiers sdf et de produire les descripteurs.
  i, ii, j, jj, k, kk: integer;
  IO:     TIOterm; // gestion de la ligne de commande
  HFile:  TextFile;
  dtmp:   double;
  tmplst: TList;
  pkk: PRAtmFrg;
  xmlfrg: TXMLFrg;
  Nfrg: integer;
  SVMSL,ARFFSL,SData,Stuple: TStringList;
begin
  Nfrg:=0;
  //Init variables
  IO := TIOterm.Create;
  IO.Display;
  //Manage the XML file
  case IO.SGRWSFrgXML of
    0: ;
    1: begin //Write
         xmlfrg:=TXMLFrg.Create(IO.SGXFileName, IO.SGOFileName, False);
         xmlfrg.WriteFragmentations();
         for i:=0 to IO.SGFrgType.Count-1 do
         begin
           pitem:=IO.SGFrgType[i];
           xmlfrg.WriteFragmentation(pitem);
         end;
         xmlfrg.Finalize;
         FreeAndNil(xmlfrg);
       end;
    2: begin //Read (add to command line interpreter)
         xmlfrg:=TXMLFrg.Create(IO.SGXFileName, True);
         xmlfrg.ReadFragmentations();
         while (not xmlfrg.stop) do
         begin
           pitem:=new(PDefFrg);
           xmlfrg.ReadFragmentation(pitem);
           IO.SGFrgType.Add(pitem);
           if (IO.SGHFileName='') then IO.SGHFileName:=xmlfrg.lof; //set the header file name, command line takes priority
           xmlfrg.NextFragmentation;
         end;
         FreeAndNil(xmlfrg);
       end;
    3: begin //Read then write
         xmlfrg:=TXMLFrg.Create(IO.SGXFileName, True);
         xmlfrg.ReadFragmentations();
         while (not xmlfrg.stop) do
         begin
           pitem:=new(PDefFrg);
           xmlfrg.ReadFragmentation(pitem);
           IO.SGFrgType.Add(pitem);
           if (IO.SGHFileName='') then IO.SGHFileName:=xmlfrg.lof; //set the header file name, command line takes priority
           xmlfrg.NextFragmentation;
         end;
         FreeAndNil(xmlfrg);
         //
         xmlfrg:=TXMLFrg.Create(IO.SGXFileName, IO.SGOFileName, False);
         xmlfrg.WriteFragmentations();
         for i:=0 to IO.SGFrgType.Count-1 do
         begin
           pitem:=IO.SGFrgType[i];
           xmlfrg.WriteFragmentation(pitem);
         end;
         xmlfrg.Finalize;
         FreeAndNil(xmlfrg);
    end
    else
  end;
  //
  nme := IO.SGIFileName;
  //writeln(nme);
  try
    h    := SDFmini.Create(nme);
  except
    writeln('ERROR: file '+nme+' not found');
    halt;
  end;
  frag := TFragment.Create;
  SLStdo:=TStringList.Create;
  //Set fragmentors
  for i := 0 to IO.SGFrgType.Count - 1 do
  begin
    pitem := IO.SGFrgType[i];
    if pitem^.GetAtomFrg=True then
       frag.GetFrgPerAtom:=True;
    frag.AddFragmentor(pitem^);
  end;
  //set SDF property field and eventually a fragmentation header
  sfield   := IO.SGifield;
  h.prptag := sfield; writeln('h.prptag '+sfield);
  hdrlst   := TStringList.Create;
  if (IO.SGHFileName <> '') then
  begin
    try
      AssignFile(HFile, IO.SGHFileName);
      Reset(HFile);
    except
      writeln('ERROR: Impossible to open header file: ' + IO.SGHFileName);
      writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      writeln('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      raise TIOtermException.Create(UseMsg);
      Halt;
    end;
    repeat
      readln(HFile, dtmp, stmp);
      stmp := DelSpace(stmp);
      hdrlst.Add(stmp);
      writeln(stmp);
    until (EOF(HFile));
    frag.SetFrgLst(hdrlst);
  end;
  //Start fragmentation compound per compound
  if (FileExists(IO.SGOFileName+'.svm') and (not IO.SGbPipe)) then
  begin
    writeln('WARNING: file '+IO.SGOFileName+'.svm overwritten');
    writeln('WARNING: file '+IO.SGOFileName+'.hdr overwritten');
    DeleteFile(IO.SGOFileName+'.svm');
    DeleteFile(IO.SGOFileName+'.hdr');
  end;
  if (FileExists(IO.SGOFileName+'.atm') and (not IO.SGbPipe)) then
  begin
    writeln('WARNING: file '+IO.SGOFileName+'.atm overwritten');
    DeleteFile(IO.SGOFileName+'.atm');
  end;
  Write('/');
  k := 0;
  repeat     //lis la molécule
    tmol := h.NextTMol;
    Inc(k);
    if tmol.Count < 4 then
      continue;
    Write('.' + IntToStr(k));
    frag.FrgReset;    // purge le fragmentor des fragments générés pour la molécule précédente.
    frag.MolToFrgLst(tmol);   //prend le fichier sdf et génère les fragments.
    totfraglst := frag.GetFrgLst;  // Dictionnaire mot qui décrit le frg et pointeur vers une valeur
    if (frag.IsAnyNewFrg) then write('!');
    if frag.GetFrgPerAtom then
    begin
      IO.WriteAtmFrg(frag.FrgPerAtom);
      //Debug Coloration
      {writeln;
      for i:=1 to frag.FrgPerAtom.Count-1 do
      begin
        tmplst:=frag.FrgPerAtom.Items[i] as TList;
        write('Atom: '+IntToStr(i)+' Fragments:');
        for j:=0 to tmplst.Count-1 do begin
          pkk:=PRAtmFrg(tmplst.Items[j]);
          write(' '+IntToStr(pkk^.idx+1)+'('+IntToStr(pkk^.len)+')|'+totfraglst[pkk^.idx]);
        end;
        writeln();
      end;}
      //Debug Coloration
    end;
    svalue := h.prpval;
    svalue:=Trim(svalue);
    if (IO.SGoformat = 'SVM') or (IO.SGoformat= 'ARFF') then
    begin
      IO.WriteSVM(svalue, totfraglst);
    end
    else if (IO.SGoformat = 'SMF') then
    begin
      IO.WriteSMF(totfraglst);
      IO.WriteProp(svalue);
    end
    else if (IO.SGoformat = 'CSV') then
    begin
      IO.WriteSParse(svalue,totfraglst);
    end
    else if (IO.SGoformat = 'STDO') then
    begin
      stmp:=svalue;
      for i := 0 to totfraglst.Count - 1 do
      begin
        j := PInteger(totfraglst.Objects[i])^;
        if (j <> 0) then
        begin
          stmp:=stmp+';' + IntToStr(j);
        end else stmp:=stmp+';0';
      end;
      SLStdo.Add(stmp);
    end;
    if (IO.SGiRS>0) and (k>1) and ((k-1) mod IO.SGiRS = 0) then
    begin
      if (IO.SGoformat = 'SMF') or (IO.SGoformat = 'SVM') or (IO.SGoformat= 'CSV') or (IO.SGoformat = 'ARFF') then begin
         IO.WriteHeader(totfraglst);
         writeln('|');
      end;
    end;
  until (h.fend);
  writeln('/');
  if (IO.SGoformat = 'CSV') then IO.CleanSParse(totfraglst.Count);
  if (IO.SGoformat = 'SMF') or (IO.SGoformat = 'SVM') or (IO.SGoformat= 'CSV') or (IO.SGoformat = 'ARFF') then begin
    IO.WriteHeader(totfraglst)
  end else if (IO.SGoformat = 'STDO') then
  begin
    writeln('-----');
    write(QuotedStr('prop'));
    for i := 0 to totfraglst.Count - 1 do
      write(';'+QuotedStr(totfraglst[i]));
      writeln;
    write(SLStdo.Text);
  end else
    for i := 0 to totfraglst.Count - 1 do
      write(format('%6s', [IntToStr(i + 1)]) + '. ' + totfraglst[i]);
  writeln('-----');
  if (IO.SGoformat='SVM') then
  begin
     Nfrg:=totfraglst.Count;
     sNFrg:=' '+IntToStr(Nfrg)+':0';//last fragment on an SVM line
     SVMSL:=TStringList.Create;
     SData:=TStringList.Create;
     SData.Delimiter:=' ';
     Stuple:=TStringList.Create;
     Stuple.Delimiter:=':';
     SVMSL.LoadFromFile(IO.SGOFileName+'.svm');
     for i:=0 to SVMSL.Count-1 do
     begin
       stmp:=SVMSL[i];
       SData.DelimitedText:=stmp;
       if SData.Count>1 then
       begin
         Stuple.DelimitedText:=SData[SData.Count-1];
         kk:=StrToInt(Stuple[0]);
         if (kk<Nfrg) then
            stmp:=stmp+sNFrg;
       end else stmp:=stmp+sNFrg; //Add the last fragment even if the list of fragments is empty
       SVMSL[i]:=stmp;
     end;
     SVMSL.SaveToFile(IO.SGOFileName+'.svm');
     FreeAndNil(Stuple);
     FreeAndNil(SData);
     FreeAndNil(SVMSL);
  end;
  if (IO.SGoformat='ARFF') then
  begin
    // ARFF conversion
    Nfrg := totfraglst.Count;
    SVMSL := TStringList.Create;
    ARFFSL := TStringList.Create;
    //Head of ARFF
    ARFFSL.Add('@RELATION "'+IO.SGIFileName+'"');
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
    SVMSL.LoadFromFile(IO.SGOFileName+'.svm');
    for ii := 0 to SVMSL.Count - 1 do
    begin
      stmp := SVMSL[ii];
      SData.DelimitedText := stmp;
      //ARFF output
      stmp:='{';
      for jj:=1 to SData.Count-1 do
      begin
          STuple.DelimitedText:=SData[jj];
          kk:=StrToInt(STuple[0])-1;
          stmp:=stmp+IntToStr(kk)+' '+STuple[1]+', ';
      end;
      //In ARFF no need to add a 0 for the last descriptor
      stmp:=stmp+IntToStr(Nfrg)+' '+SData[0]+'}';
      ARFFSL.Add(stmp);
    end;
    FreeAndNil(Stuple);
    FreeAndNil(SData);
    FreeAndNil(SVMSL);
    DeleteFile(IO.SGOFileName+'.svm');//Remove the temporary SVM file
    ARFFSL.SaveToFile(IO.SGOFileName+'.arff');
    FreeAndNil(ARFFSL);
  end;
  FreeAndNil(frag);
  FreeAndNil(hdrlst);
  FreeAndNil(h);
  FreeAndNil(IO);
  FreeAndNil(SLStdo);
end.

