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
unit UnitMixture;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, contnrs, BasicSDF, UnitFragment, ISIDA_descriptor, math;

type

  EnumComboMode=(CONCAT,DIFF,ADIFF,SUM,PROD,UNK);

  { TMixture }

  TMixture = class(TObject)
  private
    fNumComponent: integer;
    fNumMixture: integer;
    ffrgs: TFragment;
    faMFrgLst: TStringList; //List of mixture fragments
    faMHdrLst: TStringList; //List of headers for strict fragmentation
    fLfrgs:   TObjectList; //List of fragmentations
    fComboMode: EnumComboMode;
    faMol: LRDsc; // Store the result of a combination of descriptors for a mixture
    fLWeights: TList; //List of PDouble
    fBasicSDFList: TObjectList; //List of iterators on SDF/ fbasicSDFList own SDFmini objects and SDFmini objects own SDF containing TStringList
    fbeom: Boolean;// End of file indicator / True if any of the component of the mixture reach the end of file during mixture fragmentation
    fstrict: Boolean;
    function GetWeight(id: integer): Double;
    procedure SetFrgs(id: integer; afrgs: TFragment);
    function GetFrgs(id: integer): TFragment;
    procedure SetWeight(id: integer; AValue: Double);
  public
    constructor Create;
    destructor Destroy; override;
    property NumComponent: integer Read fNumComponent;
    property NumMixture: integer Read fNumMixture;
    property ComboMode: EnumComboMode Read fComboMode Write fComboMode;
    property aMFrgLst: TStringList Read faMFrgLst;
    property aMHdrLst: TStringList Read faMHdrLst;
    property aMol: LRDsc Read faMol;
    property Lfrgs[id: integer]: TFragment Read GetFrgs Write SetFrgs;
    property LWeight[id: integer]: Double Read GetWeight Write SetWeight;
    property beom: Boolean Read fbeom;
    property strict: Boolean Read fstrict;
    procedure AddStrict;
    procedure AddFrgs;
    procedure AddFrgs(afrgs: TFragment);
    procedure SetSDFList(flenms: TStringList); //Setup the list of connectors to SDF
    procedure ReplaceSDF(flenms: TStringList);
    procedure ResetSDF; //Reset all SDF connectors to their original state
    procedure ResetFrgs; //Reset all Fragments
    procedure FragmentEachComponent;//Enumerate fragments from each component of the mixture -> to do first
    procedure FragmentNextMixture;//Fragmentations of the mixture as they are read from the fBasicSDFList
    procedure FragmentMixture(aSDFLst: TObjectList);//Fragmentations of the components of the mixture by their own fragmentor
    procedure Combine;
    procedure Concatenate;
    procedure Difference;
    procedure AbsDifference;
    procedure Summation;
    procedure Product;
    procedure WriteComponentHeader(id: integer; outfle: string);
    procedure WriteComponentHeaders(outfle: string);
    procedure WriteMixtureHeader(outfle:string);
  protected

  end;

implementation

{ TModelBase }

procedure TMixture.SetFrgs(id: integer; afrgs: TFragment);
begin
  try
    fLfrgs[id]:=TObject(afrgs);
  except
    raise Exception.Create('ERROR UnitMixture.SetFrgs: Impossible to access to fragments of mixture component '+IntToStr(id));
    halt;
  end;
end;

function TMixture.GetWeight(id: integer): Double;
var
  pWeight: PDouble;
begin
  pWeight:=fLWeights[id];
  if pWeight<>nil then Result:=pWeight^ else Result:=NaN;
end;

function TMixture.GetFrgs(id: integer): TFragment;
begin
  try
    Result:=fLfrgs[id] as TFragment;
  except
    raise Exception.Create('ERROR UnitMixture.GetFrgs: Impossible to access to fragments of mixture component '+IntToStr(id));
    halt;
  end;
end;

procedure TMixture.SetWeight(id: integer; AValue: Double);
var
  pWeight: PDouble;
begin
  pWeight:=fLWeights[id];
  if pWeight<>nil then pWeight^:=AValue else
  begin
    raise Exception.Create('ERROR UnitMixture.SetWeight: No weight for mixture component '+IntToStr(id));
    halt;
  end;
end;

constructor TMixture.Create;
begin
  faMFrgLst:=TStringList.Create;
  fBasicSDFList:=TObjectList.Create;
  fBasicSDFList.OwnsObjects:=True;
  faMol  := LRDsc.Create;
  fLfrgs := TObjectList.Create;
  fLfrgs.OwnsObjects:=True;
  fComboMode:=UNK;
  fNumComponent:=0;
  fLWeights:=TList.Create;
  fbeom:=False;
end;

destructor TMixture.Destroy;
var
  i: integer;
  pWeight: PDouble;
begin
  FreeAndNil(faMFrgLst);
  FreeAndNil(fBasicSDFList);
  FreeAndNil(faMol);
  FreeAndNil(fLfrgs);
  for i:=0 to fLWeights.Count-1 do
  begin
    pWeight:=fLWeights[i];
    if pWeight<>nil then dispose(pWeight);
  end;
  FreeAndNil(fLWeights);
  inherited Destroy;
end;

procedure TMixture.AddFrgs;
var
  afrgs: TFragment;
begin
  afrgs:=TFragment.Create;
  AddFrgs(afrgs);
end;

procedure TMixture.AddStrict;
begin
  fstrict:=True;
end;

procedure TMixture.AddFrgs(afrgs: TFragment);
var
  pWeight: PDouble;
begin
  fLfrgs.Add(afrgs);
  if fLfrgs.Count>fNumComponent then
  begin
    raise Exception.Create('ERROR TMixture.AddFrgs: too many fragmentation in the list compared to the number of mixture components ('+IntToStr(fLfrgs.Count)+'>'+IntToStr(fNumComponent)+')');
    halt;
  end;
  //
  new(pWeight);
  pWeight^:=1.0;//Default weight is 1 -> no modification to descriptors value
  fLWeights.Add(pWeight);
end;

procedure TMixture.SetSDFList(flenms: TStringList);
var
  i: integer;
  aBasicSDF: SDFmini;
begin
  if fNumComponent=0 then
    fNumComponent:=flenms.Count
  else if fNumComponent<>flenms.Count then
  begin
    raise Exception.Create('ERROR TMixture.SetSDFList: fragmentation list and mixture component mismatch ('+IntToStr(fNumComponent)+'<>'+IntToStr(flenms.Count)+')');
    halt;
  end;
  for i:=0 to flenms.Count-1 do
  begin
    aBasicSDF:=SDFmini.Create(flenms[i]);
    fBasicSDFList.Add(aBasicSDF);
  end;
end;

procedure TMixture.ResetSDF;
var
  i: integer;
  aSDFmini: SDFmini;
begin
  for i:=0 to fBasicSDFList.Count-1 do
  begin
    aSDFmini:=(fBasicSDFList[i] as SDFmini);
    aSDFmini.Clear;
//    fBasicSDFList.Remove(aSDFmini);
  end;
end;

procedure TMixture.ReplaceSDF(flenms: TStringList);
var
  i: integer;
  aSDFmini,aBasicSDF: SDFmini;
begin
  for i:=0 to fBasicSDFList.Count-1 do
  begin
    aSDFmini:=(fBasicSDFList[i] as SDFmini);
    aSDFmini.Clear;
    aBasicSDF:=SDFmini.Create(flenms[i]);
    fBasicSDFList[i]:=aBasicSDF;
  end;
end;


procedure TMixture.ResetFrgs;
var
  i: integer;
begin
  for i:=0 to NumComponent-1 do
  begin
    Lfrgs[i].FrgReset;
  end;
end;


procedure TMixture.FragmentEachComponent;
//Enumerates the list of fragments of each component of the mixture
var
  i,k,kold,j: integer;
  aBasicSDF: SDFmini;
  aSDF: TStringList;
  stmp: string;
begin
  kold:=0;
  for i:=0 to NumComponent-1 do
  begin
    aBasicSDF:=(fBasicSDFList[i] as SDFmini);
    k:=1;
    repeat
      aSDF:=aBasicSDF.NextTMol;

      //for j:=0 to aSDF.Count-1 do begin writeln('--------------'+aSDF[j]); end;

      Lfrgs[i].FrgReset;
      Lfrgs[i].MolToFrgLst(aSDF);
      stmp:=LFrgs[i].GetFrgLst.Text;
      Inc(k);
    until aBasicSDF.fend;
    if kold=0 then
      kold:=k
    else if kold<>k then
    begin
      {$IfDef CGIMODE}
      writeln('<font color="red"><b>ERROR: incompatible number of entries in input SDFs ('+IntToStr(k)+'<>'+IntToStr(kold)+')</b></font>');
      writeln('<br><br><img src="http://infochim.u-strasbg.fr/IMG/png/back_35x35.png" border="0" height="35"/><a href="predictor_mixtures.cgi" style="text-decoration: none;">Back to Main Menu</a>');
      {$EndIf}
      raise Exception.Create('ERROR TMixture.FragmentEachComponent: incompatible number of entries in input SDFs ('+IntToStr(k)+'<>'+IntToStr(kold)+')');
      halt;
    end;
    aBasicSDF.Clear;
  end;
end;

procedure TMixture.FragmentNextMixture;
//Fragmentation of mixtures from the SDF files / each TFragment in Lfrgs must have a properly setup TotFragLst
var
  i,k,kold: integer;
  aBasicSDF: SDFmini;
  aSDF: TStringList;
  aSDFLst: TObjectList;
begin
  if (not fbeom) then
  begin
    aSDFLst:=TObjectList.Create;
    aSDFLst.OwnsObjects:=False;
    for i:=0 to fNumComponent-1 do
    begin
      aBasicSDF:=(fBasicSDFList[i] as SDFmini);
      aSDF:=aBasicSDF.NextTMol;
      aSDFLst.Add(aSDF);
      //
      if (not fbeom) and (aBasicSDF.fend) then fbeom:=True; //a mixture component has reached the end of file
    end;
    FragmentMixture(aSDFLst);
    Combine;
    //
    FreeAndNil(aSDFLst);
  end;
end;

procedure TMixture.FragmentMixture(aSDFLst: TObjectList);
//Fragmentation of the components of a mixture / each TFragment in Lfrgs must have a properly setup TotFragLst
var
  i: integer;
  aSDF: TStringList;
  afrgs: TFragment;
begin
  if aSDFLst.Count <> fLfrgs.Count then
  begin
    raise Exception.Create('ERROR Mixture.FragmentComponents: fragmentation list and molecule list mismatch');
  end;
  for i:=0 to fNumComponent-1 do
  begin
    afrgs:=fLfrgs.Items[i] as TFragment;
    aSDF:=aSDFLst.Items[i] as TStringList;
    afrgs.FrgReset;
    afrgs.MolToFrgLst(aSDF);
  end;
end;

procedure TMixture.Combine;
// Common gate for combination of descriptors in a mixture
// The LRDsc is managed by this object.
// The calling subprogramm shall create a deep copy if one if needed
begin
  case fComboMode of
    CONCAT: Concatenate;
    DIFF: Difference;
    ADIFF: AbsDifference;
    SUM: Summation;
    PROD: Product;
    else begin
      raise Exception.Create('ERROR UnitMixture.Combine: unknown or undefined descriptors mixture combination');
      halt;
    end;
  end;
end;

procedure TMixture.Concatenate;
//Descriptors of mixture component are concatenated
var
  i,j,IdCnt: integer;
  afrgs: TFragment;
  aFrgLst: TStringList;
  aRdsc: RDsc;
begin
  faMol.Clear;
  faMFrgLst.Clear;
  IdCnt:=0;
  for i:=0 to fLfrgs.Count-1 do
  begin
    afrgs:=fLfrgs.Items[i] as TFragment;
    aFrgLst:=afrgs.GetFrgLst;
    for j:=0 to afrgs.TotFrgList.Count-1 do
    begin
      aRdsc.Idx:=IdCnt;
      aRdsc.Nme:=aFrgLst.Strings[j]+'_'+IntToStr(i);
      aRdsc.Cnt:=LWeight[i]*PInteger(aFrgLst.Objects[j])^;
      faMFrgLst.Add(aRdsc.Nme);
      if aRDsc.Cnt<>0 then faMol.Add(aRdsc);//Sparse storage
      Inc(IdCnt);
    end;
  end;
end;

procedure TMixture.Difference;
var
  afrgs1,afrgs2: TFragment;
  aFrgLst1, aFrgLst2: TStringList;
  j,nFrgs: integer;
  aRdsc: RDsc;
  w1, w2: Double;
begin
  faMol.Clear;
  faMFrgLst.Clear;
  if fLfrgs.Count<>2 then
  begin
    raise Exception.Create('ERROR UnitMixture.Diffrence: number of mixture components <> 2');
    halt;
  end;
  afrgs1:=fLfrgs.Items[0] as TFragment;
  afrgs2:=fLfrgs.Items[1] as TFragment;
  aFrgLst1:=afrgs1.GetFrgLst;
  aFrgLst2:=afrgs2.GetFrgLst;
  w1:=LWeight[0];
  w2:=LWeight[1];
  if aFrgLst1.Count <> aFrgLst2.Count then
  begin
    raise Exception.Create('ERROR UnitMixture.Difference: mixture components not sharing the same fragment list');
    halt;
  end;
  nFrgs:=aFrgLst1.Count;
  for j:=0 to nFrgs-1 do
  begin
    aRdsc.Idx:=j+1;
    aRdsc.Nme:=aFrgLst1.Strings[j]+'_1-'+aFrgLst2.Strings[j]+'_2';
    aRdsc.Cnt:=w1*PInteger(aFrgLst1.Objects[j])^-w2*PInteger(aFrgLst2.Objects[j])^;
    faMFrgLst.Add(aRdsc.Nme);
    if aRDsc.Cnt<>0 then faMol.Add(aRdsc);//Sparse storage
  end;
end;

procedure TMixture.AbsDifference;
var
  i: integer;
  aRDsc: RDsc;
begin
  Difference;
  for i:=0 to faMol.Count-1 do
  begin
    aRDsc:=faMol.Objects[i];
    aRDsc.Cnt:=abs(aRDsc.Cnt);
    aRDsc.Nme:='|'+aRDsc.Nme+'|';
  end;
end;

procedure TMixture.Summation;
var
  i,j: integer;
  afrgs: TFragment;
  aFrgLst,aFrgLst0: TStringList;
  aRdsc: RDsc;
  aPRdsc: PRdsc;
begin
  faMol.Clear;
  faMFrgLst.Clear;
  aFrgLst0:=TFragment(fLfrgs.Items[0]).GetFrgLst;
  for j:=0 to aFrgLst0.Count-1 do
    faMFrgLst.Add(aFrgLst0[j]+'_1');
  for i:=0 to fLfrgs.Count-1 do
  begin
    afrgs:=fLfrgs.Items[i] as TFragment;
    aFrgLst:=afrgs.GetFrgLst;
    for j:=0 to afrgs.TotFrgList.Count-1 do
    begin
      if (i>0) then faMFrgLst[j]:=faMFrgLst[j]+'+'+aFrgLst.Strings[j]+'_'+IntToStr(i+1);
      aRdsc.Idx:=j+1;
      aRdsc.Nme:=aFrgLst.Strings[j];
      aRdsc.Cnt:=LWeight[i]*PInteger(aFrgLst.Objects[j])^;
      faMFrgLst[j]:=aRdsc.Nme;
      if aRDsc.Cnt<>0 then
      begin
        aPRdsc:=aMol.FindRDscNme(aRdsc.Nme);//Search if descriptor already exists
        if aPRdsc<>nil then
        begin
          //aPRdsc^.Nme:=aPRdsc^.Nme+'+'+aRdsc.Cnt;
          aPRdsc^.Cnt:=aPRdsc^.Cnt+aRdsc.Cnt//add descriptor value to existing value
        end else
        begin
          //aPRdsc^.Nme:=aRdsc.Nme;
          aMol.Add(aRdsc);//setup descriptor value
        end;
      end;
    end;
  end;
end;

procedure TMixture.Product;
var
  i,j: integer;
  afrgs: TFragment;
  aFrgLst,aFrgLst0: TStringList;
  aRdsc: RDsc;
  aPRdsc: PRdsc;
begin
  faMol.Clear;
  faMFrgLst.Clear;
  aFrgLst0:=TFragment(fLfrgs.Items[0]).GetFrgLst;
  for j:=0 to aFrgLst0.Count-1 do
    faMFrgLst.Add(aFrgLst0[j]+'_1');
  for i:=0 to fLfrgs.Count-1 do
  begin
    afrgs:=fLfrgs.Items[i] as TFragment;
    aFrgLst:=afrgs.GetFrgLst;
    for j:=0 to afrgs.TotFrgList.Count-1 do
    begin
      if (i>0) then faMFrgLst[j]:=faMFrgLst[j]+'*'+aFrgLst.Strings[j]+'_'+IntToStr(i+1);
      aRdsc.Idx:=j+1;
      aRdsc.Nme:=aFrgLst.Strings[j];
      aRdsc.Cnt:=Power(PInteger(aFrgLst.Objects[j])^,LWeight[i]);
      faMFrgLst[j]:=aRdsc.Nme;
      aPRdsc:=aMol.FindRDscNme(aRdsc.Nme);//Search if descriptor already exists
      if aPRdsc<>nil then
      begin
        aPRdsc^.Cnt:=aPRdsc^.Cnt*aRdsc.Cnt//multiply descriptor value to existing value
      end else // No sparse storage: if one of the component has a zero value, the mixture has a zero value
      begin
        aMol.Add(aRdsc);//setup descriptor value
      end;
    end;
  end;
end;

procedure TMixture.WriteComponentHeader(id: integer; outfle: string);
var
  Outfile: Text;
  i: integer;
  frglst: TStringList;
begin
  AssignFile(OutFile, ChangeFileExt(Outfle,'_'+IntToStr(id)+'.hdr'));
  ReWrite(OutFile);
  frglst:=Lfrgs[id].GetFrgLst;
  for i := 0 to frglst.Count - 1 do
  begin
    writeln(Outfile, format('%6s', [IntToStr(i + 1)]) + '. ' +
      format('%120s', [frglst.Strings[i]]));
  end;
  Close(Outfile);
end;

procedure TMixture.WriteComponentHeaders(outfle: string);
var
  i: integer;
begin
  for i:=0 to fNumComponent-1 do
    WriteComponentHeader(i,outfle);
end;

procedure TMixture.WriteMixtureHeader(outfle: string);
var
  Outfile: Text;
  i: integer;
  frglst: TStringList;
begin
  AssignFile(OutFile, ChangeFileExt(Outfle,'_M.hdr'));
  ReWrite(OutFile);
  for i := 0 to faMFrgLst.Count - 1 do
  begin
    writeln(Outfile, format('%6s', [IntToStr(i + 1)]) + '. ' +
      format('%120s', [faMFrgLst.Strings[i]]));
  end;
  Close(Outfile);
end;

end.
