program multifragmentor;

{$mode objfpc}{$H+}

uses {$IFDEF UNIX} {$IFDEF UseCThreads}
  cthreads, {$ENDIF} {$ENDIF}
  Classes,
  SysUtils,
  CustApp { you can add units after this },
  Contnrs,
  TypInfo,
  StrUtils,
  BasicSDF,
  UnitFragment,
  UnitMixture;

const
  ShtArgs='hi:o:m:f:';
  LngArgs='help input: output: mode: fragmentor: StrictFrg hdr:';

type

  { TMultiFragmentor }

  TMultiFragmentor = class(TCustomApplication)
  protected
    procedure DoRun; override;
  public
    constructor Create(TheOwner: TComponent); override;
    destructor Destroy; override;
    procedure WriteHelp; virtual;
    function FragmentOptInterpretor(sopt: string): TFragment;
  end;

  { TMultiFragmentor }

  procedure TMultiFragmentor.DoRun;
  var
    ErrorMsg: string;
    sinput,soutput,smode: string;
    sw: string;
    afrgt: TFragment;
    adeffrgt: DefFrg;
    padeffrgt: PDefFrg;
    amixture: TMixture;
    aSVM: TStringList;
    i, j, k: integer;
    sdfles: TStringList;
    sfrgnts: TStringList;
    strictfrg: boolean; //if Strict fragmentation
    SGHFileNames, hdrlst: TStringList; //list of headers for Strict fragmentation
    aofhdrlst:Array of TStringList;
    SGHFileName: string;
    HFile: TextFile;
    dtmp: double;
    //
    stmp: string;
  begin
    sdfles:=TStringList.Create;
    SGHFileNames:=TStringList.Create;
    SGHFileNames.Delimiter:=',';
    SGHFileNames.StrictDelimiter:=True;
    amixture := TMixture.Create;
    // quick check parameters
    ErrorMsg := CheckOptions(ShtArgs, LngArgs);
    if ErrorMsg <> '' then
    begin
      ShowException(Exception.Create(ErrorMsg));
      Terminate;
      Exit;
    end;

    // parse parameters
    if HasOption('h', 'help') then
    begin
      WriteHelp;
      FreeAndNil(sdfles);
      FreeAndNil(amixture);
      FreeAndNil(SGHFileNames);
      Terminate;
      Exit;
    end;

    //input
    if HasOption('i', 'input') then
    begin
      sinput:=GetOptionValue('i','input');
      sdfles.CommaText:=sinput;
    end else
    begin
      raise Exception.Create('ERROR MultiFragmentor: input files are mandatory');
      FreeAndNil(sdfles);
      FreeAndNil(amixture);
      FreeAndNil(SGHFileNames);
      Terminate;
      Exit;
    end;
    //File connectors
    amixture.SetSDFList(sdfles);

    //output
    if HasOption('o', 'output') then
    begin
      soutput:=GetOptionValue('o','output');
    end;

    //combination mode
    if HasOption('m', 'mode') then
    begin
      smode:=GetOptionValue('m','mode');
      try
        amixture.ComboMode:=EnumComboMode(GetEnumValue(Typeinfo(EnumComboMode), smode));
      Except
        raise Exception.Create('ERROR MultiFragmentor: unknown combination of descriptors ( '+smode+' )');
        FreeAndNil(sdfles);
        FreeAndNil(amixture);
        FreeAndNil(SGHFileNames);
        Terminate;
        Exit;
      end;
    end;

    //StrictFrg
    if HasOption('StrictFrg') then
    begin
      strictfrg:=True;
      if (not HasOption('hdr')) then
      begin
        raise Exception.Create('ERROR MultiFragmentor: Missing header file for strict fragmentation (option --hdr)');
        FreeAndNil(sdfles);
        FreeAndNil(amixture);
        FreeAndNil(SGHFileNames);
        strictfrg:=False;
        Terminate;
        Exit;
      end;
    end;

    //Headers
    if HasOption('hdr') then
    begin
      if (not HasOption('StrictFrg')) then
      begin
        raise Exception.Create('ERROR MultiFragmentor: Missing --StrictFrg option for strict fragmentation');
        FreeAndNil(sdfles);
        FreeAndNil(amixture);
        FreeAndNil(SGHFileNames);
        Terminate;
        Exit;
      end else begin
        strictfrg:=True;
        stmp:=GetOptionValue('hdr');
        SGHFileNames.DelimitedText:=stmp;//Components are comma separated
        for i:=0 to SGHFileNames.Count-1 do
        begin
          SGHFileName:=SGHFileNames[i];
          if not (FileExists(SGHFileName)) then begin
            raise Exception.Create('ERROR MultiFragmentor: Impossible to open header file: ' + SGHFileName);
            FreeAndNil(sdfles);
            FreeAndNil(amixture);
            FreeAndNil(SGHFileNames);
            Terminate;
            Exit;
          end else begin
            hdrlst := TStringList.Create;
            SetLength(aofhdrlst,i+1);
            if (SGHFileName <> '') then
            begin
              AssignFile(HFile, SGHFileName);
              Reset(HFile);
              repeat
                readln(HFile, dtmp, stmp);
                stmp := DelSpace(stmp);
                hdrlst.Add(stmp);
                writeln(stmp);
              until (EOF(HFile));
              //frag.SetFrgLst(hdrlst);
              aofhdrlst[i]:=hdrlst;
            end;
          end;
        end;
      end;
    end;

    //fragmentors of the mixture components
    sfrgnts:=TStringList.Create;
    sfrgnts.Delimiter:=',';
    sfrgnts.StrictDelimiter:=True;
    if HasOption('f','fragmentor') then
    begin
      stmp:=GetOptionValue('f','fragmentor');
      sfrgnts.DelimitedText:=stmp;//Components are comma separated
      for i:=0 to sfrgnts.Count-1 do
      begin
        stmp:=sfrgnts[i];
        afrgt:=FragmentOptInterpretor(stmp);
        if(strictfrg) then begin
          afrgt.SetFrgLst(aofhdrlst[i]);
        end;
        amixture.AddFrgs(afrgt);
      end;
    end;
    FreeAndNil(sfrgnts);

    { add your program here }
    //Setup the header for all components; the fragments of mixtures can be computed only when the fragments of each component are defined
    amixture.FragmentEachComponent;
    amixture.WriteComponentHeaders(soutput);
    // If in need to redefine the fragments
    //for i:=0 to amixture.NumComponent-1 do
    //  (amixture.Lfrgs[i] as TFragment).SetFrgLst(Lhdr[i] as TStringList);
    aSVM := TStringList.Create;
    repeat
      amixture.FragmentNextMixture;
      //output
      aSVM.Add(amixture.aMol.WriteSVM());
      //
    until amixture.beom;
    amixture.WriteMixtureHeader(soutput);
    aSVM.SaveToFile(ChangeFileExt(soutput,'.svm'));
    Writeln(aSVM.Text);
    FreeAndNil(aSVM);
    // stop program loop
    FreeAndNil(sdfles);
    FreeAndNil(amixture);
    FreeAndNil(SGHFileNames);
    Terminate;
  end;

  constructor TMultiFragmentor.Create(TheOwner: TComponent);
  begin
    inherited Create(TheOwner);
    StopOnException := True;
  end;

  destructor TMultiFragmentor.Destroy;
  begin
    inherited Destroy;
  end;

  procedure TMultiFragmentor.WriteHelp;
  begin
    { add your help code here }
    writeln('Usage: ', ExeName, ' -h -i <component1.sdf>,<component2.sdf>,...'+
    '-o output -m [CONCAT,DIF,ADIF,ADD,PROD]'+
    '-f "t <int> l <int> u <int> ...","t <int> l <int> u <int> ...",...'+
    '--StrictFrg --hdr <header1.hdr>,<header2.hdr>,...');
    writeln('-h, --help: help message');
    writeln('-i, --input: input SDF file names of the components, seperated by commas');
    writeln('-o, --output: base name for the output');
    writeln('-m, --mode: how components descriptors are combined');
    writeln('CONCAT: concatenated');
    writeln('DIF: difference');
    writeln('ADIF: absolute valude of the difference');
    writeln('ADD: summation');
    writeln('PROD: product');
    writeln('-f, --fragmentor: definition of a fragmentors, comma separated; fragmentors are applied to SDF in the same order they are listed on the -i option');
    writeln('--StrictFrg: can only be used with the --hdr option to indicate the header file (.hdr). The outputed svm will be limited to the descriptors indicated in that header file and keeping the same order.');
    writeln('--hdr: Name of a header file. If present, the fragmentation will reproduce the list of fragments the header contains. The output header file will match this input concatenated with new fragments discovered at the end.');
  end;

  function TMultiFragmentor.FragmentOptInterpretor(sopt: string): TFragment;
  var
    sfrgdef: TStringList;
    afrgt: TFragment;
    adeffrgt: DefFrg;
    bStore: boolean;
    i,j: integer;
    stmp: string;
  begin
    sfrgdef:=TStringList.Create;
    sfrgdef.Delimiter:=' ';
    sfrgdef.StrictDelimiter:=True;
    //
    sfrgdef.DelimitedText:=sopt;//Options are space separated / a new fragmentor starts with 't'
    afrgt:=TFragment.Create;//Destruction of this object must be managed by a calling subporgram
    bStore:=False;
    for j:=0 to sfrgdef.Count -1 do
    begin
      stmp:=sfrgdef[j];
      if (Pos('t',stmp)=1) then //type
      begin
        if (bStore) then afrgt.AddFragmentor(adeffrgt);
        stmp:=stmp.Split(':')[1];//Option:Value
        InitDefFrg(adeffrgt);
        adeffrgt.FrgType:=StrToInt(stmp);
        bStore:=True;
      end else
      if (Pos('l',stmp)=1) then //len min
      begin
        stmp:=stmp.Split(':')[1];//Option:Value
        adeffrgt.lmin:=StrToInt(stmp);
      end else
      if (Pos('u',stmp)=1) then //len max
      begin
        stmp:=stmp.Split(':')[1];//Option:Value
        adeffrgt.lmax:=StrToInt(stmp);
      end;
    end;
    if (bStore) then afrgt.AddFragmentor(adeffrgt);
    FreeAndNil(sfrgdef);
    //
    Result:=afrgt;
  end;

var
  Application: TMultiFragmentor;
begin
  Application := TMultiFragmentor.Create(nil);
  Application.Title:='MultiFragmentor';
  Application.Run;
  Application.Free;
end.
