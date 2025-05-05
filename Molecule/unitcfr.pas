unit UnitCFR;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, UnitModelBase, ISIDA_descriptor, TypePredictor;

type

    { TCFR }

    TCFR = class(TModelBase)
    private
           fcoef: LRDsc;
    public
          property coef: LRDsc read fcoef;
          constructor Create;
          destructor Destroy; override;
          procedure ReadCFR(afile: string);
          function NmeDesSel(Mol: LRDsc):LRDsc;
          function Predict(Mol: LRDsc): double; override;

    end;

implementation

{ TCFR }

constructor TCFR.Create;
begin
     modtyp:=MLR;
     fcoef:=LRDsc.Create;
     inherited Create;
end;

destructor TCFR.Destroy;
begin
     FreeAndNil(fcoef);
     inherited Destroy;
end;

procedure TCFR.ReadCFR(afile: string);
var
   CFRfle, TSData: TStringList;
   i: integer;
   cf: RDsc;
begin
     CFRfle:=TStringList.Create;
     TSData:=TStringList.Create;
     CFRfle.LoadFromFile(afile);
     for i:=0 to CFRfle.Count-1 do begin
         TSData.Clear;
         TSData.DelimitedText:=CFRfle[i];
         //Coefficient is the first element
         //Standard deviation is the second element
         //Name of the descriptor is the third element
         cf.Idx:=i;
         cf.Cnt:=double(StrToFloat(TSData[0]));
         cf.Nme:=TSData[2];
         fcoef.Add(cf);
     end;
     FreeAndNil(TSData);
     FreeAndNil(CFRfle);
end;

function TCFR.NmeDesSel(Mol: LRDsc): LRDsc;
var
   i,j: integer;
   pitem: PRDsc;
   new: RDsc;
   stmp: string;
begin
     //Initialise
     Result:=LRDsc.Create;
     //Scan the model list of descriptors
     for i:=0 to fcoef.Count-1 do
     begin
          pitem:=Mol.FindRDscNme(fcoef.Objects[i].Nme);
          new:=pitem^;
          Result.Add(new);
     end;
end;

function TCFR.Predict(Mol: LRDsc): double;
var
   i,lo: integer;
   pitem: PRDsc;
begin
     //Initialise considering the case of the constant
     if (fcoef.Objects[0].Nme='a0') then
     begin
          Result:=fcoef.Objects[0].Cnt;
          lo:=1;
     end else
     begin
          Result:=0;
          lo:=0;
     end;
     //Add to the prediction the contribution of each desciptor
     for i:=lo to fcoef.Count-1 do
     begin
          pitem:=Mol.FindRDscNme(fcoef.Objects[i].Nme);
          if (pitem<>nil) then
             Result:=Result+fcoef.Objects[i].Cnt*pitem^.Cnt;
     end;
end;

end.

