unit UnitModelBase;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, contnrs, UnitFragment, ISIDA_descriptor, TypePredictor;

type

    { TModelBase }

    TModelBase = class(TObject)
    private
           ffrgs: TFragment;
           fADs: TObjectList;
           fmodtyp: EnumModelType;
    public
          constructor Create;
          destructor Destroy; override;
          property frgs: TFragment read ffrgs write ffrgs;
          property modtyp: EnumModelType read fmodtyp write fmodtyp;
          property ADs: TObjectList read fADs write fADs;
          function Predict(Mol: LRDsc): double; virtual; abstract;
    protected

    end;
implementation

{ TModelBase }

constructor TModelBase.Create;
begin
     ffrgs:=TFragment.Create;
     fADs:=TObjectList.create;
end;

destructor TModelBase.Destroy;
begin
     FreeAndNil(ffrgs);
     FreeAndNil(fADs);
     inherited Destroy;
end;

end.

