program xFragmentor;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Interfaces, // this includes the LCL widgetset
  Forms, unitxfragmentor, UnitFragment, UnitXMLFrg, BasicSDF, TypePredictor
  { you can add units after this };

{$R *.res}

begin
  RequireDerivedFormResource := True;
  Application.Initialize;
  Application.CreateForm(TFrmxFragmentor, FrmxFragmentor);
  Application.Run;
end.

