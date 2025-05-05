{ Fragmentor of the ISIDA Project

  Copyright (C) 2024 Laboratoire de Chemoinformatique, UMR 7140 CNRS (https://complex-matter.unistra.fr/equipes-de-recherche/laboratoire-de-chemoinformatique/home/)
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
unit BasicSDF;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, strutils;
  
type
    FileException=Class(Exception);

    { SDFmini }

    SDFmini=Class(TObject)
    private
           ffle: TextFile;
           fprptag: string;
           fprpval: string;
           fnme: string;
           ffend: boolean;
           fmol: TStringList;
           fmolString: string;
           fMolIndex: integer;
           fNb_mol: integer;
    public
          constructor Create;
          constructor Create(flenme:string);
          destructor Destroy; override;
          procedure Clear;
          procedure SetFile(flenme:string);
          function NextTMol: TStringList;
          function NextTMolFirstString: string;
          procedure countMol;
          property nme:string read fnme write fnme;
          property prptag: string read fprptag write fprptag;
          property prpval: string read fprpval write fprpval;
          property fend: boolean read ffend write ffend;
          property Nb_mol: integer read fNb_mol write  fNb_mol;
          property MolIndex: integer read fMolIndex write fMolIndex;
    end;
    
implementation

{ SDFmini }

constructor SDFmini.Create;
begin
     inherited Create;
     ffend:=false;
     fmol:=TStringList.Create;
     fMolIndex:=-1;
end;

constructor SDFmini.Create(flenme:string);
begin
     inherited Create;
     fmol:=TStringList.Create;
     SetFile(flenme);
{     fnme:=flenme;
     AssignFile(ffle,fnme);
     if not (FileExists(flenme)) then Raise FileException.Create('ERROR: File '+fnme+' not found')
     else Reset(ffle);
     ffend:=false;
     fMolIndex:=0;
     countMol;}
end;

destructor SDFmini.Destroy;
begin
     close(ffle);
     FreeAndNil(fmol);
     inherited Destroy;
end;

procedure SDFmini.Clear;
begin
     reset(ffle);
     fMolIndex:=0;
     ffend:=false;
end;

procedure SDFmini.SetFile(flenme: string);
begin
     fnme:=flenme;
     AssignFile(ffle,fnme);
     if not (FileExists(flenme)) then Raise FileException.Create('ERROR: File '+fnme+' not found')
     else Reset(ffle);
     ffend:=False;
     fMolIndex:=0;
     countMol;
end;

function SDFmini.NextTMol: TStringList;
var
   line: string;
   getprp: boolean;
begin
     fmol.Clear;
     getprp:=false;
     fprpval:='?';
     repeat
           readln(ffle,line);
           fmol.Add(line);
           if getprp then
           begin
                fprpval:=line;
                getprp:=false;
           end;
           if (AnsiContainsStr(line,fprptag)) then getprp:=true;
     until(eof(ffle) or (line='$$$$'));
     inc(fMolIndex);
     if (eof(ffle)) then ffend:=true;
     Result:=fmol;
end;

//Reads an sdf and puts the first record into a String
function SDFmini.NextTMolFirstString: string;
var
   line: string;
   getprp: boolean;
begin
     fmolString:='';
     getprp:=false;
     fprpval:='?';
     repeat
           readln(ffle,line);
           fmolString:=fmolString+line;
           if getprp then
           begin
                fprpval:=line;
                getprp:=false;
           end;
           if (AnsiContainsStr(line,fprptag)) then getprp:=true;
     until(eof(ffle) or (line='$$$$'));
     inc(fMolIndex);
     if (eof(ffle)) then ffend:=true;
     Result:=fmolString;
end;

procedure SDFmini.countMol;
var
   i,cnt : integer;
begin
     cnt:=fMolIndex;
     Clear;
     repeat
           NextTMol;
     until (ffend);
     fNb_Mol:=fMolIndex;
     //
     Clear;
     for i:=1 to cnt do NextTMol;
end;

end.

