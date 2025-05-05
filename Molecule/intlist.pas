 {*****************************************************************//}
 {                                                                 //}
 {  TIntList Class                                                 //}
 {  Copyright? BrandsPatch LLC                                     //}
 {  http://www.explainth.at                                        //}
 {                                                                 //}
 {  All Rights Reserved                                            //}
 {                                                                 //}
 {  Permission is granted to use, modify and redistribute          //}
 {  the code in this Delphi unit on the condition that this        //}
 {  notice is retained unchanged.                                  //}
 {                                                                 //}
 {  BrandsPatch  declines all responsibility for any losses,       //}
 {  direct or indirect, that may arise  as a result of using       //}
 {  this code.                                                     //}
 {                                                                 //}
 {*****************************************************************//}
unit IntList;

{$mode objfpc}{$H+}{$R+}

interface

uses SysUtils, Classes;

type
  WordInt = smallint;

type
  TIntList = class(TObject)
  private
    FList: TList;
    function GetCount: integer;
    function GetInteger(Index: integer): integer;
    procedure SetInteger(Index, Value: integer);
  public
    property IntCount: integer Read GetCount;
    property Integers[Index: integer]: integer Read GetInteger Write SetInteger; default;
    constructor Create;
    destructor Destroy; override;
    procedure ReadFromStream(AStream: TStream);
    procedure WriteToStream(AStream: TStream);
    function Add(Value: integer): integer;
    procedure ClearEx(ACount: integer);
    function Decrement(Index: integer): boolean;
    procedure Delete(Index: integer);
    function Discard(Value: integer): integer;
    procedure Exchange(Index1, Index2: integer);
    function Find(Value: integer): integer;
    procedure Increment(Index: integer);
    function Insert(Index, Value: integer): integer;
    procedure Pack;
    function Remove(Index: integer): boolean;
    procedure Replicate(AList: TIntList);
    procedure Merge(AList: TIntList);
    procedure Sort;
  end;

implementation

constructor TIntList.Create;
begin
  inherited Create;
  FList := TList.Create;
end;

destructor TIntList.Destroy;
begin
  FList.Free;
  inherited;
end;


procedure TIntList.Merge(AList: TIntList);
var
  i: integer;
begin
  for i := 0 to AList.IntCount - 1 do
    if (Find(Alist[i]) = -1) then
      Add(AList[i]);
end; { TIntList }


procedure TIntList.ReadFromStream(AStream: TStream);
var
  ACount: WordInt;
  i, j:   integer;
begin
  with FList, AStream do
  begin
    Clear; {empty the current contents of the list}
    ReadBuffer(ACount, sizeof(ACount));
    {first read the number of entries in the list}
    Capacity := ACount;
    for i := 1 to ACount do
    begin
      ReadBuffer(j, sizeof(integer));
      Add(Pointer(j));
      {read in each entry and add it to the list}
    end;
  end;
end;

procedure TIntList.WriteToStream(AStream: TStream);
var
  i: WordInt;
  j: integer;
begin
  i := FList.Count;
  with AStream, FList do
  begin
    i := Count;
    WriteBuffer(i, sizeof(i));
    {first write the number of entries in the list}
    for i := 0 to Count - 1 do
    begin
      j := integer(Items[i]);
      WriteBuffer(j, sizeof(j));
      {then write out each entry}
    end;
  end;
end;

function TIntList.Add(Value: integer): integer;
begin
  Result := FList.Add(Pointer(Value));
end;

procedure TIntList.ClearEx(ACount: integer);
begin
  with FList do
  begin
    Clear;
    Capacity := ACount;
  end;
end;

function TIntList.Decrement(Index: integer): boolean;
begin
  with FList do
    if ((Index >= 0) and (Index < Count)) then
    begin
      Items[Index] := Pointer(integer(Items[Index]) - 1);
      Result := (integer(FList[Index]) = 0);
    end
    else
      raise EListError.Create('Invalid List Index');
end;

procedure TIntList.Delete(Index: integer);
begin
  with FList do
  begin
    if (Index < 0) or (Index >= Count) then
      exit;
    Delete(Index);
  end;
end;

function TIntList.Discard(Value: integer): integer;
begin
  Result := Find(Value);
  if (Result >= 0) then
    Delete(Result);
end;

procedure TIntList.Exchange(Index1, Index2: integer);
begin
  FList.Exchange(Index1, Index2);
end;

function TIntList.Find(Value: integer): integer;
begin
  with FList do
    for Result := 0 to Count - 1 do
      if (integer(Items[Result]) = Value) then
        exit;
  Result := -1;
end;

procedure TIntList.Increment(Index: integer);
begin
  with FList do
    if ((Index >= 0) and (Index < Count)) then
      Items[Index] := Pointer(integer(Items[Index]) + 1);
end;

function TIntList.Insert(Index, Value: integer): integer;
begin
  with FList do
    if ((Index >= 0) and (Index < Count)) then
    begin
      Result := Index;
      Insert(Index, Pointer(Value));
    end
    else
      Result := Add(Pointer(Value));
end;

procedure TIntList.Pack;
begin
  FList.Pack;
end;

function TIntList.Remove(Index: integer): boolean;
begin
  with FList do
    if ((Index >= 0) and (Index < Count)) then
    begin
      Delete(Index);
      Result := True;
    end
    else
      Result := False;
end;

procedure TIntList.Replicate(AList: TIntList);
var
  AStream: TStream;
begin
  AStream := TMemoryStream.Create;
  with AStream do
    try
      AList.WriteToStream(AStream);
      Position := 0;
      ReadFromStream(AStream);
    finally
      Free
    end;
  {replicate the contents of AList in the current list}
end;


function SortIntegers(P1, P2: Pointer): integer;
var
  i: integer absolute P1;
  j: integer absolute P2;
begin
  Result := (i - j); {return -ve, 0 or +ve}
end;

procedure TIntList.Sort;
begin
  FList.Sort(@SortIntegers);
end;

function TIntList.GetCount: integer;
begin
  Result := FList.Count;
end;

function TIntList.GetInteger(Index: integer): integer;
begin
  with FList do
    if ((Index >= 0) and (Index < Count)) then
      Result := integer(Items[Index])
    else
      raise EListError.Create('Invalid list index.');
end;

procedure TIntList.SetInteger(Index, Value: integer);
begin
  with FList do
    if ((Index >= 0) and (Index < Count)) then
      Items[Index] := Pointer(Value)
    else
      raise EListError.Create('Invalid list index.');
end;

end.
