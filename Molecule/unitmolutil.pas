unit unitmolutil;

{$mode objfpc}{$H+}

interface

uses {$IFDEF UNIX} {$IFDEF UseCThreads}
  cthreads, {$ENDIF} {$ENDIF}
  Classes,
  SysUtils,
  CustApp,
  Math,
  { you can add units after this }
  U_BASE_GRAPHES,
  U_GRAPHES,
  unitatomandbondtype,
  unitmoleculebase,
  U_TYPE_GRAPHES;

type

  TTreeNode = record
    Depth, Degree: integer;
    Data: array of integer;
  end;

  //AoI = array of integer;
  //AoII = array of array of integer;
  AoNode = array of TTreeNode;

  { TTree }
  TTreeException = class(Exception);

  TTree = class(TObject)
  private
    fParent: iarray;
    fNodes: AoNode;
    fMxSze: integer;
    fMxData: integer;
    finode: integer;
    procedure SetMxSze(n: integer);
  public
    constructor Create;
    constructor Create(MxData: integer);
    destructor Destroy; override;
    procedure Clear;
    function AddChild(parent: integer): integer;
    function AddChild(parent: integer; MxD: integer): integer;
    property Parent: iarray read fParent write fParent;
    property Nodes: AoNode read fNodes write fNodes;
    property MxSze: integer read fMxSze write SetMxSze;
    property MxData: integer read fMxData write fMxData;
    property inode: integer read finode write finode;
  end;

  procedure SubstrucSearch1(molA: TMoleculeBase; molB: TMoleculeBase;
  var AtmIdx: iarray; out bSS: boolean);
  procedure SubstrucSearch2(molA: TMoleculeBase; molB: TMoleculeBase;
  var AtmIdx: iarray; out bSS: boolean);

implementation

  procedure SubstrucSearch1(molA: TMoleculeBase; molB: TMoleculeBase;
  var AtmIdx: iarray; out bSS: boolean);
  //Procedure of Ullman, J Assoc Comput Mach 1976, 23:31–42.
  //A, substructure to search for in B.
  //Array AtmIdx must have length nA.
  //AtmIdx is indexed at position 1.
  var
    M, Md: aiarray;
    F: iarray;
    nA, nB: integer;
    i, j, k, l: integer;
    k1, TBlast: integer;
    at1, bt1: integer;
    AofTB: iarray;
    PBdA, PBdB: PRBond;
    Aa, Bb: TNodeInfo;
    lastA, lastB: Node;
    bEnded, bExit: boolean;
    cpt: integer;
  const
    MxSzeC = 5000000;

    procedure displayM(M: aiarray; nA: integer; nB: integer);
    var
      i, j: integer;
    begin
      writeln('*******M*******');
      for i := 1 to nA do
      begin
        for j := 1 to nB do
          Write(IntToStr(M[i, j]) + ' ');
        writeln;
      end;
      writeln('***************');
    end;

  begin
    bSS:= False;
    nA := molA.p_NX;
    nB := molB.p_NX;
    //Don't set this array here, the same array can be used multiple times.
    //SetLength(AtmIdx,nA+1);//Set size of output list of atom match
    SetLength(M, nA + 1, nB + 1);//Compatibility matrix
    SetLength(Md, nA + 1, nB + 1);
    SetLength(AofTB, molB.p_MaxOutDeg);//Labels of output nodes
    SetLength(F, nA + 1);//
    //The outer degree of the query cannot be larger than the one of the molecule
    for i := 1 to nA do
      for j := 1 to nB do
        if (molA.p_OutDeg[i] <= molB.p_OutDeg[j]) then
          M[i, j] := 1
        else
          M[i, j] := 0;
    //Test chemistry
    for i := 1 to nA do
    begin
      for j := 1 to nB do
      begin//Test nature of atoms
        if M[i, j] = 1 then
        begin
          if molA.AtmSet[i]^.S <> molB.AtmSet[j]^.S then
            M[i, j] := 0;
          //Test nature of bonds
          l := 0;
          for k := molB.p_HEAD[j] to molB.p_HEAD[j + 1] - 1 do
          begin
            AofTB[l] := molB.BndSet[k]^.B;
            Inc(l);
          end;
          TBlast := l - 1;
          for k := molA.p_HEAD[i] to molA.p_HEAD[i + 1] - 1 do
          begin
            l := 0;
            while (l <= TBlast) and (molA.BndSet[k]^.B <> AofTB[l]) do
              Inc(l);
            if molA.BndSet[k]^.B = AofTB[l] then
            begin
              for k1 := l to TBlast - 1 do
                AofTB[k1] := AofTB[k1 + 1];
              Dec(TBlast);
            end
            else
              M[i, j] := 0;
          end;
        end;
      end;
    end;
    //displayM(M, nA, nB);
    for i := 1 to nA do
      for j := 1 to nB do
        Md[i, j] := M[i, j];
    for i := 0 to nA do
      F[i] := 0;

    at1 := 1;
    cpt := 0;
    bExit := False;
    while (at1 > 0) and (not bExit) do
    begin
      Inc(cpt);
      //writeln('Level: ' + IntToStr(at1) + ' Iter: ' + IntToStr(cpt));
      bt1 := F[at1] + 1;
      bEnded := True;
      if bt1 <= nB then
      begin
        while (Md[at1, bt1] = 0) and (bt1 < nB) do
          Inc(bt1);
        if (Md[at1, bt1] = 1) and (at1 <= nA) then
        begin
          bEnded := False;
          //Keep only one 1. Record its position
          for j := 1 to nB do
            if (j <> bt1) then
              Md[at1, j] := 0;
          //Go down one level in the tree
          F[at1] := bt1;
          for i := 1 to nA do
            if (i <> at1) then
              Md[i, bt1] := 0;

          molA.GetSucc(at1, Aa, lastA);
          molB.GetSucc(bt1, Bb, lastB);
          i := 1;
          bSS := True;
          while (i <= lastA) and bSS do
          begin
            bSS := False;
            for j:=1 to lastB do
            begin
              if Md[Aa[i], Bb[j]] = 1 then
              begin
                PBdA := molA.FindBond(at1, Aa[i]);
                PBdB := molB.FindBond(bt1, Bb[j]);
                if (PBdA <> nil) and (PBdB <> nil) then
                  if PBdA^.B = PBdB^.B then
                    bSS := True;
              end;
            end;
            {j := 1;
            while (j < lastB) and (Md[Aa[i], Bb[j]] = 0) do
              Inc(j);
            if Md[Aa[i], Bb[j]] = 1 then
            begin
              PBdA := molA.FindBond(at1, Aa[i]);
              PBdB := molB.FindBond(bt1, Bb[j]);
              if (PBdA <> nil) and (PBdB <> nil) then
                if PBdA^.B = PBdB^.B then
                  bSS := True;}
              {if bSS = True then
              begin
                for k := j to LastB do
                  Bb[k] := Bb[k + 1];
                LastB := LastB - 1;
              end;}
            //end;
            {if bSS = False then
            begin
              writeln('%%%Not SS%%%');
              displayM(Md, nA, nB);
            end;}
            Inc(i);
          end;
          if bSS = False then
            bEnded := True;

          Inc(at1);
          {writeln('%%%down%%%');
          displayM(M, nA, nB);
          displayM(Md, nA, nB);}
        end;
        if (at1 > nA) and (not bEnded) then
        begin
          {writeln('%%%leaf%%%');
          displayM(Md, nA, nB);}
          for j := 1 to nB do
          begin
            k := 1;
            while (F[k] <> j) and (k < at1) do
              Inc(k);
            if (F[k] = j) and (k < at1) then
              Md[nA, j] := 0
            else
              Md[nA, j] := M[nA, j];
          end;
          at1 := nA;
          if bSS = True then
          begin
            bExit := True;
            //writeln('///Stop here\\\');
          end;
        end;
      end;
      if (bEnded) then
      begin
        //Expoloration of this node is finished; go up
        for i := at1 to nA do //reinitialize the cursor
          F[at1] := 0;
        Dec(at1);
        for i := at1 to nA do //reinit Md using M and already positioned 1.
          for j := 1 to nB do
          begin
            k := 1;
            while (F[k] <> j) and (k < at1) do
              Inc(k);
            if (F[k] = j) and (k < at1) then
              Md[i, j] := 0
            else
              Md[i, j] := M[i, j];
          end;
        {writeln('%%%up%%%');
        displayM(M, nA, nB);
        displayM(Md, nA, nB);}
      end;
    end;
    for i:=1 to nA do
        AtmIdx[i]:=F[i];
  end;

  procedure SubstrucSearch2(molA: TMoleculeBase; molB: TMoleculeBase;
  var AtmIdx: iarray; out bSS: boolean);
  //Procedure of Ullman, J Assoc Comput Mach 1976, 23:31–42.
  //A, substructure to search for in B.
  //Array AtmIdx must have length nA.
  var
    M, Md: aiarray;
    F: iarray;
    nA, nB: integer;
    i, j, k, l: integer;
    k1, TBlast: integer;
    at1, bt1: integer;
    AofTB: iarray;
    PBdA, PBdB: PRBond;
    Aa, Bb: TNodeInfo;
    lastA, lastB: Node;
    bEnded, bExit: boolean;
    cpt: integer;
  const
    MxSzeC = 5000000;

    procedure displayM(M: aiarray; nA: integer; nB: integer);
    var
      i, j: integer;
    begin
      writeln('*******M*******');
      for i := 1 to nA do
      begin
        for j := 1 to nB do
          Write(IntToStr(M[i, j]) + ' ');
        writeln;
      end;
      writeln('***************');
    end;

  begin
    bSS:= False;
    nA := molA.p_NX;
    nB := molB.p_NX;
    SetLength(M, nA + 1, nB + 1);
    SetLength(Md, nA + 1, nB + 1);
    SetLength(AofTB, molB.p_MaxOutDeg);
    SetLength(F, nA + 1);
    //The outer degree of the query cannot be larger than the one of the molecule
    for i := 1 to nA do
      for j := 1 to nB do
        if (molA.p_OutDeg[i] <= molB.p_OutDeg[j]) then
          M[i, j] := 1
        else
          M[i, j] := 0;
    //Test chemistry
    for i := 1 to nA do
    begin
      for j := 1 to nB do
      begin//Test nature of atoms
        if M[i, j] = 1 then
        begin
          if MolA.AtmSet[i]^.S='A' then
             M[i,j]:=1
          else if MolA.AtmSet[i]^.S='X'then
               if MolB.AtmSet[j]^.S='F'then
                  M[i,j]:=1
               else if MolB.AtmSet[j]^.S='Cl' then
                  M[i,j]:=1
               else if MolB.AtmSet[j]^.S='Br' then
                  M[i,j]:=1
               else if MolB.AtmSet[j]^.S='I' then
                  M[i,j]:=1
               else if MolB.AtmSet[j]^.S='X' then
                  M[i,j]:=1
               else
                  M[i,j]:=0
          else if molA.AtmSet[i]^.S <> molB.AtmSet[j]^.S then
            M[i, j] := 0;
          //Test nature of bonds
          l := 0;
          for k := molB.p_HEAD[j] to molB.p_HEAD[j + 1] - 1 do
          begin
            AofTB[l] := molB.BndSet[k]^.B;
            Inc(l);
          end;
          TBlast := l - 1;
          for k := molA.p_HEAD[i] to molA.p_HEAD[i + 1] - 1 do
          begin
            l := 0;
            while (l <= TBlast) and (molA.BndSet[k]^.B <> AofTB[l]) do
              Inc(l);
            if molA.BndSet[k]^.B = AofTB[l] then
            begin
              for k1 := l to TBlast - 1 do
                AofTB[k1] := AofTB[k1 + 1];
              Dec(TBlast);
            end
            else
              M[i, j] := 0;
          end;
        end;
      end;
    end;
    //displayM(M, nA, nB);
    for i := 1 to nA do
      for j := 1 to nB do
        Md[i, j] := M[i, j];
    for i := 0 to nA do
      F[i] := 0;

    at1 := 1;
    cpt := 0;
    bExit := False;
    while (at1 > 0) and (not bExit) do
    begin
      Inc(cpt);
      //writeln('Level: ' + IntToStr(at1) + ' Iter: ' + IntToStr(cpt));
      bt1 := F[at1] + 1;
      bEnded := True;
      if bt1 <= nB then
      begin
        while (Md[at1, bt1] = 0) and (bt1 < nB) do
          Inc(bt1);
        if (Md[at1, bt1] = 1) and (at1 <= nA) then
        begin
          bEnded := False;
          //Keep only one 1. Record its position
          for j := 1 to nB do
            if (j <> bt1) then
              Md[at1, j] := 0;
          //Go down one level in the tree
          F[at1] := bt1;
          for i := 1 to nA do
            if (i <> at1) then
              Md[i, bt1] := 0;

          molA.GetSucc(at1, Aa, lastA);
          molB.GetSucc(bt1, Bb, lastB);
          i := 1;
          bSS := True;
          while (i <= lastA) and bSS do
          begin
            bSS := False;
            for j:=1 to lastB do
            begin
              if Md[Aa[i], Bb[j]] = 1 then
              begin
                PBdA := molA.FindBond(at1, Aa[i]);
                PBdB := molB.FindBond(bt1, Bb[j]);
                if (PBdA <> nil) and (PBdB <> nil) then
                  if PBdA^.B = PBdB^.B then
                    bSS := True;
              end;
            end;
            {j := 1;
            while (j < lastB) and (Md[Aa[i], Bb[j]] = 0) do
              Inc(j);
            if Md[Aa[i], Bb[j]] = 1 then
            begin
              PBdA := molA.FindBond(at1, Aa[i]);
              PBdB := molB.FindBond(bt1, Bb[j]);
              if (PBdA <> nil) and (PBdB <> nil) then
                if PBdA^.B = PBdB^.B then
                  bSS := True;}
              {if bSS = True then
              begin
                for k := j to LastB do
                  Bb[k] := Bb[k + 1];
                LastB := LastB - 1;
              end;}
            //end;
            {if bSS = False then
            begin
              writeln('%%%Not SS%%%');
              displayM(Md, nA, nB);
            end;}
            Inc(i);
          end;
          if bSS = False then
            bEnded := True;

          Inc(at1);
          {writeln('%%%down%%%');
          displayM(M, nA, nB);
          displayM(Md, nA, nB);}
        end;
        if (at1 > nA) and (not bEnded) then
        begin
          {writeln('%%%leaf%%%');
          displayM(Md, nA, nB);}
          for j := 1 to nB do
          begin
            k := 1;
            while (F[k] <> j) and (k < at1) do
              Inc(k);
            if (F[k] = j) and (k < at1) then
              Md[nA, j] := 0
            else
              Md[nA, j] := M[nA, j];
          end;
          at1 := nA;
          if bSS = True then
          begin
            bExit := True;
            //writeln('///Stop here\\\');
          end;
        end;
      end;
      if (bEnded) then
      begin
        //Expoloration of this node is finished; go up
        for i := at1 to nA do //reinitialize the cursor
          F[at1] := 0;
        Dec(at1);
        for i := at1 to nA do //reinit Md using M and already positioned 1.
          for j := 1 to nB do
          begin
            k := 1;
            while (F[k] <> j) and (k < at1) do
              Inc(k);
            if (F[k] = j) and (k < at1) then
              Md[i, j] := 0
            else
              Md[i, j] := M[i, j];
          end;
        {writeln('%%%up%%%');
        displayM(M, nA, nB);
        displayM(Md, nA, nB);}
      end;
    end;
    for i:=1 to nA do
        AtmIdx[i-1]:=F[i];
  end;
  { TTree }

  procedure TTree.SetMxSze(n: integer);
  var
    i: integer;
  begin
    fMxSze := n;
    SetLength(fParent, fMxSze);
    SetLength(fNodes, fMxSze);
    for i := Low(fParent) to High(fParent) do
    begin
      fParent[i] := 0;
      fNodes[i].Depth := 0;
      fNodes[i].Degree := 0;
    end;
    if fMxSze >= 1 then
    begin
      SetLength(fNodes[0].Data, fMxData);
    end;
  end;

  constructor TTree.Create;
  begin
    fMxSze := 0;
    fMxData := 0;
    finode := 0;
    inherited Create;
  end;

  constructor TTree.Create(MxData: integer);
  begin
    fMxSze := 0;
    finode := 0;
    fMxData := MxData;
    inherited Create;
  end;

  destructor TTree.Destroy;
  var
    i: integer;
  begin
    for i := Low(fNodes) to High(fNodes) do
      SetLength(fNodes[i].Data, 0);
    SetMxSze(0);
    inherited Destroy;
  end;

  procedure TTree.Clear;
  var
    i: integer;
  begin
    fMxSze := 0;
    finode := 0;
    for i := Low(fNodes) to High(fNodes) do
      SetLength(fNodes[i].Data, 0);
    SetMxSze(0);
  end;

  function TTree.AddChild(parent: integer): integer;
  begin
    Inc(finode);
    if (inode > fMxSze) then
    begin
      raise TTreeException.Create('BUG: Too many combinations to scan!' +
        ' The program will crash.');
      Halt;
    end;
    Inc(fNodes[parent].Degree);
    fParent[finode] := parent;
    fNodes[finode].Depth := fNodes[parent].Depth + 1;
    fNodes[finode].Degree := 0;
    SetLength(fNodes[finode].Data, fMxData);
    Result := finode;
  end;

  function TTree.AddChild(parent: integer; MxD: integer): integer;
  begin
    Inc(finode);
    if (inode > fMxSze) then
    begin
      raise TTreeException.Create('BUG: Too many combinations to scan!' +
        ' The program will crash.');
      Halt;
    end;
    Inc(fNodes[parent].Degree);
    fParent[finode] := parent;
    fNodes[finode].Depth := fNodes[parent].Depth + 1;
    fNodes[finode].Degree := 0;
    SetLength(fNodes[finode].Data, MxD);
    Result := finode;
  end;

end.

