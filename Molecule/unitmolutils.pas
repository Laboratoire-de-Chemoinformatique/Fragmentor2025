unit UnitMolUtils;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, U_GRAPHES, U_TYPE_GRAPHES, U_BASE_GRAPHES, matrix;


type
   iarray = Array of Integer;
   aiarray = Array of iarray;
   //Array of matrices
   marray = Array of aiarray;


    { MolUtils /!\ EN CONSTRUCTION}

    MolUtils=Class(TObject)

    public
          constructor Create;
          function MultiplyMat(M1,M2: aiarray) : aiarray;
          function TransposeMat(M1: aiarray) : aiarray;
          function Isomorphism(ga,gb:T_GRAPHE_MATRICIEL; M1: aiarray) : Integer;
          function LineOk(M1: aiarray; x: Integer) : Integer;
          procedure SubstrucSearch(ga,gb:T_GRAPHE_MATRICIEL);
    end;

implementation

{ MolUtils }

constructor MolUtils.Create;
begin
     inherited Create;
end;

//Multiply 2 matrices
function MultiplyMat(M1,M2: aiarray) : aiarray;
var
     MR : aiarray;
     r1,r2,c1,c2 : integer;
     s : String;
     i,j,k,res : Integer;
begin
     //Sizes M1
     r1:= High(M1);
     c1:= High(M1[r1]);

     //Sizes M2
     r2:= High(M2);
     c2:= High(M2[r2]);

     SetLength(MR, r1+1, c2+1);

     for i:=0 to r1 do begin
      for j:=0 to c2 do begin
       res:=0;
       for k:=0 to c1 do begin
        res:= res+(M1[i][k]*M2[k][j]);
       end;
       MR[i][j]:=res;
      end;
     end;

     Result:=MR;
end;

//Transpose a matrix
function TransposeMat(M1: aiarray) : aiarray;
var
     MR : aiarray;
     r1,c1 : integer;
     s : String;
     i,j,k,res : Integer;
begin
     //Sizes M1
     r1:= High(M1);
     c1:= High(M1[r1]);

     SetLength(MR, c1+1, r1+1);

     for i:=0 to r1 do begin
      for j:=0 to c1 do begin
       MR[j][i]:=M1[i][j];
      end;
     end;

     Result:=MR;
end;

//Isomorphism
function Isomorphism(ga,gb: T_GRAPHE_MATRICIEL; M1: aiarray) : Integer;
var
     //Matrice de resultats
     Iso : Integer;
     C,MR,MR2 : aiarray;
     ra,ca : integer;
     s : String;
begin
     Iso:=1;

     //Operations
     MR:=gb.MultiplyWithMat1(M1);
     MR2:=TransposeMat(MR,M);
     C:=MultiplyMat(M1,MR2,M);

     //Check if its ok
     Iso:=ga.CompareWithMat(C,M);

     Result:=Iso;
end;

//Checking a specific line in a matrix
function LineOk(M1: aiarray; x: Integer) : Integer;
var
     LineFree,i : Integer;
     s : String;
begin
     LineFree:=1;

     for i:=0 to High(M1) do begin
      if (M1[i][x]=1) then LineFree:=0;
     end;

     Result:=LineFree;
end;

//Ullmann
Procedure SubstrucSearch(ga,gb:T_GRAPHE_MATRICIEL);
var
     M0,Mat,MV : aiarray;
     s : String;
     i,j,d,pa,pb,x,k,num,Iso : Integer;
     Back, Found: Boolean;
     F,H: Array of integer;
     //ga, gb: T_GRAPHE_MATRICIEL;
begin
     //Initialisations
     d:=0;
     num:=0;
     Back:=True;
     Found:=False;

     //Read the adjacency matrices
     {ga:= T_GRAPHE_MATRICIEL.CREATE;
     gb:= T_GRAPHE_MATRICIEL.CREATE;
     ga.ReadMatrix('DATA/GA2.MAT');
     gb.ReadMatrix('DATA/GB2.MAT');
     ga.AFFMatrix(M,'Matrice A pour Galpha',78,99);
     gb.AFFMatrix(M,'Matrice B pour Gbeta',78,99);}

     //Look for the adjacency matrices sizes
     pa:= ga.p_LastCol;
     pb:= gb.p_LastCol;

     M0:= CreateMat0(ga,gb,pa,pb,M);

     //Initialisation of Mat, MV, H and F
     SetLength(F,pb);
     for i:=0 to High(F) do F[i]:=0;

     SetLength(MV,pa,pb);
     for i:=0 to High(MV) do begin
      for j:=0 to High(MV[i]) do begin
       MV[i][j]:=0;
      end;
     end;

     SetLength(H,pa);
     for i:=0 to High(H) do H[i]:=-1;

     SetLength(Mat,pa,pb);
     //Mat:=M0;
     for i:=0 to High(M0) do begin
      for j:=0 to High(M0[i]) do begin
       Mat[i][j]:=M0[i][j];
      end;
     end;

     //------------------BEGINNING--------------------------------
    while d<>-1 do begin
     Back:=False;
     Found:=False;
     M.Lines.Add('--------d:'+IntTOStr(d)+'----------');

     //If we are checking the last col
     if d=pa-1 then begin
      for x:=0 to pb-1 do begin
       //Check every col that has not already been checked
       //If an x is found, keep it and put the corresponding line of Mat to zero
       //if (MV[d][x]<>1) then begin
       if (LineOk(MV, x, M)=1) then begin
         if (Mat[d][x]<>0) and (F[x]<>1) then begin
         Found:=True;
         Back:=False;
         H[d]:=x;
         F[x]:=1;
         for i:=0 to High(Mat[d]) do begin
          if i<>x then Mat[d][i]:=0
          else Mat[d][i]:=1;
         end;
         MV[d][x]:=1;

         //Check condition and keep H if it's ok
         //------------------------------------
         Iso:=isomorphism(ga,gb,Mat,M);
         M.Lines.Add('Iso : '+IntTOStr(Iso));

         if Iso=1 then begin
          //By now just keep H
          for i:=0 to High(H) do begin
           s:= 'H AT THE END['+IntTOStr(i)+']='+IntTOStr(H[i]);
           M.Lines.Add(s);
          end;
          M.Lines.Add('\o/');
          M.Lines.Add('Mat at the end');
          for i:=0 to High(Mat) do begin
           for j:=0 to High(Mat[i]) do begin
            s:= 'Mat['+IntTOStr(i)+','+IntTOStr(j)+']='+IntTOStr(Mat[i,j]);
            M.Lines.Add(s);
          end;
         end;
        end;
        //------------------------------------

         //Put the line back for further investigations (Mat[d]:=M0[d])
         for i:=0 to High(Mat[d]) do begin
          Mat[d][i]:=M0[d][i];
         end;
        end;
       end;
      end;
     Back:=True;
     end

     //We are not at the end
     else begin
      for x:=0 to pb-1 do begin
       //Check every col that has not already been checked
       //If an x is found, keep it and put the corresponding line of Mat to zero
       //if MV[d][x]<>1 then begin
       if (LineOk(MV, x, M)=1) then begin
        if (Mat[d][x]<>0) and (F[x]<>1) then begin
         Found:=True;
         Back:=False;
         H[d]:=x;
         F[x]:=1;
         for i:=0 to High(Mat[d]) do begin
          if i<>x then Mat[d][i]:=0
          else Mat[d][i]:=1;
         end;
         MV[d][x]:=1;
         Break;
        end;
       end;
      end;
     end;

     //Stop the search if we don't find anything
     //Go back if all cols have been searched
     if Found=False then begin
      M.Lines.Add('rien trouve');
      Back:=True;
      //If it is the first line we are searching on
      if d=0 then Break;
     end;

     //Check if we go down or back up
     if Back=False then begin
      M.Lines.Add('Not Back');
      Found:=False;
      d:=d+1;
     end

     else begin
      M.Lines.Add('Back');
      Found:=False;

      //Put everything back for the rollback
      //MV
      for i:=d to High(MV) do begin
       for j:=0 to High(MV[i]) do begin
        MV[i][j]:=0;
       end;
      end;

      //F
      for i:=0 to High(F) do begin
       F[i]:=0;
      end;
      if (d-2)>=0 then begin
       k:= High(MV[d-2]);
       while k<>-1 do begin
        if MV[d-2][k]=1 then begin
         F[k]:=1;
         Break;
        end
        else
         k:=k-1;
       end;
      end;

      //M0 and Mat
      for i:=d-1 to High(M0) do begin
       for j:=0 to High(M0[i]) do begin
        Mat[i][j]:=M0[i][j];
       end;
      end;

      //H
      for i:=d-1 to High(M0) do begin
       H[i]:=-1;
      end;

      //Go back
      d:=d-1;
     end;

     num:=num+1;
     end;

     FreeAndNil(ga);
     FreeAndNil(gb);
end;




end.

