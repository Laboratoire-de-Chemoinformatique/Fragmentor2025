// ------------------------------------
// Version 1.0 : Algorithmes de graphes
// ------------------------------------

// Copyright (C) 2003 (Lacomme, Prins, Sevaux)

// Ce programme (ainsi que la biblioth�que objets livr�e avec le livre)
// est libre, vous pouvez le redistribuer et/ou
// le modifier selon les termes de la
// Licence Publique G�n�rale GNU publi�e par la
// Free Software Foundation (version 2 ou bien toute autre version
// ult�rieure choisie par vous).
// Ce programme est distribu� car potentiellement utile,
// mais SANS AUCUNE GARANTIE, ni explicite ni implicite,
// y compris les garanties de commercialisation ou
// d'adaptation dans un but sp�cifique. Reportez-vous � la
// Licence Publique G�n�rale GNU pour plus de d�tails.
// Vous devez avoir re�u une copie de la Licence
// Publique G�n�rale GNU en m�me temps que
// ce programme ; si ce n'est pas le cas, �crivez � la
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
// 02111-1307, �tats-Unis.

// Diff�rents sites donnent des copies officielles ou non de cette licence :
// http://www.linux-france.org/article/these/gpl.html
// http://www.gnu.org/copyleft/gpl.html


unit U_TYPE_GRAPHES;

interface

CONST

MaxNode      =   333;                    {Nb max de noeuds possibles}

MaxNodeEns   =   250;  // utilise dans le type StackNode dans
                       // la m�thode BackTrack

MaxArcNum    =  200*MaxNode;//*MaxNode+MaxNode;{Nb max d'arcs utilisables}
//Chaque noeud peut �tre connect� � 100 autres par des arcs ou � 50 autres par des c�t�s
LoopLess     =  True;                    {Constante "sans boucles"}
Simple       =  True;                    {Constante "simple"}
Layered      =  True;                    {Constante "en couches"=sans circuit}
MaxCost      =  maxsmallint;//MaxLongInt-1;              {Valeur max permise pour les co�ts}

MaxNbCosts = 3;     {Nb max de couts par graphe geres automatiquement}
// utilis� dans les buckets
UMax      = 1023;
Trace: Boolean=false;
    {Trace pour Local,SimAn,Tabu,Backtrack}
    {a n'utiliser que si le code est recompile en mode console}
    { car le mode trace utilise les anciennes instructions writeln}

// utilis� dans les chemins
MaxPathLength = 20;
MaxPathNumber = 20;

TYPE

// type permettant d'eviter la liaison avec les types dependants du compilateur

Cost         = LongInt;                  {Type cout}

Node         = 0..MaxNode+1;             {Type noeud (sommet)}

NodeEns      = 0..MaxNodeEns+1;
                        // utilise dans le type StackNode dans
                       // la m�thode BackTrack

ArcNum       = 0..MaxArcNum+1;          {Type numero d'arc (listes d'adjacence)}


// type generaux de manipulation

TNodeCost    = Array[Node] of Cost;      {Type tableau de co�ts sur noeuds}
TNodeInfo    = Array[Node] of Node;      {Type pour colorations de noeuds,etc}
PTNodeInfo   = ^TNodeInfo;
TNodeBool    = Array[Node] of Boolean;   {Type tableau de bool�ens sur noeuds}

TArcCost     = Array[ArcNum] of Cost;    {Type tableau de valuations d'arcs}
PTArcCost    = ^TArcCost;                {Type pointeur sur TArcCost}
TInverse     = Array[ArcNum] of ArcNum;  {Type tableau arc --> arc inverse}
PTInverse    = ^TInverse;                {Type pointeur sur TInverse}
TArcBool     = Array[ArcNum] of Boolean; {Type tableau de boolen sur arcs}
PTArcBool    = ^TArcBool;

// types utilises pour coder le graphe par liste de successeurs

THead        = Array[Node] of ArcNum;    {Type tableau de tetes de listes}
TSucc        = Array[ArcNum] of Node;    {Type tableau de listes d'adjacence}
PTSucc       = ^TSucc;                   {Type pointeur sur TSucc}

// types utilises pour coder le graphe sous forme de matrice

CostMatrix   = Array[Node,Node] of Cost; {Type matrice de couts}
NodeMatrix   = Array[Node] of TNodeInfo; {Type matrice de noeuds (pour Floyd)}

// utilises dans l'unite graphes
CostNb     = 0..MaxNbCosts;              {N� de valuation (de tableau cout)}
TPTArcCost = Array[CostNb] of PTArcCost; {Tableau de pointeurs sur valuations}

// utilises dans l'unit� U_BASE_GRAPHE pour la classe Buckets
BuckNo    = 0..UMax;

PathMatrix   = Array[1..MaxPathNumber,1..MaxPathLength] of Node; {Type matrice de noeuds (pour stocker des chemins courts)}

implementation

end.
