unit main;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, Menus, ExtCtrls, Spin,
  StdCtrls, Math;

type

  { TfrmMain }

  TfrmMain = class(TForm)
    cbbBasisLijn: TComboBox;
    cbbModel: TComboBox;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    MainMenu1: TMainMenu;
    miWolk: TMenuItem;
    miMandel: TMenuItem;
    miJuliab: TMenuItem;
    miMira: TMenuItem;
    miHenon: TMenuItem;
    miCollet: TMenuItem;
    miBrownl: TMenuItem;
    miPythbs: TMenuItem;
    miMondriaan: TMenuItem;
    miStof: TMenuItem;
    miStofa: TMenuItem;
    miStofbt: TMenuItem;
    miStofb: TMenuItem;
    miSterfractal: TMenuItem;
    miBoomm: TMenuItem;
    miPytht: TMenuItem;
    miPythb3: TMenuItem;
    miPythb2: TMenuItem;
    miPythb1: TMenuItem;
    miWervel: TMenuItem;
    miBolspira: TMenuItem;
    miLogspira: TMenuItem;
    miArchi: TMenuItem;
    miWikkel: TMenuItem;
    miDraak1: TMenuItem;
    miDraak0: TMenuItem;
    miDraak: TMenuItem;
    miLevy: TMenuItem;
    miKronkel: TMenuItem;
    miMink: TMenuItem;
    miKoch: TMenuItem;
    miKam: TMenuItem;
    miSier: TMenuItem;
    miBoom2: TMenuItem;
    miBoom3: TMenuItem;
    miBoomH2: TMenuItem;
    miBoomH1: TMenuItem;
    miFractals: TMenuItem;
    pbMain: TPaintBox;
    SpinEdit1: TSpinEdit;
    SpinEdit2: TSpinEdit;
    procedure FormCreate(Sender: TObject);
    procedure FormResize(Sender: TObject);
    procedure miArchiClick(Sender: TObject);
    procedure miBolspiraClick(Sender: TObject);
    procedure miBoom2Click(Sender: TObject);
    procedure miBoom3Click(Sender: TObject);
    procedure miBoomH1Click(Sender: TObject);
    procedure miBoomH2Click(Sender: TObject);
    procedure miBoommClick(Sender: TObject);
    procedure miBrownlClick(Sender: TObject);
    procedure miColletClick(Sender: TObject);
    procedure miDraak0Click(Sender: TObject);
    procedure miDraak1Click(Sender: TObject);
    procedure miDraakClick(Sender: TObject);
    procedure miHenonClick(Sender: TObject);
    procedure miJuliabClick(Sender: TObject);
    procedure miKamClick(Sender: TObject);
    procedure miKochClick(Sender: TObject);
    procedure miKronkelClick(Sender: TObject);
    procedure miLevyClick(Sender: TObject);
    procedure miLogspiraClick(Sender: TObject);
    procedure miMandelClick(Sender: TObject);
    procedure miMinkClick(Sender: TObject);
    procedure miMiraClick(Sender: TObject);
    procedure miMondriaanClick(Sender: TObject);
    procedure miPythb1Click(Sender: TObject);
    procedure miPythb2Click(Sender: TObject);
    procedure miPythb3Click(Sender: TObject);
    procedure miPythbsClick(Sender: TObject);
    procedure miPythtClick(Sender: TObject);
    procedure miSierClick(Sender: TObject);
    procedure miSterfractalClick(Sender: TObject);
    procedure miStofaClick(Sender: TObject);
    procedure miStofbClick(Sender: TObject);
    procedure miStofbtClick(Sender: TObject);
    procedure miStofClick(Sender: TObject);
    procedure miWervelClick(Sender: TObject);
    procedure miWikkelClick(Sender: TObject);
    procedure miWolkClick(Sender: TObject);
    procedure SpinEdit1Change(Sender: TObject);
    procedure SpinEdit2Change(Sender: TObject);
  private
    procedure Teken;
    procedure Clear;
    procedure MinMaxPercNegYfToSmallestFactShift(xlo,xhi,ylo,yhi,perc : double; NegYf : boolean);
    procedure MinMaxPercXPercYToFactShift(xlo, xhi, ylo, yhi, PercX,PercY: double; NegYf: boolean);
    procedure MinMaxPercXPercYNegXfNegYfToFactShift(xlo, xhi, ylo, yhi, PercX,PercY: double; NegXf,NegYf : boolean);
  public

  end;

var
  frmMain: TfrmMain;
  a, b, c, n, x, y, xf, yf, xs, ys: Double;
  j, m, p, s: Integer;
  x1, x2 , x3, x4, y1, y2, y3, y4: Array[0..4] of Double;
  EgaColor : array[0..15] of TColor =
    (TColor($000000),TColor($AA0000),TColor($00AA00),TColor($AAAA00),
     TColor($0000AA),TColor($AA00AA),TColor($0055AA),TColor($AAAAAA),
     TColor($555555),TColor($FF5555),TColor($55FF55),TColor($FFFF55),
     TColor($5555FF),TColor($FF55FF),TColor($55FFFF),TColor($FFFFFF));

implementation

uses invoer;

{$R *.lfm}

{ TfrmMain }

procedure TfrmMain.MinMaxPercNegYfToSmallestFactShift(xlo,xhi,ylo,yhi,perc : double; NegYf : boolean);
begin
  // xf * xlo + xs =      perc  * ImgW
  // xf * xhi + xs = (1-  perc) * ImgW
  // xf *(xhi-xlo) = (1-2*perc) * ImgW           // subtract to eliminate xs
  xf:=(1-2*perc)*pbMain.Width/(xhi-xlo);

  // yf * ylo + ys =      perc  * ImgH
  // yf * yhi + ys = (1 - perc) * ImgH
  // yf *(yhi-ylo) = (1-2*perc) * ImgH
  yf:=(1-2*perc)*pbMain.Height/(yhi-ylo);        // subtract to eliminate ys

  if xf<yf then yf:=xf else xf:=yf;              // use smallest factor

  if NegYf then yf:=-yf;                         // is sometimes necessary

  xs:=pbMain.Width /2-(xlo+(xhi-xlo)/2)*xf;      // middle value -> midscreen horz
  ys:=pbMain.Height/2-(ylo+(yhi-ylo)/2)*yf;      // middle value -> midscreen vert
end;

procedure TfrmMain.MinMaxPercXPercYToFactShift(xlo, xhi, ylo, yhi, PercX,
  PercY: double; NegYf: boolean);
begin
  // xf * xlo + xs =      percx  * ImgW
  // xf * xhi + xs = (1-  percx) * ImgW
  // xf *(xhi-xlo) = (1-2*percx) * ImgW
  xf:=(1-2*percx)*pbMain.Width/(xhi-xlo);        // subtract to eliminate xs

  // yf * ylo + ys =      percy  * ImgH
  // yf * yhi + ys = (1 - percy) * ImgH
  // yf *(yhi-ylo) = (1-2*percy) * ImgH
  yf:=(1-2*percy)*pbMain.Height/(yhi-ylo);       // subtract to eliminate ys
  if NegYf then yf:=-yf;
  xs:=percx*pbMain.Width - xf*xlo;               // from first x equation
  ys:=percy*pbMain.Height- yf*ylo;               // from first y equation
end;

procedure TfrmMain.MinMaxPercXPercYNegXfNegYfToFactShift(xlo, xhi, ylo, yhi,
  PercX, PercY: double; NegXf, NegYf: boolean);
begin
  // xf * xlo + xs =      percx  * ImgW
  // xf * xhi + xs = (1-  percx) * ImgW
  // xf *(xhi-xlo) = (1-2*percx) * ImgW
  xf:=(1-2*percx)*pbMain.Canvas.Width/(xhi-xlo);        // subtract to eliminate xs

  if NegXf then xf:=-xf;

  // yf * ylo + ys =      percy  * ImgH
  // yf * yhi + ys = (1 - percy) * ImgH
  // yf *(yhi-ylo) = (1-2*percy) * ImgH
  yf:=(1-2*percy)*pbMain.Canvas.Height/(yhi-ylo);       // subtract to eliminate ys

  if NegYf then yf:=-yf;

  //xs:=percx*Image1.Width - xf*xlo;                 // from first x equation, no longer valid !!!
  xs:=(pbMain.Canvas.Width-xf*(xhi+xlo))/2;               // from adding the 2 x-equations
  //ys:=percy*Image1.Height- yf*ylo;                 // from first y equation, no longer valid !!!
  ys:=(pbMain.Canvas.Height-yf*(yhi+ylo))/2;              // from adding the 2 y-equations
end;

procedure TfrmMain.miBoomH1Click(Sender: TObject);
var
  m, n, n2, p, s: Integer;
  a, am: Double;
  x, y: Array[1..8192] of Double;
begin
  Clear;
  Label4.Caption := 'Boomfractal op letter H';
  Label1.Visible := True;
  Label1.Caption := 'Orde:';
  SpinEdit1.Visible := True;
  SpinEdit1.MinValue := 1;
  SpinEdit1.MaxValue := 10;
  SpinEdit1.Increment := 1;
  xf := pbMain.Width/4;
  yf := pbMain.Height/2.82;
  xs := pbMain.Width div 2;
  ys := pbMain.Height div 2;
  p := SpinEdit1.Value;
  a := Sqrt(1/2);
  x[1] := 0;
  y[1] := 0;
  for m:=0 to p do
  begin
    s := m mod 2;
    for n:=Trunc(Power(2,m)) to Trunc(Power(2,(m+1))-1) do
    begin
      am := Power(a,m);
      n2 := n+n;
      if s=1 then
      begin
        x[n2] := x[n];
        y[n2] := y[n] + am;
        x[n2+1] := x[n];
        y[n2+1] := y[n] - am;
      end
      else
      begin
        x[n2] := x[n] + am;
        y[n2] := y[n];
        x[n2+1] := x[n] - am;
        y[n2+1] := y[n];
      end;
    end;
  end;
  n := 2;
  while n <= Trunc(Power(2,(p+2))-2) do
  begin
    pbMain.Canvas.Line(Trunc(xs+xf*x[n]),Trunc(ys+yf*y[n]),Trunc(xs+xf*x[n+1]),Trunc(ys+yf*y[n+1]));
    n := n + 2;
  end;
end;

procedure TfrmMain.miArchiClick(Sender: TObject);
var
  a, r, t, x, y: Double;
begin
  Clear;
  Label4.Caption := 'Spiraal van Archimedes';
  MinMaxPercNegYfToSmallestFactShift(-12,12,-9,9,0.05,True);
  a := 0.1;
  pbMain.Canvas.MoveTo(Round(xs),Round(ys));
  t := 0;
  while Round(t)<=Round(16*pi) do
  begin
    r := a*t;
    x := r*Cos(t);
    y := r*Sin(t);
    pbMain.Canvas.LineTo(Round(xs+xf*x),Round(ys+yf*y));
    t := t+0.1;
  end;
end;

procedure TfrmMain.miBolspiraClick(Sender: TObject);
var
  a, c, p, q, s, t, u, v, x, y, z: Double;
  n: Integer;
begin
  Clear;
  Label4.Caption := 'Bolspiraal';
  MinMaxPercNegYfToSmallestFactShift(-2,2,-1.5,1.5,0.05,True);
  a := 0.2;
  c := 0.9;
  p := 1/Sqrt(2);
  q := p*Sqrt(1-c*c);
  for n:=-500 to 500 do
  begin
    s := n*pi/50;
    t := ArcTan(a*s);
    x := Cos(s)*Cos(t);
    y := Sin(s)*Cos(t);
    z := -Sin(t);
    u := p*(y-x);
    v := c*z-q*(x+y);
    if n=-500 then
      pbMain.Canvas.MoveTo(Round(xs+xf*u),Round(ys+yf*v))
    else
      pbMain.Canvas.LineTo(Round(xs+xf*u),Round(ys+yf*v));
  end;
end;

procedure TfrmMain.miBoom2Click(Sender: TObject);
var
  h, x, y: Double;
  k, l: Integer;
begin
  Clear;
  Label4.Caption := 'Boomstructuur in tweetallig talstelsel';
  xf := pbMain.Width/4;
  yf := -pbMain.Height/2;
  xs := pbMain.Width div 2;
  ys := pbMain.Height;
  //MinMaxPercXPercYToFactShift(-2,6,0,2,0,0.05,False);
  pbMain.Canvas.Line(Trunc(xs+xf*0),Trunc(ys+yf*0),Trunc(xs+xf*0),Trunc(ys+yf*1));
  for k := 1 to 7 do
  begin
    h := Power(2,(-k+1));
    for l := 1 to Trunc(Power(2,k)) do
    begin
      x := -2+(4*l-2)*h;
      y := 2-h;
      with pbMain.Canvas do
      begin
        Line(Trunc(xs+xf*(x-h)),Trunc(ys+yf*(y+h/2)),Trunc(xs+xf*(x-h)),Trunc(ys+yf*y));
        Line(Trunc(xs+xf*(x-h)),Trunc(ys+yf*y),Trunc(xs+xf*(x+h)),Trunc(ys+yf*y));
        LineTo(Trunc(xs+xf*(x+h)),Trunc(ys+yf*(y+h/2)));
      end;
    end;
  end;
end;

procedure TfrmMain.miBoom3Click(Sender: TObject);
var
  k, l, m, n, n1, p: Integer;
  a, f, x, y: Double;
  t: Array[1..5] of Integer;
begin
  Clear;
  Label4.Caption := 'Boomstructuur in drietallig talstelsel';
  {xf := pbMain.Width/4;
  yf := pbMain.Height/4;
  xs := pbMain.Width div 2;
  ys := pbMain.Height div 2;}
  MinMaxPercNegYfToSmallestFactShift(-0.48,0.96,-0.831,0.831,0.05,False);
  //MinMaxPercXPercYToFactShift(-1.2,1.2,-0.9,0.9,0.05,0.05,False);
  p := 5;
  a := 0.5;
  for m := 0 to p do
  begin
    for n := 0 to Trunc(Power(3,m))-1 do
    begin
      n1 := n;
      for l := 1 to m do
      begin
        t[l] := n1 mod 3 ;
        n1 := Trunc(n1/3);
      end;
      x := 0;
      y := 0;
      for k := 1 to m do
      begin
        f := (2*t[k]*pi)/3;
        x := x+Cos(f)*Power(a,k);
        y := y+Sin(f)*Power(a,k);
      end;
      k := m + 1;
      with pbMain.Canvas do
      begin
        Line(Trunc(xs+xf*x),Trunc(ys+yf*y),Trunc(xs+xf*(x+Power(a,k))),Trunc(ys+yf*y));
        Line(Trunc(xs+xf*x),Trunc(ys+yf*y),Trunc(xs+xf*(x-0.5*Power(a,k))),Trunc(ys+yf*(y+Sqrt(3)/2*Power(a,k))));
        Line(Trunc(xs+xf*x),Trunc(ys+yf*y),Trunc(xs+xf*(x-0.5*Power(a,k))),Trunc(ys+yf*(y-Sqrt(3)/2*Power(a,k))));
      end;
    end;
  end;
end;

procedure TfrmMain.FormResize(Sender: TObject);
var
  i: Integer;
  IsChecked: Boolean = False;
begin
  for i := 0 to miFractals.Count-1 do
  begin
    if miFractals.Items[i].Checked then
    Begin
      IsChecked := True;
      Break;
    end;
  end;
  if IsChecked then
  case i of
    0: miBoomH1Click(Self);
    1: miBoomH2Click(Self);
    2: miBoom2Click(Self);
    3: miBoom3Click(Self);
    4: miSierClick(Self);
    5: miKamClick(Self);
    6: miKochClick(Self);
    7: miMinkClick(Self);
    8: miKronkelClick(Self);
    9: miLevyClick(Self);
    10: miDraakClick(Self);
    11: miDraak0Click(Self);
    12: miDraak1Click(Self);
    13: miWikkelClick(Self);
    14: miArchiClick(Self);
    15: miLogspiraClick(Self);
    16: miBolspiraClick(Self);
    17: miWervelClick(Self);
    18: miPythb1Click(Self);
    19: miPythb2Click(Self);
    20: miPythb3Click(Self);
    21: miPythtClick(Self);
    22: miBoommClick(Self);
    23: miSterfractalClick(Self);
    24: miStofbClick(Self);
    25: miStofbtClick(Self);
    26: miStofaClick(Self);
    27: miStofClick(Self);
    28: miMondriaanClick(Self);
    29: miPythbsClick(Self);
    30: miBrownlClick(Self);
    31: miColletClick(Self);
    32: miHenonClick(Self);
    33: miMiraClick(Self);
    34: miJuliabClick(Self);
    35: miMandelClick(self);
    36: miWolkClick(Self);
  end;
end;

procedure TfrmMain.FormCreate(Sender: TObject);
begin
  Height := 800;
  Width := 1000;
  Clear;
  Randomize;
end;

procedure TfrmMain.miBoomH2Click(Sender: TObject);
begin
  Clear;
  Label4.Caption := 'Boomfractal op letter H, backtrack methode';
  xf := pbMain.Width/2;
  yf := pbMain.Height/1.5;
  xs := pbMain.Width div 2;
  ys := pbMain.Height div 2;
  p := 5;
  a := 0.5;
  x1[0] := 0;
  y1[0] := 0;
  s := 1;
  Teken;
  for m := 1 to Trunc(Power(4,(p-1)) - 1) do
  begin
    n := m;
    s := p;
    while n mod 4 = 0 do
    begin
      n := n/4;
      s := s - 1;
      // goto 76?
    end;
    x1[s-1] := x2[s-1];
    x2[s-1] := x3[s-1];
    x3[s-1] := x4[s-1];
    y1[s-1] := y2[s-1];
    y2[s-1] := y3[s-1];
    y3[s-1] := y4[s-1];
    Teken;
  end;
end;

procedure TfrmMain.miBoommClick(Sender: TObject);
var
  m,n,p,s: Integer;
  x1,y1,x2,y2,u1,v1,u2,v2: Array[0..11] of Double;
  r1,r2,a,b,a1,a2,b1,b2,e1,e2,f1,f2,c1,c2,d1,d2,h: Double;

  procedure Gosub210;
  var
    j: Integer;
    X,Y,U,V,X3,Y3,U3,U4,V3,V4 : double;
  begin
    FOR J:=S TO P do
    begin
      X:=X1[J-1];
      Y:=Y1[J-1];
      U:=U1[J-1];
      V:=V1[J-1];
      X3:=U-X;
      Y3:=V-Y;
      X1[J]:=X+A1*X3-A2*Y3;
      Y1[J]:=Y+A2*X3-A1*Y3;
      U1[J]:=X+B1*X3-B2*Y3;
      V1[J]:=Y+B2*X3+B1*Y3;
      X2[J]:=X+E1*X3-E2*Y3;
      Y2[J]:=Y+E2*X3+E1*Y3;
      U2[J]:=X+F1*X3-F2*Y3;
      V2[J]:=Y+F2*X3+F1*Y3;
      U3:=X+C1*X3-C2*Y3;
      V3:=Y+C2*X3+C1*Y3;
      U4:=X+D1*X3-D2*Y3;
      V4:=Y+D2*X3+D1*Y3;
      IF J=S THEN
      begin
        H :=A2;
        A2:=F2;
        F2:= H;
        H :=B2;
        B2:=E2;
        E2:= H;
        H :=C2;
        C2:=D2;
        D2:= H;
      end;
      with pbMain.Canvas do
      begin
        Line(Trunc(xs+xf*X),Trunc(ys+yf*Y),Trunc(xs+xf*X1[J]),Trunc(ys+yf*Y1[J]));  // left side of stem (repeating in 1 continuous line)
        Line(Trunc(xs+xf*U1[J]),Trunc(ys+yf*V1[J]),Trunc(xs+xf*U3),Trunc(ys+yf*V3)); // lower left horz part of Z
        LineTo(Trunc(xs+xf*U4),Trunc(ys+yf*V4)); // middle vert part of Z
        LineTo(Trunc(xs+xf*x2[j]),Trunc(ys+yf*y2[j]));
        Line(Trunc(xs+xf*u2[j]),Trunc(ys+yf*v2[j]),Trunc(xs+xf*u),Trunc(ys+yf*v));
      end;
    end;
  end;

begin
  Clear;
  Label4.Caption := 'Boom van Mandelbrot, backtrackmethode';
  pbMain.Canvas.Pen.Width := 2;
  MinMaxPercNegYfToSmallestFactShift(-9.5,10.5,-3,12,0.05,True);
  p := 11;
  r1 := 0.72;
  r2 := 0.67;
  a := 3.98;
  b := 4.38;
  a1 := 0;
  a2 := a;
  b1 := 0;
  b2 := a+r1;
  e1 := 1;
  e2 := b+r2;
  f1 := 1;
  f2 := b;
  c1 := 0.5;
  c2 := b2;
  d1 := 0.5;
  d2 := e2;
  x1[0] := 0;
  y1[0] := 0;
  u1[0] := 1;
  v1[0] := 0;
  pbMain.Canvas.Line(Trunc(xs),Trunc(ys),Trunc(xs+xf),Trunc(ys));
  s := 1;
  Gosub210;
  for m := 1 to Trunc(Power(2,p-1)-1) do
  begin
    s := p;
    n := m;
    while n mod 2 = 0 do
    begin
      n := n div 2;
      s := s-1;
    end;
    h := a2;
    a2 := f2;
    f2 := h;
    h := b2;
    b2 := e2;
    e2 := h;
    h := c2;
    c2 := d2;
    d2 := h;
    x1[s-1] := x2[s-1];
    y1[s-1] := y2[s-1];
    u1[s-1] := u2[s-1];
    v1[s-1] := v2[s-1];
    Gosub210;
  end;
  pbMain.Canvas.Brush.Color := clGreen;
  pbMain.Canvas.FloodFill(Trunc(xs+xf*0.5),Trunc(ys+yf*0.5),clBlack,fsBorder);
end;

procedure TfrmMain.miBrownlClick(Sender: TObject);
var
  w,x,y: Double;
  k: Integer;
begin
  Clear;
  Label4.Caption := 'Brownse lijn';
  MinMaxPercNegYfToSmallestFactShift(-1.2,1.2,-0.9,0.9,0.05,True);
  w := 40;
  y := 0;
  pbMain.Canvas.Line(Trunc(xs+xf),Trunc(ys),Trunc(xs-xf),Trunc(ys));
  for k := 0 to 2000 do
  begin
    x := -1+k/1000;
    y := y+w*(Random-0.5)/2000;
    pbMain.Canvas.LineTo(Trunc(xs+xf*x),Trunc(ys+yf*y));
  end;
end;

procedure TfrmMain.miColletClick(Sender: TObject);
var
  k,n,w: Integer;
  a,x: Double;
begin
  Clear;
  Label4.Caption := 'Bifurcatiediagram x:=ax(1-x)';
  //MinMaxPercNegYfToSmallestFactShift(0,1,0,1,0.05,True);
  w := pbMain.Canvas.Width;
  MinMaxPercXPercYNegXfNegYfToFactShift(0,w,0,1,0.05,0.005,False,True);
  //h := pbMain.Canvas.Height;
  for n := 0 to w-1 do
  begin
    a := 2.8+1.2*n/(w-1);
    x := 0.7;
    for k := 1 to 400 do
    begin
      x := a*x*(1-x);
      if k>100 then pbMain.Canvas.Pixels[Trunc(xs+xf*n),Trunc(ys+yf*x)] := clBlack;
    end;
  end;
end;

procedure TfrmMain.miDraak0Click(Sender: TObject);
var
  d, m, n, p, s: Integer;
  f, h, x, y: Double;
begin
  Clear;
  Label4.Caption := 'Draakkromme';
  MinMaxPercNegYfToSmallestFactShift(-0.7,1.7,-1,0.8,0.05,True);
  p := 12;
  h := Power(2,-p/2);
  s := 0;
  x := h*Cos(p*pi/4);
  y := h*Sin(p*pi/4);
  pbMain.Canvas.Line(Round(xs+xf*0),Round(ys+yf*0),Round(xs+xf*x),Round(ys+yf*y));
  for n := 1 to Round(Power(2,p)-1) do
  begin
    m := n;
    repeat
      if m mod 2 = 0 then m := Round(m/2);
    until m mod 2 <> 0;
    if m mod 4 = 1 then d := 1 else d := -1;
    s := s+d;
    f := (s-p/2)*pi/2;
    x := x+h*Cos(f);
    y := y+h*Sin(f);
    pbMain.Canvas.LineTo(Round(xs+xf*x),Round(ys+yf*y));
  end;
end;

procedure TfrmMain.miDraak1Click(Sender: TObject);
var
  d, m, n, p, s: Integer;
  h, x1, x2, xa, xb, y1, y2, ya, yb: Double;
begin
  Clear;
  Label4.Caption := 'Draakkromme met afgeronde hoeken';
  MinMaxPercNegYfToSmallestFactShift(-0.7,1.7,-1.1,0.7,0.05,True);
  p := 10;
  h := Power(2,-p/2);
  s := 0;
  x1 := h*Cos(p*pi/4);
  y1 := h*Sin(p*pi/4);
  pbMain.Canvas.Line(Round(xs+xf*0),Round(ys+yf*0),Round(xs+xf*(0.75*x1)),Round(ys+yf*(0.75*y1)));
  for n:=1 to Round(Power(2,p)-1) do
  begin
    m := n;
    while m mod 2 = 0 do
      m := Round(m / 2);
    if m mod 4 = 1 then d := 1 else d := -1;
    s := (s+d) mod 4;
    x2 := x1+h*Cos((s-p/2)*pi/2);
    y2 := y1+h*Sin((s-p/2)*pi/2);
    xa := (3*x1+x2)/4;
    ya := (3*y1+y2)/4;
    xb := (x1+3*x2)/4;
    yb := (y1+3*y2)/4;
    pbMain.Canvas.LineTo(Round(xs+xf*xa),Round(ys+yf*ya));
    pbMain.Canvas.LineTo(Round(xs+xf*xb),Round(ys+yf*yb));
    x1 := x2;
    y1 := y2;
  end;
  pbMain.Canvas.LineTo(Round(xs+xf*1),Round(ys+yf*0));
end;

procedure TfrmMain.miDraakClick(Sender: TObject);
var
  d, m, n, p, s: Integer;
  a, b, h, x, y: Double;
begin
  Clear;
  Label4.Caption := 'Draakkromme met gegeven hoek';
  Label1.Caption := 'Orde';
  Label1.Visible := True;
  SpinEdit1.Visible := True;
  SpinEdit1.MinValue := 0;
  SpinEdit1.MaxValue := 9;
  SpinEdit1.Increment := 1;
  MinMaxPercNegYfToSmallestFactShift(-2.5,1.5,-2,1,0.05,False);
  p := SpinEdit1.Value;
  h := Power(2,-p/2);
  a := 1.7453;
  b := pi - a;
  x := h;
  y := 0;
  pbMain.Canvas.Line(Round(xs+xf*0),Round(ys+yf*0),Round(xs+xf*h),Round(ys+yf*0));
  s := 0;
  for n := 1 to Round(Power(2,p)-1) do
  begin
    m := n;
    while m mod 2 = 0 do m := Round(m/2);
    if m mod 4 = 1 then d := 1 else d := -1;
    s := s + d;
    x := x + h * Cos(s*b);
    y := y + h * Sin(s*b);
    pbMain.Canvas.LineTo(Round(xs+xf*x),Round(ys+yf*y));
  end;
end;

procedure TfrmMain.miHenonClick(Sender: TObject);
var
  a,b,x,y,z: Double;
  n,p: Integer;

begin
  Clear;
  Label4.Caption := 'Banen van Henons Quadratisch Systeem';
  MinMaxPercNegYfToSmallestFactShift(-1.2,1.2,-0.9,0.9,0.05,True);
  a := 0.24;
  b := Sqrt(1-Power(a,2));
  repeat
    frmInvoer.ShowModal;
    x := StrToFloat(frmInvoer.edX.Text);
    y := StrToFloat(frmInvoer.edY.Text);
    p := StrToInt(frmInvoer.edP.Text);
    //if (x=0) and (y=0) and (p=0) then Exit;
    n := 1;
    while n <= p do
    begin
      pbMain.Canvas.Pixels[Trunc(xs+xf*x),Trunc(ys+yf*y)] := clBlack;
      z := x;
      x := x*a-(y-Power(x,2))*b;
      y := z*b+(y-Power(z,2))*a;
      if Abs(x)+Abs(y)>10 then Break;
      n := n+1;
    end;
  until p=0;
end;

procedure TfrmMain.miJuliabClick(Sender: TObject);
var
  n,p: Integer;
  r,t,x,y: Double;
  x1,x2,y1,y2: Array[0..16] of Double;

  procedure gosub130;
  var
    j: Integer;
  begin
    for j := s to p do
    begin
      x := x1[j-1];
      y := y1[j-1];
      r := Sqrt(Power(x-a,2)+Power(y-b,2))/2;
      t := (x-a)/2;
      x1[j] := Sqrt(r+t);
      x2[j] := -x1[j];
      y1[j] := Sqrt(r-t)*Sign(y-b);
      y2[j] := -y1[j];
      pbMain.Canvas.Pixels[Trunc(xs+xf*x1[j]),Trunc(ys+yf*y1[j])] := clBlack;
      pbMain.Canvas.Pixels[Trunc(xs+xf*x2[j]),Trunc(ys+yf*y2[j])] := clBlack;
    end;
  end;

begin
  Clear;
  Label4.Caption := 'Julia fractal van Z:=Z^2+C';
  MinMaxPercNegYfToSmallestFactShift(-2,2,-1.5,1.5,0.05,True);
  p := 16;
  a := 0;
  b := 1;
  x1[0] := -1.3002;
  y1[0] := 0.6248;
  pbMain.Canvas.Pixels[Trunc(xs+xf*x1[0]),Trunc(ys+yf*y1[0])] := clBlack;
  s := 1;
  gosub130;
  for m := 1 to Trunc(Power(2,p-1)-1) do
  begin
    s := p;
    n := m;
    while n mod 2 = 0 do
    begin
      n := n div 2;
      s := s-1;
      x1[s-1] := x2[s-1];
      y1[s-1] := y2[s-1];
      gosub130;
    end;
  end;
end;

procedure TfrmMain.miKamClick(Sender: TObject);
var
  b1: Double;
  i, j, k, p: Integer;
  a, b: Array[0..729] of Double;
begin
  Clear;
  Label4.Caption := 'Kam van Cantor';
  xf := pbMain.Width/2;
  yf := pbMain.Height/2;
  xs := pbMain.Width div 4;
  ys := pbMain.Height div 4;
  a[0] := 0;
  a[1] := 1;
  b1 := 0.35;
  with pbMain.Canvas do
  begin
    Line(Trunc(xs+xf*0),Trunc(ys+yf*0),Trunc(xs+xf*1),Trunc(ys+yf*0));
    LineTo(Trunc(xs+xf*1),Trunc(ys+yf*-b1));
    LineTo(Trunc(xs+xf*0),Trunc(ys+yf*-b1));
    LineTo(Trunc(xs+xf*0),Trunc(ys+yf*0));
    for p := 1 to 6 do
    begin
      for i := 0 to Trunc(Power(2,p))-1 do
      begin
        b[i] := a[i]/3;
        b[i+Trunc(Power(2,p))] := 1-(1-a[i])/3;
      end;
      for j := 1 to Trunc(Power(2,p+1))-1 do
        a[j] := b[j];
      k := 0;
      while k <= Trunc(Power(2,p+1))-1 do
      begin
        Line(Trunc(xs+xf*a[k]),Trunc(ys+yf*(b1*p)),Trunc(xs+xf*a[k+1]),Trunc(ys+yf*(b1*p)));
        Line(Trunc(xs+xf*a[k]),Trunc(ys+yf*(b1*p)),Trunc(xs+xf*a[k]),Trunc(ys+yf*(b1*p-b1)));
        Line(Trunc(xs+xf*a[k+1]),Trunc(ys+yf*(b1*p)),Trunc(xs+xf*a[k+1]),Trunc(ys+yf*(b1*p-b1)));
        k := k + 2;
      end;
    end;
  end;
end;

procedure TfrmMain.miKochClick(Sender: TObject);
var
  k, l, m, n, p, s: Integer;
  h, x, y: Double;
  t: Array of Integer;
begin
  Clear;
  Label4.Caption := 'Kromme van Helge von Koch';
  xf := pbMain.Width/2;
  yf := pbMain.Height/2;
  xs := pbMain.Width div 3;
  ys := pbMain.Height div 2;
  MinMaxPercNegYfToSmallestFactShift(0,1,0,0.289,0.05,True);
  p := 4;         // keuze orde
  SetLength(t,p);
  h := Power(3,-p);
  x := 0;
  y := 0;
  pbMain.Canvas.MoveTo(Trunc(xs),Trunc(ys));
  for n := 0 to Trunc(Power(4,p) - 1) do  // notatie in viertallig stelsel
  begin
    m := n;
    for l := 0 to p - 1 do
    begin
      t[l] := m mod 4;
      m := m div 4;
    end;
    // bepaling richting lijnstuk n
    s := 0;
    for k := 0 to p - 1 do
    begin
      s := (s + (t[k] + 1) mod 3) - 1;
    end;
    // tekening lijnstuk n
    x := x + Cos(pi * s / 3) * h;
    y := y + Sin(pi * s / 3) * h;
    pbMain.Canvas.LineTo(Trunc(xs+xf*x),Trunc(ys+yf*y));
  end;
end;

procedure TfrmMain.miKronkelClick(Sender: TObject);
var
  x, y: Array [0..4096] of Double;
  a, b, c, d: Array [0..7] of Double;
  i, j, k, l, m, m1, m2, n, Pmax, u, v: Integer;
  aa, bb, xx, yy, x1, y1: Double;
begin
  Clear;
  Label4.Caption := 'Fractal kromme met basislijn en motief, Kochkruis';
  Label2.Caption := 'Basislijn';
  Label3.Caption := 'Model';
  Label2.Visible := True;
  Label3.Visible := True;
  cbbBasisLijn.Visible := True;
  cbbModel.Visible := True;
  if cbbBasislijn.ItemIndex = 1 then
    MinMaxPercNegYfToSmallestFactShift(-1.4,1.4,-1.5,1.5,0.05,False)
  else
    MinMaxPercNegYfToSmallestFactShift(-1,1,-1,1,0.05,False);
  if cbbBasislijn.ItemIndex=0 then
  begin
  // kruis anti-clockwise
    A[0]:= 1;                  //     .  x
    B[0]:= 1;                  //     .  .

    A[1]:= -1;                 //     x  .
    B[1]:= 1;                  //     .  .

    A[2]:= -1;                 //     .  .
    B[2]:= -1;                 //     x  .

    A[3]:= 1;                  //     .  .
    B[3]:= -1;                 //     .  x

    A[4]:= 1;                  //     .  x
    B[4]:= 1;                  //     .  .

    U:=4;                      // elementen basislijn
  end;

  if cbbBasislijn.ItemIndex=1 then
  begin
  // vierkant clockwise
    A[0]:= 1;                  //     .  x
    B[0]:= 1;                  //     .  .

    A[1]:= 1;                  //     .  .
    B[1]:= -1;                 //     .  x

    A[2]:= -1;                 //     .  .
    B[2]:= -1;                 //     x  .

    A[3]:= -1;                 //     x  .
    B[3]:= 1;                  //     .  .

    A[4]:= 1;                  //     .  x
    B[4]:= 1;                  //     .  .   }

    U:=4;                      // elementen basislijn
  end;

  if cbbBasislijn.ItemIndex=2 then
  begin
  // ruit anticlockwise
    A[0]:= -1;                 //       .
    B[0]:= 0;                  //    x     .
                               //       .

    A[1]:= 0;                  //       .
    B[1]:= -1;                 //    .     .
                               //       x

    A[2]:=  1;                 //       .
    B[2]:=  0;                 //    .     x
                               //       .

    A[3]:= 0;                  //       x
    B[3]:= 1;                  //    .     .
                               //       .

    A[4]:= -1;                 //      .
    B[4]:= 0;                  //   x.    .
                               //      .
    U:=4;                      // elementen basislijn
  end;

  if cbbBasislijn.ItemIndex=3 then
  begin
    // ruit clockwise
    A[0]:= -1;                 //       .
    B[0]:= 0;                  //    x     .
                               //       .

    A[1]:= 0;                  //       x
    B[1]:= 1;                  //    .     .
                               //       .

    A[2]:=  1;                 //       .
    B[2]:=  0;                 //    .     x
                               //       .

    A[3]:= 0;                  //       .
    B[3]:= -1;                 //    .     .
                               //       x

    A[4]:= -1;                 //      .
    B[4]:= 0;                  //   x.    .
                               //      .
    U:=4;                      // elementen basislijn
  end;

  if cbbBasislijn.ItemIndex=4 then
  begin
  // driehoek clockwise
    A[0]:= -0.5;               //    *
    B[0]:= 0.866;              //         .
                               //    .

    A[1]:= 1;                  //    .
    B[1]:= 0;                  //         *
                               //    .

    A[2]:= -0.5;               //    .
    B[2]:= -0.866;             //         .
                               //    *

    A[3]:= -0.5;               //    *
    B[3]:=  0.866;             //         .
                               //    .

    U:=3;                      // elementen basislijn
  end;

  if cbbBasislijn.ItemIndex=5 then
  begin
  // driehoek anti-clockwise
    A[0]:= -0.5;               //    *
    B[0]:= 0.866;              //         .
                               //    .

    A[1]:= -0.5;               //    .
    B[1]:= -0.866;             //         .
                               //    *

    A[2]:= 1;                  //    .
    B[2]:= 0;                  //         *
                               //    .


    A[3]:= -0.5;               //    *
    B[3]:=  0.866;             //         .
                               //    .

    U:=3;                      // elementen basislijn
  end;

  if cbbModel.ItemIndex=0 then
    begin
      // motief 0 naar buiten
      C[0]:= 0;
      D[0]:= 0;

      C[1]:= 0.33333;
      D[1]:= 0;

      C[2]:= 0.5;
      D[2]:= 0.2887;

      C[3]:= 0.66667;
      D[3]:= 0;

      V:=4;
      Pmax:=9;
    end;

  // motief 0 naar binnen
 { C[0]:= 0;
  D[0]:= 0;

  C[1]:= 0.33333;
  D[1]:= 0;

  C[2]:= 0.5;
  D[2]:= -0.2887;

  C[3]:= 0.66667;
  D[3]:= 0; }

  if cbbModel.ItemIndex=1 then
    begin
      // motief 1 naar buiten
      C[0]:= 0;
      D[0]:= 0;

      C[1]:= 0.4;
      D[1]:= 0.2;

      C[2]:= 0.6;
      D[2]:=-0.2;

      C[3]:= 0;
      D[3]:= 0;

      V:=3;
      Pmax:=9;
    end;

  if cbbModel.ItemIndex=2 then
    begin
      // motief 2 naar buiten  not OK
      C[0]:= 0;
      D[0]:= 0;

      C[1]:= 0.5;
      D[1]:= 0.25;

      C[2]:= 0.5;
      D[2]:=-0.25;

      C[3]:= 0;
      D[3]:= 0;

      V:=3;
      Pmax:=9;
    end;

  if cbbModel.ItemIndex=3 then
    begin
      // motief 3 naar buiten
      C[0]:= 0;
      D[0]:= 0;

      C[1]:= 0.2;
      D[1]:= 0.25;

      C[2]:= 0.8;
      D[2]:=-0.25;

      C[3]:= 0;
      D[3]:= 0;

      V:=3;
      Pmax:=9;
    end;

  if cbbModel.ItemIndex=4 then
    begin
      // motief 4 naar buiten
      C[0]:= 0;
      D[0]:= 0;

      C[1]:= 0.47;
      D[1]:= 0;

      C[2]:= 0.5;
      D[2]:= 0.288;   // was 0.47       // Driehoek+ and 0.288 interesting

      C[3]:= 0.53;
      D[3]:= 0;

      C[4]:= 0;
      D[4]:= 0;

      V:=4;
      Pmax:=9;
    end;

  if cbbModel.ItemIndex=5 then
    begin
      // motief 5 naar buiten
      C[0]:= 0;
      D[0]:= 0;

      C[1]:= 0.5;
      D[1]:= 0;

      C[2]:= 0.5;
      D[2]:= 0.33;

      C[3]:= 0.5;
      D[3]:= 0;

      C[4]:= 0;
      D[4]:= 0;

      V:=4;
      Pmax:=9;
    end;

  if cbbModel.ItemIndex=6 then
    begin
      // motief 6 naar buiten
      C[0]:= 0;
      D[0]:= 0;

      C[1]:= 0.5;
      D[1]:= 0.5;

      C[2]:= 0.5;
      D[2]:= 0.5;

      C[3]:= 0.0;
      D[3]:= 0;

      C[4]:= 0;
      D[4]:= 0;

      V:=3;
      Pmax:=9;
    end;

  if cbbModel.ItemIndex=7 then
    begin
      // motief 7 naar buiten
      C[0]:= 0;
      D[0]:= 0;

      C[1]:= 0.25;
      D[1]:= 0;

      C[2]:= 0.25;
      D[2]:= 0.25;

      C[3]:= 0.5;
      D[3]:= 0.25;

      C[4]:= 0.5;
      D[4]:= -0.25;

      C[5]:= 0.75;
      D[5]:= -0.25;

      C[6]:= 0.75;
      D[6]:= 0.0;

      C[7]:= 0;
      D[7]:= 0;

      V:=4;
      Pmax:=9;
    end;
  p := 5;
  x[0] := 0;
  y[0] := 0;
  x[Trunc(Power(v,p))] := 1;
  y[Trunc(Power(v,p))] := 0;
  for i := 0 to p-1 do
  begin
    j := 0;
    while j <= Trunc(Power(v,p)-1) do
    begin
      m1 := j + Trunc(Power(v,p-i));
      x1 := x[m1] - x[j];
      y1 := y[m1] - y[j];
      for k := 1 to v-1 do
      begin
        m2 := j + k*Trunc(Power(v,p-i-1));
        x[m2] := x1*c[k]-y1*d[k]+x[j];
        y[m2] := y1*c[k]+x1*d[k]+y[j];
      end;
      j := j + Trunc(Power(v,p-i));
    end;
  end;
  pbMain.Canvas.MoveTo(Trunc(xs+xf*a[0]),Trunc(ys+yf*b[0]));
  for m := 0 to u-1 do
  begin
    aa := a[m+1]-a[m];
    bb := b[m+1]-b[m];
    for n := 0 to Trunc(Power(v,p)) do
    begin
      xx := aa*x[n]-bb*y[n]+a[m];
      yy := bb*x[n]+aa*y[n]+b[m];
      pbMain.Canvas.LineTo(Trunc(xs+xf*xx),Trunc(ys+yf*yy));
    end;
  end;
end;

procedure TfrmMain.miLevyClick(Sender: TObject);
var
  m, n, p, s: Integer;
  a, b, h, x, y: Double;
begin
  Clear;
  Label4.Caption := 'Pythagoraslijn';
  MinMaxPercNegYfToSmallestFactShift(-0.484,1.484,-1.000,0.250,0.05,True);
  p := 12;
  h := Power(2,-(p/2));
  a := h*Cos(p*pi/4);
  b := h*Sin(p*pi/4);
  with pbMain.Canvas do
  begin
    Line(Round(xs+xf*0),Round(ys+yf*0),Round(xs+xf*a),Round(ys+yf*-b));
    LineTo(Round(xs+xf*(a+b)),Round(ys+yf*(a-b)));
  end;
  x := 1;
  y := 1;
  for n := 2 to Round(Power(2,p)-1) do
  begin
    m := n;
    s := 1;
    while m > 1 do
    begin
     if m mod 2 = 1 then s := s + 1;
     m := m div 2;
    end;
    if s mod 4 = 0 then x := x + 1;
    if s mod 4 = 1 then y := y + 1;
    if s mod 4 = 2 then x := x - 1;
    if s mod 4 = 3 then y := y - 1;
    pbMain.Canvas.LineTo(Round(xs+xf*(a*x+b*y)),Round(ys+yf*(a*y-b*x)));
  end;
end;

procedure TfrmMain.miLogspiraClick(Sender: TObject);
var
  a, b, r, t, x, y: Double;
begin
  Clear;
  Label4.Caption := 'Logaritmische spiraal';
  MinMaxPercNegYfToSmallestFactShift(-4,4,-3,3,0.05,True);
  a := 0.1;
  b := 0.1;
  pbMain.Canvas.MoveTo(Round(xs+xf*a),Round(ys));
  t := 0;
  While Round(t)<=35 do
  begin
    r := a*Exp(b*t);
    x := r*Cos(t);
    y := r*Sin(t);
    pbMain.Canvas.LineTo(Round(xs+xf*x),Round(ys+yf*y));
    t := t+0.1;
  end;
end;

procedure TfrmMain.miMandelClick(Sender: TObject);
var
  a,b,p1,p2,p3,p4,x,y,z: Double;
  i,j,k,n1,n2: Integer;
begin
  Clear;
  Label4.Caption := 'Mandelbrot set';
  p1 := -2.5;
  p2 := -1.5;
  p3 := 1.5;
  p4 := 1.5;
  n1 := pbMain.Canvas.Width div 2;
  n2 := Trunc(0.833*n1*(p4-p2)/(p3-p1));
  for i := -n1 to n1 do
  begin
    a := ((n1-i)*p1+(n1+i)*p3)/(2*n1);
    for j := 0 to n2 do
    begin
      b := ((n2-j)*p2+(n2+j)*p4)/(2*n2);
      x := a;
      y := b;
      for k := 1 to 50 do
      begin
        z := x;
        x := x*x-y*y+a;
        y := 2*y*z+b;
        if x*x+y*y > 16 then break;
      end;
      pbMain.Canvas.Pixels[pbMain.Canvas.Width div 2+i,pbMain.Canvas.Height div 2-50-j] := EgaColor[k mod 16];
      pbMain.Canvas.Pixels[pbMain.Canvas.Width div 2+i,pbMain.Canvas.Height div 2-50+j] := EgaColor[k mod 16];
    end;
  end;
end;

procedure TfrmMain.miMinkClick(Sender: TObject);
var
  k, l, m, n, p, s: Integer;
  h, x, y: Double;
  a: Array of Integer;
  t: Array of Integer;
begin
  Clear;
  Label4.Caption := 'De worst van Minkowski';
  MinMaxPercNegYfToSmallestFactShift(0,1,-0.333,0.333,0.05,False);
  SetLength(a,7);
  a[0] := 0;
  a[1] := 1;
  a[2] := 0;
  a[3] := 3;
  a[4] := 3;
  a[5] := 0;
  a[6] := 1;
  a[7] := 0;
  p := 3;
  SetLength(t,p);
  h := Power(4,-p);
  x := 0;
  y := 0;
  pbMain.Canvas.MoveTo(Round(xs),Round(ys));
  for n := 0 to Trunc(Power(8,p)-1) do
  begin
    m := n;
    for l := 0 to p-1 do
    begin
      t[l] := m mod 8;
      m := m div 8;
    end;
    s := 0;
    for k := 0 to p-1 do
    begin
      s := s + a[t[k]];
      s := s mod 4;
    end;
    case s of
      0: x := x+h;
      1: y := y+h;
      2: x := x-h;
      3: y := y-h;
    end;
    pbMain.Canvas.LineTo(Trunc(xs+xf*x),Trunc(ys+yf*y));
  end;
end;

procedure TfrmMain.miMiraClick(Sender: TObject);
var
  n: Integer;
  a,b,c,w,x,y,z: Double;
begin
  Clear;
  Label4.Caption := 'Banen van Myra ''s algemeen systeem';
  MinMaxPercNegYfToSmallestFactShift(-40,40,-30,30,0.05,True);
  a := 0.31;
  b := 1;
  c := 2-2*a;
  x := 12;
  y := 0;
  w := a*x+c*x*x/(1+x*x);
  for n := 0 to 10000 do
  begin
    pbMain.Canvas.Pixels[Trunc(xs+xf*x),Trunc(ys+yf*y)] := clBlack;
    z := x;
    x := b*y+w;
    w := a*x+c*x*x/(1+x*x);
    y := w-z;
  end;
end;

procedure TfrmMain.miMondriaanClick(Sender: TObject);
var
  h,k,x,y: Double;
  n: Integer;
begin
  Clear;
  Label4.Caption := 'Moderne kunst';
  MinMaxPercNegYfToSmallestFactShift(-0.3,1.3,-0.1,1.1,0.05,True);
  Randomize;
  h := 0.1;
  pbMain.Canvas.Line(Trunc(xs),Trunc(ys),Trunc(xs+xf),Trunc(ys));
  pbMain.Canvas.LineTo(Trunc(xs+xf),Trunc(ys+yf));
  pbMain.Canvas.LineTo(Trunc(xs),Trunc(ys+yf));
  pbMain.Canvas.LineTo(Trunc(xs),Trunc(ys));
  for n := 1 to 100 do
  begin
    x := 0.8*Round(Random(100))/100+0.1;
    y := 0.8*Round(Random(100))/100+0.1;
    k := h*(1-Sqrt(Random));
    if Random > 0.5 then
      pbMain.Canvas.Line(Trunc(xs+xf*(x-k)),Trunc(ys+yf*y),Trunc(xs+xf*(x+k)),Trunc(ys+yf*y))
    else
      pbMain.Canvas.Line(Trunc(xs+xf*x),Trunc(ys+yf*(y-k)),Trunc(xs+xf*x),Trunc(ys+yf*(y+k)));
  end;
end;

procedure TfrmMain.miPythb1Click(Sender: TObject);
var
  j, k, l, m, n, p: Integer;
  c, f, h, u, v, x, y: Double;
  a: Array[0..8] of Integer;
begin
  Clear;
  Label4.Caption := 'Boom van Pythagoras';
  MinMaxPercNegYfToSmallestFactShift(-8,8,-4,8,0.05,True);
  p := 8;
  x := 0;
  y := 0;
  u := 1;
  v := 1;
  c := 1/Sqrt(2);
  for m := 0 to p do
  begin
    for n := Round(Power(2,m)) to Round(Power(2,(m+1))-1) do
    begin
      l := n;
      h := 1;
      x := 0;
      y := 0;
      f := 0;
      for k := 0 to m-1 do
      begin
        a[m-k] := l mod 2;
        l := Trunc(l/2)
      end;
      x := 0;
      y := 0;
      for j := 1 to m do
      begin
        if a[j] = 0 then
        begin
          x := x-h*(Cos(f)+2*Sin(f));
          y := y+h*(2*Cos(f)-Sin(f));
          f := f+pi/4;
          h := c*h;
        end
        else
        begin
          x := x+h*(Cos(f)-2*Sin(f));
          y := y+h*(2*Cos(f)+Sin(f));
          f := f-pi/4;
          h := c*h;
        end;
      end;
      u := h*(Cos(f)+Sin(f));
      v := h*(Cos(f)-Sin(f));
      with pbMain.Canvas do
      begin
        Line(Round(xs+xf*(x-v)),Round(ys+yf*(y-u)),Round(xs+xf*(x+u)),Round(ys+yf*(y-v)));
        LineTo(Round(xs+xf*(x+v)),Round(ys+yf*(y+u)));
        LineTo(Round(xs+xf*(x-u)),Round(ys+yf*(y+v)));
        LineTo(Round(xs+xf*(x-v)),Round(ys+yf*(y-u)));
      end;
    end;
  end;
end;

procedure TfrmMain.miPythb2Click(Sender: TObject);
var
  x, y: Array[0..4096] of Double;
  a1, a2, b1, b2, c, c1, c2, d1, d2, f, s, x0, x1, xa, xb, xc, xd, y0, y1, ya, yb, yc, yd: Double;
  j, m: Integer;
begin
  Clear;
  Label4.Caption := 'Scheve Pythagorasboom';
  MinMaxPercNegYfToSmallestFactShift(-1.7,4.8,-1.1,3.33,0.05,True);
  f := pi/3;
  c := Cos(f);
  s := Sin(f);
  a1 := -c*s;
  a2 := Power(c,2);
  b1 := a1+a2;
  b2 := -a1+a2;
  c1 := b2;
  c2 := 1-b1;
  d1 := 1-a1;
  d2 := 1-a2;
  x[2] := 0;
  y[2] := 0;
  x[3] := 1;
  y[3] := 0;
  with pbMain.Canvas do
  begin
    Line(Round(xs),Round(ys),Round(xs+xf*1),Round(ys));
    LineTo(Round(xs+xf*1),Round(ys+yf*-1));
    LineTo(Round(xs),Round(ys+yf*-1));
    LineTo(Round(xs),Round(ys+yf*-1));
    LineTo(Round(xs),Round(ys));
  end;
  for m:=1 to 10 do
  begin
    for j:=0 to Round(Power(2,m-1)-1) do
    begin
      x0 := x[Round(Power(2,m)+2*j)];
      y0 := y[Round(Power(2,m)+2*j)];
      x1 := x[Round(Power(2,m)+2*j+1)];
      y1 := y[Round(Power(2,m)+2*j+1)];
      xa := x0+a1*(x1-x0)-a2*(y1-y0);
      ya := y0+a2*(x1-x0)+a1*(y1-y0);
      xb := x0+b1*(x1-x0)-b2*(y1-y0);
      yb := y0+b2*(x1-x0)+b1*(y1-y0);
      xc := x0+c1*(x1-x0)-c2*(y1-y0);
      yc := y0+c2*(x1-x0)+c1*(y1-y0);
      xd := x0+d1*(x1-x0)-d2*(y1-y0);
      yd := y0+d2*(x1-x0)+d1*(y1-y0);
      x[Round(Power(2,m+1)+4*j)] := xa;
      y[Round(Power(2,m+1)+4*j)] := ya;
      x[Round(Power(2,m+1)+4*j+1)] := xb;
      y[Round(Power(2,m+1)+4*j+1)] := yb;
      x[Round(Power(2,m+1)+4*j+2)] := xc;
      y[Round(Power(2,m+1)+4*j+2)] := yc;
      x[Round(Power(2,m+1)+4*j+3)] := xd;
      y[Round(Power(2,m+1)+4*j+3)] := yd;
      with pbMain.Canvas do
      begin
        Line(Round(xs+xf*x0),Round(ys+yf*y0),Round(xs+xf*xa),Round(ys+yf*ya));
        LineTo(Round(xs+xf*xb),Round(ys+yf*yb));
        LineTo(Round(xs+xf*x1),Round(ys+yf*y1));
        LineTo(Round(xs+xf*xd),Round(ys+yf*yd));
        LineTo(Round(xs+xf*xc),Round(ys+yf*yc));
        LineTo(Round(xs+xf*x0),Round(ys+yf*y0));
      end;
    end;
  end;

end;

procedure TfrmMain.miPythb3Click(Sender: TObject);
var
  x1,y1,x2,y2,u1,v1,u2,v2: Array [0..12] of Double;
  p, j, n ,m,s1: Integer;
  f,c,s,a1,a2,b1,b2,c1,c2,d1,d2,x,y,u,v,x3,y3: Double;
  procedure Tekenboom;
  var
    j: Integer;
  begin
    for j := s1 to p do
    begin
      x := x1[j-1];
      y := y1[j-1];
      u := u1[j-1];
      v := v1[j-1];
      x3 := u-x;
      y3 := v-y;
      x1[j] := x+a1*x3-a2*y3;
      y1[j] := y+a2*x3+a1*y3;
      u1[j] := x+b1*x3-b2*y3;
      v1[j] := y+b2*x3+b1*y3;
      x2[j] := x+c1*x3-c2*y3;
      y2[j] := y+c2*x3+c1*y3;
      u2[j] := x+d1*x3-d2*y3;
      v2[j] := y+d2*x3+d1*y3;
      with pbMain.Canvas do
      begin
        Line(Trunc(xs+xf*x),Trunc(ys+yf*y),Trunc(xs+xf*x1[j]),Trunc(ys+yf*y1[j]));
        LineTo(Trunc(xs+xf*u1[j]),Trunc(ys+yf*v1[j]));
        LineTo(Trunc(xs+xf*u),Trunc(ys+yf*v));
        LineTo(Trunc(xs+xf*x),Trunc(ys+yf*y));
        LineTo(Trunc(xs+xf*x2[j]),Trunc(ys+yf*y2[j]));
        LineTo(Trunc(xs+xf*u2[j]),Trunc(ys+yf*v2[j]));
        LineTo(Trunc(xs+xf*u),Trunc(ys+yf*v));
      end;
    end;
  end;

begin
  Clear;
  Label4.Caption := 'Pythagorasboom, backtrackmethode';
  MinMaxPercNegYfToSmallestFactShift(-3.2,3,-1,3.3,0.05,True);
  //MinMaxPercXPercYToFactShift(-5,5,-3,4.5,0.05,0.05,True);
  p := 12;
  f := pi/5;
  c := Cos(f);
  s := Sin(f);
  a1 := -c*s;
  a2 := Power(c,2);
  b1 := a1+a2;
  b2 := -a1+a2;
  c1 := b2;
  c2 := 1-b1;
  d1 := 1-a1;
  d2 := 1-a2;
  x1[0] := 0;
  y1[0] := 0;
  u1[0] := 1;
  v1[0] := 0;
  with pbMain.Canvas do
  begin
    Line(Trunc(xs),Trunc(ys),Trunc(xs),Trunc(ys+yf*-1));
    LineTo(Trunc(xs+xf*1),Trunc(ys+yf*-1));
    LineTo(Trunc(xs+xf*1),Trunc(ys));
  end;
  s1 := 1;
  Tekenboom;
  for m:=1 to Round(Power(2,p-1)-1) do
  begin
    s := p;
    n := m;
    while n mod 2 = 0 do
    begin
      n := n div 2;
      s := s-1;
    end;
    s1 := Trunc(s);
    x1[s1-1] := x2[s1-1];
    y1[s1-1] := y2[s1-1];
    u1[s1-1] := u2[s1-1];
    v1[s1-1] := v2[s1-1];
    Tekenboom;
  end;
end;

procedure TfrmMain.miPythbsClick(Sender: TObject);
var
  m,n,p,s: Integer;
  a1,a2,b1,b2,c1,c2,w,u,v,x,x3,y,y3: Double;
  x1,y1,x2,y2,u1,v1,u2,v2: Array [0..12] of Double;

  procedure Gosub260;
  begin
    a2 := a2*(1+(Random-0.5)*w);
    c2 := c2*(1+(Random-0.5)*w);
    b2 := (a2+c2)/2+0.5;
  end;

  procedure Gosub160;
  var
    j: Integer;
    //: Double;
  begin
    for j := s to p do
    begin
      x := x1[j-1];
      y := y1[j-1];
      u := u1[j-1];
      v := v1[j-1];
      x3 := u-x;
      y3 := v-y;
      x1[j] := x+a1*x3-a2*y3;
      y1[j] := y+a2*x3+a1*y3;
      u1[j] := x+b1*x3-b2*y3;
      v1[j] := y+b2*x3+b1*y3;
      x2[j] := u1[j];
      y2[j] := v1[j];
      u2[j] := x+c1*x3-c2*y3;
      v2[j] := y+c2*x3+c1*y3;
      with pbMain.Canvas do
      begin
        Line(Trunc(xs+xf*x),Trunc(ys+yf*y),Trunc(xs+xf*x1[j]),Trunc(ys+yf*y1[j]));
        Line(Trunc(xs+xf*u1[j]),Trunc(ys+yf*v1[j]),Trunc(xs+xf*x2[j]),Trunc(ys+yf*y2[j]));
        Line(Trunc(xs+xf*u2[j]),Trunc(ys+yf*v2[j]),Trunc(xs+xf*u),Trunc(ys+yf*v));
      end;
    end;
  end;

  procedure Gosub140;
  begin
    x1[s-1] := x2[s-1];
    y1[s-1] := y2[s-1];
    u1[s-1] := u2[s-1];
    v1[s-1] := v2[s-1];
    Gosub160;
  end;

begin
  Clear;
  Label4.Caption := 'Pythagorasboom, backtrackmethode';
  MinMaxPercNegYfToSmallestFactShift(-9.5,10.5,-3,12,0.05,True);
  w := 0.15;
  p := 11;
  a1 := 0;
  a2 := 3;
  b1 := 0.5;
  b2 := 3.5;
  c1 := 1;
  c2 := 3;
  x1[0] := 0;
  y1[0] := 0;
  u1[0] := 1;
  v1[0] := 0;
  pbMain.Canvas.Line(Trunc(xs),Trunc(ys),Trunc(xs+xf),Trunc(ys));
  s := 1;
  Gosub260;
  Gosub160;
  for m := 1 to Trunc(Power(2,(p-1))-1) do
  begin
    s := p;
    n := m;
    if s<5 then Gosub260;
    while n mod 2 = 0 do
    begin
      n := n div 2;
      s := s-1;
    end;
    Gosub140;
  end;
end;

procedure TfrmMain.miPythtClick(Sender: TObject);
var
  a1, a2, b1, b2, c, c1, c2, d1, d2, f, s, u, v, x0, x1, xa, xb, xc, xd, y0, y1, ya, yb, yc, yd: Double;
  x, y: Array[0..4096] of Double;
  j, m: Integer;
begin
  Clear;
  Label4.Caption := 'Kale pythagorasboom';
  MinMaxPercNegYfToSmallestFactShift(-3.5,4.5,-2,4,0.05,True);
  f := pi/4;
  c := Cos(f);
  s := Sin(f);
  a1 := -c*s;
  a2 := Power(c,2);
  b1 := a1+a2;
  b2 := -a1+a2;
  c1 := b2;
  c2 := 1-b1;
  d1 := 1-a1;
  d2 := 1-a2;
  x[2] := 0;
  y[2] := 0;
  x[3] := 1;
  y[3] := 0;
  pbMain.Canvas.Line(Trunc(xs+xf*0.5),Trunc(ys+yf*-1),Trunc(xs+xf*0.5),Trunc(ys));
  for m := 1 to 10 do
  begin
    for j := 0 to Round(Power(2,m-1)-1) do
    begin
      x0 := x[Round(Power(2,m)+2*j)];
      y0 := y[Round(Power(2,m)+2*j)];
      x1 := x[Round(Power(2,m)+2*j+1)];
      y1 := y[Round(Power(2,m)+2*j+1)];
      u := x1-x0;
      v := y1-y0;
      xa := x0+a1*u-a2*v;
      ya := y0+a2*u+a1*v;
      xb := x0+b1*u-b2*v;
      yb := y0+b2*u+b1*v;
      xc := x0+c1*u-c2*v;
      yc := y0+c2*u+c1*v;
      xd := x0+d1*u-d2*v;
      yd := y0+d2*u+d1*v;
      x[Round(Power(2,m+1)+4*j)] := xa;
      y[Round(Power(2,m+1)+4*j)] := ya;
      x[Round(Power(2,m+1)+4*j+1)] := xb;
      y[Round(Power(2,m+1)+4*j+1)] := yb;
      x[Round(Power(2,m+1)+4*j+2)] := xc;
      y[Round(Power(2,m+1)+4*j+2)] := yc;
      x[Round(Power(2,m+1)+4*j+3)] := xd;
      y[Round(Power(2,m+1)+4*j+3)] := yd;
      with pbMain.Canvas do
      begin
        Line(Trunc(xs+xf*((x0+x1)/2)),Trunc(ys+yf*((y0+y1)/2)),Trunc(xs+xf*((xa+xb)/2)),Trunc(ys+yf*((ya+yb)/2)));
        Line(Trunc(xs+xf*((x0+x1)/2)),Trunc(ys+yf*((y0+y1)/2)),Trunc(xs+xf*((xc+xd)/2)),Trunc(ys+yf*((yc+yd)/2)));
      end;
    end;
  end;
end;

procedure TfrmMain.miSierClick(Sender: TObject);
var
  a, n1, u1, u2, v1, v2, x, y: Double;
  k, l, m, n, p: Integer;
  t: Array[0..5] of Double;
begin
  Clear;
  Label4.Caption := 'Zeef van Sierpinski';
  {xf := pbMain.Width/4;
  yf := pbMain.Height/4;
  xs := pbMain.Width div 2;
  ys := pbMain.Height div 2;}
  MinMaxPercNegYfToSmallestFactShift(-1.723,1.723,-1.973,1.01,0.05,False);
  p := 5;
  a := Sqrt(3);
  for m := 0 to p do
  begin
    for n := 0 to Round(Power(3,m)-1) do
    begin
      n1 := n;
      for l := 0 to m-1 do
      begin
        t[l] := n1 mod 3;
        n1 := Trunc(n1/3);
      end;
      x := 0;
      y := 0;
      for k := 0 to m-1 do
      begin
        x := x+Cos((4*t[k]+1)*pi/6)/Power(2,k);
        y := y+Sin((4*t[k]+1)*pi/6)/Power(2,k);
      end;
      u1 := x+a/Power(2,m+1);
      u2 := x-a/Power(2,m+1);
      v1 := y-1/Power(2,m+1);
      v2 := y+1/Power(2,m);
      with pbMain.Canvas do
      begin
        Line(Trunc(xs+xf*u1),Trunc(ys+yf*v1),Trunc(xs+xf*x),Trunc(ys+yf*v2));
        LineTo(Trunc(xs+xf*u2),Trunc(ys+yf*v1));
        LineTo(Trunc(xs+xf*u1),Trunc(ys+yf*v1));
      end;
    end;
  end;
end;

procedure TfrmMain.miSterfractalClick(Sender: TObject);
var
  f,m,n,p,v: Integer;
  a,b,r,x,y: Double;
begin
  Clear;
  Label4.Caption := 'Sterfractal';
  MinMaxPercNegYfToSmallestFactShift(-0.5,1.5,-0.8,0.7,0.05,True);
  p := 5;
  v := 4;
  a := 144;
  r := 0.35;
  a := a*pi/180;
  pbMain.Canvas.MoveTo(Trunc(xs),Trunc(ys));
  x := 0;
  y := 0;
  for n := 0 to Trunc((v+1)*Power(v,p-1)-1) do
  begin
    m := n;
    b := n*a;
    f := 0;
    while (m mod v=0) and (f<p-1) do
    begin
      f := f+1;
      m := m div v;
    end;
    x := x+Power(r,(p-f-1))*Cos(b);
    y := y+Power(r,(p-f-1))*Sin(b);
    pbMain.Canvas.LineTo(Trunc(xs+xf*x),Trunc(ys+yf*y));
  end;
end;

procedure TfrmMain.miStofaClick(Sender: TObject);
var
  a,b,c,r,x,y,z: Double;
  k: Integer;
begin
  Clear;
  Label4.Caption := 'Stoffractal, varia, Monte Carlo methode';
  MinMaxPercNegYfToSmallestFactShift(-0.5,1,-0.866,0.866,0.05,False);
  Randomize;
  r := 1;
  a := r*Cos(2*pi/3);
  b := r*Sin(2*pi/3);
  c := 2.75;
  x := 1;
  y := 0;
  for k := 1 to 10000 do
  begin
    if Random < 0.5 then
    begin
      z := x;
      x := a*x-b*y;
      y := b*z+a*y;
      {z := x;
      x := x/2+y;
      y := -z-y/2;}
      {z := -x;
      x := y/2;
      y := z+y;}
    end
    else
    begin
      z := x;
      x := (x*x-y*y+c-1)/c;
      y := 2*z*y/c;
    end;
    pbMain.Canvas.Pixels[Trunc(xs+xf*x),Trunc(ys+yf*y)] := clBlue;
  end;
end;

procedure TfrmMain.miStofbClick(Sender: TObject);
var
  a,b,c,d,x,y: Double;
  n,m,p, s: Integer;
  x1,y1,x2,y2: Array [0..12] of Double;

  procedure Gosub150;
  var
    j: Integer;
  begin
    for j := s to p do
    begin
      x := x1[j-1];
      y := y1[j-1];
      x1[j] := a*x-b*y;
      y1[j] := b*x+a*y; // rotatie
      x2[j] := c*x-d*y+1-c;
      y2[j] := d*x+c*y-d; // rotatie
      pbMain.Canvas.Pixels[Trunc(xs+xf*x1[j]),Trunc(ys+yf*y1[j])] := clBlue;
      pbMain.Canvas.Pixels[Trunc(xs+xf*x2[j]),Trunc(ys+yf*y2[j])] := clBlue;
    end;
  end;

begin
  Clear;
  Label4.Caption := 'Stoffractal backtrack methode';
  MinMaxPercNegYfToSmallestFactShift(-1.5,2.5,-1.3,1.7,0.05,True);
  a := 0;
  b := 0.7071;
  c := 0.5;
  d := -0.5;
  p := 12;
  pbMain.Canvas.Pixels[Trunc(xs),Trunc(ys)] := clBlue;
  pbMain.Canvas.Pixels[Trunc(xs+xf*1),Trunc(ys)] := clBlue;
  x1[0] := a;
  y1[0] := b;
  pbMain.Canvas.Pixels[Trunc(xs+xf*x1[0]),Trunc(ys+yf*y1[0])] := clBlue;
  s := 1;
  Gosub150;
  for m := 1 to Trunc(Power(2,p-1)-1) do
  begin
    s := p;
    n := m;
    while n mod 2 = 0 do
    begin
      n := n div 2;
      s := s-1;
      x1[s-1] := x2[s-1];
      y1[s-1] := y2[s-1];
      Gosub150;
    end;
  end;
end;

procedure TfrmMain.miStofbtClick(Sender: TObject);
var
  m,n,p,s: Integer;
  x1,x2,x3,y1,y2,y3: Array [0..7] of Double;
  a,b,c,d,e,f,g,h,t1,t2: Double;

  procedure Goto150;
  var
    j: Integer;
  begin
    for j := s to p do
    begin
      x := x1[j-1];
      y := y1[j-1];
      x1[j] := a*x-b*y;
      y1[j] := b*x+a*y;
      x2[j] := c*x-d*y+1-c;
      y2[j] := d*x+c*y-d;
      x3[j] := e*x-f*y+g;
      y3[j] := f*x+e*y+h;
      pbMain.Canvas.Pixels[Trunc(xs+xf*x1[j]),Trunc(ys+yf*y1[j])] := clBlue;
      pbMain.Canvas.Pixels[Trunc(xs+xf*x2[j]),Trunc(ys+yf*y2[j])] := clBlue;
      pbMain.Canvas.Pixels[Trunc(xs+xf*x3[j]),Trunc(ys+yf*y3[j])] := clBlue;
    end;
  end;

begin
  Clear;
  Label4.Caption := 'Stoffractal backtrack methode, drietallig';
  MinMaxPercNegYfToSmallestFactShift(-0.8,1.6,-0.6,1.2,0.05,True);
  p := 7;
  t1 := 0.5;
  t2 := 0.866; // positie top
  a := 0.43;
  b := 0.3;
  c := a;
  d := b;
  e := a;
  f := b;
  g := t1*(1-e)+t2*f;
  h := t1*f+t2*(1-e);
  x1[0] := 0.5;
  y1[0] := 0.289;
  with pbMain.Canvas do
  begin
    Pixels[Trunc(xs),Trunc(ys)] := clBlue;
    Pixels[Trunc(xs+xf*1),Trunc(ys)] := clBlue;
    Pixels[Trunc(xs+xf*t1),Trunc(ys+yf*t2)] := clBlue;
    Pixels[Trunc(xs+xf*x1[0]),Trunc(ys+yf*y1[0])] := clBlue;
  end;
  for m := 0 to Trunc(Power(3,p-1)-1) do
  begin
    s := p;
    n := m;
    if m = 0 then
    begin
      s := 1;
      Goto150;
      Continue;
    end;
    while n mod 3 = 0 do
    begin
      n := n div 3;
      s := s-1;
    end;
    x1[s-1] := x2[s-1];
    y1[s-1] := y2[s-1];
    x2[s-1] := x3[s-1];
    y2[s-1] := y3[s-1];
    Goto150;
  end;
end;

procedure TfrmMain.miStofClick(Sender: TObject);
var
  a,b,c,d,r1,r2,x,y,z: Double;
  k: Integer;
begin
  Clear;
  Label4.Caption := 'Stoffractal, Monte Carlo methode';
  MinMaxPercNegYfToSmallestFactShift(-1.1,2.1,-1.2,1.2,0.05,True);
  Randomize;
  r1 := 0.6;
  r2 := 0.6;
  a := r1*Cos(2*pi/3);
  b := r1*Sin(2*pi/3);
  c := r2*Cos(2*pi/3);
  d := -r2*Sin(2*pi/3);
  x := a;
  y := b;
  for k := 1 to 10000 do
  begin
    if Random < 0.5 then
    begin
      z := x;
      x := a*x-b*y;
      y := b*z+a*y;
    end
    else
    begin
      z := x;
      x := c*x-d*y+1-c;
      y := d*z+a*y-d;
    end;
    pbMain.Canvas.Pixels[Trunc(xs+xf*x),Trunc(ys+yf*y)] := clBlue;;
  end;
end;

procedure TfrmMain.miWervelClick(Sender: TObject);
var
  b, c, t, z: Double;
  k, l, m, n: Integer;
  x, y: Array[0..4] of Double;
begin
  Clear;
  Label4.Caption := 'Draaiend en krimpend vierkant';
  MinMaxPercNegYfToSmallestFactShift(-4/3,4/3,-1,1,0.05,True);
  b := 0.05;
  c := 1/(Sin(b)+Cos(b));
  for k := 0 To 4 do
  begin
    t := (2*k+1)*pi/4;
    x[k] := Sin(t);
    y[k] := Cos(t);
  end;
  for n := 1 to 64 do
  begin
    pbMain.Canvas.MoveTo(Round(xs+xf*x[0]),Round(ys+yf*y[0]));
    for l := 1 to 4 do
      pbMain.Canvas.LineTo(Round(xs+xf*x[l]),Round(ys+yf*y[l]));
    for m := 0 to 4 do
    begin
      z := x[m];
      x[m] := (x[m]*Cos(b)-y[m]*Sin(b))*c;
      y[m] := (z*Sin(b)+y[m]*Cos(b))*c;
    end;
  end;
end;

procedure TfrmMain.miWikkelClick(Sender: TObject);
var
  a, n: Integer;
  t, x, y: Double;
begin
  Clear;
  Label4.Caption := 'Wikkellijn van cirkel';
  MinMaxPercNegYfToSmallestFactShift(-4.72,7.86,-6.29,3.15,0.05,True);
  a := 1;
  pbMain.Canvas.Brush.Color := clBlue;
  pbMain.Canvas.EllipseC(Trunc(xs),Trunc(ys),Round(xf*a),Round(yf*a));
  pbMain.Canvas.MoveTo(Round(xs+xf*a),Round(ys));
  for n := 0 to 100 do
  begin
    t := 2*pi*n/80;
    x := a*(Cos(t)+t*Sin(t));
    y := a*(Sin(t)-t*Cos(t));
    pbMain.Canvas.LineTo(Round(xs+xf*x),Round(ys+yf*y));
    if n mod 10 = 0 then
      pbMain.Canvas.Line(Round(xs+xf*(Cos(t))),Round(ys+yf*(Sin(t))),Round(xs+xf*x),Round(ys+yf*y));
  end;
end;

procedure TfrmMain.miWolkClick(Sender: TObject);
var
  a,b,w,x,y,z: Double;
  n,p: Integer;

   procedure Gosub170;
  begin
    if x > 1 then
      w := a*x+b*(x-1)
    else
      if x < -1 then
        w := a*x+b*(x+1)
      else
      begin
        case SpinEdit2.Value of
          110: w := a*x+b*Sin(x);
          120: w := a*x+b*Cos(x);
          130: w := a+b*Sin(x);
          140: w := a+b*Cos(x);
          150:
            begin
              if Abs(x) < 1 then
                w := a*x
              else
                w := b*x+(a-b)/x;
            end;
          else
            w := a*x;
        end;   // end case
      end;  //end if
  end;

begin
  Clear;
  Label5.Caption := 'Functies';
  Label5.Visible := True;
  SpinEdit2.MinValue := 110;
  SpinEdit2.MaxValue := 150;
  SpinEdit2.Increment := 10;
  SpinEdit2.Visible := True;
  Label4.Caption := 'Banen van dynamisch systeem';
  MinMaxPercNegYfToSmallestFactShift(-200,200,-150,150,0.05,True);
  a := 3.5;
  b := -3;
  w := 0;
  Screen.Cursor := crHourGlass;
  for p := 1 to 1500 do
  begin
    x := 200*Random-200;
    y := 150*Random-150;
    Gosub170;
    for n := 0 to p do
    begin
      pbMain.Canvas.Pixels[Trunc(xs+xf*x),Trunc(ys+yf*y)] := n*4000*p; //clBlack;
      z := x;
      x := y+w;
      Gosub170;
      y := w-z;
    end;
  end;
  Screen.Cursor := crDefault;
end;

procedure TfrmMain.SpinEdit1Change(Sender: TObject);
var
  i: Integer;
begin
  for i := 0 to miFractals.Count-1 do
  begin
    if miFractals.Items[i].Checked then Break;
  end;
  //frmMain.Caption := 'i=' + i.ToString + ' count=' + miFractals.Count.ToString;
  case i of
    0: miBoomH1Click(Self);
    10: miDraakClick(Self);
  end;
end;

procedure TfrmMain.SpinEdit2Change(Sender: TObject);
begin
  miWolkClick(Self);
end;

procedure TfrmMain.Teken;
var
  j: Integer;
begin
  for j:= s to p do
  begin
    x := x1[j-1];
    y := y1[j-1];
    b := Power(a,j);
    c := a * b * 1.5;
    x1[j] := x + b;
    y1[j] := y + c;
    x2[j] := x + b;
    y2[j] := y - c;
    x3[j] := x - b;
    y3[j] := y + c;
    x4[j] := x - b;
    y4[j] := y - c;
    with pbMain.Canvas do
    begin
      Line(Trunc(xs+xf*(x-b)),Trunc(ys+yf*y),Trunc(xs+xf*(x+b)),Trunc(ys+yf*y));
      Line(Trunc(xs+xf*x1[j]),Trunc(ys+yf*y1[j]),Trunc(xs+xf*x2[j]),Trunc(ys+yf*y2[j]));
      Line(Trunc(xs+xf*x3[j]),Trunc(ys+yf*y3[j]),Trunc(xs+xf*x4[j]),Trunc(ys+yf*y4[j]));
    end;
  end;
end;

procedure TfrmMain.Clear;
begin
  Label1.Visible := False;
  Label2.Visible := False;
  Label3.Visible := False;
  Label5.Visible := False;
  SpinEdit1.Visible := False;
  SpinEdit2.Visible := False;
  cbbBasisLijn.Visible := False;
  cbbModel.Visible := False;
  pbMain.Canvas.Brush.Color := clWhite;
  pbMain.Canvas.FillRect(0,0,pbMain.Width,pbMain.Height);
end;

end.

