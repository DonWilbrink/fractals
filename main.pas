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
    MainMenu1: TMainMenu;
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
    procedure FormCreate(Sender: TObject);
    procedure FormResize(Sender: TObject);
    procedure miArchiClick(Sender: TObject);
    procedure miBolspiraClick(Sender: TObject);
    procedure miBoom2Click(Sender: TObject);
    procedure miBoom3Click(Sender: TObject);
    procedure miBoomH1Click(Sender: TObject);
    procedure miBoomH2Click(Sender: TObject);
    procedure miDraak0Click(Sender: TObject);
    procedure miDraak1Click(Sender: TObject);
    procedure miDraakClick(Sender: TObject);
    procedure miKamClick(Sender: TObject);
    procedure miKochClick(Sender: TObject);
    procedure miKronkelClick(Sender: TObject);
    procedure miLevyClick(Sender: TObject);
    procedure miLogspiraClick(Sender: TObject);
    procedure miMinkClick(Sender: TObject);
    procedure miPythb1Click(Sender: TObject);
    procedure miPythb2Click(Sender: TObject);
    procedure miPythb3Click(Sender: TObject);
    procedure miSierClick(Sender: TObject);
    procedure miWervelClick(Sender: TObject);
    procedure miWikkelClick(Sender: TObject);
    procedure SpinEdit1Change(Sender: TObject);
  private
    procedure Teken;
    procedure Clear;
    procedure MinMaxPercNegYfToSmallestFactShift(xlo,xhi,ylo,yhi,perc : double; NegYf : boolean);
    procedure MinMaxPercXPercYToFactShift(xlo, xhi, ylo, yhi, PercX,PercY: double; NegYf: boolean);

  public

  end;

var
  frmMain: TfrmMain;
  a, b, c, n, x, y, xf, yf, xs, ys: Double;
  j, m, p, s: Integer;
  x1, x2 , x3, x4, y1, y2, y3, y4: Array[0..4] of Double;

implementation

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

procedure TfrmMain.MinMaxPercXPercYToFactShift(xlo, xhi, ylo, yhi, PercX,PercY: double; NegYf: Boolean);
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
  end;
end;

procedure TfrmMain.FormCreate(Sender: TObject);
begin
  Height := 800;
  Width := 900;
  Clear;
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
  SpinEdit1.Visible := False;
  cbbBasisLijn.Visible := False;
  cbbModel.Visible := False;
  pbMain.Canvas.Brush.Color := clWhite;
  pbMain.Canvas.FillRect(0,0,pbMain.Width,pbMain.Height);
end;

end.

