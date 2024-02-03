unit invoer;

{$mode ObjFPC}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, StdCtrls, Buttons;

type

  { TfrmInvoer }

  TfrmInvoer = class(TForm)
    BitBtn1: TBitBtn;
    edX: TEdit;
    edY: TEdit;
    edP: TEdit;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
  private

  public

  end;

var
  frmInvoer: TfrmInvoer;

implementation

{$R *.lfm}

end.

