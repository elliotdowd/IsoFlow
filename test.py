import wx
from wx.lib.delayedresult import startWorker
from subprocess import Popen, PIPE

class MainWindow(wx.Frame):
    def __init__(self, *args, **kwargs):
        wx.Frame.__init__(self, *args, **kwargs)

        self.panel = wx.Panel(self)
        self.button = wx.Button(self.panel, label="Run!")
        self.command = wx.TextCtrl(self.panel)
        self.result = wx.TextCtrl(self.panel, style=wx.TE_MULTILINE)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.command, 0, wx.EXPAND)
        self.sizer.Add(self.button, 0, wx.EXPAND)
        self.sizer.Add(self.result, 1, wx.EXPAND)

        self.command.SetValue("dir")
        self.button.Bind(wx.EVT_BUTTON, self.CallCommand)

        self.panel.SetSizerAndFit(self.sizer) 
        self.Show()

    def CallCommand(self, e):
        startWorker(self.WorkCommandDone, self.WorkCommand)

    def WorkCommand(self):
        self.button.Disable()
        p = Popen(self.command.GetValue(), shell=True,
                  stdin=PIPE, stdout=PIPE, stderr=PIPE)
        while True:
            line = p.stdout.readline()
            if line != '':
                wx.CallAfter(self.result.AppendText, line)
            else:
                break

    def WorkCommandDone(self, result):
        self.button.Enable()

app = wx.App(False)
win = MainWindow(None)
app.MainLoop()