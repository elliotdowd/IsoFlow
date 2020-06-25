import wx
import webbrowser

class MyApp(wx.App):
    def __init__(self):
        super().__init__(clearSigInt=True)

        # init frame
        self.InitFrame()

    def InitFrame(self):
        frame = MyFrame(parent=None, title="My Frame", pos=(100, 100))
        frame.Show()


class MyFrame(wx.Frame):
    def __init__(self, parent, title, pos):
        super(MyFrame, self).__init__(parent=parent, title=title, pos=pos, size=(300,200))
        self.OnInit()

    def OnInit(self):
        panel = MyPanel(parent=self)


class MyPanel(wx.Panel):
    def __init__(self, parent):
        super(MyPanel, self).__init__(parent)

        gridSizer = wx.GridSizer(4, 4, 5, 5)

        for i in range(1, 17):
            btn = "Button" + str(i)

            gridSizer.Add(wx.Button(self,label=btn), 0, wx.EXPAND)
            self.SetSizer(gridSizer)

        # # box sizer
        # vbox = wx.BoxSizer(wx.VERTICAL)
        # hbox = wx.BoxSizer(wx.HORIZONTAL)

        # # add message to panel
        # text = wx.StaticText(self, id=wx.ID_ANY, label = 'yeeet', style = wx.ALIGN_CENTER_HORIZONTAL)
        # # ID_ANY means we don't care about the ID

        # vbox.Add(text, 0, wx.EXPAND)
        # self.SetSizer(vbox)

        # text2 = wx.StaticText(self, id=wx.ID_ANY, label = 'yeeet2', style = wx.ALIGN_CENTER_HORIZONTAL)
        # hbox.Add(text2, 0, wx.EXPAND)
        # vbox.Add(hbox)
        # self.SetSizer(vbox)

        # # add button here
        # button = wx.Button(parent=self, label="Clicc", pos = (20, 80))
        # button.Bind(event=wx.EVT_BUTTON, handler=self.onSubmit)


    def onSubmit(self, event):
        # add button action
        webbrowser.open('https://google.com')

if __name__== "__main__":
    app = MyApp()
    app.MainLoop()