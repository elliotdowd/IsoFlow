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
        super().__init__(parent=parent, title=title, pos=pos)
        self.OnInit()

    def OnInit(self):
        panel = MyPanel(parent=self)


class MyPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent=parent)

        # add message to panel
        text = wx.StaticText(self, id=wx.ID_ANY, label = 'yeeet', pos=(20,20))
        # ID_ANY means we don't care about the ID

        # add button here
        button = wx.Button(parent=self, label="Clicc", pos = (20, 80))
        button.Bind(event=wx.EVT_BUTTON, handler=self.onSubmit)

    def onSubmit(self, event):
        # add button action
        webbrowser.open('https://google.com')

if __name__== "__main__":
    app = MyApp()
    app.MainLoop()