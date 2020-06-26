import wx
import wx.grid as grid
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

        # geometry table
        domainGrid = grid.Grid(self)
        domainGrid.CreateGrid(6, 1)

        # set row and column values in domain table
        #domainGrid.SetCellFont(0, 0, wx.Font(12, wx.ROMAN, wx.ITALIC, wx.NORMAL))
        domainGrid.SetColLabelValue(0, "Input")
        domainGrid.SetCellValue(0, 0, "1.5")
        domainGrid.SetCellValue(1, 0, "1.3")
        domainGrid.SetCellValue(2, 0, "0.5")
        domainGrid.SetCellValue(3, 0, "20")
        domainGrid.SetCellValue(4, 0, "30")
        domainGrid.SetCellValue(5, 0, "26")

        # set cell editors for M and N
        #domainGrid.SetCellEditor(4, 0, grid.GridCellNumberEditor(4, 0))
        #domainGrid.SetCellEditor(5, 0, grid.GridCellNumberEditor(5, 0))

        # set rownames
        domainGrid.SetRowLabelValue(0, "Length")
        domainGrid.SetRowLabelValue(1, "Height")
        domainGrid.SetRowLabelValue(2, "Wedge Start")
        domainGrid.SetRowLabelValue(3, "Wedge Angle")
        domainGrid.SetRowLabelValue(4, "M")
        domainGrid.SetRowLabelValue(5, "N")





        domainSizer = wx.BoxSizer(wx.VERTICAL)
        domainSizer.Add(domainGrid, 1, wx.EXPAND)
        self.SetSizer(domainSizer)


    #     # add text
    #     sizer = wx.BoxSizer(wx.HORIZONTAL)
    #     self.label = wx.StaticText(self, label = "Hello")
    #     sizer.Add(self.label, 1, wx.EXPAND)


    #     # add button
    #     self.btn = wx.Button(self, label = "Click Here")
    #     self.btn.Bind(wx.EVT_BUTTON, self.onClick)
    #     sizer.Add(self.btn, 1)

    #     self.SetSizer(sizer)


    #     # add checkbox
    #     self.check1 = wx.CheckBox(self, label = 'Checc')
    #     self.Bind(wx.EVT_CHECKBOX, self.onCheck)
    #     sizer.Add(self.check1)


    # def onCheck(self, event):
    #     cb = event.GetEventObject()
    #     self.label.SetLabelText("Selected" + cb.GetLabel())
        
        

    # def onClick(self, event):
    #     self.label.SetLabelText("Text has been changed")


        # gridSizer = wx.GridSizer(4, 4, 5, 5)

        # for i in range(1, 17):
        #     btn = "Button" + str(i)

        #     gridSizer.Add(wx.Button(self,label=btn), 0, wx.EXPAND)
        #     self.SetSizer(gridSizer)

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


    # def onSubmit(self, event):
    #     # add button action
    #     webbrowser.open('https://google.com')

if __name__== "__main__":
    app = MyApp()
    app.MainLoop()