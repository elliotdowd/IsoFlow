import wx
  
# import the newly created GUI file 
import gui1
class MainFrame(gui1.MainFrame): 
   def __init__(self,parent): 
      gui1.MainFrame.__init__(self,parent)
		
        
app = wx.App(False)
frame = MainFrame(None) 
frame.Show(True) 
#start the applications 
app.MainLoop() 