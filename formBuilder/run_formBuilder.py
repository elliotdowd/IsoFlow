import wx 
  
#import the newly created GUI file 
import gui
class CalcFrame(gui.MainFrame): 
   def __init__(self,parent): 
      gui.MainFrame.__init__(self,parent)  
		
        
app = wx.App(False) 
frame = CalcFrame(None) 
frame.Show(True) 
#start the applications 
app.MainLoop() 