import wx 
  
#import the newly created GUI file 
import import_formBuilder
class CalcFrame(import_formBuilder.MainFrame): 
   def __init__(self,parent): 
      import_formBuilder.MainFrame.__init__(self,parent)  
		
        
app = wx.App(False) 
frame = CalcFrame(None) 
frame.Show(True) 
#start the applications 
app.MainLoop() 