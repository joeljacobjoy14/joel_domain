import rhinoscriptsyntax as rs
import os 


def plotcsv(filepath, xloc = 0):
    with open(filepath, 'r') as _file:
        linelist = _file.readlines()
        pointslist = [(float(x.strip()), float(y.strip()),float(z.strip()) )\
            for(x,y,z) in [ v.strip().split(",")  for v  in linelist ]]
                
            
        rs.AddInterpCurve(pointslist, degree=3)
        #for i,pt in enumerate(pointslist):
           # rs.AddTextDot(str(i), pt)

if __name__== "__main__":
    folderpath = r"C:\Users\Admaren02\Desktop\duppoints"
    for f in os.listdir(folderpath):
        fullpath = os.path.join(folderpath, f)
        plotcsv(fullpath)
    