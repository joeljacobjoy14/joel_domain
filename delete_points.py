import csv
import string, math
import os
def pltpointlist(ptlist):
    import matplotlib.pyplot as plt
    import numpy as np
    
    pointarr = np.array(ptlist)
    plt.plot(pointarr[:,1], pointarr[:,2], 'o-')
    plt.show()
def genPointList(csvpath):
    
    
    with open(csvpath, 'r') as _mod:
        linelist = _mod.readlines()
        pointslist = [(float(x.strip()), float(y.strip()),float(z.strip()) )\
           for(x,y,z) in [ v.strip().split(",")  for v  in linelist ]]
        
    return pointslist



folderpath = r"C:\Users\Admaren02\Desktop\shipmodel"
for f in os.listdir(folderpath):
    
    
    
    filepath = os.path.join(folderpath, f)
    
    with open(filepath, 'r') as _file:
        
        linelist = _file.readlines()
        pointslist = [(float(x.strip()), float(y.strip()),float(z.strip()) )\
           for(x,y,z) in [ v.strip().split(",")  for v  in linelist ]]
        f_p = pointslist[0][0]
        out_file = open(r"C:\Users\Admaren02\Desktop\shipmodel\modxlocation%0.1f.csv"%f_p, 'w')
#        print("the previous section")
#        print(pointslist)
#        pltpointlist(pointslist)
        pt = pointslist[0]
                                                                                
       # out_file.write("%0.3f, %0.3f , %0.3f\n" %(pt.X, pt.Y, pt.Z))
      #  pltpointlist(pointslist)
        for ptno, t in enumerate(pointslist[1:]):
    
            y1 = pt[1]
            z1 = pt[2]
            y2 = t[1]
            z2 = t[2]
            
            dist = math.sqrt((y1-y2)**2 + (z1-z2)**2)
            if dist  <.0001:
                pass
            else:
                out_file.write("{0:.3f},{1:.3f},{2:.3f}\n".format(pt[0], pt[1], pt[2]))
                
            pt = t
            if ptno == len(pointslist) - 2:
                out_file.write("{0:.3f},{1:.3f},{2:.3f}\n".format(pt[0], pt[1], pt[2]))
                out_file.close()
                
#        print("new section is ")
#        modpoints = genPointList(r"C:\Users\Admaren02\Desktop\ship sections\frwd\modxlocation%0.1f.csv"%f_p)
#        print(modpoints)
#    #    print("the modified section is")
#        pltpointlist(modpoints)

    
_file.close()

        
                

            
            
                
            
#        rs.AddInterpCurve(pointslist, degree=1)
#        for i,pt in enumerate(pointslist):
#            rs.AddTextDot(str(i), pt)

if __name__== "__main__":
    
    
    
    
    filepath =  r"C:\Users\Admaren02\Desktop\shipmodel"
    