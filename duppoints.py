import csv
import string



def duppoints(filepath, xloc = 0):
    with open(filepath, 'r') as _file:
        
        x = filepath
        linelist = _file.readlines()
        pointslist = [(float(x.strip()), float(y.strip()),float(z.strip()) )\
           for(x,y,z) in [ v.strip().split(",")  for v  in linelist ]]
        out_file = open('%smodified'%x, 'w')
        pt = pointslist[0]
        out_file.write("%0.3f, %0.3f , %0.3f\n" %(pt.X, pt.Y, pt.Z))
        
        for t in pointslist:
            for i in range(pointslist):
                
                y1 = temp.y
                z1 = temp.z
                y2 = t[i+1].y
                z2 = t[i+1].z
                if((math.sqrt((y1-y2)**2 + (z1-z2)**2) <= 10**-3)):
                    print("duplicate point")
                else:               
                    out_file.write(t[i+1])
        out_file.close()
    _file.close()
    
            
                
            
            
            
            
                
            
#        rs.AddInterpCurve(pointslist, degree=1)
#        for i,pt in enumerate(pointslist):
#            rs.AddTextDot(str(i), pt)

if __name__== "__main__":
    
    duppoints(r"C:\Users\Admaren02\Desktop\ship sections\aft\xlocation900.00.csv")