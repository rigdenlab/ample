'''
Created on 21 Feb 2019

@author: jmht
'''


resolution = 2.4
coiledcoil = True

if resolution <= 1.0:
    argd = {
        'a' : 8.0,
        'e' : 1.0,
        'I' : 200,
        'm' : 200,
        's' : 0.25,
        't' : 10,
        'v' : 0.5
        }
    argd_final = argd.copy()
    argd_final.update({ 'a' : 1,
                       's' : 0.2 })
elif resolution <= 1.3:
    argd = {
        'a' : 8.0,
        'e' : 1.0,
        'I' : 100,
        'm' : 200,
        's' : 0.35,
        't' : 10,
        'v' : 0.25
        }
    argd_final = argd.copy()
    argd_final.update({ 'a' : 1,
                        'm' : 100,
                        's' : 0.3 })
elif resolution <= 1.5:
    argd = {
        'a' : 8.0,
        'e' : 1.0,
        'I' : 150,
        'm' : 50,
        's' : 0.45,
        't' : 10,
        'v' : 0.1
        }
    argd_final = argd.copy()
    argd_final.update({ 'a' : 1,
                       's' : 0.4 })
elif resolution <= 2.0:
    argd = {
        'a' : 8.0,
        'e' : 1.0,
        'I' : 200,
        'm' : 15,
        's' : 0.5,
        't' : 10,
        'v' : 0.0
        }
    argd_final = argd.copy()
    argd_final.update({ 'a' : 1,
                        'e' : resolution - 0.5,
                        's' : 0.45 })
elif resolution <= 2.5 and coiledcoil: 
        argd = {
            'a' : 8.0,
            'e' : resolution - 0.3,
            'I' : 10,
            'm' : 10,
            's' : 0.6,
            't' : 10,
            'v' : 0.0
            }
        argd_final = argd.copy()
        argd_final.update({ 'a' : 1,
                            'e' : resolution - 0.5,
                            's' : 0.55 })
elif  resolution <= 3.0 and coiledcoil:
        argd = {
            'a' : 8.0,
            'e' : resolution - 0.3,
            'I' : 5,
            'm' : 5,
            's' : 0.6,
            't' : 10,
            'v' : 0.0
            }
        argd_final = argd.copy()
        argd_final.update({ 'a' : 1,
                            'e' : resolution - 0.5,
                            's' : 0.55 })
else:  # resolution > 2.0 and not coiled-coil
    argd = {
        'a' : 8.0,
        'e' : 1.0,
        'I' : 15,
        'm' : 10,
        's' : 0.6,
        't' : 10,
        'v' : 0.0
        }
    argd_final = argd.copy()
    argd_final.update({ 'a' : 1,
                        'e' : resolution - 0.5,
                        'm' : 15,
                        's' : 0.55 })


fstr = "-q -m{m} -a{a} -s{s} -v{v} -t{t} -y{resolution} -e{e}"
linesh = fstr.format(resolution=resolution, **argd)
linesh_last = fstr.format(resolution=resolution, **argd_final)
if coiledcoil:
    linesh += " -I{I} -Q".format(**argd)
    linesh_last += " -I{I} -Q".format(**argd_final)


print "GOT  ",linesh
print "GOT2 ",linesh_last