import os

three2one = {
    'ALA' : 'A',    
    'ARG' : 'R',    
    'ASN' : 'N',    
    'ASP' : 'D',    
    'CYS' : 'C',    
    'GLU' : 'E',    
    'GLN' : 'Q',    
    'GLY' : 'G',    
    'HIS' : 'H',    
    'ILE' : 'I',    
    'LEU' : 'L',    
    'LYS' : 'K',    
    'MET' : 'M',    
    'PHE' : 'F',    
    'PRO' : 'P',    
    'SER' : 'S',    
    'THR' : 'T',    
    'TRP' : 'W',    
    'TYR' : 'Y',   
    'VAL' : 'V'
}

one2three =  dict((v, k) for (k, v) in three2one.items())



fn = "0remc30.pdb"
fn = "9remc35.pdb"
fasta = "/media/data/shared/quark_sequences/testset/1MIX_1.fasta"

headers = []
coord_list = []

with open(fn) as f:
    while True:
        # read header
        header = f.readline().strip()
        if not header: break
        headers.append(header)
        fields = header.strip().split()
        nca = int(fields[0])
        ca = []
        for _ in range(nca):
            x, y, z = f.readline().strip().split()
            ca.append((float(x),float(y),float(z)))
        coord_list.append(ca)
        
assert len(headers) == len(coord_list)


with open(fasta) as f:
    f.readline()
    seq=""
    while True:
        s = f.readline().strip()
        if not s: break
        seq += s

seq3 = [ one2three[s] for s in seq ]
assert len(seq3) == len(coord_list[0])

bname = os.path.splitext(os.path.basename(fn))[0]

for i in range(len(headers)):
    fname = "{0}_{1}.pdb".format(bname,i)
    with open(fname,'w') as f:
        f.write("REMARK   3 QUARK MODEL {0} FROM FILE: {1}\n".format(i,fn))
        f.write("REMARK   3 {0}\n".format(headers[i]))
        f.write("MODEL 1\n")
        for j, coord in enumerate(coord_list[i]):
            x,y,z = coord
            f.write("ATOM  {0:5d}  CA  {1} A {2:3d}    {3:8.3F}{4:8.3F}{5:8.3F}  0.50 30.00           C\n".format(int(j),seq3[j],int(j),x,y,z))
        #           "ATOM      1  N   VAL A   2      35.075  18.239 -14.019  0.50 41.75           N"
        f.write("TER\n")
        f.write("ENDMDL\n\n")

        
