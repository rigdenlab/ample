'''
Various miscellaneous functions.
Might end up somewhere else at somepoint.
'''

# get a program test for existsnce
def which(program):
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
#########
def check_for_exe(exename, varname):
    exepath = ''
    print 'looking for', exename
    if not varname:
        print 'no ' + exename + ' given on the command line, looking in the PATH'
        print which(exename)
        if not which(exename):
            print 'You need to give the path for ' + exename
            sys.exit()
        else:

            exepath = which(exename)

    else:
        exepath = varname
        print 'using here', exepath




    if not os.path.exists(exepath):
        print 'You need to give the path for ' + exename + ', executable in the PATH dosnt exist'
        sys.exit()
    else:
        return exepath
########
def get_psipred_prediction(psipred):
    string = ''
    for line in open(psipred):
        get_stat1 = re.compile('^Pred\:\s*(\w*)')
        result_stat1 = get_stat1.match(line)
        if result_stat1:
            stat1_get = re.split(get_stat1, line)
            # print stat1_get[1]
            string = string + stat1_get[1]


    C = 0
    H = 0
    E = 0
    length = len(string)

    for c in string:
        if c == 'C':
            C = C + 1
        if c == 'H':
            H = H + 1
        if c == 'E':
            E = E + 1


    H_percent = float(H) / length * 100
    E_percent = float(E) / length * 100

    if H > 0 and E > 0:
        print  'Your protein is predicted to be mixed alpha beta, your chances of success are intermediate'
    if H == 0 and E > 0:
        print  'Your protein is predicted to be all beta, your chances of success are low'
    if H > 0 and E == 0:
        print  'Your protein is predicted to be all alpha, your chances of success are high'
    if  H == 0 and E == 0:
        print  'Your protein is has no predicted secondary structure, your chances of success are low'
