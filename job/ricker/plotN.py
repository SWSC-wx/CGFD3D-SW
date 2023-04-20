import os
import sys
def run(argv):
    n = int(argv[1])
    var = argv[2]
    aix = argv[3]
    print(n)
    isk = 20
    for i in range(20, n, isk):
        if aix == 'x':
            cmd = 'python3 plotYZ.py %d %s' %(i, var)
        if aix == 'y':
            cmd = 'python3 plotXZ.py %d %s' %(i, var)
        if aix == 'z':
            cmd = 'python3 plotXY.py %d %s' %(i, var)
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    sys.exit(run(sys.argv))
