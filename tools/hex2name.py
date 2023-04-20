import sys
import os
import subprocess

timing = ''
bld = ''

def runcmd(cmd):
    subp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, encoding="utf-8")
    try:
        outs, errs = subp.communicate(timeout=15)
        return(outs, errs)
    except Exception:
        subp.kill()
        outs, errs = subp.communicate()
        print("RUN COMMAND ERROR! %s"%errs)

def getName(pc):
    global bld
    cmd = "swaddr2line -Cfe %s %s" % (bld, pc)
    (outs, err) = runcmd(cmd)
    str = outs.split("\n")[0]
    if '(' in str:
        name = str.split('(')[0]
        return(name)
    else:
        return(str.strip('\n'))

def hex2name(timing):
    with open(timing, 'r') as fr:
        with open("r.%s" % (timing), 'w') as fw:
            lines = fr.readlines()
            for line in lines:
                strs = line.split(' ')
                outline = ''
                for str in strs:
                    if '[' in str and ']' in str and '0x4ff' in str:
                        pc = str.strip('\n')[1:-1]
                        name = getName(pc)
                        str = name + ' '*(len(pc) - len(name) + 2)
                    outline += ' ' + str
                fw.write(outline)

def main(argv):
    global timing
    global bld
    timing = argv[1]
    bld = argv[2]
    hex2name(timing)

if __name__ == "__main__":
        sys.exit(main(sys.argv))
