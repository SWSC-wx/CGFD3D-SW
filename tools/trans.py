import sys
import os
import subprocess

def runcmd(cmd):
    subp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, encoding="utf-8")
    try:
        outs, errs = subp.communicate(timeout=15)
        return(outs, errs)
    except subprocess.TimeoutExpired:
        subp.kill()
        outs, errs = subp.communicate()
        print("RUN COMMAND ERROR! %s"%errs)


def main(argv):
    cmd = "ls *.c"
    (out, err) = runcmd(cmd)
    for i in out.split("\n"):
        #cmd = "mv %s %s"%(i, i.split(".")[0]+".cpp")
        cmd = "rm %s"%(i.split(".")[0]+".cpp")
        os.system(cmd)

if __name__ == "__main__":
        sys.exit(main(sys.argv))
