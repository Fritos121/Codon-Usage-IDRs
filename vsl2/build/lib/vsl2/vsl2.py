# Copyright Christopher Oldfield

import tempfile
import subprocess
import os
import re

check_aa_re = re.compile(r'[^ACDEFGHIKLMNPQRSTVWY]')

try:
    import pkg_resources
    vsl2_path = pkg_resources.resource_filename(
        "vsl2.data", 'VSL2.jar'
    )
except:
    vsl2_path = os.path.join(os.path.dirname(__file__), 'data', 'VSL2.jar')


def index_by_prefix(lst, prefix):
    for i, line in enumerate(lst):
        if line.startswith(prefix):
            return i

def vsl2b(seq):
    non_standard_aa = check_aa_re.findall(seq)
    if non_standard_aa:
        raise ValueError('Non-standard AAs in sequence: "{}"'.format("".join(non_standard_aa)))
    #create temporary file
    tmpFH = tempfile.NamedTemporaryFile(mode='w', delete=False)
    tmp_name = tmpFH.name
    tmpFH.write(seq)
    tmpFH.close()

    #run
    p = subprocess.Popen(
        ['java', '-jar', vsl2_path, '-s:{}'.format(tmp_name)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    stdout, stderr = p.communicate()
    if p.returncode:
        os.remove(tmp_name)
        raise RuntimeError("error running disopred '{0}'".format(stderr))
    lines = stdout.decode('utf-8').splitlines()

    #trim leading and lagging
    lines = lines[index_by_prefix(lines, "----------")+1:]
    lines = lines[:index_by_prefix(lines, "==========")]

    #parse
    num, res, pred, clss = zip(*[line.strip().split() for line in lines])

    #clean up
    os.remove(tmp_name)

    return num, res, list(map(float,pred)), clss

