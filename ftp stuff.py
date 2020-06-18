import ftplib
import os


def retrieve_ftp_file(file_name):
    file = open(file_name, "wb")
    ftp.retrbinary("RETR " + file_name, file.write)
    file.close()


def unzip_gz(file_name):
    import gzip
    import shutil
    with gzip.open(file_name, 'rb') as f_in:
        with open(file_name[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


# go to ncbi ftp site
ftp_site = "ftp.ncbi.nlm.nih.gov"
ftp = ftplib.FTP(ftp_site)
ftp.login()
fnbase = "_cds_from_genomic.fna.gz"

# go to current assembly refseq location (change ftp_pattern in try.py)
ftp.cwd("genomes/all/GCF/000/005/575/GCF_000005575.2_AgamP3")
# dirs = ftp.nlst()
#print(dirs)

# where orgs are being saved
base_dir = "D:/Orthologs/Codon_Dist/"
# create dir for specific org using tax_id and make it current working directory
try:
    os.chdir(base_dir + "7165")
except OSError:
    os.mkdir(base_dir + "7165")
    os.chdir(base_dir + "7165")

# download archived fna file for all CDS in org
# make sure file doesn't already exist
file_name = 'GCF_000005575.2_AgamP3' + fnbase
if file_name not in os.listdir(os.getcwd()):
    retrieve_ftp_file(file_name)
    unzip_gz(file_name)  # unzip CDS fna file

# go through each dir and create codon dist dictionaries as txt files, delete uncompressed
