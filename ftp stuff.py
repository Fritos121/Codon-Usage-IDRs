import ftplib
import os


def retrieve_ftp_file(local_target_dir, file_name, ftp_obj):
    file = open(os.path.join(local_target_dir, file_name), "wb")
    ftp_obj.retrbinary("RETR " + file_name, file.write)
    file.close()


def unzip_gz(file_name):
    import gzip
    import shutil
    print("Unzipping \"" + file_name + "\"")
    with gzip.open(file_name, 'rb') as f_in:
        with open(file_name[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


# go to ncbi ftp site
ftp_site = "ftp.ncbi.nlm.nih.gov"
ftp = ftplib.FTP(ftp_site)
ftp.login()
fnbase = "_cds_from_genomic.fna.gz"

# where orgs are being saved
base_dir = os.path.abspath("D:/Orthologs/Codon_Dist")

# go to current assembly refseq location
assembly_links = [('/genomes/all/GCF/000/005/575/GCF_000005575.2_AgamP3', '7165'),
                  ('/genomes/all/GCF/000/209/225/GCF_000209225.1_ASM20922v1', '45351')]

for assembly, tax_id in assembly_links:
    # create dir for specific org using tax_id
    dir_name = os.path.join(base_dir, tax_id)
    try:
        os.mkdir(dir_name)
    except OSError:
        pass

    # download archived fna file for all CDS in org
    # make sure file doesn't already exist
    file_name = assembly.split("/")[-1] + fnbase
    if file_name not in os.listdir(dir_name):
        ftp.cwd(assembly)  # make it current working directory
        # dirs = ftp.nlst()
        # print(dirs)
        print("Downloading \"" + file_name + "\" from FTP server...")
        retrieve_ftp_file(dir_name, file_name, ftp)

        file_name = os.path.join(dir_name, file_name)
        unzip_gz(file_name)  # unzip CDS fna file

# go through each dir and create codon dist dictionaries as txt files, delete uncompressed
