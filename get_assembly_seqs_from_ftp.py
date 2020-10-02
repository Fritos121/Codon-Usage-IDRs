import ftplib
import os
import sys
import glob


def retrieve_ftp_file(local_target_dir, ftp_file_name, ftp_obj):
    """
    Download file into directory using ftplib object
    :param local_target_dir: string; where the file is to be stored locally
    :param ftp_file_name: string; name of file to be downloaded from FTP (local filename will be the same as this)
    :param ftp_obj: ftp object from ftplib
    :return: N/A
    """
    print("Downloading \"" + ftp_file_name + "\" from FTP server...")
    file = open(os.path.join(local_target_dir, ftp_file_name), "wb")
    ftp_obj.retrbinary("RETR " + ftp_file_name, file.write)
    file.close()


def unzip_gz(file_name):
    """
    Unzip GZ archive file
    :param file_name: string; name of file to unzip
    :return: N/A
    """
    import gzip
    import shutil
    print("Unzipping \"" + file_name + "\"")
    with gzip.open(file_name, 'rb') as f_in:
        with open(file_name[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


if __name__ == "__main__":

    # ensure proper command line arguments are passed.
    if len(sys.argv) != 3:
        exit("Required positional arguments: {} <infile> <base_directory>".format(sys.argv[0]))

    infile = sys.argv[1]
    in_fh = open(infile)

    # where orgs are being saved
    base_dir = sys.argv[2]

    # expects csv file
    assembly_links = [tuple(line.strip().split(", ")) for line in in_fh.readlines()]

    # go to ncbi ftp site
    ftp_site = "ftp.ncbi.nlm.nih.gov"
    ftp = ftplib.FTP(ftp_site)
    ftp.login()
    fnbase = "_cds_from_genomic.fna.gz"


    # base_dir = r"D:/Orthologs/Source_Org_Codon_Dist"

    # assembly_links = [('/genomes/all/GCF/000/005/575/GCF_000005575.2_AgamP3', '7165'),
                      # ('/genomes/all/GCF/000/209/225/GCF_000209225.1_ASM20922v1', '45351')]

    for assembly, tax_id in assembly_links:
        # create dir for specific org using tax_id; pass if dir already exists
        dir_name = os.path.join(base_dir, tax_id)
        try:
            os.mkdir(dir_name)
        except OSError:
            pass

        # assembly name always at end of ftp link
        file_name = assembly.split("/")[-1] + fnbase

        # download archived fna file for all CDS in org
        # same org can have multiple assemblies... if assembly has already been downloaded, don't download new assembly
        fna_list = glob.glob(os.path.join(dir_name, '') + '*.fna')
        gz_list = glob.glob(os.path.join(dir_name, '') + '*.fna.gz')
        if len(fna_list) > 0:
            continue

        # unzip assembly file already present
        elif len(fna_list) == 0 and len(gz_list) > 0:
            for file in gz_list:
                unzip_gz(file)
            continue

        # make sure file doesn't already exist
        if file_name not in os.listdir(dir_name):
            ftp.cwd(assembly)  # make it current working directory
            retrieve_ftp_file(dir_name, file_name, ftp)

            file_name = os.path.join(dir_name, file_name)
            unzip_gz(file_name)  # unzip CDS fna file

