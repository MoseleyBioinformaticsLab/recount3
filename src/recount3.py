"""
recount3old.py
"""

BASE_URL = "http://duffel.rail.bio/recount3/"  # "http://idies.jhu.edu/"
# http://duffel.rail.bio/recount3/human/homes_index this is base of list of projects
import urllib.request
import ssl
import os
import urllib.parse
import typing as t
import zipfile
import traceback
import gzip
import re
import pandas as pd
import csv
import numpy as np

ex_rest_params_1 = {"type": "annotations", "organism": "human", "genomic_unit": "gene", "file_extension": "G026"}
ex_rest_params_2 = {"type": "count_files_gene_or_exon", "organism": "human", "junction_type": "metadata",
                    "genomic_unit": "gene", "data_source": "sra", "project": "SRP107565", "file_extension": "G026"}
ex_rest_params_3 = {"type": "count_files_junctions", "organism": "human", "data_source": "sra",
                    "junction_type": "metadata", "junction_file_extension": "recount_qc", "project": "SRP107565"}
ex_rest_params_4 = {"type": "metadata_files", "organism": "human", "data_source": "sra", "project": "SRP096765",
                    "table_name": "recount_pred"}
ex_rest_params_5 = {"type": "bigwig_files", "organism": "mouse", "data_source": "sra", "project": "DRP001299",
                    "sample": "DRR014697"}
ex_rest_params_6 = {"type": "data_sources", "organism": "human"}
ex_rest_params_7 = {"type": "data_source_metadata", "organism": "human", "data_source": "sra",
                    "junction_type": "metadata"}


def test_download():
    try:
        print(GenericRecount3URL(ex_rest_params_1))
        download(GenericRecount3URL(ex_rest_params_1))
    except Exception:
        print(traceback.format_exc())
    try:
        print(GenericRecount3URL(ex_rest_params_2))
        download(GenericRecount3URL(ex_rest_params_2))
    except Exception:
        print(traceback.format_exc())
    try:
        print(GenericRecount3URL(ex_rest_params_3))
        download(GenericRecount3URL(ex_rest_params_3))
    except Exception:
        print(traceback.format_exc())
    try:
        print(GenericRecount3URL(ex_rest_params_4))
        download(GenericRecount3URL(ex_rest_params_4))
    except Exception:
        print(traceback.format_exc())
    try:
        print(GenericRecount3URL(ex_rest_params_5))
        download(GenericRecount3URL(ex_rest_params_5))
    except Exception:
        print(traceback.format_exc())
    try:
        print(GenericRecount3URL(ex_rest_params_6))
        download(GenericRecount3URL(ex_rest_params_6))
    except Exception:
        print(traceback.format_exc())
    try:
        print(GenericRecount3URL(ex_rest_params_7))
        download(GenericRecount3URL(ex_rest_params_7))
    except Exception:
        print(traceback.format_exc())


class GenericRecount3URL:

    def __init__(self, rest_params: dict[str, str], base_url: str = BASE_URL):
        """
        :param rest_params: Key-value pairs that describe parts of the URL.
        :param base_url: Recount3 repository URL.
        """
        self.rest_params = rest_params
        self.base_url = base_url
        self._validate()
        self.url = self._create_url()

    def __repr__(self):
        return self.url

    expected_params = {
        "data_sources": ["organism"],
        "data_source_metadata": ["organism", "data_source", "junction_type"],
        "annotations": ["organism", "genomic_unit", "file_extension"],
        "count_files_gene_or_exon": ["organism", "data_source", "genomic_unit", "project", "file_extension"],
        "count_files_junctions": ["organism", "data_source", "junction_type", "junction_file_extension", "project"],
        "metadata_files": ["organism", "data_source", "project", "table_name"],
        "bigwig_files": ["organism", "data_source", "project", "sample"]
    }

    def _validate(self) -> None:
        """Tests whether the correct key-value pairs are present to generate a URL.
        :raises:
            KeyError: if self.rest_params has wrong or missing key-value pairs.
        """
        if self.rest_params["type"] not in self.expected_params:
            raise KeyError()

        for params_type, params in self.expected_params.items():
            if self.rest_params["type"] == params_type and any(key not in self.rest_params for key in params):
                raise KeyError()

    def _create_url(self) -> t.Optional[str]:
        """Returns URL based on rest_params["type"].

        :return: URL
        """
        if self.rest_params["type"] == "data_sources":
            return self.base_url + "{organism}/homes_index".format(**self.rest_params)
        elif self.rest_params["type"] == "data_source_metadata":
            return self.base_url + "{organism}/data_sources/{data_source}/{junction_type}/{data_source}.recount_project.MD.gz".format(
                **self.rest_params)
        elif self.rest_params["type"] == "annotations":
            return self.base_url + "{organism}/annotations/{genomic_unit}_sums/{organism}.{genomic_unit}_sums.{file_extension}.gtf.gz".format(
                **self.rest_params)
        elif self.rest_params["type"] == "count_files_gene_or_exon":
            return self.base_url + "{organism}/data_sources/{data_source}/{junction_type}/{project_2_char}/{project}/{data_source}.recount_project.{project}.MD.gz".format(
                project_2_char=self.rest_params["project"][-2:], **self.rest_params)
        elif self.rest_params["type"] == "count_files_junctions":
            return self.base_url + "{organism}/data_sources/{data_source}/{junction_type}/{project_2_char}/{project}/{data_source}.{junction_file_extension}.{project}.MD.gz".format(
                project_2_char=self.rest_params["project"][-2:], **self.rest_params)
        elif self.rest_params["type"] == "metadata_files":
            return self.base_url + "{organism}/data_sources/{data_source}/metadata/{project_2_char}/{project}/{data_source}.{table_name}.{project}.MD.gz".format(
                project_2_char=self.rest_params["project"][-2:], **self.rest_params)
        elif self.rest_params["type"] == "bigwig_files":
            return self.base_url + "{organism}/data_sources/{data_source}/base_sums/{project_2_char}/{project}/{sample_2_char}/{data_source}.base_sums.{project}_{sample}.ALL.bw".format(
                project_2_char=self.rest_params["project"][-2:], sample_2_char=self.rest_params["sample"][-2:],
                **self.rest_params)
        else:
            return None

    @staticmethod
    def open_url(url: str) -> t.BinaryIO:
        """ Open and return URL file object.

        :param url:
        :return:
        """
        ssl_context = ssl.create_default_context();
        ssl_context.check_hostname = False
        ssl_context.verify_mode = ssl.CERT_NONE
        return urllib.request.urlopen(url, data=None, cafile=None, capath=None, cadefault=False, context=ssl_context)

    def open(self) -> t.BinaryIO:
        return self.open_url(self.url)


def download(url: t.Union[str, GenericRecount3URL], path: str = "", mode: str = "bw") -> None:
    """Downloads a URL via a URL proxy object and saves it locally to a file.

    :param url:
    :param path: File destination.
    :param mode: Mode of opening the file.
    """
    if type(url) != str:
        url = str(url)
    parsed_url = urllib.parse.urlparse(url)
    filename = os.path.basename(parsed_url.path)
    if path.endswith(".zip"):
        with GenericRecount3URL.open_url(url) as url_file, zipfile.ZipFile(path, mode=mode) as zip_file:
            with zip_file.open(filename, mode=mode) as out_file:
                out_file.write(url_file.read())
    else:
        file_path = os.path.join(path, filename)
        with GenericRecount3URL.open_url(url) as url_file, open(file_path, mode=mode) as out_file:
            out_file.write(url_file.read())


def create_sample_project_lists(organism: str = "") -> t.Tuple[list, list]:
    organisms = [organism] if organism else ["mouse", "human"]
    data_sources_rest_params = []
    regex = re.compile(r'/|\n')
    for organism in organisms:
        url = str(GenericRecount3URL({"type": "data_sources", "organism": organism}))
        with GenericRecount3URL.open_url(url) as url_file:
            decoded = url_file.read().decode()
            location = regex.split(decoded)
            locations = list(filter(lambda val: val != 'data_sources' and val != '', location))
            for repo in locations:
                data_sources_rest_params.append(
                    {"type": "data_source_metadata", "organism": organism, "junction_type": "metadata",
                     "data_source": repo})
    projects = set()
    samples = set()

    for rest_params in data_sources_rest_params:
        url = str(GenericRecount3URL(rest_params))
        with GenericRecount3URL.open_url(url) as url_file:
            decoded = gzip.decompress(url_file.read()).decode()
            split = decoded.split("\n")
            for line in split:
                columns = line.split("\t")
                if len(columns) > 3 and columns[2] == columns[3] and columns[2][:3] in ["SRP", "ERP", "DRP"]:
                    projects.add(columns[2])
                    samples.add(columns[1])

    return list(samples), list(projects)


def load_gene_sums_table(filepath: str = "") -> pd.DataFrame:

    """Converts a compressed gene_sums file to pandas dataframe representation

    :param filepath: location of gzip compressed gene_sums file
    :return: pandas dataframe representation of gene counts
    """

    if type(filepath) != str:
        filepath = str(filepath)
    if filepath.endswith(".gz"):
        with gzip.open(filename=filepath, mode='rb') as decompressedGeneSumsFile:
            gene_sums_table = pd.read_csv(decompressedGeneSumsFile, sep="\t", skiprows=2, index_col=0)
    else:
        with open(file=filepath, mode='rb') as GeneSumsFile:
            gene_sums_table = pd.read_csv(GeneSumsFile, sep="\t", skiprows=2, index_col=0)
    return gene_sums_table


def load_gene_sums_matrix(filepath: str = "") -> np.array:

    """Converts a compressed gene_sums file to numpy matrix representation

    :param filepath:
    :return: numpy 2d numpy array containing gene counts
    """
    gene_sums_list = []
    if type(filepath) != str:
        filepath = str(filepath)
    if filepath.endswith(".gz"):
        with gzip.open(filename=filepath, mode='rt', encoding="UTF-8") as decompressedGeneSumsFile:
            file_reader = csv.reader(decompressedGeneSumsFile, delimiter="\t")
            for row in file_reader:
                gene_sums_list.append([float(expression) for expression in row])
    gene_sums_matrix = np.array(gene_sums_list)

    return gene_sums_matrix

def load_annotations(filepath: str = "", annotation_columns: tuple = ("seqname", "source", "feature", "start", "end",
                                                                      "score", "strand", "frame",
                                                                      "attribute")) -> pd.DataFrame:
    """

    :param filepath: location of annotation file
    :param annotation_columns: names of all the columns in the annotation file
    :return: a pandas DataFrame containing all annotations
    """
    if type(filepath) != str:
        filepath = str(filepath)
    if filepath.endswith(".gz"):
        with gzip.open(filename=filepath, mode='rb') as decompressedGeneSumsFile:
            gene_sums_table = pd.read_csv(decompressedGeneSumsFile, sep="\t", columns=annotation_columns, skiprows=2,
                                          index_col=0)
    else:
        with open(file=filepath, mode='rb') as geneSumsFile:
            gene_sums_table = pd.read_csv(geneSumsFile, sep="\t", columns=annotation_columns, skiprows=2, index_col=0)
    return gene_sums_table


if __name__ == "__main__":
    samplist, projlist = create_sample_project_lists("human")
    with open("samplist.txt", "w") as f:
        f.write('\n'.join(samplist))
    with open("projlist.txt", "w") as f:
        f.write('\n'.join(projlist))
    # test_download()
