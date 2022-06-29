"""
recount3old.py
"""

BASE_URL = "http://duffel.rail.bio/recount3/"  #"http://idies.jhu.edu/"
import urllib.request
import ssl
import os
import urllib.parse
import typing as t
import zipfile
import traceback

ex_rest_params_1 = {"type": "annotations", "organism": "human", "genomic_unit": "gene", "file_extension": "G026"}
ex_rest_params_2 = {"type": "count_files_gene_or_exon", "organism": "human", "junction_type": "metadata", "genomic_unit": "gene", "data_source": "sra", "project": "SRP107565", "file_extension": "G026"}
ex_rest_params_3 = {"type": "count_files_junctions", "organism": "human", "data_source": "sra", "junction_type": "metadata", "junction_file_extension": "recount_qc", "project":"SRP107565"}
ex_rest_params_4 = {"type": "metadata_files", "organism": "human", "data_source": "sra", "project": "SRP096765", "table_name": "recount_pred"}
ex_rest_params_5 = {"type": "bigwig_files", "organism": "mouse", "data_source": "sra", "project": "DRP001299", "sample": "DRR014697"}

def test_download():
    try:
        download(GenericRecount3URL(ex_rest_params_1))
    except Exception as error:
        print(traceback.format_exc())
    try:
        download(GenericRecount3URL(ex_rest_params_2))
    except Exception as error:
        print(traceback.format_exc())
    try:
        download(GenericRecount3URL(ex_rest_params_3))
    except Exception as error:
        print(traceback.format_exc())
    try:
        download(GenericRecount3URL(ex_rest_params_4))
    except Exception as error:
        print(traceback.format_exc())
    try:
        download(GenericRecount3URL(ex_rest_params_5))
    except Exception as error:
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

    expected_params = {
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

        for params_type,params in self.expected_params.items():
            if self.rest_params["type"] == params_type and any(key not in self.rest_params for key in params):
                raise KeyError()

    def _create_url(self) -> t.Optional[str]:
        """Returns URL based on rest_params["type"].

        :return: URL
        """
        if self.rest_params["type"] == "annotations":
            return self.base_url + "{organism}/annotations/{genomic_unit}_sums/{organism}.{genomic_unit}_sums.{file_extension}.gtf.gz".format(**self.rest_params)
        elif self.rest_params["type"] == "count_files_gene_or_exon":
            return self.base_url + "{organism}/data_sources/{data_source}/{junction_type}/{project_2_char}/{project}/{data_source}.recount_project.{project}.MD.gz".format(project_2_char=self.rest_params["project"][-2:], **self.rest_params)
        elif self.rest_params["type"] == "count_files_junctions":
            return self.base_url + "{organism}/data_sources/{data_source}/{junction_type}/{project_2_char}/{project}/{data_source}.{junction_file_extension}.{project}.MD.gz".format(project_2_char=self.rest_params["project"][-2:], **self.rest_params)
        elif self.rest_params["type"] == "metadata_files":
            return self.base_url + "{organism}/data_sources/{data_source}/metadata/{project_2_char}/{project}/{data_source}.{table_name}.{project}.MD.gz".format(project_2_char=self.rest_params["project"][-2:], **self.rest_params)
        elif self.rest_params["type"] == "bigwig_files":
            return self.base_url + "{organism}/data_sources/{data_source}/base_sums/{project_2_char}/{project}/{sample_2_char}/{data_source}.base_sums.{project}_{sample}.ALL.bw".format(project_2_char=self.rest_params["project"][-2:], sample_2_char=self.rest_params["sample"][-2:], **self.rest_params)
        else:
            return None

    def open(self):
        """

        :return: Opened URL
        """
        ssl_context = ssl.create_default_context();
        ssl_context.check_hostname = False
        ssl_context.verify_mode = ssl.CERT_NONE
        return urllib.request.urlopen(self.url, data=None, cafile=None, capath=None, cadefault=False, context=ssl_context)

def download(url_object: GenericRecount3URL, path: str = "", mode: str = "bw") -> None:
    """Downloads a URL via a URL proxy object and saves it locally to a file.

    :param url_object: Object of GenericRecount3URL.
    :param path: File destination.
    :param mode: Mode of opening the file.
    """
    parsed_url = urllib.parse.urlparse(url_object.url)
    filename = os.path.basename(parsed_url.path)
    if path.endswith(".zip"):
        with url_object.open() as url_file, zipfile.ZipFile(path, mode=mode) as zip_file:
            with zip_file.open(filename, mode=mode) as out_file:
                out_file.write(url_file.read())
    else:
        file_path = os.path.join(path, filename)
        with url_object.open() as url_file, open(file_path, mode=mode) as out_file:
            out_file.write(url_file.read())

if __name__ == "__main__":
    test_download()