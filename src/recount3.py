"""
recount3old.py
"""

BASE_URL = "http://idies.jhu.edu/"
import urllib.request
import os
import urllib.parse
import typing as t

ex_rest_params_1 = {"type": "annotations", "organism": "human", "genomic_unit": "gene", "file_extension": "G026"}
ex_rest_params_2 = {"type": "count_files_gene_or_exon", "organism": "human", "genomic_unit": "gene", "data_source": "recount_project", "project": "SRP107565", "file_extension": "G026"}
ex_rest_params_3 = {"type": "count_files_junctions", "organism": "human", "data_source": "recount_project", "junction_type": "all", "junction_file_extension": "RR"}
ex_rest_params_4 = {"type": "metadata_files", "organism": "human", "data_source": "recount_qc", "project": "SRP096765", "table_name": "recount_pred"}
ex_rest_params_5 = {"type": "bigwig_files", "organism": "human", "data_source": "sra", "project": "SRP009615", "sample": "SRR387777"}

def test_download():
    download(GenericRecount3URL(ex_rest_params_1))
    download(GenericRecount3URL(ex_rest_params_2))
    download(GenericRecount3URL(ex_rest_params_3))
    download(GenericRecount3URL(ex_rest_params_4))
    download(GenericRecount3URL(ex_rest_params_5))

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
        "count_files_junctions": ["organism", "data_source", "junction_type", "junction_file_extension"],
        "metadata_files": ["organism", "data_source", "project", "table_name"],
        "bigwig_files": ["organism", "data_source", "project", "sample"]
    }

    def _validate(self) -> None:
        """Tests whether the correct key-value pairs are present to generate a URL.
        :raises:
            KeyError: if self.rest_params has wrong or missing key-value pairs.
        """
        if self.rest_params["type"] not in expected_params:
            raise KeyError()

        for params_type,params in expected_params.item():
            if self.rest_params["type"] == params_type and any(key not in self.rest_params for key in params):
                raise KeyError()

    def _create_url(self) -> t.Optional[str]:
        """Returns URL based on rest_params["type"].

        :return: URL
        """
        if self.rest_params["type"] == "annotations":
            return self.base_url + "{organism}/annotations/{genomic_unit}_sums/{organism}.{genomic_unit}_sums.{file_extension}.gtf.gz".format(**self.rest_params)
        elif self.rest_params["type"] == "count_files_gene_or_exon":
            return self.base_url + "{organism}/data_sources/{data_source}/{genomic_unit}_sums/{project_2_char}/{project}/{data_source}.gene_sums.{project}.{file_extension}.gz".format(project_2_char=self.rest_params["project"][-2:], **self.rest_params)
        elif self.rest_params["type"] == "count_files_junctions":
            return self.base_url + "{organism}/data_sources/{data_source}/{junction_type}.{junction_file_extension}.gz".format(**self.rest_params)
        elif self.rest_params["type"] == "metadata_files":
            return self.base_url + "{organism}/data_sources/{data_source}/metadata/{project_2_char}/{project}/{data_source}.{table_name}.{project}.MD.gz".format(project_2_char=self.rest_params["project"][-2:], **self.rest_params)
        elif self.rest_params["type"] == "bigwig_files":
            return self.base_url + "{organism}/data_sources/{data_source}/base_sums/{project_2_char}/{project}/{sample_2_char}/{data_source}.base_sums.{project}_{sample}.ALL.bw".format(project_2_char=self.rest_params["project"][-2:], sample_2_char=self.rest_params["sample"][-2], **self.rest_params)
        else:
            return None

    def open(self):
        """

        :return: Opened URL
        """
        return urllib.request.urlopen(self.url, data=None, cafile=None, capath=None, cadefault=False, context=None)

def download(url_object: GenericRecount3URL, path: str = "", mode: str = "w") -> None:
    """Downloads a URL via a URL proxy object and saves it locally to a file.

    :param url_object: Object of GenericRecount3URL.
    :param path: File destination.
    :param mode: Mode of opening the file.
    """
    parsed_url = urllib.parse.urlparse(url_object.url)
    file_path = path + os.path.basename(parsed_url.path)
    with url_object.open() as url_file, open(file_path, mode=mode) as out_file:
        out_file.write(url_file.read())

class GenericFilePath(object):
    """`GenericFilePath` class knows how to open local files or files over URL."""

    def __init__(self, path):
        """Initialize path.
        :param str path: String representing a path to local file(s) or valid URL address of file(s).
        """
        self.path = path

    def open(self):
        """Generator that opens and yields filehandles using appropriate facilities:
        test if path represents a local file or file over URL, if file is compressed
        or not.
        :return: Filehandle to be processed into an instance.
        """
        is_url = self.is_url(self.path)
        compression_type = self.is_compressed(self.path)

        if not compression_type:
            if is_url:
                filehandle = urlopen(self.path)
            else:
                filehandle = open(self.path, "r", encoding="utf-8")
            source = self.path
            yield filehandle, source
            filehandle.close()

        elif compression_type:
            if is_url:
                response = urlopen(self.path)
                path = response.read()
                response.close()
            else:
                path = self.path

            if compression_type == "zip":
                ziparchive = zipfile.ZipFile(io.BytesIO(path), "r") if is_url else zipfile.ZipFile(path)
                for name in ziparchive.infolist():
                    if not name.filename.endswith("/"):
                        filehandle = ziparchive.open(name)
                        source = self.path + "/" + name.filename
                        yield filehandle, source
                        filehandle.close()

            elif compression_type in ("tar", "tar.bz2", "tar.gz"):
                tararchive = tarfile.open(fileobj=io.BytesIO(path)) if is_url else tarfile.open(path)
                for name in tararchive:
                    if name.isfile():
                        filehandle = tararchive.extractfile(name)
                        source = self.path + "/" + name.name
                        yield filehandle, source
                        filehandle.close()

            elif compression_type == "bz2":
                filehandle = bz2.BZ2File(io.BytesIO(path)) if is_url else bz2.BZ2File(path)
                source = self.path
                yield filehandle, source
                filehandle.close()

            elif compression_type == "gz":
                filehandle = gzip.open(io.BytesIO(path)) if is_url else gzip.open(path)
                source = self.path
                yield filehandle, source
                filehandle.close()

    @staticmethod
    def is_compressed(path):
        """Test if path represents compressed file(s).
        :param str path: Path to file(s).
        :return: String specifying compression type if compressed, "" otherwise.
        :rtype: :py:class:`str`
        """
        if path.endswith(".zip"):
            return "zip"
        elif path.endswith(".tar.gz"):
            return "tar.gz"
        elif path.endswith(".tar.bz2"):
            return "tar.bz2"
        elif path.endswith(".gz"):
            return "gz"
        elif path.endswith(".bz2"):
            return "bz2"
        elif path.endswith(".tar"):
            return "tar"
        return ""

    @staticmethod
    def is_url(path):
        """Test if path represents a valid URL.
        :param str path: Path to file.
        :return: True if path is valid url string, False otherwise.
        :rtype: :py:obj:`True` or :py:obj:`False`
        """
        try:
            parse_result = urlparse(path)
            return all((parse_result.scheme, parse_result.netloc, parse_result.path))
        except ValueError:
            return False