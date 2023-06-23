import sys
import traceback
from pathlib import Path
recount_path = (Path.cwd() / '../src/recount').as_posix()
sys.path.insert(1, recount_path)
import recount3 as r3

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
        print(r3.GenericRecount3URL(ex_rest_params_1))
        r3.download(r3.GenericRecount3URL(ex_rest_params_1))
    except Exception:
        print(traceback.format_exc())
    try:
        print(r3.GenericRecount3URL(ex_rest_params_2))
        r3.download(r3.GenericRecount3URL(ex_rest_params_2))
    except Exception:
        print(traceback.format_exc())
    try:
        print(r3.GenericRecount3URL(ex_rest_params_3))
        r3.download(r3.GenericRecount3URL(ex_rest_params_3))
    except Exception:
        print(traceback.format_exc())
    try:
        print(r3.GenericRecount3URL(ex_rest_params_4))
        r3.download(r3.GenericRecount3URL(ex_rest_params_4))
    except Exception:
        print(traceback.format_exc())
    try:
        print(r3.GenericRecount3URL(ex_rest_params_5))
        r3.download(r3.GenericRecount3URL(ex_rest_params_5))
    except Exception:
        print(traceback.format_exc())
    try:
        print(r3.GenericRecount3URL(ex_rest_params_6))
        r3.download(r3.GenericRecount3URL(ex_rest_params_6))
    except Exception:
        print(traceback.format_exc())
    try:
        print(r3.GenericRecount3URL(ex_rest_params_7))
        r3.download(r3.GenericRecount3URL(ex_rest_params_7))
    except Exception:
        print(traceback.format_exc())


if __name__ == "__main__":
    test_download()
