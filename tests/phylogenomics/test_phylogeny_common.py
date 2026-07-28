from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parents[2] / "scripts" / "06_phylogenomics"
sys.path.insert(0, str(SCRIPT_DIR))

from phylogeny_common import accession_from_text, normalize, parse_attributes, read_fasta, read_tsv


def test_accession_normalization():
    assert accession_from_text("bacteraemia__GCA_050469755.1") == "GCA_050469755.1"
    assert accession_from_text("no_accession") is None


def test_gff_attributes_decode_and_ignore_order():
    attrs = parse_attributes("product=Chaperonin%20GroEL;ID=abc;locus_tag=EF_1")
    assert attrs == {"product": "Chaperonin GroEL", "ID": "abc", "locus_tag": "EF_1"}


def test_normalize_is_exact_match_friendly():
    assert normalize("Phosphate-binding protein PstS 1") == "phosphate binding protein psts 1"


def test_read_fasta_rejects_duplicate_identifiers(tmp_path):
    path = tmp_path / "duplicate.fasta"
    path.write_text(">A\nACG\n>A\nTTT\n")
    try:
        read_fasta(path)
    except ValueError as error:
        assert "Duplicate FASTA identifier" in str(error)
    else:
        raise AssertionError("duplicate identifiers were accepted")


def test_read_tsv_normalizes_embedded_windows_carriage_return(tmp_path):
    path = tmp_path / "windows.tsv"
    path.write_bytes(b"accession\tgroup\nGCA_1.1\r\tReproductive\n")
    assert read_tsv(path) == [{"accession": "GCA_1.1", "group": "Reproductive"}]
