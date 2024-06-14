from io import StringIO

from Bio import SeqIO
from django.test import SimpleTestCase


class TestDownloadSequencesView(SimpleTestCase):

    prot_seq = 'METTKPSFQDVLEFVRLYRRKNKLQREIQDVEKKIRDNQKRVLLLDNLSDYIKPGMSVEAIQGIIASMKSDYEDRVDDYIIKNAELSKERRDISKKLKVMGEAKVEG'

    def test_get_is_forbidden(self):
        resp = self.client.get("/download_sequences")
        self.assertEqual(405,  resp.status_code)
        self.assertEqual('Method Not Allowed', resp.reason_phrase)

    def test_download_prot_sequence(self):
        resp = self.client.post("/download_sequences",
                                {"loci": ["CHUV_00025"]},
                                content_type="application/json")
        seq = SeqIO.read(StringIO(resp.content.decode("utf8")), "fasta")
        self.assertEqual("CHUV_00025", seq.id)
        self.assertEqual(self.prot_seq, str(seq.seq))

    def test_download_prot_sequences(self):
        resp = self.client.post("/download_sequences",
                                {"loci": ["CHUV_00025", "CHUV_00010"]},
                                content_type="application/json")
        sequences = SeqIO.parse(StringIO(resp.content.decode("utf8")), "fasta")
        self.assertEqual(["CHUV_00025", "CHUV_00010"],
                         [seq.id for seq in sequences])

    def test_download_dna_sequence(self):
        resp = self.client.post("/download_sequences",
                                {"loci": ["CHUV_00025"], "dna": True},
                                content_type="application/json")
        seq = SeqIO.read(StringIO(resp.content.decode("utf8")), "fasta")
        self.assertEqual("CHUV_00025", seq.id)
        # Drop last character of sequence after translation as it's a *
        self.assertEqual(self.prot_seq, str(seq.translate().seq)[:-1])
