import subprocess
from Bio import SeqIO
from io import StringIO


class SequenceExtractor:
    def extract(self, reads_to_extract, database_fasta_file, output_file):
        """Extract the reads_to_extract from the database_fasta_file and put them in
        output_file.

        Parameters
        ----------
        reads_to_extract: Iterable of str
            IDs of reads to be extracted
        database_fasta_file: str
            path the fasta file that containing the reads
        output_file: str
            path to the file where they are put

        Returns
        -------
        Nothing"""
        cmd = (
            "mfqe --fasta-read-name-lists /dev/stdin "
            "--input-fasta-files {0} "
            "--output-fasta-files {1} "
            "--output-uncompressed"
        ).format(database_fasta_file, output_file)

        # früher: extern.run(cmd, stdin='\n'.join(reads_to_extract))
        subprocess.run(
            cmd,
            input="\n".join(reads_to_extract),  # geht an /dev/stdin von mfqe
            text=True,  # input als String, nicht als Bytes
            shell=True,  # weil cmd ein String ist
            check=True,  # Fehler, wenn Exitcode != 0
        )

    def extract_forward_and_reverse_complement(
        self,
        forward_reads_to_extract,
        reverse_reads_to_extract,
        database_fasta_file,
        output_file,
    ):
        """As per extract except also reverse complement the sequences."""
        self.extract(forward_reads_to_extract, database_fasta_file, output_file)
        cmd_rev = (
            "mfqe --fasta-read-name-lists /dev/stdin --input-fasta-files {0} "
            "--output-fasta-files /dev/stdout --output-uncompressed".format(
                database_fasta_file
            )
        )

        # extern.run(cmd_rev, stdin='...\n...')
        proc = subprocess.run(
            cmd_rev,
            input="\n".join(reverse_reads_to_extract),  # entspricht stdin=...
            text=True,  # stdin/stdout als str statt bytes
            shell=True,  # weil cmd_rev ein String ist
            check=True,  # Fehler, falls Exitcode != 0
            capture_output=True,  # stdout/stderr abfangen
        )
        output = proc.stdout

        with open(output_file, "a") as f:
            for record in SeqIO.parse(StringIO(output), "fasta"):
                record.seq = record.reverse_complement().seq
                SeqIO.write(record, f, "fasta")
