from __future__ import annotations

import csv
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Tuple


FASTA_EXTENSIONS = {".fa", ".fna", ".fasta"}


@dataclass(frozen=True)
class Allele:
    locus: str
    allele_id: int
    seq: str


@dataclass(frozen=True)
class Hit:
    qname: str
    locus: str
    allele_id: int
    tname: str
    tstart: int
    tend: int
    strand: str
    identity: float


def iter_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    header = None
    chunks: List[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def parse_allele_header(
    qname: str, allow_locus_only: bool, locus_hint: Optional[str]
) -> Tuple[str, Optional[int]]:
    if "_allele_" in qname:
        locus, allele_str = qname.rsplit("_allele_", 1)
        if not locus or not allele_str.isdigit():
            raise ValueError(
                f"Allele header '{qname}' must end with numeric allele ID (e.g. ABC_0001_allele_12)."
            )
        return locus, int(allele_str)
    match = re.match(r"^(.+)_([0-9]+)$", qname)
    if match:
        locus, allele_str = match.groups()
        if locus:
            return locus, int(allele_str)
    if qname.isdigit() and locus_hint:
        return locus_hint, int(qname)
    if allow_locus_only:
        if not qname:
            raise ValueError("Allele header is empty; expected a locus name.")
        return qname, None
    raise ValueError(
        (
            f"Allele header '{qname}' must be 'locus_allele_ID' or 'locus_ID'. "
            "For directory inputs, numeric-only headers are accepted when the filename is the locus. "
            "Use --allow-locus-only for locus-only headers."
        )
    )


def collect_allele_files(alleles_path: Path) -> List[Path]:
    if alleles_path.is_file():
        return [alleles_path]
    if not alleles_path.is_dir():
        raise FileNotFoundError(f"Allele input '{alleles_path}' was not found.")
    allele_files = sorted(
        [p for p in alleles_path.iterdir() if p.suffix.lower() in FASTA_EXTENSIONS]
    )
    if not allele_files:
        raise ValueError(f"No FASTA files found in allele directory '{alleles_path}'.")
    return allele_files


def parse_alleles(
    allele_files: Iterable[Path],
    allow_locus_only: bool,
    locus_from_filename: bool,
) -> Tuple[Dict[str, Allele], List[str], Dict[str, int], List[str]]:
    qname_to_allele: Dict[str, Allele] = {}
    locus_order: List[str] = []
    qname_order: List[str] = []
    max_allele_id: Dict[str, int] = {}
    seen_ids: Dict[str, set] = {}

    for allele_file in allele_files:
        locus_hint = allele_file.stem if locus_from_filename else None
        for header, seq in iter_fasta(allele_file):
            qname = header.split()[0]
            locus, allele_id = parse_allele_header(qname, allow_locus_only, locus_hint)
            if locus not in seen_ids:
                seen_ids[locus] = set()
            if allele_id is None:
                allele_id = max_allele_id.get(locus, 0) + 1
            qname = f"{locus}_allele_{allele_id}"
            if qname in qname_to_allele:
                raise ValueError(f"Duplicate allele ID '{qname}' in input files.")
            if allele_id in seen_ids[locus]:
                raise ValueError(
                    f"Duplicate allele number {allele_id} for locus '{locus}'."
                )
            seen_ids[locus].add(allele_id)
            if locus not in locus_order:
                locus_order.append(locus)
            max_allele_id[locus] = max(max_allele_id.get(locus, 0), allele_id)
            qname_order.append(qname)
            qname_to_allele[qname] = Allele(locus=locus, allele_id=allele_id, seq=seq.upper())

    if not qname_to_allele:
        raise ValueError("No alleles were parsed from input FASTA files.")

    return qname_to_allele, locus_order, max_allele_id, qname_order


def write_combined_fasta(
    output_path: Path, qname_order: Iterable[str], qname_to_allele: Dict[str, Allele]
) -> None:
    with output_path.open("w") as handle:
        for qname in qname_order:
            allele = qname_to_allele[qname]
            handle.write(f">{qname}\n")
            seq = allele.seq
            for i in range(0, len(seq), 60):
                handle.write(seq[i : i + 60] + "\n")


def collect_genomes(genome_paths: Iterable[str], genome_dir: Optional[str]) -> List[Path]:
    genomes: List[Path] = []
    for genome in genome_paths:
        path = Path(genome)
        if not path.exists():
            raise FileNotFoundError(f"Genome file '{path}' was not found.")
        genomes.append(path)
    if genome_dir:
        gdir = Path(genome_dir)
        if not gdir.is_dir():
            raise FileNotFoundError(f"Genome directory '{gdir}' was not found.")
        genomes.extend(sorted([p for p in gdir.iterdir() if p.suffix.lower() in FASTA_EXTENSIONS]))
    unique: List[Path] = []
    seen: set = set()
    for genome in genomes:
        if genome in seen:
            continue
        seen.add(genome)
        unique.append(genome)
    if not unique:
        raise ValueError("No genome FASTA files provided.")
    return unique


def require_minimap2() -> None:
    if shutil.which("minimap2") is None:
        raise RuntimeError("minimap2 was not found in PATH. Please install minimap2.")


def run_minimap2_index(reference_fasta: Path, index_path: Path) -> None:
    subprocess.run(["minimap2", "-d", str(index_path), str(reference_fasta)], check=True)


def run_minimap2_map(index_path: Path, query_fasta: Path, paf_path: Path, threads: int) -> None:
    with paf_path.open("w") as handle:
        subprocess.run(
            [
                "minimap2",
                "-x",
                "asm10",
                "-c",
                "-t",
                str(threads),
                str(index_path),
                str(query_fasta),
            ],
            stdout=handle,
            check=True,
        )


def extract_cigar(tags: List[str]) -> Optional[str]:
    for tag in tags:
        if tag.startswith("cg:Z:"):
            return tag[5:]
    return None


def cigar_has_indel_or_clip(cigar: str) -> bool:
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op not in {"M", "=", "X"}:
            return True
    return False


def parse_paf(
    paf_path: Path,
    qname_to_allele: Dict[str, Allele],
    min_identity: float,
    strict: bool,
) -> Dict[str, List[Hit]]:
    hits_by_locus: Dict[str, List[Hit]] = {}
    with paf_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            qname = fields[0]
            allele = qname_to_allele.get(qname)
            if allele is None:
                continue
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            strand = fields[4]
            tname = fields[5]
            tstart = int(fields[7])
            tend = int(fields[8])
            nmatch = int(fields[9])
            alnlen = int(fields[10])

            if qlen != len(allele.seq):
                continue
            if qstart != 0 or qend != qlen:
                continue
            if alnlen <= 0:
                continue
            identity = nmatch / alnlen
            if identity < min_identity:
                continue
            if strict:
                cigar = extract_cigar(fields[12:])
                if cigar is None:
                    continue
                if cigar_has_indel_or_clip(cigar):
                    continue

            hit = Hit(
                qname=qname,
                locus=allele.locus,
                allele_id=allele.allele_id,
                tname=tname,
                tstart=tstart,
                tend=tend,
                strand=strand,
                identity=identity,
            )
            hits_by_locus.setdefault(allele.locus, []).append(hit)

    return hits_by_locus


def read_genome(genome_path: Path) -> Tuple[Dict[str, str], int]:
    sequences: Dict[str, str] = {}
    total_length = 0
    for header, seq in iter_fasta(genome_path):
        name = header.split()[0]
        if name in sequences:
            raise ValueError(f"Duplicate contig name '{name}' in {genome_path}.")
        seq_upper = seq.upper()
        sequences[name] = seq_upper
        total_length += len(seq_upper)
    if not sequences:
        raise ValueError(f"No sequences found in genome FASTA '{genome_path}'.")
    return sequences, total_length


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def call_alleles_for_genome(
    locus_order: List[str],
    hits_by_locus: Dict[str, List[Hit]],
    qname_to_allele: Dict[str, Allele],
    genome_seqs: Dict[str, str],
    max_allele_id: Dict[str, int],
    novel_seq_by_locus: Dict[str, Dict[str, int]],
    missing_value: str,
    novel_records: Dict[Tuple[str, int], str],
) -> Tuple[List[str], int]:
    calls: List[str] = []
    good_targets = 0
    identity_tolerance = 1e-6

    for locus in locus_order:
        hits = hits_by_locus.get(locus, [])
        if not hits:
            calls.append(missing_value)
            continue
        hits_sorted = sorted(hits, key=lambda h: h.identity, reverse=True)
        best = hits_sorted[0]
        tied = [
            hit
            for hit in hits_sorted
            if abs(hit.identity - best.identity) <= identity_tolerance
        ]
        if len(tied) != 1:
            calls.append(missing_value)
            continue
        hit = best
        allele = qname_to_allele[hit.qname]
        contig_seq = genome_seqs.get(hit.tname)
        if contig_seq is None:
            calls.append(missing_value)
            continue

        genome_seq = contig_seq[hit.tstart : hit.tend]
        if hit.strand == "-":
            genome_seq = reverse_complement(genome_seq)
        genome_seq = genome_seq.upper()
        allele_seq = allele.seq.upper()

        if len(genome_seq) != len(allele_seq):
            calls.append(missing_value)
            continue

        if genome_seq == allele_seq:
            calls.append(str(allele.allele_id))
            good_targets += 1
            continue

        locus_novels = novel_seq_by_locus.setdefault(locus, {})
        allele_id = locus_novels.get(genome_seq)
        if allele_id is None:
            max_allele_id[locus] += 1
            allele_id = max_allele_id[locus]
            locus_novels[genome_seq] = allele_id
            novel_records[(locus, allele_id)] = genome_seq
        calls.append(str(allele_id))
        good_targets += 1

    return calls, good_targets


def format_percent(value: float) -> str:
    return f"{value:.1f}"


def format_genome_size(total_length: int) -> str:
    return f"{total_length / 1_000_000:.1f}"


def build_output_header(locus_order: List[str]) -> List[str]:
    meta = [
        "Avg. Coverage (Assembled)",
        "Approximated Genome Size (Mbases)",
        "Genome positions N called",
        "Top Species Match Identity",
        "Top Species Match",
        "Sample ID",
        "Epi Info",
        "Cluster/Outbreak",
        "Collection Date",
        "Country of Isolation",
        "City of Isolation",
        "ZIP of Isolation",
        "Lat/Long of Isolation ",
        "Lat/Long Resolution",
        "Perc. Good Targets",
    ]
    return meta + locus_order


def write_output_csv(output_path: Path, header: List[str], rows: List[List[str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", quoting=csv.QUOTE_ALL, lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def write_novel_fasta(novel_path: Path, novel_records: Dict[Tuple[str, int], str]) -> None:
    if not novel_records:
        return
    with novel_path.open("w") as handle:
        for (locus, allele_id), seq in sorted(novel_records.items()):
            handle.write(f">{locus}_allele_{allele_id} novel\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i : i + 60] + "\n")


def run_pipeline(
    alleles_path: str,
    genome_paths: Iterable[str],
    genome_dir: Optional[str],
    output_path: str,
    min_identity: float,
    strict: bool,
    threads: int,
    missing_value: str,
    novel_fasta: Optional[str],
    tmp_dir: Optional[str],
    keep_temp: bool,
    allow_locus_only: bool,
    max_genomes: Optional[int],
    progress_callback: Optional[Callable[[int, int, str], None]],
) -> None:
    require_minimap2()

    allele_input = Path(alleles_path)
    allele_files = collect_allele_files(allele_input)
    qname_to_allele, locus_order, max_allele_id, qname_order = parse_alleles(
        allele_files, allow_locus_only, allele_input.is_dir()
    )

    genomes = collect_genomes(genome_paths, genome_dir)
    if max_genomes is not None:
        genomes = genomes[:max_genomes]
    total_genomes = len(genomes)

    if keep_temp:
        work_dir = Path(tmp_dir) if tmp_dir else Path(tempfile.mkdtemp(prefix="cgmlst_"))
        work_dir.mkdir(parents=True, exist_ok=True)
        cleanup = False
    else:
        work_dir_obj = tempfile.TemporaryDirectory(prefix="cgmlst_", dir=tmp_dir)
        work_dir = Path(work_dir_obj.name)
        cleanup = True

    try:
        allele_fasta = work_dir / "cgmlst_alleles.fasta"
        write_combined_fasta(allele_fasta, qname_order, qname_to_allele)

        header = build_output_header(locus_order)
        rows: List[List[str]] = []
        novel_seq_by_locus: Dict[str, Dict[str, int]] = {}
        novel_records: Dict[Tuple[str, int], str] = {}

        for idx, genome_path in enumerate(genomes, start=1):
            genome_path = Path(genome_path)
            sample_id = genome_path.stem
            safe_id = re.sub(r"[^A-Za-z0-9._-]+", "_", sample_id)
            paf_path = work_dir / f"{safe_id}.paf"
            genome_index = work_dir / f"{safe_id}.mmi"

            run_minimap2_index(genome_path, genome_index)
            run_minimap2_map(genome_index, allele_fasta, paf_path, threads)
            hits_by_locus = parse_paf(paf_path, qname_to_allele, min_identity, strict)
            genome_seqs, genome_length = read_genome(genome_path)

            calls, good_targets = call_alleles_for_genome(
                locus_order,
                hits_by_locus,
                qname_to_allele,
                genome_seqs,
                max_allele_id,
                novel_seq_by_locus,
                missing_value,
                novel_records,
            )

            perc_good = 0.0
            if locus_order:
                perc_good = good_targets / len(locus_order) * 100

            row = [
                "",  # Avg. Coverage (Assembled)
                format_genome_size(genome_length),
                "",  # Genome positions N called
                "",  # Top Species Match Identity
                "",  # Top Species Match
                sample_id,
                "",  # Epi Info
                "",  # Cluster/Outbreak
                "",  # Collection Date
                "",  # Country of Isolation
                "",  # City of Isolation
                "",  # ZIP of Isolation
                "",  # Lat/Long of Isolation
                "",  # Lat/Long Resolution
                format_percent(perc_good),
            ] + calls
            rows.append(row)
            if progress_callback is not None:
                progress_callback(idx, total_genomes, sample_id)

        write_output_csv(Path(output_path), header, rows)

        if novel_fasta is None:
            base = Path(output_path)
            novel_fasta = str(base.with_suffix("")) + "_novel_alleles.fasta"
        if novel_records:
            write_novel_fasta(Path(novel_fasta), novel_records)

    finally:
        if cleanup:
            work_dir_obj.cleanup()

    if keep_temp:
        print(f"Temporary files kept in {work_dir}", file=sys.stderr)
