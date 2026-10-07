"""Create an NPR1c alignment diagram showing three aligned NbL genes, mapped sgRNAs, and optional VIGS annotation.

This plot renders the three aligned NbL gene tracks from a gapped FASTA. The
reference gene NbL05g11340 is used for the coordinate axis, and sgRNA sites
are mapped to the ungapped reference sequence.
"""
from pathlib import Path
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from PRp27.PRp27_timecourse_202606.figure_styles import ANALYSIS_RCPARAMS, FONT_STACK, TEXT_COLOR


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def read_fasta(path: Path):
    seqs = {}
    current = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            seqs[current] = ""
        else:
            seqs[current] += line.strip()
    return seqs


def find_subseq_positions(ungapped: str, query: str):
    """Return list of (start,end) 1-based positions for exact matches in ungapped sequence."""
    q = query.upper()
    hits = []
    i = 0
    while True:
        idx = ungapped.find(q, i)
        if idx == -1:
            break
        hits.append((idx + 1, idx + len(q)))
        i = idx + 1
    return hits


def map_ungapped_to_gapped(ref_seq: str):
    mapping = {}
    ungapped_coord = 0
    for gapped_pos, base in enumerate(ref_seq, start=1):
        if base != '-':
            ungapped_coord += 1
            mapping[ungapped_coord] = gapped_pos
    return mapping


def gapped_blocks(seq: str):
    blocks = []
    start = None
    for pos, base in enumerate(seq, start=1):
        if base != '-' and start is None:
            start = pos
        elif base == '-' and start is not None:
            blocks.append((start, pos - 1))
            start = None
    if start is not None:
        blocks.append((start, len(seq)))
    return blocks


def plot_alignment_schematic(gene_seqs, ref_name, sg_list, vigs_region, outdir: Path):
    ref_seq = gene_seqs[ref_name]
    ref_ungapped = ref_seq.replace('-', '')
    length = len(ref_seq)
    coord_map = map_ungapped_to_gapped(ref_seq)

    plt.rcParams.update(ANALYSIS_RCPARAMS)
    fig, ax = plt.subplots(figsize=(10, 4.0))
    ax.set_xlim(0.5, length + 0.5)
    ax.set_ylim(0, 5.0)
    ax.axis('off')

    row_height = 0.40
    track_gap = 0.92
    y_start = 4.0
    label_x = -0.06
    gene_names = list(gene_seqs.keys())

    for idx, gene_name in enumerate(gene_names):
        y = y_start - idx * track_gap
        for start, end in gapped_blocks(gene_seqs[gene_name]):
            rect = patches.Rectangle((start - 0.5, y - row_height / 2), end - start + 1, row_height,
                                     facecolor='#222222', edgecolor='none')
            ax.add_patch(rect)
        ax.text(label_x, y, gene_name, ha='right', va='center', fontsize=11, fontfamily=FONT_STACK,
                color='#222222', transform=ax.transAxes)

    if vigs_region is not None:
        vs, ve = vigs_region
        rect = patches.Rectangle((vs - 0.5, y_start - row_height / 2 - 0.22), ve - vs + 1,
                                 track_gap * len(gene_names) - 0.22,
                                 facecolor='#D8E1E8', alpha=0.34, edgecolor='none')
        ax.add_patch(rect)
        ax.text((vs + ve) / 2, y_start + 0.58, 'VIGS fragment', ha='center', va='bottom', fontsize=10,
                fontfamily=FONT_STACK, color='#222222')

    sg_colors = {'sg1': '#4A4A4A', 'sg2': '#6C6C6C', 'sg3': '#8A6F58'}
    marker_y = y_start + 0.20
    for name, hits in sg_list.items():
        for st, en, strand in hits:
            pos_start = coord_map.get(st)
            pos_end = coord_map.get(en)
            if pos_start is None or pos_end is None:
                continue
            x = (pos_start + pos_end) / 2
            ax.scatter([x], [marker_y], marker='v', color=sg_colors.get(name), s=100, zorder=6)
            ax.text(x, marker_y + 0.16, name, ha='center', va='bottom', fontsize=10, fontfamily=FONT_STACK,
                    color='#333333')

    axis_y = y_start + 0.55
    ax.plot([0.5, length + 0.5], [axis_y, axis_y], color='#444444', linewidth=0.8)

    tick_coords = list(range(1, len(ref_ungapped) + 1, max(1, len(ref_ungapped) // 6)))
    if tick_coords[-1] != len(ref_ungapped):
        tick_coords.append(len(ref_ungapped))
    for coord in tick_coords:
        pos = coord_map.get(coord)
        if pos is None:
            continue
        ax.plot([pos, pos], [axis_y - 0.08, axis_y + 0.08], color='#444444', linewidth=0.8)
        ax.text(pos, axis_y + 0.14, str(coord), ha='center', va='bottom', fontsize=9,
                fontfamily=FONT_STACK, color='#3C3C3C')

    ax.text(length / 2, axis_y + 0.36, 'Reference coordinates (NbL05g11340)',
            ha='center', va='bottom', fontsize=10, fontfamily=FONT_STACK, color='#3C3C3C')
    ax.text(0.01, 0.98, 'A', transform=ax.transAxes, fontsize=20, fontweight='bold', va='top',
            fontfamily=FONT_STACK, color='#222222')
    ax.text(0.01, 0.90, 'Aligned NbL genes with reference-mapped sgRNAs',
            transform=ax.transAxes, ha='left', va='top', fontsize=11, fontfamily=FONT_STACK,
            color='#222222')

    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / '202607_npr1c_alignment_schematic.png'
    pdf = outdir / '202607_npr1c_alignment_schematic.pdf'
    fig.savefig(png, bbox_inches='tight', dpi=300)
    fig.savefig(pdf, bbox_inches='tight')
    plt.close(fig)
    return png, pdf


def main():
    parser = argparse.ArgumentParser(description='Create a simplified NPR1c alignment diagram.')
    parser.add_argument('--vigs', type=Path, help='Path to a text file containing VIGS fragment sequence.')
    parser.add_argument('--vigs-seq', type=str, help='VIGS fragment sequence as plain text.')
    parser.add_argument('--fasta', type=Path, default=Path('NPR1/npr1c/alignment.fasta'))
    args = parser.parse_args()

    fasta = args.fasta
    if not fasta.exists():
        raise SystemExit(f'Fasta not found: {fasta}')

    seqs = read_fasta(fasta)
    genes = ['NbL05g11340', 'NbL05g11410', 'NbL18g20120']
    for gene in genes:
        if gene not in seqs:
            raise SystemExit(f'Required gene {gene} not found in {fasta}')

    gene_seqs = {gene: seqs[gene] for gene in genes}
    ref_seq = gene_seqs['NbL05g11340']
    ref_ungapped = ref_seq.replace('-', '').upper()

    sg_seqs = {
        'sg1': 'TTGATAAGGCCTTGCCTCAT',
        'sg2': 'GAGCAGAACTTGGTCTACAA',
        'sg3': 'AACATGTTAAGAGGATACAT',
    }
    sg_list = {}
    for name, sg in sg_seqs.items():
        hits = []
        for st, en in find_subseq_positions(ref_ungapped, sg):
            hits.append((st, en, '+'))
        for st, en in find_subseq_positions(ref_ungapped, revcomp(sg)):
            hits.append((st, en, '-'))
        sg_list[name] = hits

    vigs_region = None
    vigs_seq = None
    if args.vigs:
        vigs_seq = args.vigs.read_text()
    elif args.vigs_seq:
        vigs_seq = args.vigs_seq
    if vigs_seq:
        seqtxt = vigs_seq.upper().replace('\n', '').replace(' ', '')
        if seqtxt:
            v_hits = find_subseq_positions(ref_ungapped, seqtxt)
            if v_hits:
                vigs_region = v_hits[0]
            else:
                print('No exact VIGS match found in reference; continuing without VIGS highlight.')

    outdir = Path('NPR1/npr1c/figures')
    png, pdf = plot_alignment_schematic(gene_seqs, 'NbL05g11340', sg_list, vigs_region, outdir)

    methods = outdir.parent / '202607_npr1c_alignment_methods_note.txt'
    methods.write_text(
        'Plot renders three aligned NbL genes from a gapped FASTA as horizontal tracks.\n'
        'sgRNA sites are mapped to exact matches on the ungapped NbL05g11340 reference.\n'
        'The VIGS fragment is highlighted if an exact reference match is provided.\n'
    )
    print('Generated:', png, pdf, methods)


if __name__ == '__main__':
    main()
