"""Reproducible protein CA evaluation, keeping map placement and alignment separate.

All matching is one-to-one at residue level. The sequence-registered metric uses
deposited construct sequences, not author residue numbers in predicted models.
It measures residue placement, not peptide connectivity or whole-chain integrity.
"""
from pathlib import Path
import argparse, hashlib, json, re, subprocess
import gemmi
import numpy as np
from Bio.Align import PairwiseAligner
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching


def aa(name):
    info = gemmi.find_tabulated_residue(name)
    code = info.one_letter_code.upper()
    return code if len(code) == 1 and code.isalpha() else 'X'


def read_ca(path, reference=False):
    structure = gemmi.read_structure(str(path))
    if not len(structure):
        raise ValueError('No model in structure')
    constructs = {}
    if reference:
        for entity in structure.entities:
            if entity.full_sequence:
                constructs[entity.name] = ''.join(aa(x.split(',')[0]) for x in entity.full_sequence)
        if not constructs:
            raise ValueError('Reference requires deposited construct sequences in mmCIF')
    chains = []
    records = []
    for chain in structure[0]:
        groups = {}
        for residue in chain:
            if not gemmi.find_tabulated_residue(residue.name).is_amino_acid():
                continue
            atoms = [a for a in residue if a.name == 'CA']
            if not atoms:
                continue
            atom = sorted(atoms, key=lambda a: (-a.occ, a.altloc not in ['\x00', 'A'], a.altloc))[0]
            if not np.isfinite([atom.pos.x, atom.pos.y, atom.pos.z]).all():
                raise ValueError('Nonfinite CA coordinate')
            key = (chain.name, residue.subchain)
            rec = dict(aa=aa(residue.name), xyz=[atom.pos.x, atom.pos.y, atom.pos.z],
                       entity=residue.entity_id, position=residue.label_seq,
                       auth_residue=str(residue.seqid), chain=chain.name,
                       subchain=residue.subchain, index=len(records))
            if reference:
                seq = constructs.get(rec['entity'])
                if not seq or not rec['position'] or not 1 <= rec['position'] <= len(seq):
                    raise ValueError('Missing/invalid reference entity or label_seq_id')
                expected = seq[rec['position']-1]
                if rec['aa'] != expected and 'X' not in (rec['aa'], expected):
                    raise ValueError('Reference coordinate sequence disagrees with construct')
            records.append(rec)
            groups.setdefault(key, []).append(rec)
        for key, values in groups.items():
            chains.append(dict(id='|'.join(key), records=values, sequence=''.join(r['aa'] for r in values)))
    if not records:
        raise ValueError('No protein CA coordinates')
    return dict(records=records, chains=chains, constructs=constructs,
                xyz=np.asarray([r['xyz'] for r in records]), source=str(path))


def write_ca_cif(data, path):
    """Canonicalize one CA per protein residue, removing ligand/alternate ambiguity."""
    structure = gemmi.Structure()
    structure.name = 'protein_CA_evaluation'
    model = gemmi.Model('1')
    for ci, original in enumerate(data['chains']):
        chain = gemmi.Chain('C'+str(ci+1))
        for ri, row in enumerate(original['records']):
            residue = gemmi.Residue()
            residue.name = gemmi.expand_one_letter(row['aa'], gemmi.ResidueKind.AA) if row['aa'] != 'X' else 'UNK'
            residue.seqid = gemmi.SeqId(ri+1, ' ')
            residue.het_flag = 'A'
            atom = gemmi.Atom(); atom.name = 'CA'; atom.element = gemmi.Element('C')
            atom.pos = gemmi.Position(*row['xyz']); atom.occ = 1; atom.b_iso = 0
            residue.add_atom(atom); chain.add_residue(residue)
        model.add_chain(chain)
    structure.add_model(model)
    structure.setup_entities()
    structure.make_mmcif_document().write_file(str(path))


def maximum_match(edges, n_pred, n_ref):
    rows, cols = [], []
    for i, partners in enumerate(edges):
        rows.extend([i]*len(partners)); cols.extend(partners)
    graph = csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n_pred, n_ref))
    assignment = maximum_bipartite_matching(graph, perm_type='column')
    return int(np.count_nonzero(assignment >= 0))


def fractions(count, n_pred, n_ref):
    precision, recall = (count/n_pred if n_pred else 0.0), count/n_ref
    return dict(matched=count, precision=precision, recall=recall,
                f1=2*precision*recall/(precision+recall) if count else 0.0)


def register_prediction(predicted, reference):
    aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1,
                              open_gap_score=-3, extend_gap_score=-0.1)
    aligner.end_gap_score = 0
    allowed = [set() for _ in predicted['records']]
    provenance = []
    for chain in predicted['chains']:
        choices = []
        for entity, sequence in sorted(reference['constructs'].items()):
            alignments = aligner.align(sequence, chain['sequence'])
            first = alignments[0]
            try:
                alignments[1]; ambiguous = True
            except IndexError:
                ambiguous = False
            choices.append((first.score, entity, first, ambiguous))
        best_score = max(x[0] for x in choices)
        chosen = [x for x in choices if abs(x[0]-best_score) < 1e-8]
        for score, entity, alignment, ambiguous in chosen:
            for (rs, re_), (ps, pe) in zip(*alignment.aligned):
                assert re_-rs == pe-ps
                for rp, pp in zip(range(rs, re_), range(ps, pe)):
                    allowed[chain['records'][pp]['index']].add((entity, rp+1))
        provenance.append(dict(chain=chain['id'], best_alignment_score=float(best_score),
                               tied_entities=[x[1] for x in chosen],
                               alternative_optimal_alignment=any(x[3] for x in chosen)))
    return allowed, provenance


def fixed_coordinate_scores(predicted, reference, threshold=3.0):
    n_pred, n_ref = len(predicted['records']), len(reference['records'])
    neighbors = cKDTree(reference['xyz']).query_ball_point(predicted['xyz'], threshold)
    geometry = maximum_match(neighbors, n_pred, n_ref)
    typed = [[j for j in js if predicted['records'][i]['aa'] == reference['records'][j]['aa']
              and predicted['records'][i]['aa'] != 'X'] for i, js in enumerate(neighbors)]
    aa_typed = maximum_match(typed, n_pred, n_ref)
    allowed, provenance = register_prediction(predicted, reference)
    registered = [[j for j in js if (reference['records'][j]['entity'], reference['records'][j]['position'])
                   in allowed[i]] for i, js in enumerate(typed)]
    seq_correct = maximum_match(registered, n_pred, n_ref)
    observed_positions = {(r['entity'], r['position']) for r in reference['records']}
    evaluable_pred = sum(bool(x & observed_positions) for x in allowed)
    outside_reference = sum(bool(x) and not bool(x & observed_positions) for x in allowed)
    unregistered = sum(not bool(x) for x in allowed)
    return dict(threshold_angstrom=threshold, superposition=False, predicted_ca=n_pred,
                reference_ca=n_ref, geometry_one_to_one=fractions(geometry, n_pred, n_ref),
                aa_type_one_to_one=fractions(aa_typed, n_pred, n_ref),
                construct_registered_one_to_one=fractions(seq_correct, evaluable_pred, n_ref),
                sequence_evaluable_predicted_ca=evaluable_pred,
                predicted_at_unobserved_construct_positions=outside_reference,
                unregistered_predicted_ca=unregistered,
                all_prediction_reference_supported_fraction=seq_correct/n_pred,
                registration=provenance,
                definition='Maximum-cardinality residue matching within3A, same amino-acid type and construct position. Recall uses observed reference residues; sequence precision uses predictions registered to reference-observed positions. Unobserved construct positions and unregistered predictions are counted separately. Geometry/AA-type precision use all predictions. Homomer copies may interchange. Peptide connectivity is not scored. First optimal sequence alignment is deterministic; ties are reported.')


def usalign_scores(executable, predicted_cif, reference_cif, output_dir, n_pred, n_ref):
    command = [str(executable), str(predicted_cif), str(reference_cif), '-mm', '1',
               '-ter', '1', '-split', '2', '-mol', 'prot', '-outfmt', '2']
    result = subprocess.run(command, capture_output=True, text=True, timeout=600)
    (output_dir/'usalign.log').write_text(result.stdout+result.stderr)
    result.check_returncode()
    lines = result.stdout.splitlines()
    header = next(line for line in lines if line.startswith('#PDBchain1')).lstrip('#').split('\t')
    values = [line for line in lines if line and not line.startswith('#')][-1].split('\t')
    row = dict(zip(header, values))
    # Native USalign excludes chains shorter than three residues. Preserve its
    # prediction-normalized score and label its actual denominator; require the
    # primary reference-normalized score to include every reference CA.
    st=gemmi.read_structure(str(predicted_cif))
    excluded=sum(len(c) for c in st[0] if len(c)<3)
    assert int(row['L1']) == n_pred-excluded and int(row['L2']) == n_ref, 'Unexpected USalign CA denominator mismatch'

    return dict(command=command, tm_reference=float(row['TM2']), tm_prediction=float(row['TM1']),
                aligned_ca=int(row['Lali']), aligned_sequence_identity=float(row['IDali']),
                aligned_rmsd=float(row['RMSD']), alignment='sequence-independent, one global assembly transform',
                raw_table=row, prediction_ca_total=n_pred, prediction_ca_omitted_short_chains=excluded, tm_prediction_normalization_length=int(row['L1']), tm_reference_normalization_length=n_ref, short_chain_policy='Native USalign excludes <3-residue prediction chains from structural alignment; all predictions retained in fixed-coordinate metrics. Reference-normalized primary TM score requires full reference denominator.')


def phenix_scores(activation, output_dir, n_pred, n_ref):
    command=['bash','-c','source "$1"; shift; exec phenix.chain_comparison "$@"','--',str(activation),
             str(output_dir/'reference_ca.cif'),str(output_dir/'prediction_ca.cif'),
             'unique_query_only=False','unique_part_of_target_only=False',
             'test_unique_part_of_target_only=False','use_crystal_symmetry=False','max_dist=3.0']
    run=subprocess.run(command,cwd=output_dir,capture_output=True,text=True,timeout=180)
    (output_dir/'phenix.log').write_text(run.stdout+run.stderr);run.check_returncode()
    line=next((x for x in run.stdout.splitlines() if x.startswith('All residues near target:')),None)
    if line is None:
        return dict(status='requires_parser_review',command=command)
    near=re.search(r'All residues near target:\s*(\d+)',line)
    seq=re.search(r'Seq match \(%\):\s*([0-9.]+)',line)
    count=int(near[1]);precision=count/n_pred
    return dict(status='parsed',command=command,matched_ca=count,ca_match=precision,
                ca_quality=precision*n_pred/n_ref,sequence_match=float(seq[1])/100 if seq else None,
                raw_summary=line,threshold_angstrom=3.0,superposition=False)


def evaluate(predicted_path, reference_path, output_dir, usalign=None, phenix=None):
    output_dir = Path(output_dir); output_dir.mkdir(parents=True, exist_ok=True)
    if (output_dir/'metrics.json').exists():
        raise FileExistsError('Preserve prior scores; use a new evaluation directory')
    pred, ref = read_ca(predicted_path), read_ca(reference_path, reference=True)
    pc, rc = output_dir/'prediction_ca.cif', output_dir/'reference_ca.cif'
    write_ca_cif(pred, pc); write_ca_cif(ref, rc)
    result = dict(prediction=str(predicted_path), reference=str(reference_path),
                  prediction_sha256=hashlib.sha256(Path(predicted_path).read_bytes()).hexdigest(),
                  reference_sha256=hashlib.sha256(Path(reference_path).read_bytes()).hexdigest(),
                  fixed_map_coordinates=fixed_coordinate_scores(pred, ref))
    if usalign:
        result['usalign'] = usalign_scores(usalign, pc, rc, output_dir, len(pred['records']), len(ref['records']))
    if phenix:
        result['phenix'] = phenix_scores(phenix, output_dir, len(pred['records']), len(ref['records']))
    (output_dir/'metrics.json').write_text(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--prediction', required=True); parser.add_argument('--reference', required=True)
    parser.add_argument('--output', required=True); parser.add_argument('--usalign'); parser.add_argument('--phenix')
    args = parser.parse_args()
    print(json.dumps(evaluate(args.prediction, args.reference, args.output, args.usalign,args.phenix), indent=2))
