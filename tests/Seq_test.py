def test_is_protein(protein_sequence, dna_sequence):
    assert protein_sequence.is_protein() == True
    assert dna_sequence.is_protein() == False


def test_re_search(protein_sequence, dna_sequence):
    assert {m.group(): m.start() for m in dna_sequence.re_search(r'TAG')} == {'TAG': 2}
    assert {m.group(): m.start() for m in dna_sequence.re_search(r'CG{2}')} == {'CGG': 9}
    assert {m.group(): m.start() for m in dna_sequence.re_search(r'(GC){2,3}')} == {'GCGCGC': 4}
    assert {m.group(): m.start() for m in dna_sequence.re_search(r'[GC]{2,10}')} == {'GCGCGCGGC': 4}
    assert {m.group(): m.start() for m in dna_sequence.re_search(r'[^GC]{2,}')} == {'ATA': 1, 'TAAA': 13}
    assert {m.group(): m.start() for m in protein_sequence.re_search(r'THERE')} == {'THERE': 2}
    assert {m.group(): m.start() for m in protein_sequence.re_search(r'(RE[FD]){2,3}')} == {'REFRED': 5}
    assert {m.group(): m.start() for m in protein_sequence.re_search(r'[^GC]{2,}')} == {'HITHEREFREDAND': 0, 'RE': 15}


def test_search(protein_sequence, dna_sequence):
    assert {m.group(): m.start() for m in dna_sequence.search('TAG')} == {'TAG': 2}
    assert {m.group(): m.start() for m in dna_sequence.search('CG(2)')} == {'CGG': 9}
    assert {m.group(): m.start() for m in dna_sequence.search('(GC)(2,3)')} == {'GCGCGC': 4}
    assert {m.group(): m.start() for m in dna_sequence.search('[GC](2,10)')} == {'GCGCGCGGC': 4}
    assert {m.group(): m.start() for m in dna_sequence.search('{GC}(2,)')} == {'ATA': 1, 'TAAA': 13}
    assert {m.group(): m.start() for m in protein_sequence.search('THERE')} == {'THERE': 2}
    assert {m.group(): m.start() for m in protein_sequence.search('(RE[FD])(2,3)')} == {'REFRED': 5}
    assert {m.group(): m.start() for m in protein_sequence.search('{GC}(2,)')} == {'HITHEREFREDAND': 0, 'RE': 15}
    assert protein_sequence.search('RRRRRRR') == []
