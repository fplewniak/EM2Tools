def test_is_protein(protein_rec, dna_rec):
    assert protein_rec.seq.is_protein() == True
    assert dna_rec.seq.is_protein() == False


def test_re_search(protein_rec, dna_rec):
    assert {m.group(): m.start() for m in dna_rec.seq.re_search(r'GAG')} == {'GAG': 2}
    assert {m.group(): m.start() for m in dna_rec.seq.re_search(r'CG{2}')} == {'CGG': 6}
    assert {m.group(): m.start() for m in dna_rec.seq.re_search(r'(GC){2,3}')} == {'GCGCGC': 20}
    assert {m.group(): m.start() for m in dna_rec.seq.re_search(r'[GC]{4,10}')} == {'GCGCGCCGC': 20}
    assert {m.group(): m.start() for m in dna_rec.seq.re_search(r'[^GC]{2,}')} == {'AT': 18, 'TAAA': 9}
    assert [(m.group(), m.start()) for m in dna_rec.seq.re_search(r'[^GC]{2,}')] == [('AT', 0), ('TAAA', 9),
                                                                                     ('AT', 14), ('AT', 18)]
    assert {m.group(): m.start() for m in protein_rec.seq.re_search(r'THERE')} == {'THERE': 2}
    assert {m.group(): m.start() for m in protein_rec.seq.re_search(r'(RE[FD]){2,3}')} == {'REFRED': 5}
    assert {m.group(): m.start() for m in protein_rec.seq.re_search(r'[^GC]{2,}')} == {'HITHEREFREDAND': 0, 'RE': 15}


def test_search(protein_rec, dna_rec):
    assert {m.group(): m.start() for m in dna_rec.seq.search('GAG')} == {'GAG': 2}
    assert {m.group(): m.start() for m in dna_rec.seq.search('CG(2)')} == {'CGG': 6}
    assert {m.group(): m.start() for m in dna_rec.seq.search('(GC)(2,3)')} == {'GCGCGC': 20}
    assert {m.group(): m.start() for m in dna_rec.seq.search('[GC](4,10)')} == {'GCGCGCCGC': 20}
    assert {m.group(): m.start() for m in dna_rec.seq.search('{GC}(2,)')} == {'AT': 18, 'TAAA': 9}
    assert [(m.group(), m.start()) for m in dna_rec.seq.search('{GC}(2,)')] == [('AT', 0), ('TAAA', 9),
                                                                                 ('AT', 14), ('AT', 18)]
    assert {m.group(): m.start() for m in protein_rec.seq.search('THERE')} == {'THERE': 2}
    assert {m.group(): m.start() for m in protein_rec.seq.search('(RE[FD])(2,3)')} == {'REFRED': 5}
    assert {m.group(): m.start() for m in protein_rec.seq.search('{GC}(2,)')} == {'HITHEREFREDAND': 0, 'RE': 15}
    assert protein_rec.seq.search('RRRRRRR') == []


def test_get_orfs(dna_rec1_overlap_rec2):
    assert [(orf.start, orf.end) for orf in dna_rec1_overlap_rec2.seq.get_orfs()] ==\
           [(0, 9), (2, 29), (12, 36), (32, 44), (4, 58), (47, 59), (39, 60)]
    assert [(orf.start, orf.end) for orf in dna_rec1_overlap_rec2.seq.get_orfs(start=None)] ==\
           [(0, 9), (14, 29), (18, 36), (35, 44), (22, 58), (53, 59)]