"""
primers.py -- sequence handling and primer design checks (no Streamlit in here)

Everything the web app needs besides the melting-curve maths:
  * reading pasted sequences / FASTA files (including files with several records)
  * simple primer quality checks: GC%, GC clamp, runs, self-dimer, hairpin, 3' dimer
  * scanning a template for candidate primers
  * pairing forward and reverse primers by similar Tm

The dimer and hairpin checks are quick screening rules (longest stretch of
complementary bases), not free-energy calculations like IDT OligoAnalyzer.
"""

from bisect import bisect_left, bisect_right

from thermo import calc_tm

_COMPLEMENT = str.maketrans("ACGT", "TGCA")

# Screening thresholds (bases). Change here if your course guidelines differ.
MAX_RUN = 3          # a run of 4 identical bases is rejected
MAX_SELF_DIMER = 4   # a stretch of 5 or more complementary bases is flagged
MAX_HAIRPIN = 4      # a hairpin stem of 5 or more bases is flagged
MIN_LOOP = 3         # smallest loop a hairpin can close
DIMER_3P = 4         # how many 3'-end bases are tested for self-binding
MAX_TEMPLATE = 10000


def revcomp(seq):
    """Reverse complement: the other strand, read 5' to 3'."""
    return seq.translate(_COMPLEMENT)[::-1]


def gc_percent(seq):
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


# ---------------------------------------------------------------- reading input

def parse_sequences(text):
    """Turn pasted text or a FASTA file into (records, errors).

    records: list of (name, sequence). Sequences are upper-cased; spaces and
    digits (e.g. GenBank line numbers) are ignored.
    * If the text has FASTA headers ('>name'), each header starts a new record
      and the lines under it are joined.
    * Otherwise every non-empty line is one sequence (a list of primers).
    errors: list of readable messages for records that cannot be used.
    """
    records, errors = [], []
    lines = [ln.strip() for ln in text.replace("\r", "\n").split("\n") if ln.strip()]
    if not lines:
        return records, errors

    def add(name, raw):
        seq = "".join(ch for ch in raw.upper() if not ch.isspace() and not ch.isdigit())
        bad = sorted(set(seq) - set("ACGT"))
        if not seq:
            errors.append("%s: no sequence found." % name)
        elif bad:
            errors.append("%s: only A, C, G, T are allowed (found %s)." % (name, ", ".join(bad)))
        elif len(seq) < 4:
            errors.append("%s: sequence is too short (%d bases)." % (name, len(seq)))
        else:
            records.append((name, seq))

    if any(ln.startswith(">") for ln in lines):
        name, chunks = None, []
        for ln in lines:
            if ln.startswith(">"):
                if name is not None:
                    add(name, "".join(chunks))
                elif chunks:
                    add("sequence", "".join(chunks))   # text that came before the first header
                name, chunks = ln[1:].strip() or "record %d" % (len(records) + len(errors) + 1), []
            else:
                chunks.append(ln)
        if name is not None:
            add(name, "".join(chunks))
        else:
            add("sequence", "".join(chunks))
    else:
        for i, ln in enumerate(lines, 1):
            add("Seq %d" % i, ln)
    return records, errors


# ------------------------------------------------------------------ the checks

def longest_complementary(a, b):
    """Longest run of bases in `a` that can pair, antiparallel, with a run in `b`.

    Equivalent to the longest common substring of `a` and the reverse complement of `b`.
    """
    rb = revcomp(b)
    best = 0
    prev = [0] * (len(rb) + 1)
    for i in range(1, len(a) + 1):
        cur = [0] * (len(rb) + 1)
        for j in range(1, len(rb) + 1):
            if a[i - 1] == rb[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best:
                    best = cur[j]
        prev = cur
    return best


def hairpin_stem(seq, min_loop=MIN_LOOP):
    """Longest stem where an earlier stretch pairs with a later stretch of the same primer,
    leaving a loop of at least `min_loop` bases."""
    n, best = len(seq), 0
    for k in range(2, n // 2 + 1):
        found = False
        for i in range(n - 2 * k - min_loop + 1):
            target = revcomp(seq[i:i + k])
            if seq.find(target, i + k + min_loop) != -1:
                found = True
                break
        if found:
            best = k
        else:
            break  # if no stem of length k exists, none longer can
    return best


def three_prime_self_match(seq, k=DIMER_3P):
    """True if the last k bases at the 3' end can pair with another copy of the primer
    (the dangerous kind of primer-dimer, because the polymerase extends from the 3' end)."""
    return revcomp(seq[-k:]) in seq


def has_run(seq, n=MAX_RUN + 1):
    return any(base * n in seq for base in "ACGT")


def gc_clamp_ok(seq):
    """3' end ends in G or C, and the last 5 bases hold no more than 3 G/C."""
    tail = seq[-5:]
    return seq[-1] in "GC" and (tail.count("G") + tail.count("C")) <= 3


def check_primer(seq):
    """Run every quality check on one primer (5'->3'). Returns a dict of results."""
    out = {
        "gc": gc_percent(seq),
        "clamp": gc_clamp_ok(seq),
        "run": has_run(seq),
        "self_dimer": longest_complementary(seq, seq),
        "hairpin": hairpin_stem(seq),
        "dimer_3p": three_prime_self_match(seq),
    }
    return out


def problems(chk, gc_min=40.0, gc_max=60.0):
    """List the reasons a primer fails (empty list = passes everything)."""
    p = []
    if not gc_min <= chk["gc"] <= gc_max:
        p.append("GC%% outside %g-%g" % (gc_min, gc_max))
    if chk["run"]:
        p.append("run of 4 identical bases")
    if not chk["clamp"]:
        p.append("weak GC clamp")
    if chk["self_dimer"] > MAX_SELF_DIMER:
        p.append("self-dimer (%d bp)" % chk["self_dimer"])
    if chk["hairpin"] > MAX_HAIRPIN:
        p.append("hairpin (%d bp stem)" % chk["hairpin"])
    if chk["dimer_3p"]:
        p.append("3' end pairs with itself")
    return p


# ---------------------------------------------------------- scanning a template

def scan_template(template, length, gc_min=40.0, gc_max=60.0, tm_min=0.0, tm_max=100.0,
                  na_mM=50.0, mg_mM=0.0, oligo_nM=250.0, dntp_mM=0.0):
    """Slide a window along the template; return one row per candidate that passes every check.

    Forward primers read the template as written.  Reverse primers are the reverse
    complement of a window, so they bind the other strand.
    """
    rows = []
    for strand in ("fwd", "rev"):
        for i in range(len(template) - length + 1):
            window = template[i:i + length]
            seq = window if strand == "fwd" else revcomp(window)
            chk = check_primer(seq)
            if problems(chk, gc_min, gc_max):
                continue
            tm = calc_tm(seq, na_mM, oligo_nM, mg_mM, dntp_mM)["tm"]
            if not tm_min <= tm <= tm_max:
                continue
            rows.append({"strand": strand, "start": i, "end": i + length, "seq": seq,
                         "tm": tm, "chk": chk})
    return rows


def find_pairs(candidates, product_min=100, product_max=500, max_dtm=2.0, limit=25):
    """Pair each forward primer with a reverse primer downstream of it.

    A pair is kept if the product length is in range, the primers' Tm differ by at most
    `max_dtm`, and the two primers do not bind each other (cross-dimer).
    Best pairs (smallest Tm difference, then shorter product) come first.
    """
    fwd = [c for c in candidates if c["strand"] == "fwd"]
    rev = sorted((c for c in candidates if c["strand"] == "rev"), key=lambda c: c["start"])
    rev_starts = [c["start"] for c in rev]
    pairs = []
    for f in fwd:
        lo = bisect_left(rev_starts, f["end"])                      # reverse primer must start after f ends
        hi = bisect_right(rev_starts, f["start"] + product_max)     # product cannot be longer than product_max
        for r in rev[lo:hi]:
            product = r["end"] - f["start"]
            if not product_min <= product <= product_max:
                continue
            dtm = abs(f["tm"] - r["tm"])
            if dtm > max_dtm:
                continue
            cross = longest_complementary(f["seq"], r["seq"])
            if cross > MAX_SELF_DIMER:
                continue
            pairs.append({"fwd": f, "rev": r, "product": product, "dtm": dtm, "cross": cross})
    pairs.sort(key=lambda p: (p["dtm"], p["product"]))
    return pairs[:limit]
