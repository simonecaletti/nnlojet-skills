# The full-subtraction workflow as a graph

Implementing the subtraction for a whole process is usually described as
a line:

> understand what is missing → write a channel → spike test → fix it →
> next channel

That description is not wrong, but it is lossy in three places, and each
loss has a measured cost:

1. **"spike test" is not one node.** The BUILD is process-wide (one
   `maple makefort<RR|RV> -Diprocess=N` regenerates every `auto*.f` of
   the iprocess), while the TEST is channel-local (`./check<N>to<M> <CH>`
   selects the (ME, subtraction) pair at runtime, and the mode-range
   arguments narrow it to individual limits). Conflating them makes the
   build look like a per-channel cost when it is a shared, batchable one.
2. **"fix it" is not one edge.** It is at least seven, landing on four
   different nodes, with cost ratios around 100:1 between the cheapest
   and the dearest. A procedure with only one repair edge routes every
   failure to "edit the block in front of me" — which is exactly how the
   `C1g0ZepemS` rebuild spent ~20 build cycles repairing blocks
   downstream of a wrong `S,a`.
3. **"next channel" is not a plain loop.** Facts learned in channel N
   propagate backwards to channels already finished, and at least one
   channel must stay untouched as a control.

This file writes the procedure out as a graph so the routing is explicit.
It contains no process-specific answers and is safe to load during a
clean-room rebuild.

---

## Node table

Each node has a precondition, an artifact, and an exit test. A node is
done when its artifact exists and its exit test passes — not when it
feels finished. `∥` marks nodes that may run concurrently across
channels; `SERIAL` marks the critical section.

| # | node | precondition | artifact | exit test |
|---|---|---|---|---|
| ① | **Process intake** | a process and layer are named | iprocess number, `FLAVlist<name>.map` entry per channel, channel list, test program + its argument convention | every channel has a FLAVlist entry, or a derivation of one by the reduced-Born rule |
| ② | **Conventions contract** | ① | a written note: absorption convention, Full-vs-split policy, `colflag` ME ordering, argument-order convention | the contract names a choice for each, with a reason |
| ③ | **Derive skeleton** `∥` | ② | `modes.json` (`genuine_modes.py`), block prediction (`predict_blocks.py`), `master.map` with `# block:` markers and axes | `predict_blocks` runs clean; every genuine mode has a predicted block |
| ④ | **Measure antennae** `∥` | ③ named an antenna not in the datasheet | new datasheet entries: pole graph, `requires_split`, split identity, cluster semantics | `antenna_datasheet.py show <name>` returns the entry with provenance |
| ⑤ | **Compose term** `∥` | ③ ④ | `<TERM>.map` + `<TERM>.spec.json` | file parses; `aN` gap-free (guaranteed by `compose_blocks.py`) |
| ⑥ | **Static gate** | ⑤ | `audit_blocks` + `pole_ledger` + `pairing_balance` output | `audit_blocks` and `pole_ledger` both exit zero |
| ⑦ | **Build** `SERIAL` | ⑥ green for every term being promoted | rebuilt check binary | `regen_rebuild.sh` exits zero |
| ⑧ | **Preflight + run** | ⑦ | `preflight.py` clean, run log per channel | preflight exits zero; every genuine mode produced output |
| ⑨ | **Verdict** | ⑧ | the limit table, in chat | every genuine mode's MEDIAN sits at 1 across the x scan |

Nodes ①② are per PROCESS and run once. ③④⑤ are per CHANNEL and are the
fan-out. ⑥ is per channel in principle. ⑦ is shared. ⑧⑨ are per channel.

### Why ⑦ is the critical section

Three couplings, none of them physics:

- `maple makefort<RR|RV>` regenerates **every** `auto*.f` of the
  iprocess, not just the term you edited;
- `regen_rebuild.sh` runs the pole ledger over **every** `.map` in
  `maple/process/<DIR>/` and aborts the whole build if any one fails —
  so a sibling's work-in-progress blocks your build, with an error
  naming the sibling's file;
- there is **one binary** in `test/process/<PROC>/`, and the script
  deletes it before rebuilding.

Concurrent writers therefore deadlock each other and, worse, can run a
binary built from someone else's `.map`. Stage terms outside the live
`maple/process/<DIR>/` and let one owner promote them in batches.

### Why the shared binary is an asset, not just a constraint

Because the test is channel-local and the binary serves every channel, a
**validated sibling channel is a control that runs on the very same
binary** as the channel under test. Any artefact common to both is
provably harness, not term. That is a stronger control than a separate
build could give — which is the reason for the sequencing rule below.

---

## The repair edges

The single arrow labelled "fix it" resolves into these. Diagnose in this
order: the cheap edges are also the ones most often mistaken for the
expensive ones.

| # | signature | reading | goes to |
|---|---|---|---|
| E1 | `preflight.py` exits non-zero | stale build chain, or a collinear cluster only partly rotated | **⑦** — rebuild, no `.map` edit |
| E2 | outliers shrink monotonically as `n_azim` goes 4 → 8 → 16, median tightens onto 1 | the azimuthal average, not the term | **pass** — no edit at all |
| E3 | `audit_blocks` / `pole_ledger` error | a whole block missing or spurious, a split half used alone, a stale cluster, an orphaned counterterm | **⑤** — recompose |
| E4 | `wt_attribute.py` shows one line at an anomalous exponent or magnitude | wrong argument order or wrong coefficient on a named line | **⑤** |
| E5 | `fit_lines.py` returns CLOSABLE | the printed rationals ARE the term | **⑤** |
| E6 | `fit_lines.py` returns NOT CLOSABLE | the deficit is outside the span of the lines you wrote — stop searching, derive | **④** (measure), then ③ |
| E7 | half the modes exact, the mirror half broken; a coefficient change cannot move both | the `S,a` decomposition is wrong, and every later block is being asked to pair with an `S,a` no X40 matches | **③** — re-open the radiator axis, sweep |

Two properties of this table are worth stating out loud:

- **E2 is a pass, not a repair.** Reporting it as a failure is a false
  negative that costs real debugging time; it has been observed on the
  soft+collinear modes of a fully correct term.
- **E7 is the edge the linear procedure cannot express**, and it is the
  expensive one. Its signature is a PATTERN across modes, not a single
  mode's number, so it is only visible if you are looking at the whole
  limit table at once. That is why ⑨'s artifact is the table.

### Edge selection is itself a rung on the cost ladder

E1 and E2 cost no `.map` edit. E3 costs seconds. E4/E5/E6 cost one run.
E7 costs a sweep. Walk them in order; the temptation is to jump to E7's
conclusion ("the structure is wrong") from a symptom that E1 explains.

---

## Cross-channel edges

These are absent from the linear description and are where a
multi-channel effort actually goes wrong.

- **X1 — a datasheet fact learned in channel N invalidates channels
  1…N−1.** A split identity or a corrected pole graph changes what was
  legal in already-composed terms. Re-run ⑥ on every composed term; it
  is a static check and costs seconds.
- **X2 — a change to the conventions contract ② invalidates every
  composed term.** This is why ② is a node and not an afterthought:
  convention drift across channels of one process is very hard to detect
  afterwards, because each term looks locally correct.
- **X3 — a validated channel becomes infrastructure.** It is the
  `--baseline` for `pairing_balance.py` and the control channel for ⑧.
  So: **never have every channel in flight at once.** Keep at least one
  validated or untouched sibling, or you lose the ability to separate
  harness artefacts from your term.
- **X4 — sibling terms are one unit, not two channels.** An
  identical-flavour term can hold the X40 lines while its
  single-unresolved limits live in the non-identical sibling (the
  `"families"` declaration in the spec). Split those across independent
  workers and you get one term correctly declaring itself partial and
  nobody supplying the remainder. The same applies to the RR and RV
  layers of one channel, linked by the arrow structure
  `run-layer-check` verifies.

---

## Channel ordering

The loop edge "next channel" has a preferred order:

1. the structurally **simplest complete** channel first — it validates
   the toolchain and fixes the conventions contract against something
   real;
2. then channels **sharing the most antenna content** with it, so ④ is
   already paid;
3. sibling pairs (X4) together, by one worker;
4. keep one channel back as the control (X3) until the end.

---

## Delegating this graph

The node table is a delegation contract. If you hand work to a
sub-worker, the boundaries fall out of the graph rather than being
imposed on it:

- **safe to fan out**: ③④⑤ — static, no build, each producing a named
  artifact. ④ in particular parallelises over antennae rather than
  channels, and its artifact (the datasheet) is shared, so it must be
  written back through one owner to avoid four workers measuring the
  same antenna four times and disagreeing.
- **must not fan out**: ⑥⑦⑧ as a unit, for the reasons above.
- **worth fanning out for the wrong-looking reason**: an independent
  re-derivation at ③ that has NOT seen the incumbent `.map`. The value
  is context isolation, not throughput — a worker with no stake in the
  current term is the one that can see E7.

**Do not fan out by BLOCK.** S,a / S,b1 / S,b2 / S,c are a dependency
chain, not siblings: S,b2 is indexed by the (S,a line, X40 line) pairs
that share a pole, an S,b2 line's absorption is fixed by the X40 half it
cancels, and S,c is defined as a difference against the iterated product.
The minimal buildable and testable unit is S,a + S,b1 + S,b2 for one
flavour sector, so a per-block worker has no local success signal and
will report "my block looks right" for a term that fails. S,d is the one
genuinely separable block. A per-CHANNEL worker, by contrast, does have a
local signal — its own ME, its own genuine modes — which is the whole
difference between the two decompositions.

---

## Instrumenting it

Log which edge fires on each iteration. If E3/E4/E5 fire fifteen times
and E7 never does, you are watching a foundation error happen in
telemetry instead of discovering it twenty builds later. The graph's main
practical benefit over the linear procedure is that this question becomes
askable at all.

One dependency: every edge is selected by node ⑨'s verdict, so the
verdict has to be trustworthy. The median-plateau rule (and the
demotion of outlier counts to a tiebreaker) is what makes it so — see
run-spike-test, THE SCORING RULE in `scan_blocks.py`.
