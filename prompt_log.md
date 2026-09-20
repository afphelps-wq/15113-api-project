# Prompt Log

## AI tools used

**Claude Code** (Anthropic) in VS Code — Claude Sonnet 5 for the planning and setup phase, then
Claude Opus 5 for the implementation. Claude Code runs commands and edits files directly, so most of
the work was it writing code, running it, showing me the output, and fixing what broke.

Below are the prompts that actually changed the project, in order, with what came out of each. It is
not the whole conversation.

## Finding the project

**1. The original idea, which did not survive contact with reality.**
> "I need to integrate an API key, I want something that I could to create some type of tool for
> researchers that is related to spatial transcriptomics that takes visium data and does some type
> of analysis. Would that be do complicated or is there an API key I could use for that?"

The useful answer was that Visium analysis is not an API at all — it is a local Python library. That
reframing is why the project ended up combining local analysis with genuinely external APIs.

**2. Changing direction to something achievable.**
> "Would it be more simple to create a bulk rna sequencing pipeline and then wrap it would a user
> friendly interface so anyone can input the details about the samples and then the analysis is done
> and results are matched and NCBI utilties can be used for top genes in a region and open AI api
> for prompting"

Two corrections came out of this: bulk data has no spatial "regions", so the unit of analysis became
a two-group comparison; and the app takes a count matrix rather than raw FASTQ files, because read
alignment needs far more compute than a web form.

**3. Verifying the dataset instead of trusting the abstract.**
> "Also search and verify the data accession to find the metadata and make sure it is available then
> copy it into this test_dat folder" (with the paper abstract pasted in)

The abstract had no accession number. Claude searched NCBI GEO, identified **GSE205154**, and
confirmed all 289 sample IDs matched the count matrix columns. It also found two things I would have
hit later: the counts are non-integer estimates that must be rounded before DESeq2, and the paper's
liver-versus-lung comparison is not in the public metadata, so the demo had to use primary versus
metastasis instead.

**4. Planning before building — the highest-leverage prompt in the project.**
> "before this I want to plan. ask me as many questions that you can think of to clarify the this
> project as best possible and save my answers to a txt file in this folder that is part of the
> context"

Six rounds of questions settled the audience, the input formats, the must-have versus stretch
features, the privacy rule and the deadline. The answers are in
[`project_context.txt`](project_context.txt) and the plan in [`build_plan.md`](build_plan.md).
Writing them to a file meant later work stayed consistent with decisions made early, and made the
scope cuts obvious when time got short.

## Building it

**5. Grounding the AI summary.** Rather than one prompt, this was a decision recorded during
planning: the model may use *only* the statistics and abstracts supplied to it, must cite PMIDs, and
must say "insufficient evidence" where the abstracts do not support a claim. The system prompt in
[`rnaseq/llm.py`](rnaseq/llm.py) also warns that automated PubMed searches return irrelevant papers
for short gene symbols, and that bulk tissue comparisons across organs largely reflect tissue
composition.

That last instruction proved its worth. The model wrote, unprompted:

> "The pattern is consistent with a tissue-composition shift between the two groups... rather than a
> simple tumour-cell intrinsic transcriptional program."

It also declined to interpret three genes whose retrieved abstracts were irrelevant. A check of the
output confirmed all 10 cited PMIDs were real papers from the supplied set, with no fabrications.

**6. "keep going" / "do documentation first".** The build ran in stages — analysis core, Flask app,
PubMed, OpenAI, README — each tested against real data rather than assumed to work.

**7. "add pathway enrichment."** Added g:Profiler (GO, KEGG, Reactome), which needs no API key. Two
decisions mattered more than the code: send the *tested* genes as the statistical background rather
than the whole genome, and drop terms larger than 500 genes, because with thousands of significant
genes the top hit was otherwise "multicellular organismal process" (533 of 6,617 genes) — unarguable
and useless.

**8. "give options for pathway enrichment to be displayed using dot plots, bar plots, or enrichment
network maps generated via R packages like clusterProfiler."** Claude checked and found R installed
but clusterProfiler absent, then asked whether to use R for real or match its look in Python. I chose
Python: using R would have meant graders installing R plus three Bioconductor packages, and would
have replaced g:Profiler's statistics with clusterProfiler's own. The figures follow enrichplot's
conventions but are drawn with matplotlib.

**9. "Add these options for deg analysis: volcano plot, heatmap, PCA plot, MA plot."** Volcano and
heatmap already existed, so this added MA and PCA. The PCA produced the most interesting finding in
the project: coloured by tissue, the first component separates **liver** metastases from everything
else, and the other metastases sit with the primaries — so the headline comparison is mostly a liver
signal. Two tests now assert it.

## Fixing and polishing

**10. Design prompts, with screenshots as the input.** Two rounds: a soft-card dashboard, then a
frosted-glass one, each supplied as a reference image. Prompts like *"fix this white boxes so they
are lined up properly"* and *"the pathway enrichment graphs are overlapping as the background of the
page"* were pasted screenshots of what was wrong. Before the second redesign, Claude wrote a test
that reads `app.js` and checks the page still provides every id, class and data attribute the script
depends on, so a redesign cannot silently break a button.

**11. "make a dark mode option."** Dark mode had to reach the figures, which are images drawn on the
server; otherwise they stay white rectangles. Each theme has its own colour steps, checked with a
colour-blind and contrast validator against its own background.

**12. "act as a reviewer and harshly skritinize every aspect of the page."** The most valuable prompt
after the planning one. Asking for a hostile review — with instructions to test rather than guess —
surfaced real defects that all the passing tests had missed, including a way for a crafted metadata
file to run code in the page, stale pathway results being fed to the AI summary after a second
comparison, and a macOS conflict on port 5000 that makes `localhost:5000` show a blank page.

## Bugs that only showed up by running the code

| Bug | How it was found |
|---|---|
| Collapsing duplicate gene labels silently reordered every gene | A unit test compared output order to input order |
| Volcano legend covered the data; gene labels overlapped | Claude rendered the PNG and looked at it |
| The upload warned it was ignoring `gene_symbol`, the column it was using | Reading the warnings from a real upload |
| Every AI summary came back empty | The API reported `max_output_tokens`: gpt-5-mini spends output tokens on internal reasoning first, so the 1,400-token cap was consumed before any text |
| The enrichment network rendered blank | A zero self-weight times an infinite self-distance gives NaN, which spread to every node position |
| Hidden figures were displayed anyway, stacked | A CSS `display` rule silently overrode the `hidden` attribute |
| The browser was requesting a `/favicon.ico` that did not exist | A scripted browser run flagged the console error |

## What I would tell someone doing this next

The planning prompt was worth more than any implementation prompt. Answering questions about
audience, scope and priorities before any code existed meant the build had a clear must-have list
when the deadline got tight.

The second lesson is to make the AI check its own work, and to ask for criticism explicitly. "Write
this feature" produces something plausible; "test it, look at the output, and tell me what's wrong
with it" is what caught the empty summaries, the fabricated-citation risk, the blank network figure
and the security hole.
