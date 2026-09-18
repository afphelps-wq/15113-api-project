# Prompt Log

## AI tools used

**Claude Code** (Anthropic) in VS Code, using Claude Sonnet 5 for the planning and setup phase and
Claude Opus 5 for the implementation. Claude Code can run commands and edit files directly, so much
of the work involved it running the code it had written, showing me the output, and fixing what
broke. The prompts below are the ones that actually shaped the project, in order.

## Key prompts

**1. Framing the problem.**
> "I need to integrate an API key, I want something that I could to create some type of tool for
> researchers that is related to spatial transcriptomics that takes visium data and does some type
> of analysis. Would that be do complicated or is there an API key I could use for that?"

This started as a spatial transcriptomics idea. The useful answer was that Visium analysis is not an
API at all — it is a local library — so the project needed a separate API for the parts that are
genuinely external. That reframing led to the eventual design.

**2. Changing direction to something achievable.**
> "Would it be more simple to create a bulk rna sequencing pipeline and then wrap it would a user
> friendly interface so anyone can input the details about the samples and then the analysis is done
> and results are matched and NCBI utilties can be used for top genes in a region and open AI api
> for prompting"

Bulk RNA-seq replaced spatial. Two corrections came out of this exchange: bulk data has no spatial
"regions", so the unit of analysis became a two-group comparison; and the app should take a count
matrix rather than raw FASTQ files, since read alignment needs far more compute than a web form.

**3. Verifying the dataset instead of trusting the abstract.**
> "Also search and verify the data accession to find the metadata and make sure it is available then
> copy it into this test_dat folder" (with the paper abstract pasted in)

The abstract had no accession number in it. Claude searched NCBI GEO, identified **GSE205154**, and
confirmed that all 289 sample IDs in the metadata matched the count matrix columns exactly. It also
found two things I would have hit later: the counts are non-integer estimates that must be rounded
before DESeq2, and the paper's liver-versus-lung comparison is not in the public metadata, so the
app's demo had to use primary versus metastasis instead.

**4. Planning before building.**
> "before this I want to plan. ask me as many questions that you can think of to clarify the this
> project as best possible and save my answers to a txt file in this folder that is part of the
> context"

This was the highest-leverage prompt in the project. Six rounds of questions settled the audience,
the input formats, the must-have versus stretch features, the privacy rule, and the deadline. The
answers live in [`project_context.txt`](project_context.txt) and the resulting plan in
[`build_plan.md`](build_plan.md). Writing them to a file meant later work stayed consistent with
decisions made early.

**5. Grounding the AI summary.**
Rather than one prompt, this was a decision recorded during planning: the model may use *only* the
statistics and abstracts supplied to it, must cite PMIDs, and must say "insufficient evidence" where
the abstracts do not support a claim. The system prompt in [`rnaseq/llm.py`](rnaseq/llm.py) also
tells the model that automated PubMed searches return irrelevant papers for short gene symbols, and
that bulk tissue comparisons across organs largely reflect tissue composition.

That last instruction proved its worth. In the generated summary the model wrote, unprompted:

> "The pattern is consistent with a tissue-composition shift between the two groups... rather than a
> simple tumour-cell intrinsic transcriptional program."

It also declined to interpret three genes whose retrieved abstracts were irrelevant. A check of the
output confirmed all 10 cited PMIDs were real papers from the supplied set, with no fabrications.

**6. Building, then verifying with real calls.**
> "keep going" / "do documentation first"

The implementation ran in stages, each one tested against real data rather than assumed to work.
Four bugs surfaced this way that a reading of the code would not have caught:

| Bug | How it was found |
|---|---|
| Collapsing duplicate gene labels silently reordered every gene in the file | A unit test compared the output order to the input order |
| The volcano plot legend covered the data points, and gene labels overlapped each other | Claude rendered the PNG and looked at it |
| The upload warned it was ignoring `gene_symbol`, the column it was actually using as the gene label | Running the real upload and reading the warnings |
| The AI summary came back empty every time | The API reported `status=incomplete, reason=max_output_tokens`: gpt-5-mini spends output tokens on internal reasoning before writing, so the 1400-token cap was consumed before any text was produced. Raised to 4000 with reasoning effort set to low |

## What I would tell someone doing this next

The planning prompt was worth more than any implementation prompt. Answering questions about
audience, scope and priorities before any code existed meant the build had a clear must-have list
when the deadline got tight.

The other lesson is to make the AI verify its own work. Asking for the analysis was easy; the value
came from checking that the cited PMIDs were real, that the biology matched what is known about
these tissues, and that the numbers in the README came from an actual run rather than an estimate.
