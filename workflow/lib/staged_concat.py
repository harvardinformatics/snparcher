"""Plan a staged concatenation of many files.

Round 1 gathers the base files in chunks of `batch_size`; each later round
gathers the previous round's chunks, until a single file is left. Staged
outputs are named by `stage_path(round, chunk)`.
"""


def stage_chunk_counts(num_files, batch_size, max_rounds):
    """Number of chunks in each round; empty when there is a single file."""
    if num_files < 1:
        raise ValueError("Staged concat requires at least one input file")

    counts = []
    current = num_files
    while current > 1:
        if len(counts) >= max_rounds:
            raise ValueError(
                f"Staged concat exceeded variant_calling.gatk.concat_max_rounds={max_rounds}. "
                "Increase concat_batch_size or concat_max_rounds."
            )
        current = (current + batch_size - 1) // batch_size
        counts.append(current)
    return counts


def stage_inputs(base_files, round_idx, chunk_idx, batch_size, max_rounds, stage_path):
    """Input files for one chunk of one round."""
    counts = stage_chunk_counts(len(base_files), batch_size, max_rounds)
    if not counts:
        raise ValueError("Requested staged concat inputs with only one base file")

    if round_idx < 1 or round_idx > len(counts):
        raise ValueError(
            f"Invalid staged concat round {round_idx}. Valid rounds are 1..{len(counts)}"
        )

    n_chunks = counts[round_idx - 1]
    if chunk_idx < 0 or chunk_idx >= n_chunks:
        raise ValueError(
            f"Invalid staged concat chunk {chunk_idx} for round {round_idx}. "
            f"Valid chunks are 0..{n_chunks - 1}"
        )

    if round_idx == 1:
        round_inputs = list(base_files)
    else:
        round_inputs = [stage_path(round_idx - 1, i) for i in range(counts[round_idx - 2])]

    selected = round_inputs[chunk_idx * batch_size : (chunk_idx + 1) * batch_size]
    if not selected:
        raise ValueError(
            f"No files selected for staged concat round={round_idx}, chunk={chunk_idx}"
        )
    return selected


def final_stage_file(base_files, batch_size, max_rounds, stage_path):
    """The gathered file: the last round's only chunk, or the single base file."""
    if not base_files:
        raise ValueError("No files available for concat")

    counts = stage_chunk_counts(len(base_files), batch_size, max_rounds)
    if not counts:
        return base_files[0]
    return stage_path(len(counts), 0)
