import pandas as pd
from pathlib import Path


try:
    # Snakemake 8.x.x
    from snakemake_interface_common.exceptions import WorkflowError
except ImportError:
    # Snakmake 7.x.x
    from snakemake.exceptions import WorkflowError


def parse_sample_sheet(config: dict) -> pd.DataFrame:
    samples = (
        pd.read_table(config["samples"], sep=",", dtype=str)
        .replace(" ", "_", regex=True)
        .infer_objects(copy=False)  # needed to maintain same behavior in future pandas versions
    )
    config_genomes = get_config_genomes(config, samples)
    refGenome = "refGenome" in samples.columns and samples["refGenome"].notna().any()
    refPath = "refPath" in samples.columns and samples["refPath"].notna().any()
    if not any([config_genomes, refGenome, refPath]):
        raise WorkflowError("No 'refGenome' or 'refPath' found in config or sample sheet.")
    if config_genomes is not None:
        config_refGenome, config_refPath = config_genomes
        samples["refGenome"] = config_refGenome
        samples["refPath"] = config_refPath
    if "refPath" in samples.columns and samples["refPath"].notna().any():
        check_ref_paths(samples)
    return samples


def get_config_genomes(config: dict, samples: pd.DataFrame):
    refGenome = config.get("refGenome", False)
    refPath = config.get("refPath", False)

    if refGenome and refPath:
        if "refGenome" in samples.columns and samples["refGenome"].notna().any():
            raise WorkflowError(
                "'refGenome' is set in sample sheet AND in config. These are mutually exclusive."
            )
        return refGenome, refPath
    elif refGenome and not refPath:
        raise WorkflowError(
            "'refGenome' is set in config, but 'refPath' is not. Both are required to use these settings."
        )
    elif refPath and not refGenome:
        raise WorkflowError(
            "'refPath' is set in config, but 'refGenome' is not. Both are required to use these settings."
        )
    return None


def check_ref_paths(samples: pd.DataFrame) -> None:
    """
    Checks reference paths to make sure they exist, otherwise we might try to download them based on refGenome.
    Also make sure only one refPath per refGenome.
    """
    for refname in samples["refGenome"].dropna().tolist():
        refs = samples[samples["refGenome"] == refname]["refPath"].dropna().unique().tolist()
        if len(refs) > 1:
            raise WorkflowError(
                f"refGenome '{refname}' has more than one unique 'refPath' specified: {refs}"
            )
        for ref in refs:
            if not Path(ref).exists:
                raise WorkflowError(
                    f"refPath: '{ref}' was specified in sample sheet, but could not be found."
                )
