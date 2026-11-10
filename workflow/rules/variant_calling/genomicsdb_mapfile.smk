localrules: create_db_mapfile


# The GenomicsDB sample-name map, shared by every GenomicsDB path: GATK with
# and without intervals, and Parabricks.
def create_db_mapfile_input(wc):
    inputs = {"gvcfs": get_joint_gvcf_paths()}
    if GATK_LONG_CONTIG_MODE:
        inputs.update(
            {
                "gvcf_indexes": get_joint_gvcf_indexes(),
                "archive_gvcfs": get_generated_gvcf_archives(),
                "archive_indexes": get_generated_gvcf_archive_indexes(),
            }
        )
    return inputs


rule create_db_mapfile:
    input:
        unpack(create_db_mapfile_input),
    output:
        mapfile="results/genomics_db/mapfile.txt",
    run:
        write_joint_gvcf_mapfile(output.mapfile)
