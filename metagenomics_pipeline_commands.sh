# Basic set of tools to run on a paired-end metagenomic Illumina dataset
# Things that need to be installed: bowtie2, fastqmcf, Metaphlan2, MIDAS, HUMAnN2, sickle, fastqc, parallel
# R packages used for binding a bunch of tab-delimted metaphlan results: dplyr, readr




export raw_dir="input_data/ihg-client.ucsf.edu/hsiaoe/170822_K00358_0112_AHL2YWBBXX_hsiaoe-Swpool/"
export adapter_file="adapters.fa"
export fastqmcf_dir="fastqmcf_out"
export sickle_dir="sickle_out"
export human_bt2_ref="~/work/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"
export host_removed_dir="host_removal_out"
export fastqc_on_adapter_removed_dir="fastqc_on_sickle_out"
export fastqc_on_adapter_removed_dir2="fastqc_on_fastqmcf_out"
export metaphlan2_dir="metaphlan2_out"
export midas_dir="MIDAS_out"
export MIDAS_DB="/pollard/home/slyalina/work/reference/midas_db_v1.2/"
export humann2_dir="humann2_out"


# Checking that files didn't get mangled in upload
export raw_dir
md5sum ${raw_dir}/*.fastq.gz > received_md5sums.txt
join <(sort received_md5sums.txt) <(sort ${raw_dir}/md5sums.txt)
# Assuming all files have same naming scheme, extract their prefixes to be used for parallelization later
find input_data/ihg-client.ucsf.edu/hsiaoe/170822_K00358_0112_AHL2YWBBXX_hsiaoe-Swpool/ -name "*fastq.gz" -exec basename {} \; | grep -o "^.*R" | sort -u > prefixes.txt

# Run fastq-mcf (adapter trimming, QC)
mkdir -p ${fastqmcf_dir}
run_fastqmcf_on_sample(){
	prefix=$1
	~/work/tools/ea-utils-read-only/clipper/fastq-mcf ${adapter_file} \
	-o ${fastqmcf_dir}/${prefix}1.fq.gz \
	-o ${fastqmcf_dir}/${prefix}2.fq.gz \
	<(gunzip -c ${raw_dir}/${prefix}1_001.fastq.gz) \
	<(gunzip -c ${raw_dir}/${prefix}2_001.fastq.gz) > ${fastqmcf_dir}/${prefix}_log.txt

}
export -f run_fastqmcf_on_sample
parallel -j 6 -a prefixes.txt run_fastqmcf_on_sample {}

# Run sickle (alternative method for adapter trimming, QC)
mkdir -p ${sickle_dir}
run_sickle_on_sample(){
	prefix=$1
	sickle pe \
	-t sanger \
	-s ${sickle_dir}/${prefix}_singles.fq.gz \
	-g \
	-o ${sickle_dir}/${prefix}1.fq.gz -p ${sickle_dir}/${prefix}2.fq.gz \
	-f ${raw_dir}/${prefix}1_001.fastq.gz -r ${raw_dir}/${prefix}2_001.fastq.gz > ${sickle_dir}/${prefix}_log.txt
}
export -f run_sickle_on_sample
parallel -j 6 -a prefixes.txt run_sickle_on_sample {}

# Files provided by IHG core already have fastqc report so don't need to run that
# Run FastQC on sickle output
mkdir -p ${fastqc_on_adapter_removed_dir}
run_fastqc_on_sample(){
	prefix=$1
	fastqc --noextract -t 2 \
	-o ${fastqc_on_adapter_removed_dir} \
	${sickle_dir}/${prefix}1.fq.gz ${sickle_dir}/${prefix}2.fq.gz
}
export -f run_fastqc_on_sample
parallel -j 6 -a prefixes.txt run_fastqc_on_sample {}

# Run FastQC on fastq-mcf output
mkdir -p ${fastqc_on_adapter_removed_dir2}
run_fastqc_on_sample2(){
	prefix=$1
	fastqc --noextract -t 2 \
	-o ${fastqc_on_adapter_removed_dir2} \
	${fastqmcf_dir}/${prefix}1.fq.gz ${fastqmcf_dir}/${prefix}2.fq.gz
}
export -f run_fastqc_on_sample2
parallel -j 6 -a prefixes.txt run_fastqc_on_sample2 {}


# Run host removal (Map against host genome, pick up unmapped reads)
mkdir -p ${host_removed_dir}
run_host_removal_on_sample(){
	prefix=$1
	bowtie2 -x ${human_bt2_ref} \
	-1 ${sickle_dir}/${prefix}1.fq.gz \
	-2 ${sickle_dir}/${prefix}2.fq.gz \
	-p 2 --mm \
	--un-conc-gz ${host_removed_dir}/${prefix}%.fq.gz > /dev/null 2>${host_removed_dir}/${prefix}_log.txt
}
export -f run_host_removal_on_sample
parallel -j 6 -a prefixes.txt run_host_removal_on_sample {}

# Run Metaphlan2
mkdir -p ${metaphlan2_dir}
run_metaphlan2_on_sample(){
	prefix=$1
	~/work/tools/Metaphlan2/metaphlan2.py \
	${host_removed_dir}/${prefix}1.fq.gz,${host_removed_dir}/${prefix}2.fq.gz \
	--no_map -t rel_ab_w_read_stats \
	--sample_id ${prefix} \
	--nproc 4 --input_type fastq -o ${metaphlan2_dir}/${prefix}_metaphlan2_output.txt
}
export -f run_metaphlan2_on_sample
parallel -j 3 -a prefixes.txt run_metaphlan2_on_sample {}


# Bind Metaphlan2 outputs together
mkdir -p outputs
R_bind_script="
files = list.files('%',pattern='*_metaphlan2_output.txt',full.names = TRUE);
df = dplyr::bind_rows(setNames(lapply(files, function(x)
    readr::read_tsv(
    x,
    col_names = c(
    'clade_name',
    'relative_abundance',
    'coverage',
    'average_genome_length_in_the_clade',
    'estimated_number_of_reads_from_the_clade'
    ),
    comment = '#'
    )), files), .id = 'file');
df = dplyr::mutate(df,sample=stringr::str_extract(file,'SWT.'));
readr::write_tsv(df,'outputs/combined_metaphlan2_table.txt');
"
R_bind_script="${R_bind_script//%/$metaphlan2_dir}"
R_bind_script=`echo $R_bind_script | tr '\n' ' '`
Rscript -e "$R_bind_script"

# Run MIDAS (species and genes)
mkdir -p ${midas_dir}
export PYTHONPATH=~/work/tools/MIDAS:$PYTHONPATH
export PATH=~/work/tools/MIDAS/bin/:$PATH
run_midas_on_sample_species(){
	prefix=$1
	~/work/tools/MIDAS/scripts/run_midas.py species \
	${midas_dir}/${prefix} \
	-1 ${host_removed_dir}/${prefix}1.fq.gz \
	-2 ${host_removed_dir}/${prefix}2.fq.gz \
	-t 4
}
export -f run_midas_on_sample_species
parallel -j 3 -a prefixes.txt run_midas_on_sample_species {}
python ~/work/tools/MIDAS/scripts/merge_midas.py species outputs -t dir -i ${midas_dir}

run_midas_on_sample_genes(){
	prefix=$1
	~/work/tools/MIDAS/scripts/run_midas.py genes \
	${midas_dir}/${prefix} \
	-1 ${host_removed_dir}/${prefix}1.fq.gz \
	-2 ${host_removed_dir}/${prefix}2.fq.gz \
	-t 4
}
export -f run_midas_on_sample_genes
parallel -j 3 -a prefixes.txt run_midas_on_sample_genes {}


# Run humann2 (skip metaphlan2 step and use outputs we have already)
mkdir -p ${humann2_dir}
run_humann2_on_sample(){ 
	prefix=$1
	mytmp=`mktemp -d`
	cat ${host_removed_dir}/${prefix}1.fq.gz ${host_removed_dir}/${prefix}2.fq.gz > ${mytmp}/input.fq.gz

	humann2 -v -i ${mytmp}/input.fq.gz \
	--o-log ${humann2_dir}/${prefix}_log.txt \
	--threads 8 -o ${humann2_dir}/ \
	--output-basename ${prefix} \
	--taxonomic-profile ${metaphlan2_dir}/${prefix}_metaphlan2_output.txt

	rm -rf ${humann2_dir}/${prefix}_humann2_temp/tmp*
	pigz -p 8 ${humann2_dir}/${prefix}_humann2_temp/*tsv ${humann2_dir}/${prefix}_humann2_temp/*txt
	rm -rf ${humann2_dir}/${prefix}_humann2_temp/*bt2
	rm -rf ${humann2_dir}/${prefix}_humann2_temp/*fa
	rm -rf ${humann2_dir}/${prefix}_humann2_temp/*ffn
}
export -f run_humann2_on_sample
parallel -j 2 -a prefixes.txt run_humann2_on_sample {}
run_humann2_unstrat(){
	prefix=$1
	humann2_split_stratified_table --input ${humann2_dir}/${prefix}_genefamilies.tsv --output ${humann2_dir}/split_genefamilies
	humann2_split_stratified_table --input ${humann2_dir}/${prefix}_pathabundance.tsv --output ${humann2_dir}/split_pathabundance
	humann2_split_stratified_table --input ${humann2_dir}/${prefix}_pathcoverage.tsv --output ${humann2_dir}/split_pathcoverage
}
export -f run_humann2_unstrat
parallel -j 6 -a prefixes.txt run_humann2_unstrat {}


# Make a bunch of humann2 annotated outputs
mkdir -p outputs/humann2

humann2_join_tables --input ${humann2_dir} --file_name "genefamilies" --output outputs/humann2/humann2_default_stratified_genefamilies.txt
humann2_join_tables --input ${humann2_dir} --file_name "pathabundance" --output outputs/humann2/humann2_default_stratified_pathabundance.txt
humann2_join_tables --input ${humann2_dir} --file_name "pathcoverage" --output outputs/humann2/humann2_default_stratified_pathcoverage.txt

humann2_join_tables --input ${humann2_dir}/split_genefamilies --file_name "unstratified" --output outputs/humann2/humann2_default_unstrat_genefamilies.txt
humann2_join_tables --input ${humann2_dir}/split_pathabundance --file_name "unstratified" --output outputs/humann2/humann2_default_unstrat_pathabundance.txt
humann2_join_tables --input ${humann2_dir}/split_pathcoverage --file_name "unstratified" --output outputs/humann2/humann2_default_unstrat_pathcoverage.txt


humann2_regroup_table --input outputs/humann2/humann2_default_stratified_genefamilies.txt --groups uniref90_ko --output outputs/humann2/humann2_ko_stratified_genefamilies.txt
humann2_regroup_table --input outputs/humann2/humann2_default_unstrat_genefamilies.txt --groups uniref90_ko --output outputs/humann2/humann2_ko_unstrat_genefamilies.txt

humann2_regroup_table --input outputs/humann2/humann2_default_stratified_genefamilies.txt --groups uniref90_eggnog --output outputs/humann2/humann2_eggnog_stratified_genefamilies.txt
humann2_regroup_table --input outputs/humann2/humann2_default_unstrat_genefamilies.txt --groups uniref90_eggnog --output outputs/humann2/humann2_eggnog_unstrat_genefamilies.txt

humann2_regroup_table --input outputs/humann2/humann2_default_stratified_genefamilies.txt --groups uniref90_infogo1000 --output outputs/humann2/humann2_infogo_stratified_genefamilies.txt
humann2_regroup_table --input outputs/humann2/humann2_default_unstrat_genefamilies.txt --groups uniref90_infogo1000 --output outputs/humann2/humann2_infogo_unstrat_genefamilies.txt




